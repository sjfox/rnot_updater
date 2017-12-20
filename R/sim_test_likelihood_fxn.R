library(tidyverse)
library(cowplot)
library(Rcpp)

##################################################
## Simulation functions
##################################################
sim_outbreak <- function(R0, k=.12){
  ## Simulates outbreaks assuming a negative binomial distribution and branching process
  if(R0 >1){
    stop("R0 > 1, so won't run because may cause infinite loop")
  }
  currently_infected <- 1
  cases <- 1
  while(currently_infected > 0){
    new_infected <- rnbinom(currently_infected, mu = R0, size = k)
    currently_infected <- sum(new_infected)
    cases <-  cases + currently_infected
  }
  cases
}

sim_detection <- function(cases_vec, reporting_rate){
  ## Function to simulate the detection process on a vector from the simulated outbreaks
  rbinom(n = length(cases_vec), size = cases_vec, prob = reporting_rate)
}

sim_known_import_detection <- function(cases_vec, reporting_rate){
  ## Function to simulate the detection process on a vector from the simulated outbreaks
  rbinom(n = length(cases_vec), size = cases_vec - 1, prob = reporting_rate) + 1
}

sim_known_import_detection_no_guarantee <- function(cases_vec, reporting_rate){
  ## Function to simulate the detection process on a vector from the simulated outbreaks
  outbrk_detect <- if_else(rbinom(n = length(cases_vec), size = 1, prob = reporting_rate) == 1, TRUE, FALSE)
  cases_vec[outbrk_detect] <- rbinom(n = length(cases_vec[outbrk_detect]), size = cases_vec[outbrk_detect] - 1, prob = reporting_rate) + 1
  cases_vec[!outbrk_detect] <- 0
  cases_vec
}


##################################################
## Analytical Calculation functions
##################################################
prob_offspring_perf_det <- function(R0, j, k){
  dnbinom(j-1, mu = R0, size = k)
}

prob_outbreak_perf_det <- function(R0, j, k){
  exp(lgamma(k*j+j-1)-lgamma(k*j)-lgamma(j+1)+(j-1)*log(R0/k)-(k*j+j-1)*log(1+R0/k))
}

prob_outbreak_imperf_det <- function(R0, jj, k, reporting_rate){
  # jj is the detected number of cases in a chain
  # j_true becomes the possible true number of cases in a chain that gave rise to j detected cases

  num_calc = 1e4
  if(R0 == 0){
    j = 1:num_calc
    true_chain_pdf <- vector(mode = "numeric", length = num_calc)[1]
    true_chain_pdf[1] <- 1
  } else{
    j = 1:num_calc
    # log_real_chain_pdf = lgamma(k*j+j-1)-lgamma(k*j)-lgamma(j+1)+(j-1)*log(r0/k)-(k*j+j-1)*log(1+r0/k)
    true_chain_pdf = prob_outbreak_perf_det(R0, j, k)
  }
  if(jj ==0){
    return(NA)
  }
  prob0 = sum(exp(j*log(1-reporting_rate)+log(true_chain_pdf)))
  denominator = 1 - prob0

  l = jj:num_calc
  numerator = sum(exp(log(true_chain_pdf[l]) +  dbinom(x= jj, size = l, prob = reporting_rate, log = T)))

  numerator/denominator
}

prob_outbreak_imperf_det_import <- function(R0, jj, k, reporting_rate){
  # jj is the detected number of cases in a chain
  # the j vector becomes the possible true number of cases in a chain that gave rise to j detected cases
  # browser()
  num_calc = 1e4
  if(R0 == 0){
    j = 1:num_calc
    true_chain_pdf <- vector(mode = "numeric", length = num_calc)[1]
    true_chain_pdf[1] <- 1
  } else{
    j = 1:num_calc
    # log_real_chain_pdf = lgamma(k*j+j-1)-lgamma(k*j)-lgamma(j+1)+(j-1)*log(r0/k)-(k*j+j-1)*log(1+r0/k)
    true_chain_pdf = prob_outbreak_perf_det(R0, j, k)
  }
  if(jj == 0){
    return(NA)
  }
  prob0 = 0
  # prob0 = sum(exp(j*log(1-reporting_rate)+log(true_chain_pdf)))
  denominator = 1 - prob0

  l = jj:num_calc
  numerator = sum( exp(log(true_chain_pdf[l]) + dbinom(x = jj-1, size = l-1, prob = reporting_rate, log=T) ))

  numerator/denominator
}

rnot = 0.7
j = 0:100
k=0.12
rr = 0.05

# 10000 %>% rerun(sim_outbreak(rnot, k)) %>% unlist() -> outbreak_sizes
outbreak_sizes <- replicate(10000, sim_outbreak(rnot, k))



## Calculate the analytical expectation for the outbreak size probabilities
size_predictions <- data_frame(obs_size = j) %>%
                      mutate(`Perfect` = map(.x = obs_size, .f = prob_outbreak_perf_det, R0=(rnot*rr), k=k) %>% unlist(),
                             `Imperfect` = map(.x = obs_size, .f = prob_outbreak_imperf_det, R0 =rnot, k=k, reporting_rate = rr) %>% unlist(),
                             `Imperfect Import` = map(.x = obs_size, .f = prob_outbreak_imperf_det_import, R0 =rnot, k=k, reporting_rate = rr) %>% unlist()) %>%
  gather(key, value, 2:4) %>%
  mutate(key = factor(key, levels = c("Perfect", "Imperfect", "Imperfect Import")))

data_frame(`Offspring` = outbreak_sizes,
           `Perfect` = outbreak_sizes,
           `Imperfect` = sim_detection(outbreak_sizes, rr),
           `Imperfect Import` = sim_known_import_detection(outbreak_sizes, rr)) %>%
  gather(key,value, 2:4) %>%
  filter(value !=0) %>%
  mutate(key = factor(key, levels = c("Perfect", "Imperfect", "Imperfect Import"))) %>%
  ggplot(aes(value)) +
  facet_wrap(~key) +
  geom_histogram(binwidth=1, aes(y=..density..)) +
  coord_cartesian(xlim = c(0,20)) +
  labs(x= "Outbreak Size", y = "Probability Density")+
  geom_point(data = size_predictions, aes ( x = obs_size, y = value), color = "red", size = 2) -> sim_v_analytic_plot
sim_v_analytic_plot

save_plot("ms_figs/sfigs/sfx_sim_v_analytic_plot.png", sim_v_analytic_plot, base_height = 5, base_aspect_ratio = 3)




########################################################
## Make figure for supplemental across many Ro and j
rnot = seq(0.1, 0.95, length = 5)
rr = 0.0574
k=0.12
nreps=10000



analytic_calcs <- as_data_frame(expand.grid(R0=rnot, jj=1:10)) %>%
  mutate(reporting_rate = rr) %>%
  mutate(`Imperfect Import` = pmap(., .f = prob_outbreak_imperf_det_import, k=k) %>% unlist(),
         `Imperfect` = pmap(., .f = prob_outbreak_imperf_det, k=k) %>% unlist(),
         `Perfect` = map2(.x = R0, .y = jj, .f = prob_outbreak_perf_det, k=k) %>% unlist()) %>%
  gather(sim, value, 4:6)

all_sims <- as_data_frame(expand.grid(R0=rnot, jj=1:10)) %>%
  mutate(reporting_rate = rr) %>%
  mutate(jj = map(R0, .f = function(x) replicate(nreps, sim_outbreak(x, k ))))

all_sims %>% unnest() %>%
  mutate(`Perfect` = jj,
         `Imperfect` = sim_detection(jj, unique(reporting_rate)),
         `Imperfect Import` = sim_known_import_detection(jj, unique(reporting_rate))) -> sim_data

sim_data %>%
  gather(sim, value, 4:6) %>%
  filter(value !=0) %>%
  mutate(sim = factor(sim, levels = c("Perfect", "Imperfect", "Imperfect Import"))) %>%
  ggplot(aes(value)) +
    facet_grid(R0~sim) +
    geom_histogram(binwidth=1, aes(y=..density..)) +
    coord_cartesian(xlim = c(0.5,10.5), ylim=c(0,1.05), expand=F) +
    scale_x_continuous(breaks = 1:10) +
    labs(x= "Outbreak Size", y = "Probability Density") +
    theme(strip.background = element_rect(fill = NA))+
    panel_border(colour = "black") +
    geom_point(data = analytic_calcs, aes ( x = jj, y = value), color = "red", size = 2) -> sim_v_analytic_plot2

## Beware this takes a long time to run
save_plot("ms_figs/sfigs/sfx_sim_v_analytic_plot2.png", sim_v_analytic_plot2, base_height = 10, base_aspect_ratio = 0.9)


##########################################
## Double check the calculations and comparisons
optimize_fn <- function(R0, p_val, fxn, exp_it=FALSE, ...){
  if(exp_it){
    abs(p_val - exp(fxn(R0, ...))^10)
  } else{
    abs(p_val - fxn(R0, ...)^10)
  }
}

## Find maximum R0 for p<0.05 for offspring
optimise(f = optimize_fn, interval = c(0.001,5), p_val = 0.05, fxn = prob_offspring_perf_det, j=1, k =.12)
## Perfect detection
optimise(f = optimize_fn, interval = c(0.001,5), p_val = 0.05, fxn = prob_outbreak_perf_det, j=1, k =.12)
## Imperfect outbreak detection
optimise(f = optimize_fn, interval = c(0.001,5), p_val = 0.05, fxn = prob_outbreak_imperf_det, j=1, k =.12, reporting_rate = 0.05)
## Imperfect outbreak detection no zeroes
optimise(f = optimize_fn, interval = c(0.001,5), p_val = 0.05, fxn = prob_outbreak_imperf_det_import, j=1, k =.12, reporting_rate = .05)
optimise(f = optimize_fn, interval = c(0.001,5), p_val = 0.05, fxn = offspring_size_obs_llike, sec_trans=0, ods =.12, reporting_rate = .05, exp_it=TRUE)



rnot = c(0.05, 0.8, 1.5)
j = 1:5
k=0.12
rr = 0.5




expand.grid(R0 = rnot, j = j) %>%
  mutate(`Perfect` = pmap(.l = ., .f = prob_outbreak_perf_det, k=k) %>% unlist(),
         `Imperfect` = pmap(.l=., .f = prob_outbreak_imperf_det, k=k, reporting_rate = rr) %>% unlist(),
         `Imperfect Import` = pmap(.l = ., .f = prob_outbreak_imperf_det_import, k=k, reporting_rate = rr) %>% unlist()) %>%
  gather(key, value, `Perfect`:`Imperfect Import`) %>%
  ggplot(aes(j, value, color = key)) +
    facet_wrap(~R0) +
    geom_line() +
    coord_cartesian(expand=F) +
    panel_border()





# ########################################################
# ## Make figure for supplemental across many Ro and j
# rnot = seq(0.1, 0.95, length = 10)
# rr = seq(0.05, 0.5, length = 10)
# k=0.12
# nreps=10000
#
#
#
# analytic_calcs <- expand(data_frame(R0=rnot,reporting_rate=rr, jj = 1:10), R0, reporting_rate, jj ) %>%
#   mutate(`Imperfect Import` = pmap(., .f = prob_outbreak_imperf_det_import, k=k) %>% unlist(),
#          `Imperfect` = pmap(., .f = prob_outbreak_imperf_det, k=k) %>% unlist(),
#          `Perfect` = map2(.x = R0, .y = jj, .f = prob_outbreak_perf_det, k=k) %>% unlist()) %>%
#   gather(sim, value, 4:6)
#
#
# all_sims <- expand(data_frame(R0=rnot,reporting_rate=rr), R0, reporting_rate) %>%
#   mutate(jj = map(R0, .f = function(x) replicate(nreps, sim_outbreak(x, k ))))
#
# all_sims %>% unnest() %>%
#   group_by(R0, reporting_rate, jj) %>%
#   summarize(prob = n()/nreps) %>%
#   filter(jj < 11) %>%
#   mutate(sim = "Perfect") -> perf_det
#
# ##Imperfect Detection
# all_sims %>% unnest() %>%
#   mutate(jj = sim_detection(jj, reporting_rate)) %>%
#   filter(jj !=0 ) %>%
#   group_by(R0, reporting_rate) %>%
#   nest() %>%
#   mutate(nrows = map(data, nrow) %>% unlist()) %>%
#   unnest() %>%
#   group_by(R0, reporting_rate, jj) %>%
#   summarize(prob = n()/unique(nrows)) %>%
#   filter(jj < 11) %>%
#   mutate(sim = "Imperfect")-> imperf_det
#
#
# all_sims %>% unnest() %>%
#   mutate(jj = sim_known_import_detection(jj, reporting_rate)) %>%
#   filter(jj !=0 ) %>%
#   group_by(R0, reporting_rate) %>%
#   nest() %>%
#   mutate(nrows = map(data, nrow) %>% unlist()) %>%
#   unnest() %>%
#   group_by(R0, reporting_rate, jj) %>%
#   summarize(prob = n()/unique(nrows)) %>%
#   filter(jj < 11) %>%
#   mutate(sim = "Imperfect Import")-> imperf_import_det
#
#
#
# bind_rows(list(perf_det, imperf_det, imperf_import_det)) %>%
#   left_join(analytic_calcs) %>%
#   group_by(reporting_rate, R0, sim) %>%
#   summarize(correlation = cor(prob, value)) %>%
#   ggplot(aes(R0, reporting_rate, fill = correlation)) +
#   geom_tile() +
#   facet_wrap(~sim) +
#   scale_fill_continuous(limits = c(0,1)) +
#   coord_cartesian(expand=F)
#
# bind_rows(list(perf_det, imperf_det, imperf_import_det)) %>%
#   left_join(analytic_calcs) %>%
#   group_by(reporting_rate, R0, sim) %>%
#   gather(key, values, prob, value) %>%
#   ggplot(aes(R0, values, color = key)) +
#   geom_point() +
#   facet_wrap(~sim) +
#   # scale_fill_continuous(limits = c(0,1)) +
#   coord_cartesian(expand=F)
