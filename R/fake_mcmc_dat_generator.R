###############################################
## Generating Fake data for ms Figure 1 - MCMC
###############################################
rm(list=ls())
library(MASS)
library(tidyverse)
sapply(c("R/fitting_fxns.R", "R/mcmc_sampling.R"), source)

set.seed(235443)

load("data_produced/county_r0_actual_dists.rda")


##############################################
## Panel data
## First extracts R0 information for Aug/oct from Harris county
## Then solves for posterior distributions for Rnots in August for Harris according
##    To a certain number of importations (specified by introductions variable)
## Then extracts prior R0s from all counties for August and solves for posteriors
##    from these counties as well
## Combine all
## Then calculates the secondary case probabilities for inset
##############################################

## Fxn for setting up parameters for the mcmc
get_fake_parms <- function(dist, intros){
  parms <- vector("list", length = length(intros))
  fit_gamma <- get_gamma_parms(t(dist))
  gamma_parms <- data_frame(shape = fit_gamma$estimate["shape"], rate = fit_gamma$estimate["rate"])
  ii <- 1
  for(intro in intros){
    parms[[ii]] <- subs_parms(list(rnot = NA, rnot_dist = gamma_parms, num_intros = intro, distribution="nbinom"), zika_parms())
    ii <- ii + 1
  }
  parms
}

## These data are used for the August/October panels
months <- c("Aug", "Sep")
har_rnot_dists <- county_r0_actual_dists %>% filter(county=="harris", month %in% months, year == 2016) %>%
  gather(samp, rnot, V1:V1000) %>% mutate(month = factor(month, levels=months))


## First going to solve for posterior distributions
introductions <- c(15, round(15/0.0574))
aug_parms <- har_rnot_dists %>% filter(month =="Aug") %>% select(rnot) %>% get_fake_parms(intros = introductions)
parms <- c(aug_parms)


fake_alpha_mcmc <- parms %>% purrr::map(mcmc_zika_rnot,
                                        alpha_tuning = .1,
                                        rnot_tuning = .1,
                                        rr_tuning = 0.1,
                                        burnin = 100000,
                                        N = 200000,
                                        thin=100)

fake_alpha_df <- fake_alpha_mcmc %>% transpose()
fake_alpha_df <- fake_alpha_df[[1]] %>% purrr::map(as_data_frame) %>% purrr::map(function(x) select(x, 2:4)) %>% bind_rows()
colnames(fake_alpha_df) <- c("alpha", "reporting_rate", "Aug")
fake_alpha_df <- fake_alpha_df %>% mutate(intros = rep(introductions, each = 1000))


## Functions for scaling prior Rnots and getting the correct rnot distributions from distributions
scale_fake_rnots <- function(rnots, alphas){
  n <- 1000
  sample(rnots, size=n,replace=T) * sample(x = alphas, size = n, replace = T)
}

get_rnots <- function(month_needed, rnot_dists){
  # browser()
  rnot_dists %>% filter(month==month_needed, year==2016) %>%
    select(rnot) %>% unlist()
}

fake_alpha_df <- fake_alpha_df %>%
  group_by(intros) %>%
  mutate(Sep = scale_fake_rnots(get_rnots(months[2], har_rnot_dists), alpha)) %>%
  gather(month_scaled, scaled_rnot, Aug, Sep)

original_rnots <- har_rnot_dists %>% mutate(alpha = NA, reporting_rate = NA, month_scaled = month, scaled_rnot=rnot, intros=0) %>%
  select(-samp, -county, -rnot, -month, -year)

fake_alpha_dat <- bind_rows(fake_alpha_df, original_rnots) %>% ungroup() %>% mutate(intros = factor(intros, levels = c("0", introductions)))

########################################################
## Now get priors and posteriors for all countys
########################################################
## Convert county rnots into a usable format
county_rnots <- county_r0_actual_dists %>% filter(month=="Aug", year == 2016) %>%
  unite(county_month_year, county, month, year, sep = "_") %>%
  t() %>%
  as_tibble()
colnames(county_rnots) <- county_rnots[1,]
county_rnots <- county_rnots[-1,] %>% purrr::map(as.numeric) %>% bind_cols()

county_rnots_post <- vector(mode="list", length = length(introductions)+1)
for(intro in 1:(length(introductions)+1)){
  if(intro == (length(introductions)+1)){
    county_rnots_post[[(length(introductions)+1)]] <- county_rnots
  } else{
    alpha_samp <- fake_alpha_dat %>% filter(intros==introductions[intro], month_scaled=="Aug") %>% select(alpha) %>% unlist() %>% as.numeric()
    county_rnots_post[[intro]] <- mutate_all(county_rnots, funs(scale_fake_rnots), alphas=alpha_samp)
  }
}

county_fake_rnots <- county_rnots_post %>% purrr::map(gather) %>% bind_rows() %>%
  mutate(intros = rep(c(introductions, 0), each=254000)) %>%
  separate(key, into = c("county", "month", "year"), sep = "_") %>%
  rename(rnots = value)


# ##############################################
# ## Inset panel data
# ##############################################
# avg_secondary_prob <- function(num_secondary, dist, dispersion_df){
#   mean(dnbinom(x = num_secondary, mu = dist, size = find_rnot_ods(dist, dispersion_df)), na.rm=T)
# }
# har_rnot_dists %>% group_by(month) %>% summarise(`0` = avg_secondary_prob(0, rnot, dispersion_df),
#                                                  `1` = avg_secondary_prob(1, rnot, dispersion_df),
#                                                  `2` = avg_secondary_prob(2, rnot, dispersion_df),
#                                                  `3` = avg_secondary_prob(3, rnot, dispersion_df),
#                                                  `4` = avg_secondary_prob(4, rnot, dispersion_df),
#                                                  `5` = avg_secondary_prob(5, rnot, dispersion_df)) %>%
#   gather(secondary_cases, probability, `0`:`5`) %>%
#   mutate(secondary_cases = as.numeric(secondary_cases)) -> exp_secondary_cases


save(fake_alpha_dat, county_fake_rnots, file = "data_produced/fig2_data.rda")

