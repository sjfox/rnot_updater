#################################################
## Functions for fitting R0s
#################################################
require(Rcpp)

sourceCpp("cpp/cpp_fitting_fxns.cpp")

zika_parms <- function(rnot = 1.1,
                       num_intros = 1,
                       distribution = "pois",
                       overdispersion=1,
                       date = NA,
                       rnot_dist=NA,
                       reporting_rate = 1,
                       secondary_trans = 0,
                       county_month_year = NA,
                       inform_prior=TRUE){
  return(as.list(environment()))
}

get_gamma_parms <- function(rnots){
  require(MASS)
  fit <- try(fitdistr(as.numeric(unlist(rnots)), "gamma"), silent = T)
  if(class(fit) == "try-error"){
    # browser()
    list(estimate=c(shape=mean(as.numeric(unlist(rnots))), rate=1))
  } else{
    fit
  }
}


get_alpha_parms_r0_dist_mcmc <- function(tx_data, curr_date, county_r0_dists){
  tx_data <- tx_data %>% filter( notification_date <= curr_date) %>%
    group_by(county, month, sec_trans, year) %>%
    summarize(imports = n())

  rnot_data <- left_join(tx_data, county_r0_dists, by = c("county", "month", "year")) %>%
    ungroup() %>%
    select(starts_with("V"))


  num_rnots <- nrow(rnot_data)
  gamma_parms <- data_frame(shape = rep(0, num_rnots), rate = rep(0,num_rnots))
  for(ind in 1:nrow(rnot_data)){
    fit_gamma <- get_gamma_parms(rnot_data[ind,])
    gamma_parms[ind,] <- c(fit_gamma$estimate["shape"], fit_gamma$estimate["rate"])
  }
  subs_parms(list(rnot = NA,
                  rnot_dist = gamma_parms,
                  num_intros = tx_data$imports,
                  distribution="nbinom",
                  date=curr_date,
                  secondary_trans = tx_data$sec_trans,
                  county_month_year = paste0(tx_data$county, "_", tx_data$month, "_", tx_data$year)), zika_parms())
}

get_mcmc_parm_list <- function(include_trans,  temperature, last_only=FALSE, extra_imports = FALSE){
  ## Function gives the full parm list for running mcmc on every single date of importation
  ## Can call this function with specified parms, and get a list with elements ready-to-go for mcmc
  tx_imports <- read_csv("data/Zika Disease Cases as of 09282017.csv")
  # browser()
  tx_imports <- tx_imports %>% mutate(notification_date = mdy(`First Notification Date`)) %>%
    arrange(notification_date)%>%
    mutate(month = as.character(month(notification_date, label=TRUE, abbr = T)),
           county = tolower(str_replace_all(County, " County", ""))) %>%
    select(county, notification_date, month) %>%
    filter(year(notification_date) == 2016)

  if(extra_imports){
    ## If extra imports scenario, assume you've only seen 0.0574 % of imports, so replicate each row 17 times
    tx_imports <- tx_imports[rep(1:nrow(tx_imports), round(1/.0574)),]
  }

  tx_imports$sec_trans <- 0
  include_trans <- as.numeric(include_trans)
  if(!is.na(include_trans)){
    tx_imports$sec_trans[which(tx_imports$notification_date== "2016-11-21" & tx_imports$county == "cameron")[1]] <- include_trans
    tx_imports$sec_trans[which(tx_imports$notification_date== "2016-12-12" & tx_imports$county == "cameron")[1]] <- 1
  }

  if(temperature=="historic"){
    load("data_produced/county_r0_historic_dists.rda")
    tx_data <- tx_imports  %>% mutate(month = factor(month, levels = month.abb), year = 1960)
    county_r0s <- county_r0_historic_dists
  } else if(temperature == "actual"){
    load("data_produced/county_r0_actual_dists.rda")
    tx_data <- tx_imports  %>% mutate(month = factor(month, levels = month.abb), year = year(notification_date))
    county_r0s <- county_r0_actual_dists
  } else{
    stop("incorrect specification of temperature")
  }

  if(last_only){
    unique(tx_data$notification_date)[length(unique(tx_data$notification_date))] %>%
      purrr::map(~get_alpha_parms_r0_dist_mcmc(tx_data, curr_date=.x, county_r0_dists = county_r0s))
  } else{
    unique(tx_data$notification_date) %>%
      purrr::map(~get_alpha_parms_r0_dist_mcmc(tx_data, curr_date=.x, county_r0_dists = county_r0s))
  }

}


# scaling_loglike <- function(alpha, parms, disp_df){
#   ## Returns the Negative log likelihood for set of parameters
#
#   if(!is.na(parms$rnot)){
#     rnots <- parms$rnot * alpha
#     # ods <- unlist(purrr::map(rnots, ~find_overdispersion(.x)))
#     ods <- find_rnot_ods(rnots, disp_df)
#     parms <- subs_parms(list(rnot=rnots, overdispersion=ods), parms)
#     -sum(intro_loglike(parms))
#   } else{
#     rnot_dist <- parms$rnot_dist * alpha
#     log_likes <- vector("numeric", ncol(rnot_dist))
#     ods <- matrix(find_rnot_ods(unlist(rnot_dist), disp_df), nrow(rnot_dist))
#     for(col in 1:ncol(rnot_dist)){
#       parms$rnot <- rnot_dist[,col]
#       parms$overdispersion <- ods[,col]
#       log_likes[col] <- -sum(intro_loglike(parms))
#     }
#     parms$rnot <- NA
#     return(median(log_likes))
#   }
# }
#
#
#
# get_alpha_likes <- function(parms, disp_df){
#   # Returns likelihood values for a variety of alphas, so that distributions can be calculated post-hoc
#   alphas <- seq(0, 1, length.out=1000)
#   nllikes <- unlist(purrr::map(alphas, ~scaling_loglike(., parms=parms, disp_df)))
#   likes <- exp(-nllikes)
#
#   df <- data_frame(alpha = alphas, likelihood = likes)
#   colnames(df)[2] <- as.character(parms$date)
#   df
# }
#
# ############################################################################
# ## Fitting the Rnot distribution
# ############################################################################
#
# get_secondary_above_20 <- function(rnot){
#   # Takes in an rnot value and returns the probability of seeing
#   # > 20 secondary cases from that rnot
#   p1 <- c(0.425806451612903, 0.8458765530605259)
#   p2 <- c(4.341935483870967, 3.297197366921235)
#
#   slope <- (p2[2] - p1[2]) / (p2[1] - p1[1])
#   yint <- - slope * p1[1] + p1[2]
#   if(rnot < yint){
#     # warning("R0 is low and returning zero") # Happens very often, so not worth warning
#     return(0)
#   }
#
#   prob <- (rnot - yint) / slope / 100
#   if(prob > 1){
#     return(1)
#   }
#   return(prob)
# }
#
# find_overdispersion <- function(rnot){
#   # Find the overdispersion parameter for a given R0
#   prob_above <- get_secondary_above_20(rnot)
#
#   compare_ps <- function(x, prob_above, rnot){
#     pnbinom(q = 20, mu = rnot, size = x, lower.tail = FALSE) - prob_above
#   }
#   # print(rnot)
#   if(prob_above == 0){
#     ## If Rnot is very low
#     if(rnot==0) {
#       return(1e-16)
#     }
#     # browser()
#     seq_lower <- seq(1e-16, 0.5,length.out=1000)
#     low <- -1
#     prob_above <- 1e-5
#     ps <- compare_ps(seq_lower, prob_above, rnot)
#     while(max(ps, na.rm=T) < 0){
#       prob_above <- prob_above/10
#       ps <- compare_ps(seq_lower, prob_above, rnot)
#     }
#     low <- seq_lower[which.max(ps)]
#     # print(rnot)
#     # browser()
#     overdisp <- uniroot(f = compare_ps, interval = c(low, 1),  rnot=rnot, prob_above= prob_above)
#   } else {
#     if(prob_above >= (1-1e-4) ){
#       ## If rnot is very large
#       seq_lower <- seq(0,100,length.out=1000)
#       overdisp = list(root = seq_lower[which(abs(compare_ps(seq_lower, prob_above, rnot)) <= 1e-06)[1]])
#     } else{
#       seq_lower <- seq(0,100,length.out=10000)
#       ps <- compare_ps(seq_lower, prob_above, rnot)
#
#       if(all(diff(ps) > 0)){
#         ## If Rnot is large but not very large
#         overdisp <- uniroot(f = compare_ps, interval = c(0, 100), rnot=rnot, prob_above=prob_above)
#       } else{
#         ## If Rnot isn't miniscule, but is small (~0.5-1.5)
#         ## Begin the search from the first positive difference.
#         max_p <- which.max(ps)
#         # browser()
#         overdisp <- try(uniroot(f = compare_ps, interval = c(0, seq_lower[max_p]), rnot=rnot, prob_above=prob_above), silent = TRUE)
#         if(class(overdisp) == "try-error"){
#           overdisp <- uniroot(f = compare_ps, interval = c(seq_lower[max_p], 100), rnot=rnot, prob_above=prob_above)
#         }
#       }
#     }
#
#   }
#
#   overdisp$root
# }

#################################################
## No longer used
#################################################
# normalize_vector <- function(values){
#   ## Normalizes a vector to sum to 1
#   values/sum(values,na.rm = T)
# }

# intro_like <- function(parms) {
#   ## Calculates the likelihood based on an R0 and number of introductions
#   switch(parms$distribution,
#          pois = dpois(x = 0, lambda = parms$rnot)^parms$num_intros,
#          nbinom = dnbinom(x = 0, mu = parms$rnot, size = parms$overdispersion)^parms$num_intros)
# }

# intro_loglike <- function(parms) {
#   ## Calculates the likelihood based on an R0 and number of introductions
#   switch(parms$distribution,
#          pois = dpois(x = 0, lambda = parms$rnot, log = T)*parms$num_intros,
#          nbinom = dnbinom(x = 0, mu = parms$rnot, size = parms$overdispersion, log=T)*parms$num_intros)
# }
# lprior <- function(parms){
#   dnorm(parms$rnot, mean = parms$prior_mu, sd = parms$prior_sd, log = T)
# }

# subs_parms <- function(sub_parms=NULL,
#                        ref_parms) {
#   within(ref_parms, {
#
#     for(nm in names(sub_parms)) {
#       assign(nm, sub_parms[[nm]])
#     }
#     rm(nm)
#   })
# }

# find_rnot_ods <- function(rnot, disp_df){
#   ## Takes in a vector of rnot values and returns
#   ## The dispersions for nbinom distribution
#   ## disp_df must be sorted data table, with one column as rnot and other as ods
#   ## Can be obtained by running the calc_dispersion_table.R script
#   disp_df[J(rnot), roll = "nearest"]$ods
# }

# get_rnot_ll_ci <- function(alpha, num_intros, distribution, overdispersion=1, rnots=NULL) {
#   ## Returns the median and % confidence interval for likelihood of rnot based
#   ## On number of introductions alone
#
#   if(num_intros==0){
#     warning("With 0 Introductions you have no information")
#     return(data.frame(low = NA, median = NA,  high = NA))
#   }
#   if(length(num_intros)!=1){
#     stop("Need to send single introduction number")
#   }
#
#
#   if(!is.null(rnots)){
#     n <- length(rnots)
#     ods <- overdispersion
#   } else {
#     n <- 10000
#     max_rnot <- 10
#     rnots <- seq(0, max_rnot, length.out = n)
#     ods <- unlist(purrr::map(rnots, ~find_overdispersion(.x)))
#   }
#
#   parms <- subs_parms(list(rnot=rnots, num_intros=num_intros, distribution=distribution, overdispersion=ods), zika_parms())
#   likelihoods <- intro_like(parms)
#
#   # Find index for the low, mle and high
#   low_ind <- 1
#   mle_ind <- which.max(likelihoods)
#   high_ind <- which(likelihoods < alpha)[1]
#   if(length(high_ind)==0){
#     high_ind <- n
#   }
#
#   data.frame(low = rnots[low_ind], mle = rnots[mle_ind],  high = rnots[high_ind])
# }



# get_alpha_ci <- function(parms, disp_df, sig_level=0.01){
#   # For a set of parameters, finds the possible alphas
#   alphas <- seq(0,1, length.out = 5000)
#   nllikes <- unlist(purrr::map(alphas, ~scaling_loglike(., parms=parms, disp_df)))
#   likes <- exp(-nllikes)
#
#   ## Extract the largest alpha that fulfills
#   high <- alphas[rev(which(likes > sig_level))[1]]
#
#   data_frame(mle=0, low=0, high=high)
# }


# llike_prior <- function(rnot, ref_parms){
#   parms <- subs_parms(c(rnot=rnot), ref_parms)
#   log(intro_like(parms)) + lprior(parms)
# }



# rnot_mcmc <- function(parms,
#                       rand_init=T,
#                       iters,
#                       tuning,
#                       burnin,
#                       thin = 10){
#   if(iters %% thin != 0 ){
#     stop("Thin needs to be a multiple of iters")
#   }
#
#   samples <- matrix(nrow=iters/thin, ncol=2)
#   curr_rnot <- runif(n = 1, min = 0, max = 2)
#
#   curr_lik <- llike_prior(curr_rnot, parms)
#
#   accept <- 0
#   for(ii in 1:(iters+burnin)){
#     prop_rnot <- exp(rnorm(1, mean = log(curr_rnot), sd = tuning))
#     prop_lik <- llike_prior(prop_rnot, parms)
#
#     lmh <- prop_lik - curr_lik
#
#     if ( (lmh >= 0) | (runif(n = 1,min = 0, max = 1) <= exp(lmh)) ) {
#       curr_rnot <- prop_rnot
#       accept <- accept + 1
#       curr_lik <- prop_lik
#     }
#     if(ii>burnin & (ii-burnin) %% thin==0) {
#       samples[(ii-burnin)/thin,] <- c(curr_rnot, curr_lik)
#     }
#   }
#   aratio <- accept/(iters+burnin)
#   colnames(samples) <- c("r_not", "ll")
#   samples <- as.mcmc(samples)
#   return(list(samples, aratio = aratio))
# }

# library(Rcpp)
# sourceCpp("cpp/cpp_fitting_fxns.cpp")
#
#
#
#
# find_overdispersion_R(1.5)
#
# all_equal(compare_ps(seq(0,0.5,length.out=1000), 0.01470807, 12),
#           compare_ps_R(seq(0,0.5,length.out=1000), 0.01470807, 12))



# library(Rcpp)
#
# sourceCpp("cpp/cpp_fitting_fxns.cpp")

