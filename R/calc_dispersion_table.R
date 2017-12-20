################################
## Calculate and save the overdispersion parameters for Rnots
## Note: No longer calculates and saves the table. Rather it now
## Fits a single dispersion parameter to be used for all R0s
################################
library(tidyverse)

fit_dispersion <- function(dispersion){
  R0 = seq(0, 5, length.out=100)

  p20.est = numeric(length(R0))
  for (i in 1:length(R0)) {
    p20.est[i] = ((R0[i] - 0.58) / 0.63) / 100
    if (R0[i] < 0.58) p20.est[i] = 1e-8
  }

  p20.nb = pnbinom(q = 20, mu = R0, size = dispersion, lower.tail = FALSE)
  sum((p20.est - p20.nb)^2)
}

fit_disp_param <- optimise(fit_dispersion, interval = c(0.0,5))
fit_disp_param
# Fit Dispersion parameter = 0.1202162

## Check it and looks good!
n = 100000
R0 = seq(0, 5, length.out=100)

p20.est = numeric(length(R0))
for (i in 1:length(R0)) {
  p20.est[i] = ((R0[i] - 0.58) / 0.63) / 100
  if (R0[i] < 0.58) p20.est[i] = 1e-8
}

p20.nb = numeric(length(R0))
for (i in 1:length(R0)) {
  p20.nb[i] = sum(rnbinom(n, mu=R0[i], size=fit_disp_param$minimum) > 20)/n
}


rnot_disp_estimates <- data_frame(R0, Assumed = p20.nb, Exact = p20.est)

save(file = "data_produced/rnot_disp_estimates.rda", rnot_disp_estimates)


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
#     ## If the rnot is below the y intercept,
#     ## find the overdispersion parameter for the yintercept instead
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
#   if(rnot == 0) {
#     return(1e-16)
#   }
#   if(prob_above == 0){
#     ## If Rnot is below 0.586 create relationships for rnot convert to ods:
#     ## x = rnot, y = returned od
#     ## Dispersion parameter ranges from 1e-16 (rnot=0) to 1e-5 (rnot=.586)
#     ## Those values used to make smoothish function
#     # p1 <- c(0, 1e-16)
#     # p2 <- c(0.592, 6e-5)
#     # slope <- (p2[2] - p1[2]) / (p2[1] - p1[1])
#     # yint <- - slope * p1[1] + p1[2]
#     # overdisp = list(root = slope * rnot + yint)
#     prob_above <- 1e-8
#   }
#   # } else {
#     # if(prob_above >= (1-1e-4) ){
#     #   ## If rnot is very large
#     #   seq_lower <- seq(0,100,length.out=1000)
#     #   overdisp = list(root = seq_lower[which(abs(compare_ps(seq_lower, prob_above, rnot)) <= 1e-06)[1]])
#     # } else{
#     max_int = 1
#     seq_lower <- seq(0,max_int,length.out=1000)
#     ps <- compare_ps(seq_lower, prob_above, rnot)
#
#     max_p <- which.max(ps)
#     overdisp <- try(uniroot(f = compare_ps, interval = c(seq_lower[max_p], max_int), rnot=rnot, prob_above=prob_above), silent = TRUE)
#     if(class(overdisp) == "try-error"){
#       overdisp <- try(uniroot(f = compare_ps, interval = c(0, seq_lower[max_p]), rnot=rnot, prob_above=prob_above), silent = TRUE)
#       if(class(overdisp) == "try-error"){
#         browser()
#       } else if(overdisp$root == 0){
#         browser()
#       }
#     }
#     # }
#
#   # }
#   overdisp$root
# }
#
# #################
# ## Setup rnots for analysis
# #################
# rnots <- seq(0, 20, length.out= 10000)
# ods <- unlist(purrr::map(rnots, ~find_overdispersion(.x)))
#
# dispersion_df <- data_frame(rnots = rnots, ods = ods)
# # dispersion_dt <- data.table(dispersion_df, val=rnots)
# # setattr(dispersion_dt, "sorted", "rnots")
#
# save(dispersion_df, file = "data_produced/dispersion_df.rda")

##################################
## Testing speed of various lookups
##################################
## Switched to using data_frame, since it's sorted, and runs quickly on new Rcpp version. faster than data.table
# find_nearest_basic <- function(rnot, disp_df){
#   find_nearest <- function(rnot, disp_df){
#     disp_df$ods[which.min(abs(disp_df$rnots - rnot))]
#   }
#   unlist(purrr::map(rnot, ~find_nearest(.x, disp_df)))
# }
#
# find_nearest_dt <- function(rnot,disp_df){
#   require(data.table)
#   dt <- data.table(disp_df, val = rnots)
#   setattr(dt, "sorted", "rnots")
#   dt[J(rnot), roll = "nearest"]$ods
# }
#
# find_nearest_already_dt <- function(rnot,disp_dt){
#   disp_dt[J(rnot), roll = "nearest"]$ods
# }
# dispersion_dt <- data.table(dispersion_df, val=rnots)
# setattr(dispersion_dt, "sorted", "rnots")
#
# find_nearest_basic(1:10, dispersion_df)
# find_nearest_dt(1:10, dispersion_df)
#
#
# rnot_samp <- runif(100, min=0.000001, max = 100)
# library(microbenchmark)
# mbm <- microbenchmark(find_nearest_basic(rnot_samp, dispersion_df),
#                find_nearest_dt(rnot_samp, dispersion_df),
#                find_nearest_already_dt(rnot_samp, dispersion_dt))
# mbm
