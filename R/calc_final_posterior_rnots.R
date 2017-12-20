#########################################
## Script to calculate final posterior Rnots
#########################################
rm(list=ls())

library(MASS)
library(tidyverse)
library(stringr)
library(lubridate)


base_url <- "zikaEstimatoR"
if(grepl('spencerfox', Sys.info()['login'])) setwd(file.path("~", "projects", base_url))
if(grepl('vagrant', Sys.info()['user'])) setwd( file.path("/vagrant", base_url) )
if(grepl('stampede', Sys.info()['nodename'])) setwd(file.path('/home1/02958/sjf826', base_url))
if(grepl('wrangler', Sys.info()['nodename'])) setwd(file.path('/home/02958/sjf826', base_url))

sapply(c("R/fitting_fxns.R", "R/mcmc_sampling.R"), source)

######################################################
## Define parameters for run
set.seed(1023902)

## Uncomment to generate all
# reporting_rate <- c(0.01, 0.0282, 0.0574, 0.0866, 0.1, 0.2)
include_trans <- c(NA, "1", "5")
temperature <- c("historic", "actual")
extra_import <- c(TRUE, FALSE)

# trans_ind <- 2; temp_ind <- 2; imp_ind <- 2

for(trans_ind in seq_along(include_trans)){
  for(temp_ind in seq_along(temperature)){
    for(imp_ind in seq_along(extra_import)){
      print(paste0(include_trans[trans_ind], " ", temperature[temp_ind], " ", extra_import[imp_ind]))
      daily_parms <- get_mcmc_parm_list(include_trans[trans_ind], temperature[temp_ind], last_only=TRUE, extra_imports = extra_import[imp_ind])

      ## Get posterior distribution -- only returns posterior for county/months that had importations
      posterior <- purrr::map(daily_parms, mcmc_zika_rnot,
                              alpha_tuning = .1,
                              rnot_tuning = .1,
                              rr_tuning = .1,
                              burnin = 100000,
                              N = 200000,
                              thin=10)


      est_posterior <- posterior[[1]]$samples %>% as_data_frame() %>% select(-1)

      colnames(est_posterior) <- c("alpha", "reporting_rate", daily_parms[[1]]$county_month_year)


      if(any(duplicated(colnames(est_posterior), fromLast = TRUE))){
        # cameron months may be duplicated if it's scenario with secondary transmission,
        # so remove one from each (they have the same values, so just remove either of them)
        est_posterior <- est_posterior[,-which(duplicated(colnames(est_posterior), fromLast = TRUE))]
      }


      ## Function for drawing samples and returning county R0 estimate
      sample_rnot_alpha <- function(alphas, rnots, size_vec = 10000){
        sample(alphas,size = size_vec, replace = T) * sample(rnots, size=size_vec, replace = T)
      }


      if(temperature[temp_ind]=="historic"){
        load("data_produced/county_r0_historic_dists.rda")
        county_prior_r0s <- county_r0_historic_dists
      } else{
        load("data_produced/county_r0_actual_dists.rda")
        county_prior_r0s <- county_r0_actual_dists
      }

      ## Now draw samples from the prior distribution and posterior alpha to fill in the county R0s that didn't experience importation
      all_counties <- paste0(county_prior_r0s$county, "_", county_prior_r0s$month, "_", county_prior_r0s$year)
      for(cty_mnth_yr in all_counties){
        if(!cty_mnth_yr %in% colnames(est_posterior)){
          ## For counties not already there, first extract the prior R0 samples
          county_rnots <- as.numeric(county_prior_r0s[which(cty_mnth_yr==all_counties), -c(1,2,3)])

          ## Now sample rnot and alpha from distributions
          est_posterior[cty_mnth_yr] <- sample_rnot_alpha(alphas = est_posterior$alpha, rnots = county_rnots)
        }
      }

      est_posterior <- est_posterior %>% gather(county, rnot_samp, 3:ncol(est_posterior)) %>%
        separate(col = county, c("county", "month", "year"), sep="_")


      if(grepl('spencerfox', Sys.info()['login'])) {
        save(est_posterior, file = paste0("data_produced/posterior_estimates/county_posterior_rnots_", temperature[temp_ind], "_",
                                          ifelse(is.na(include_trans[trans_ind]), 0, include_trans[trans_ind]),"_",
                                          ifelse(extra_import[imp_ind], "true", "false"), ".rda"))
      } else if(grepl('wrangler', Sys.info()['nodename'])){
        save(est_posterior, file = file.path("..","workfolder","data","ZikaEstimatoR_data", "post_rnot_est",
                                             paste0("county_posterior_rnots_", temperature[temp_ind], "_",
                                                    ifelse(is.na(include_trans[trans_ind]), 0, include_trans[trans_ind]),"_",
                                                    ifelse(extra_import[imp_ind], "true", "false"), ".rda")))
      }
    }
  }
}


# #####################################################
# ## Debugging the posterior sampling
#
# include_trans = "1" ; temperature = "actual"; extra_import = TRUE
# imp_ind <- trans_ind <- temp_ind <- 1
#
# daily_parms <- get_mcmc_parm_list(include_trans[trans_ind], temperature[temp_ind], last_only=TRUE, extra_imports = extra_import[imp_ind])
#
# ## Get posterior distribution -- only returns posterior for county/months that had importations
# set.seed(6)
# posterior <- purrr::map(daily_parms, mcmc_zika_rnot,
#                         alpha_tuning = .1,
#                         rnot_tuning = .1,
#                         rr_tuning = .1,
#                         burnin = 100000,
#                         N = 200000,
#                         thin=10)
#
# quantile(posterior[[1]]$samples[,2], probs=c(0.0275,0.5,0.975))
# plot(posterior[[1]]$samples[,2], type = "l")
# plot(posterior[[1]]$samples[,3], type = "l")

#####################################################
## Generate posterior if using 2016 cameron temperatures -- doesn't really make sense anymore
# load("data_produced/calculated_cam_county_2016_rnots.rda")
#
# ## Substitute in 2016 temperature estimates for November
# cam_2016_rnot$county <- "cameron"
# cam_county_r0_dist <- county_r0_distributions
# cam_county_r0_dist[which(cam_county_r0_dist$county == "cameron" & cam_county_r0_dist$month == "Nov"), ] <- cam_2016_rnot
#
# daily_parms <- unique(tx_data$notification_date)[length(unique(tx_data$notification_date))] %>%
#   purrr::map(~get_alpha_parms_r0_dist_mcmc(tx_data, curr_date=.x, county_r0_dists = cam_county_r0_dist, reporting_rate=as.numeric(reporting_rate)))
#
# ## Get posterior distribution -- only returns posterior for county/months that had importations
# cam_est_posterior <- purrr::map(daily_parms, mcmc_zika_rnot,
#                             alpha_tuning = .1,
#                             rnot_tuning = .1,
#                             disp_df = dispersion_df,
#                             burnin = 100000,
#                             N = 200000,
#                             thin=10)
#
# cam_est_posterior <- cam_est_posterior[[1]]$samples %>% as_data_frame() %>% select(-1)
#
# colnames(cam_est_posterior) <- c("alpha", daily_parms[[1]]$county_month)
#
# # cameron months where importations happened are duplicated. in this case, remove the last instance of them
# # Necessary to remove last instance, because first instance isn't subject to constraints of the importation
# cam_est_posterior <- cam_est_posterior[,-which(duplicated(colnames(cam_est_posterior), fromLast = TRUE))]
#
# cam_est_posterior <- cam_est_posterior %>% gather(county, rnot_samp, 2:ncol(cam_est_posterior)) %>%
#   separate(col = county, c("county", "month"), sep="_") %>%
#   filter(county == "cameron", month == "Nov") %>%
#   select(county,month, rnot_samp) %>%
#   mutate(county = "cameron2016")
#
# save(cam_est_posterior, file = paste0("data_produced/posterior_estimates/cam2016_posterior_rnots_",
#                                       ifelse(is.na(include_trans), 0, include_trans), "_", reporting_rate,".rda"))
#
