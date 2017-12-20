rm(list=ls())

#####################################################
## Reads in the arguments from command line running
## Parameters required to be specified:
##    extra_imports - should be 0 (false) or 1 (true)
##    include_trans - should be NA, 1, or 5 depending on how many cases of secondary transmission are desired
##    temperature - should be historic or actual
#####################################################
args <- (commandArgs(TRUE)) ## load arguments from R CMD BATCH

if(length(args)>0)  { ## Then cycle through each element of the list and evaluate the expressions.
  print(paste0('loading in ', args, ' from R CMD BATCH'))
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}
library(MASS)
library(tidyverse)
library(stringr)
library(lubridate)
# library(bbmle)

base_url <- "zikaEstimatoR"
if(grepl('spencerfox', Sys.info()['login'])) setwd(file.path("~", "projects", base_url))
if(grepl('vagrant', Sys.info()['user'])) setwd( file.path("/vagrant", base_url) )
if(grepl('stampede', Sys.info()['nodename'])) setwd(file.path('/home1/02958/sjf826', base_url))
if(grepl('wrangler', Sys.info()['nodename'])) setwd(file.path('/home/02958/sjf826', base_url))

sapply(c("R/fitting_fxns.R", "R/mcmc_sampling.R"), source)

## Run this line if running locally on Rstudio, without command line parameters
# include_trans = 1; reporting_rate = 0.0574; temperature = "actual"

## Get the parameter lists as specified by parm arguments
daily_parms <- get_mcmc_parm_list(include_trans, temperature, last_only=FALSE, extra_imports = ifelse(extra_imports==0, FALSE, TRUE))
# load("data_produced/dispersion_df.rda") # No longer needed

##########################################################
## Run the MCMC for each day and combine into a single data frame holding all the alpha data
##########################################################

est_alphas <- purrr::map(daily_parms, mcmc_zika_rnot,
                         alpha_tuning = .1,
                         rnot_tuning = .1,
			                   rr_tuning = .1,
                         burnin = 100000,
                         N = 200000,
                         thin=10)

est_alphas <- est_alphas %>% transpose()
est_alphas_df <- est_alphas[[1]] %>% purrr::map(as_data_frame) %>% purrr::map(function(x) select(x, 2)) %>% bind_cols()
colnames(est_alphas_df) <- daily_parms %>% purrr::map(~.$date) %>% do.call("c", .)

save(est_alphas_df, file = file.path("..","workfolder","data","ZikaEstimatoR_data", "est_rr_post",
                                     paste0("alpha_daily_mcmc_", temperature, "_",
                                     ifelse(is.na(include_trans), 0, include_trans),
                                     ifelse(extra_imports==0, "_false", "_true"),
                                     ".rda")))




