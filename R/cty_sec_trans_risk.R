######################################
## Script getting the expected number of secondary cases
## The process works as such
## 1. The importations are summarized by county, month, and year to get total importations for that specific period
## 2. The posterior distribution is also grouped as such, and these two data frames are combined
## 3. The posterior R0s are then sampled according to how many importations were in the specific period - this sample R0s
##    for the importations that occurred. From these R0s we randomly draw a sample for the expected number of secondary
##    cases for those importations. Summing all of the expected number of secondary cases gives the number of autochthonous
##    cases we would've expected from the number of importations.
## 4. The sample is
######################################
rm(list=ls())

library(tidyverse)
library(lubridate)
library(stringr)
source("R/fitting_fxns.R")


load("data_produced/posterior_estimates/county_posterior_rnots_actual_1_0.0574.rda")

tx_imports <- read_csv("data/Zika Disease Cases by Notification Date as of 030617.csv")

tx_imports <- tx_imports %>% mutate(notification_date = mdy(`First Notification Date`)) %>%
  mutate(month = as.character(month(notification_date, label=TRUE, abbr = T)),
         county = tolower(str_replace_all(County, " County", "")),
         year = year(notification_date)) %>%
  select(county, notification_date, month, year) %>%
  group_by(county,month, year) %>%
  summarise(total_imports = n()) %>%
  ungroup()

post_summary <- est_posterior %>% group_by(county, month, year) %>%
  summarise(rnot_samp = mean(rnot_samp, na.rm=T)) %>% ungroup()


group_df <- function(df, groups_used){
  ## Function that returns a dataframe that is properly grouped
  ## groups_to_use should be character vector
  if(is_null(groups_used)){
    return(df)
  }
  df %>%
    group_by_(.dots = groups_used)
}

sample_df <- function(df, nsamples, n){
  ## Generates a sample of size n expected secondary cases from a dataframe
  ## Takes in a data frame and two values
  ## nsamples - the number of R0s to randomly select each sample (used for number imports from that county-month-year)
  ## n - the number of samples to actually take (~1000)
  ## dataframe must have two columns:
  ##    rnot_samp - column containing a sample from the rnot distribution (should also work if just one number)
  ##    ods - dispersion parameter for the R0, can be identified from running calc_dispersion_table.R (~0.12)
  ## Returns a list of size n, where each element is a column of length nsamples
  ## The values in the returned column are a sample from the expected secondary distribution for the R0
  ## This can be thought of as for a given county-month-year, the number of expected secondary cases for the number of importations
  samps <- vector(mode = "list", length = n)
  for(i in 1:n){
    samps[[i]] <- df %>% sample_n(size=nsamples, replace = T) %>%
                    mutate(samps = rnbinom(n=nsamples, mu = rnot_samp, size = ods)) %>%
                    select(samps)
  }
  samps
}

rename_cols <- function(df)  {
  ## Helper function to rename columns of a data frame generated from sample_df process
  colnames(df) = 1:ncol(df)
  df
}


ctymnthyr_group_imports <- group_df(tx_imports, groups_used = c("month","county", "year")) %>%
  summarise(total_imports = sum(total_imports)) %>% mutate(year = as.character(year))

ctymnthyr_grouping <- est_posterior %>%
  group_by(month, county, year) %>%
  mutate(ods = 0.12) %>%
  nest(rnot_samp, ods) %>%
  left_join(ctymnthyr_group_imports, by=c("county","month", "year")) %>%
  filter(!is.na(total_imports)) %>%
  rowwise() %>%
  do(samples = sample_df(df = .$data, nsamples = .$total_imports, n = 1000))

expected_cases <- ctymnthyr_grouping$samples %>% purrr::map(~bind_cols(.)) %>% purrr::map(~rename_cols(.)) %>%
  bind_rows()

expected_cases <- as.numeric(colSums(expected_cases))

quantile(expected_cases, probs = c(0.025, 0.5, 0.975))
# New:
# 2.5%     50%   97.5%
# 50.000  93.000 157.025
# Old:
# 2.5%      50%    97.5%
# 37      67.10342 114

## Take expected cases, and sample from binomial with reporting rate to get distribution
tot_samps <- length(expected_cases)*1000
exp_detected_cases <- data_frame(Low = rbinom(n = tot_samps, size = expected_cases, prob = 0.0574),
                                  High = rbinom(n = tot_samps, size = round(expected_cases/0.0574), prob = 0.0574)) %>%
  gather(estimate, det_cases, Low:High)

save(exp_detected_cases, file = "data_produced/exp_detected_cases.rda")

########-- The value above is the total number of secondary cases based on the importations
#### The following puts together a dataframe with the expected secondary cases from a single importation in each
#### county-month-year
####

################# Getting the prob of secondary transmission data from r0 data
get_prob_sec_trans <- function(rnots, dispersion_df){
  ods <- 0.12
  1 - dnbinom(0, mu = rnots, size = ods) ## Probability of secondary transmission is 1 - probability of no transmission
}
est_posterior <- est_posterior %>%
                    group_by(county, month, year) %>%
                    mutate(prob_sec_trans = get_prob_sec_trans(rnot_samp, dispersion_df))

prob_sec <- est_posterior %>%
  summarise(med = quantile(prob_sec_trans, probs = 0.5, na.rm=T),
            mean_prob = mean(prob_sec_trans, na.rm=T)) %>%
  mutate(month_num = match(month,month.abb))
prob_sec
save(prob_sec, file = "data_produced/posterior_prob_sec_trans.rda")



### NOT USED

# exp_cases <- data_frame(reporting_rate = 0.0574,
#            exp_sec_cases = c(no_group_mean_rnots, samps, samps2)*0.0574,
#            granularity = rep(c("none", "month", "county_month"), each=10000))
#
# save(exp_cases, file = "data_produced/exp_sec_cases.rda")

# get_df_sample_by_group <- function(df){
#   ## Function to take samples of a dataframe based on
#   ## The  grouping of that data frame
#   ## Number samples of each group is determined by the importation dataframe
#   ## So group columns in df need to be same as import_df
#   df %>% mutate(samples = map2(data, total_imports, sample_n, replace=T)) %>%
#     unnest(samples)
# }
#
# get_n_samples_by_group <- function(df, import_df, group_vec, n){
#   grouped_imports <- import_df %>%
#     group_df(groups_used =group_vec) %>%
#     summarise(total_imports = sum(total_imports, na.rm=TRUE))
#
#   if(is_null(group_vec)){
#     df <- df %>%
#       nest(1:3) %>% mutate(total_imports = grouped_imports[1,]) %>%
#       unnest(total_imports)
#   }else{
#     df <- df %>%
#       group_df(groups_used = group_vec) %>%
#       nest() %>%
#       left_join(grouped_imports) %>%
#       filter(!is.na(total_imports))
#   }
#
#   samps <- vector(mode = "list", length = n)
#   for(i in 1:n){
#     samps[[i]] <- get_df_sample_by_group(df)
#   }
#   samps
# }
# summary_rnot_vec <- function(l_samps){
#   l_samps %>% purrr::map(~summarise(.data = .x, rnot_sum=sum(.x$rnot_samp))) %>%
#     unlist() %>% as.numeric()
# }
# no_groups <- sample_df(est_posterior, 321, 10000)
# no_group_mean_rnots <- summary_rnot_vec(no_groups)
#
# mnth_group_imports <- group_df(tx_imports, groups_used = "month") %>%
#   summarise(total_imports = sum(total_imports))
# mnth_grouping <- est_posterior %>%
#   group_by(month) %>%
#   nest(rnot_samp) %>%
#   left_join(mnth_group_imports, by="month") %>%
#   rowwise() %>%
#   do(samples = sample_df(df = .$data, nsamples = .$total_imports, n = 10000))
#
# samps <- mnth_grouping$samples %>% purrr::map(~bind_cols(.)) %>% purrr::map(~rename_cols(.)) %>%
#   bind_rows()
#
# samps <- as.numeric(colSums(samps))
