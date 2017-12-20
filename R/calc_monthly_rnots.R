## Generate R0s
library(tidyverse)
library(stringr)

sapply(c("R/r0_calc_fxns.R"), source)

########################################
## Texas county information for R0 calculation
########################################
tx_county <- read_csv(file = "data/county_r0_parameters.csv")


################# Historic Temperature calculation
########################################
## Temperature Data
########################################

tx_temps <- read_csv("data/tx_historic_temps.csv")
tx_temps <- tx_temps %>% mutate(month = factor(month, levels = month.abb)) %>%
  rename(county = subregion) %>%
  mutate(county = if_else(county=="de witt", "dewitt", county)) %>%
  spread(key = month, value = avg_temp)

########################################
## Perkins R0 functions with h=1
########################################

load("data_produced/vector_suitability/parms_fxns_r0.RData")


tx_county <- tx_county %>%
  left_join(tx_temps, by=c("county"))

tx_county <- tx_county %>% gather(key = "month", value = "avg_temperature", Jan:Dec)

## Only needed if never run before (next load fails)
county_r0_historic_dists <- rnot_calc_dist(tx_county$mosquito.abundance,
                 tx_county$gdp,
                 tx_county$avg_temperature,
                 a=a, b=b, c.r=c.r,
                 mort.fun.list=mort.fun,
                 eip.fun.list=eip.fun,
                 scam.est.list=scam.est.list)

county_r0_historic_dists <- county_r0_historic_dists %>% mutate(county = tx_county[["county"]], month=tx_county[["month"]], year = 1960) %>%
  select(county, year, month, everything())

save(county_r0_historic_dists, file = "data_produced/county_r0_historic_dists.rda")


load("data_produced/county_r0_historic_dists.rda")
tx_historic_county <- tx_county %>%
                mutate(low_r0 = apply(county_r0_historic_dists[,-(1:3)], 1, quantile, probs=c(0.025)),
                       med_r0 = apply(county_r0_historic_dists[,-(1:3)], 1, quantile, probs=c(0.5)),
                       high_r0 = apply(county_r0_historic_dists[,-(1:3)], 1, quantile, probs=c(0.975)))


tx_historic_county_rnots <- tx_historic_county %>% mutate(month = factor(month, levels = month.abb))

save(tx_historic_county_rnots, file = "data_produced/tx_county_historic_summary_rnots.rda")


#############################################################
################# Actual Temperature calculation from 2016/17 instead of historic
########################################
## Temperature Data
########################################

tx_temps <- read_csv("data/tx_actual_temps.csv")
tx_temps <- tx_temps %>% mutate(month = factor(month, levels = month.abb)) %>%
  rename(county = subregion) %>%
  mutate(county = if_else(county=="de witt", "dewitt", county)) %>%
  spread(key = month, value = avg_temp)

########################################
## Perkins R0 functions with h=1
########################################

load("data_produced/vector_suitability/parms_fxns_r0.RData")

tx_county <- read_csv(file = "data/county_r0_parameters.csv")

tx_county <- tx_county %>%
  left_join(tx_temps, by=c("county"))

tx_county <- tx_county %>% gather(key = "month", value = "avg_temperature", Jan:Dec)

## Only needed if never run before (next load fails)
county_r0_actual_dists <- rnot_calc_dist(tx_county$mosquito.abundance,
                                           tx_county$gdp,
                                           tx_county$avg_temperature,
                                           a=a, b=b, c.r=c.r,
                                           mort.fun.list=mort.fun,
                                           eip.fun.list=eip.fun,
                                           scam.est.list=scam.est.list)

county_r0_actual_dists <- county_r0_actual_dists %>% mutate(county = tx_county[["county"]], month=tx_county[["month"]], year=tx_county[["year"]]) %>%
  select(county, year, month, everything())

save(county_r0_actual_dists, file = "data_produced/county_r0_actual_dists.rda")


load("data_produced/county_r0_actual_dists.rda")
tx_actual_county <- tx_county %>%
  mutate(low_r0 = apply(county_r0_actual_dists[,-(1:3)], 1, quantile, probs=c(0.025), na.rm=T),
         med_r0 = apply(county_r0_actual_dists[,-(1:3)], 1, quantile, probs=c(0.5), na.rm=T),
         high_r0 = apply(county_r0_actual_dists[,-(1:3)], 1, quantile, probs=c(0.99), na.rm=T))


tx_actual_county_rnots <- tx_actual_county %>% mutate(month = factor(month, levels = month.abb))

save(tx_actual_county_rnots, file = "data_produced/tx_county_actual_summary_rnots.rda")











# #######################################
# ## Run for hypothetical known 2016 cameron temperature
# #######################################
# cam_2016_temp <- 22.944 # https://www.weather.gov/climate/index.php?wfo=bro
# cam_nov_data <- tx_county %>% filter(county=="cameron", month=="Nov")
# cam_nov_data$avg_temperature <- cam_2016_temp
#
# cam_2016_rnot <- rnot_calc_dist(cam_nov_data$mosquito.abundance,
#                cam_nov_data$gdp,
#                cam_nov_data$avg_temperature,
#                a=a, b=b, c.r=c.r,
#                mort.fun.list=mort.fun,
#                eip.fun.list=eip.fun,
#                scam.est.list=scam.est.list)
#
# cam_2016_rnot <- as.data.frame(t(cam_2016_rnot))
# rownames(cam_2016_rnot) <- NULL
# cam_2016_rnot <- cam_2016_rnot %>% mutate(county = "cameron2016", month="Nov") %>% select(county,month, V1:V1000)
#
# save(cam_2016_rnot, file = "data_produced/calculated_cam_county_2016_rnots.rda")


