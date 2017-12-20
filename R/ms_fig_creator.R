##############################
## Plots for the Zika Paper
## Spencer Fox Jan 23, 2017
##############################
rm(list=ls())

library(cowplot)
library(tidyverse)
library(maps)
library(lubridate)
library(stringr)
library(scales)
library(ggridges)
library(ggrepel)
###################################################################################
## Fxn for adding proper color scaling/outline to maps
source("R/plotting_fxns.R")

## Date of results to be used.
data_produced_date <- "11-27-2017"
###################################################################################

###############################
## Conceptual figure 1
###############################
load("data_produced/fig2_data.rda")

## Prior/posterior densities
updated_rnot_plot <- fake_alpha_dat %>% filter(intros %in% c(0, 15)) %>%
  filter(month_scaled=="Aug") %>%
  mutate(intros = if_else(intros==0, "August 2016", "Future August")) %>%
  ggplot(aes(scaled_rnot, fill = as.factor(intros), color = as.factor(intros))) +
  geom_density(alpha=0.6, size=1) +
  #facet_wrap(~month_scaled, nrow = 1, scales = "free_y") +
  coord_cartesian(expand=F) +
  guides(fill = guide_legend(title = ""), color = guide_legend(title = ""))+
  scale_fill_manual(values = c("black", "gray70")) +
  scale_color_manual(values = c("black", "gray70")) +
  theme(strip.background = element_rect(fill=NA),
        strip.text = element_text(size=16), legend.position=c(0.5,0.8)) +
  labs(x = expression("R"[0]), y = "Density")
updated_rnot_plot



## Setup the maps for part b
fake_county_plot_dat <- county_fake_rnots %>%
  filter(intros %in% c(0,15)) %>%
  group_by(county, month,intros) %>%
  summarise(med_rnot = median(rnots))

city_dat <- data_frame(ID = c(1L, 2L, 3L, 4L, 5L),
         lat = c(32.77666, 29.76043, 29.42412, 30.26715, 31.76188),
         lon = c(-96.79699, -95.3698, -98.49363, -97.74306, -106.485),
        Name = c("Dallas", "Houston", "San Antonio", "Austin", "El Paso"))
city_dat$intros <- "Future August"

map_data(map = "county") %>% filter(region=="texas") %>%
  mutate(subregion = if_else(subregion=="de witt", "dewitt", subregion)) %>%
  left_join(fake_county_plot_dat, by=c("subregion" = "county"))  %>%
  mutate(month = factor(month, levels=month.abb),
         intros = if_else(intros==0, "August 2016", "Future August"),
         intros = factor(intros, levels = c("August 2016", "Future August"))) %>%
  ggplot(aes(x=long, y=lat, fill = med_rnot, group = subregion)) + facet_wrap(~intros, nrow=1)+
  geom_polygon(color = "gray", size=0.1) +
  theme_void(base_size=14) +
  geom_text_repel(data = city_dat, aes(x = lon, y = lat, label = Name), size=5,
                  box.padding = unit(x = .05, units = "npc"), inherit.aes = F) +
  theme(strip.text = element_text(size=16)) -> update_map_plot

update_map_plot <- add_map_scale(update_map_plot, max_rnot = max(fake_county_plot_dat %>% ungroup() %>% select(med_rnot), na.rm=T))
update_map_plot

combined_example_fig2 <- plot_grid(updated_rnot_plot, update_map_plot, nrow = 1,labels = "AUTO", rel_widths = c(1,2))

save_plot(filename = "ms_figs/f1_scaling_example.png", plot = combined_example_fig2, base_height = 4, base_aspect_ratio = 3)



###############################
## Raw R0 and importation data Figure 2
###############################
tx_imports <- read_csv("data/Zika Disease Cases as of 09282017.csv")

import_data <- tx_imports %>% mutate(notification_date = mdy(`First Notification Date`)) %>%
  #filter(notification_date<max_date) %>%
  mutate(month = as.character(month(notification_date, label=TRUE, abbr = T)),
         county = tolower(str_replace_all(County, " County", ""))) %>%
  select(county, notification_date, month) %>%
  group_by(notification_date) %>%
  summarise(num_imports = n()) %>%
  mutate(cum_imports = cumsum(num_imports))

arrow_data <- data_frame(date = ymd(c("2016-11-22", "2016-12-12", "2017-04-25")),
                         ystart = c(7,7,7),
                         yend = c(3,5, 1))


import_plot <- import_data %>%
  ggplot(aes(notification_date, num_imports)) +
  annotate("rect", xmin=ymd("2016-01-01"), xmax=ymd("2017-01-01"), ymin=0, ymax=Inf, alpha=0.1, fill="black") +
  annotate("rect", xmin=ymd("2017-01-01"), xmax=as.Date(as.POSIXct(Inf, origin = "2016-01-01")), ymin=0, ymax=Inf, alpha=0.2, fill="green") +
  annotate("text", x=ymd("2016-02-05"), y=6.6, color = "black", label = "Training", size=6) +
  annotate("text", x=ymd("2017-01-20"), y=6.6, color = "darkgreen", label = "Test", size=6) +
  geom_bar(stat="identity", width=2) +
  labs(y="Importations", x = NULL) +
  scale_x_date(labels = date_format("%b"), breaks=date_breaks("month"))+
  scale_y_continuous(expand=c(0,0))+
  geom_segment(data = arrow_data, aes(x =date, xend=date,y=ystart, yend=yend),
               arrow = arrow(length = unit(0.3,"cm")), size=1.5, color="#5869B1")

import_plot

load("data_produced/tx_county_actual_summary_rnots.rda")

ub_r0 <- tx_actual_county_rnots %>%
  group_by(month, year) %>%
  summarize(high_r0 = max(high_r0))

ub_r0 <-  import_data %>%
  mutate(month = month(notification_date, label = T, abbr = T),
         year = year(notification_date)) %>%
  left_join(ub_r0, by = c("month","year"))

prior_rnot_plot <- import_data %>%
  mutate(month = month(notification_date, label = T, abbr = T),
         year = year(notification_date)) %>%
  left_join(tx_actual_county_rnots, by = c("month","year")) %>%
  ggplot(aes(notification_date, med_r0, group = county)) +
  geom_line(alpha=0.3) +
  geom_line(data = ub_r0, aes(x = notification_date, y=high_r0), lty=2, inherit.aes=FALSE)+
  #geom_ribbon(aes(ymax=high, ymin=low), alpha=.2) +
  scale_x_date(labels = date_format("%b"), breaks=date_breaks("month")) +
  labs(y = expression("County R"[0]), x = NULL)
prior_rnot_plot

fig1 <- plot_grid(prior_rnot_plot, import_plot, nrow = 2, align="v", labels = "AUTO")
save_plot("ms_figs/f2_priorr0_import.png", fig1, base_height = 5, base_aspect_ratio = 2)



###########################################
## Figure 3 - Posterior R0s for each month
###########################################

f3_data <- get_posterior_data("actual", 1, FALSE)
f3_posterior_maps <- plot_monthly_post_rnots(f3_data, 0.5)
f3_posterior_maps


save_plot("ms_figs/f3_posterior_maps.png", plot = f3_posterior_maps, base_height = 5, base_aspect_ratio = 1.3)

###########################################
## Figure 4 - Changing Alpha through time
###########################################

alpha_f4_data <- get_posterior_data("actual", 1, FALSE, alpha=TRUE)

f4_alpha_ridge <- alpha_plot_fxn(alpha_f4_data) +
    annotate("point", y = ymd("2016-08-10"), x = 0.042, fill = "#5869B1", color = "#5869B1", shape=25, size=5) +
    annotate("point", y = ymd("2016-09-20"), x = 0.1, fill = "#5869B1", color = "#5869B1", shape=25, size=5)
f4_alpha_ridge
save_plot("ms_figs/f4_alpha_ridge.png", plot = f4_alpha_ridge, base_height = 8, base_aspect_ratio = 1)
# alpha_plot <- alpha_plot_fxn(post_alpha_ci, 0.0574, 1)
# alpha_plot
#
# save_plot("ms_figs/f4_alpha_notridge.png", plot = alpha_plot, base_height = 4, base_aspect_ratio = 1.8)

##########################################
## Figure 5
##########################################
# f5_data <- get_posterior_data("actual", 1, FALSE)
#
#
# f5_data %>% filter(year == 2016) %>%
#   mutate(prob_det = 0.0574*(1 - dnbinom(0, mu = rnot_samp, size = 0.12))) %>%
#   group_by(county, month) %>%
#   summarize(prob_det = mean(prob_det, na.rm=T)) %>%
#   mutate(month_num = match(month,month.abb)) -> f5_sum_data
#
# f5_prob_case <- f5_sum_data %>% ggplot(aes(month_num, prob_det, group = county)) +
#     geom_line(alpha=0.5) +
#     scale_x_continuous(breaks = 1:12, labels = month.abb) +
#     scale_y_continuous(expand = c(0,0))+
#     labs(y = "Probability Detect Secondary Case", x = "")


load("data_produced/all_expected_cases.rda")
all_expected_cases  %>%
  mutate(total_expected_cases = as.factor(total_expected_cases),
         extra_import = if_else(extra_import, "Increased Importations", "Reported Importations"),
         extra_import = factor(extra_import, levels = c("Reported Importations", "Increased Importations"))) %>%
  group_by(extra_import, total_expected_cases) %>%
  summarize(total_exp_case_dens = n()/10000) %>%
  mutate(total_expected_cases = (as.numeric(total_expected_cases) - 1),
         sum_exp_cases = sum(total_expected_cases*total_exp_case_dens, na.rm = T)) ->f5_data


f5_exp_case <- f5_data %>% ggplot(aes(x = as.numeric(total_expected_cases), y = total_exp_case_dens)) +
    facet_wrap(~extra_import, nrow=2) +
    geom_bar(stat="identity") +
    theme(strip.background = element_rect(fill=NA, color = "black"))+
    labs(x = "Expected Number of Cases", y="Probability") +
    scale_x_continuous(breaks = seq(0,10)) +
    scale_y_continuous(expand=c(0,0), limits=c(0,1)) +
    panel_border(colour = "black") +
    geom_vline(aes(xintercept=as.numeric(sum_exp_cases)), size=1) +
    geom_vline(xintercept=1, color = "#5869B1", lty=2, size=1) +
    coord_cartesian(xlim=c(0,10))
f5_exp_case

# f5_prob_exp_cases <- plot_grid(f5_prob_case, f5_exp_case, labels="AUTO", rel_widths = c(1.5,1))
# save_plot("ms_figs/f5_prob_exp_cases.png", plot = f5_prob_exp_cases, base_height = 4, base_aspect_ratio = 2.3)
save_plot("ms_figs/f5_exp_case_dist.png", plot = f5_exp_case, base_height = 5, base_aspect_ratio = 0.8)
###########################################################
## Supplementary Figures
###########################################################
tx_outline <- map_data(map="state") %>% filter(region=="texas")

###########################################################
## Supp Fig 1 - comparing predicted vs assumed dispersion stuff
###########################################################
load("data_produced/rnot_disp_estimates.rda")
rnot_disp_estimates %>%
  ggplot(aes(R0, Exact)) +
  geom_line() +
  geom_point(aes(y=Assumed), size=0.5) +
  labs(x = expression("R"[0]), y = "Probability of > 20 secondary cases") -> sf1_assumed_disp_plot

save_plot(filename = "ms_figs/sfigs/sf1_assumed_disp_plot.png", plot = sf1_assumed_disp_plot, base_height = 4, base_aspect_ratio = 1.1)

###########################################################
## Supp Fig 2 - Compare Prior median and upperbound estimates for State
###########################################################
load("data_produced/county_r0_actual_dists.rda")

med_ub_prior_rnots <- county_r0_actual_dists %>% filter(year == 2016) %>%
  gather(samp, value, -county, -year, -month) %>%
  group_by(county, month, year) %>%
  summarise(med_r0 = quantile(x = value, probs=0.5, na.rm=T),
            ub_r0 = quantile(x = value, probs=0.99, na.rm=T)) %>%
  gather(rnot_level, value, med_r0:ub_r0) %>%
  ungroup() %>%
  mutate(month = factor(month, levels=month.abb),
         rnot_level = factor(if_else(rnot_level == "med_r0", "Median", "99%"), levels = c("Median", "99%")))


map_data(map = "county") %>% filter(region=="texas") %>%
  mutate(subregion = if_else(subregion=="de witt", "dewitt", subregion)) %>%
  left_join(med_ub_prior_rnots, by=c("subregion" = "county"))  %>%
  ggplot(aes(x=long, y=lat, fill = value, group = subregion)) +
  facet_wrap(~month+rnot_level, labeller = label_wrap_gen(multi_line=FALSE), ncol = 6)+
  geom_polygon(color = "gray", size=0.1) +
  theme_void(base_size=14) -> sf2_med_ub_prior_rnot

sf2_med_ub_prior_rnot <- add_map_scale(sf2_med_ub_prior_rnot, max_rnot = max(med_ub_prior_rnots$value, na.rm=T))
sf2_med_ub_prior_rnot
save_plot("ms_figs/sfigs/sf2_med_ub_prior_rnot.png", plot = sf2_med_ub_prior_rnot, base_height = 7, base_aspect_ratio = 1.3)




###########################################################
## Supp Fig 3 - comparing the prior R0 estimates for the historic/actual temperature estimated R0s
###########################################################
load("data_produced/tx_county_historic_summary_rnots.rda")
load("data_produced/tx_county_actual_summary_rnots.rda")

## Get historic rnot data and join it with the actual temperature rnot data
historic_prior_rnots <- tx_historic_county_rnots %>% mutate(hist_rnot = med_r0) %>%
  select(month, county, hist_rnot)

tx_actual_county_rnots %>% mutate(act_rnot = med_r0) %>%
  filter(year == 2016) %>%
  select(month, county, act_rnot) %>%
  left_join(historic_prior_rnots, by = c("month", "county")) %>%
  filter(act_rnot >=0.001, hist_rnot >= 0.001) %>%
  ggplot(aes(hist_rnot, act_rnot, color = month)) +
    geom_point() +
    geom_abline(intercept = 0, slope=1, lty=2) +
    scale_x_log10(limits = c(0.001, 2)) +
    scale_y_log10(limits = c(0.001, 2)) +
    labs(x = expression("R"[0]~"Estimate - Historic Temperature"), y = expression("R"[0]~" Estimate - 2016 Temperature"), color = "") +
    scale_color_manual(values = rev(c("#023FA5", "#485CA8", "#6A76B2", "#878FBD", "#A1A6C8", "#B7BBD1", "#CBCDD9", "#DADBDF", "#E2E2E2"))) +
    theme(legend.position = c(0.8, 0.3)) -> sf3_prior_r0_comparison
save_plot(filename = "ms_figs/sfigs/sf3_prior_r0_comparison.png", plot = sf3_prior_r0_comparison, base_height = 5, base_aspect_ratio = 1)


###########################################################
## Supp Fig 4 - Median + ub for posterior R0 estimates
###########################################################

sf4_post_med_ub <- plot_med_ub_post_map("actual")
sf4_post_med_ub

save_plot("ms_figs/sfigs/sf4_post_med_ub.png", plot = sf4_post_med_ub, base_height = 5, base_aspect_ratio = 1.3)


###########################################################
## Supp Fig 5 - ub estimates for historic temperature estimation and the posterior using historic temperatures
###########################################################

sf5_historic_post_med_ub <- plot_med_ub_post_map("historic")
sf5_historic_post_med_ub

save_plot("ms_figs/sfigs/sf5_historic_post_med_ub.png", plot = sf5_historic_post_med_ub, base_height = 5, base_aspect_ratio = 1.3)


###########################################################
## Comparing secondary transmission of 0, 1, and 5 on R0 posterior results
###########################################################

get_sum_post <- function(posterior){
  posterior %>% filter(year==2016) %>%
    group_by(county, month) %>%
    summarize(median = median(rnot_samp),
              ub = quantile(rnot_samp, probs = 0.99))
}

no_trans <- get_posterior_data(temperature = "actual", include_trans = 0, extra_imports = FALSE, alpha = FALSE)
no_trans <- get_sum_post(no_trans)

trans1 <- get_posterior_data(temperature = "actual", include_trans = 1, extra_imports = FALSE, alpha = FALSE)
trans1 <- get_sum_post(trans1)

trans5 <- get_posterior_data(temperature = "actual", include_trans = 5, extra_imports = FALSE, alpha = FALSE)
trans5 <- get_sum_post(trans5)

increased_imports <- get_posterior_data(temperature = "actual", include_trans = 1, extra_imports = TRUE, alpha = FALSE)
increased_imports <- get_sum_post(increased_imports)

list(no_trans, trans1, trans5, increased_imports) %>%
  map2(.y = c("no_trans", "trans1", "trans5", "increased"), ~ mutate(.x, id = .y)) %>%
  bind_rows() -> sum_trans_post


sum_trans_post %>% gather(sum_stat, rnot, median,ub) %>%
  spread(key = id, value = rnot) %>%
  gather(key, value, no_trans, trans5, increased) %>%
  mutate(key = case_when(key == "trans5" ~ "5 Cases",
                         key == "no_trans" ~ "0 Cases",
                         key == "increased" ~ "Increased Importations"),
         sum_stat = if_else(sum_stat == "median", "Median", "Upperbound")) %>%
  ggplot(aes(trans1, value, color = key)) +
  facet_wrap(~sum_stat) +
  geom_point(alpha=1) +
  geom_abline(slope = 1, intercept =  0, lty = 2, size=1) +
  scale_color_manual(values = c("#717171", "black",   "#BBBBBB")) +
  labs(x = "Baseline Estimates", y = "Sensitivity Estimates", color = "") +
  theme(legend.position = c(0.05,0.8), strip.background = element_rect(fill=NA)) +
  coord_cartesian(xlim = c(0,1)) +
  background_grid()+
  panel_border(colour = "black") -> sf6_trans_comp_plot
sf6_trans_comp_plot

save_plot("ms_figs/sfigs/sf6_trans_comp_plot.png", plot = sf6_trans_comp_plot, base_height = 4, base_aspect_ratio = 2)

sum_trans_post %>%
  filter(id=="trans5") %>%
  summarize(max_ub=max(ub)) %>%
  filter(max_ub>1) %>%
  summarize(n())
#############################################################
# Calculating values needed for manuscript
#############################################################


est_posterior = get_posterior_data(temperature = "actual", include_trans = 5, extra_imports = FALSE, alpha = FALSE)

est_posterior %>% filter(year==2016 | year == 1960) %>%
  group_by(county,month) %>%
  summarize(median = median(rnot_samp), ub = quantile(rnot_samp, probs = 0.99)) -> posterior_summary

max(posterior_summary$median)
max(posterior_summary$ub)

hist(posterior_summary$median)
hist(posterior_summary$ub)

posterior_summary %>%
  summarize(max_ub = max(ub)) %>%
  filter(max_ub>1)

# calculating confidence intervals for alpha posterior --------------------
alpha_dat <- get_posterior_data("actual", 1, FALSE, alpha=TRUE)

quantile(alpha_dat$`2016-11-14`, probs= c(0.0275, 0.5, 0.975))
quantile(alpha_dat$`2016-12-29`, probs= c(0.0275, 0.5, 0.975))

# confidence intervals for expected cases ---------------------------------
load("data_produced/all_expected_cases.rda")
all_expected_cases %>%
  group_by(extra_import) %>%
  summarize(lb = quantile(total_expected_cases, probs = 0.0275),
            mean = mean(total_expected_cases),
            ub = quantile(total_expected_cases, probs = 0.975))

f3_data %>%
  filter(year==2016) %>%
  group_by(county, month) %>%
  summarize(ub = quantile(rnot_samp, .99, na.rm = T)) %>%
  summarize(max_rnot = max(ub)) %>%
  filter(max_rnot>1) -> high_risk_counties

high_risk_counties
## counties at high risk: grimes, houston, madison, montgomery, walker, waller
## population sizes: 27,671 + 22,754 + 13,987 + 556,203 + 71,484 + 50,115
## pop texas = 27862596
## data frome: https://www.tsl.texas.gov/ref/abouttx/census.html
# (27671 + 22754 + 13987 + 556203 + 71484 + 50115 ) / 27862596

# sum_trans_post defined above -- takes a long time to read in all data though so not rewriting ehre
sum_trans_post %>% filter(id == "trans5") %>%
  group_by(county) %>%
  summarize(max_ub = max(ub)) %>%
  filter(max_ub > 1) %>%
  nrow()


## Poster plots ##################################
posterior <- get_posterior_data("actual", 1, FALSE)
post_ub <- posterior %>% select(-alpha) %>%
  group_by(county, month) %>%
  filter(month=="Aug") %>%
  summarise(med_r0 = quantile(x = rnot_samp, probs=0.5, na.rm=T),
            ub_r0 = quantile(x = rnot_samp, probs=0.99, na.rm=T)) %>%
  gather(key, value, med_r0, ub_r0) %>%
  mutate(estimate = "Posterior")

load("data_produced/county_r0_actual_dists.rda")

prior_ub <- county_r0_actual_dists %>% gather(samp, val, V1:V1000) %>%
  filter(year==2016, month == "Aug") %>%
  group_by(county, month) %>%
  summarise(med_r0 = quantile(x = val, probs=0.5, na.rm=T),
            ub_r0 = quantile(x = val, probs=0.99, na.rm=T)) %>%
  gather(key, value, med_r0, ub_r0) %>%
  mutate(estimate = "Prior")

rnot_ests <- bind_rows(post_ub, prior_ub)



post_med_ub <- map_data(map = "county") %>% filter(region=="texas") %>%
  mutate(subregion = if_else(subregion=="de witt", "dewitt", subregion)) %>%
  left_join(rnot_ests, by=c("subregion" = "county"))  %>%
  mutate(key = factor(if_else(key == "med_r0", "Median", "99%"), levels = c("Median", "99%")),
         estimate = factor(estimate, levels = c("Prior", "Posterior"))) %>%
  ggplot(aes(x=long, y=lat, fill = value, group = subregion)) +
  facet_grid(estimate~key, labeller = label_wrap_gen(multi_line=FALSE)) +
  geom_polygon(color = "gray", size=0.1) +
  theme_void(base_size = 14)
aug_ests <- add_map_scale(post_med_ub, max_rnot = max(rnot_ests$value, na.rm=T))
save_plot("figs/rnot_ests.pdf", aug_ests, base_height = 3, base_aspect_ratio = 1.3)


# ##################################
# ## Fig for MIDAS draft
# ##################################
# load("data_produced/posterior_estimates/county_posterior_rnots_actual_1_false.rda")
#
# est_posterior %>% filter(month=="Aug", year == 2016) %>%
#   group_by(county) %>%
#   summarise(med_r0 = quantile(x = rnot_samp, probs=0.5, na.rm = T)) -> aug_post_rnot
#
# map_data(map = "county") %>% filter(region=="texas") %>%
#   mutate(subregion = if_else(subregion=="de witt", "dewitt", subregion)) %>%
#   left_join(aug_post_rnot, by=c("subregion" = "county"))  %>%
#   # mutate(month = factor(month, levels=month.abb),
#          # rnot_level = factor(if_else(rnot_level == "med_r0", "Median", "97.5%"), levels = c("Median", "97.5%"))) %>%
#   ggplot(aes(x=long, y=lat, fill = med_r0, group = subregion)) +
#   #facet_wrap(~rnot_level, labeller = label_wrap_gen(multi_line=FALSE), ncol = 6) +
#   geom_polygon(color = "gray", size=0.1) +
#   theme_void(base_size = 14) +
#   theme(legend.position = c(0.0, 0.21),
#         legend.title.align = 0.5) +
#   guides(fill = guide_colorbar(title = expression(R[0]), label.position = "bottom", title.position="top", direction = "horizontal", barwidth=10)) -> midas_aug_post_rnot
# midas_aug_post_rnot <- add_map_scale(midas_aug_post_rnot,
#                                      max_rnot = max(aug_post_rnot$med_r0, na.rm=T))
# library(rtZIKVrisk)
#
# rtzikvrisk_trigger_plot <- plot_fig4(panels = c("c"))
#
# midas_plot <- plot_grid(rtzikvrisk_trigger_plot, midas_aug_post_rnot, nrow=1, labels = "AUTO")
#
# save_plot("midas_zikv_plot.png", midas_plot, base_height = 5, base_aspect_ratio = 2)
############################################################################
## NOT USED ANYMORE
############################################################################
# ########################################
# ## Figure 5 - validation step and expected number of secondary cases
# ########################################
# load("data_produced/posterior_prob_sec_trans.rda")
# prob_detect <- prob_sec %>%
#   filter(year==2016) %>%
#   ggplot(aes(month_num, mean_prob*0.0574, group = county)) +
#   geom_line(alpha=0.5) +
#   scale_x_continuous(breaks = 1:12, labels = month.abb) +
#   scale_y_continuous(expand = c(0,0))+
#   labs(y = "Probability Detect Secondary Case", x = "")
# prob_detect
#
# load("data_produced/exp_detected_cases.rda")
# text_replace <- c("All Importations Reported", "Only Fraction Reported")
# exp_detected_cases %>%
#   mutate(estimate = if_else(estimate=="Low", text_replace[1], text_replace[2])) -> exp_detected_cases
#
# detected_case_plot <- exp_detected_cases %>%
#   ggplot(aes(det_cases)) +
#   geom_histogram(data = filter(exp_detected_cases, estimate==text_replace[1]), aes(y=..density..), binwidth=1, fill="black") +
#   geom_histogram(data = filter(exp_detected_cases, estimate==text_replace[2]), aes(y=..density..),binwidth=10, fill="black") +
#   facet_wrap(~estimate, nrow=2, scales = "free_y") +
#   theme(strip.background = element_blank(), strip.text = element_text(size=16)) +
#   coord_cartesian(expand=FALSE, xlim=c(0,210)) +
#   panel_border(colour = "black") +
#   labs(y="", x = "Expected detected autochthonous cases")
# detected_case_plot
#
# f5_exp_sec_cases <- plot_grid(prob_detect, detected_case_plot, nrow = 1, labels = "AUTO")
#
# save_plot("ms_figs/f5_exp_sec_cases.png", f5_exp_sec_cases, base_height = 5, base_aspect_ratio = 2.2)

############################
# Next fig should show the Cameron county R0 distribution
# Need to take scaling function and just give raw scaled results isntead of summaries
############################
# load("data_produced/cameron_scaled_rnot_dist.rda")
# est_posterior %>% filter(county == "willacy", month=="Nov") %>%
#   ggplot(aes(rnot_samp)) + geom_histogram()
#
# cam_nov_plot <- cam_scaled_rnots %>% group_by(month(date_predicted, label=TRUE)) %>%
#   filter(date_predicted==min(date_predicted), month(date_predicted, label=TRUE) %in% month.abb[6:11]) %>%
#   mutate(county = if_else(county=="cameron", "Historic", "Actual")) %>%
#   ggplot(aes(scaled_rnot, fill=county)) +
#   facet_wrap(~month(date_predicted, label=T), scales="free_y") +
#   geom_histogram(aes(y = ..density..), alpha=0.6, position="identity") +
#   scale_fill_manual(values = c("black", "#bdbdbd")) +
#   labs(x = expression("Scaled November Cameron R"[0])) +
#   panel_border(colour="black")+
#   theme(legend.position=c(0.8,0.85), legend.title = element_blank(), strip.background = element_rect(fill="NA"))
# cam_nov_plot
#
# save_plot("ms_figs/f4_cam_scaled_nov_rnot.png", cam_nov_plot, base_height = 4, base_aspect_ratio = 1.8)


############################
# Cameron County estimation plot
############################
# statewide_alpha_rnots %>%
#   filter(county=="cameron") %>%
#   ggplot(aes(info, med_r0*high)) + geom_point()+
#     geom_errorbar(aes(ymax = high_r0*high, ymin = low_r0*high)) +
#     scale_y_log10() +
#     geom_hline(yintercept=1, lty=2)
# load("data_produced/alpha_likelihoods/alpha_like_rnot_dist_0.05.rda")
# load("data_produced/scaled_rnots.rda")
#
# county_r0_distributions %>% filter(county=="cameron", month =="Nov",) %>% as.numeric() %>% as_data_frame() %>%
#   ggplot(aes(value))+geom_histogram() + labs(x="R0",title="November Predictions")
#
#
# county_r0_distributions %>% filter(county=="cameron", month =="Nov") %>% as.numeric() %>% as_data_frame() %>%
#   mutate(value = value*est_alphas_df$`2016-11-02`) %>%
#   ggplot(aes(value))+geom_histogram() + labs(x="R0",title="November Predictions Scaled")
#
# county_r0_distributions %>% filter(county=="cameron", month =="Oct") %>% as.numeric() %>% as_data_frame() %>%
#   ggplot(aes(value))+geom_histogram() + labs(x="R0",title="2016 Nov Temp Predictions")
#
#
# county_r0_distributions %>% filter(county=="cameron", month =="Oct") %>% as.numeric() %>% as_data_frame() %>%
#   mutate(value = value*est_alphas_df$`2016-11-02`) %>%
#   ggplot(aes(value))+geom_histogram() + labs(x="R0",title="2016 Nov Temp Predictions Scaled")
#
#
#
#
# cameron_data <- statewide_alpha_rnots %>%
#   filter(county=="cameron")
#
# nov_data <- cameron_data %>% group_by(month) %>% filter(month=="Nov", info == min(info)) %>% ungroup() %>% select(low_r0:high_r0)
# oct_data <- cameron_data %>% group_by(month) %>% filter(month=="Oct", info == min(info)) %>% ungroup() %>% select(low_r0:high_r0)
#
# cam_nov_predict <- cameron_data %>% mutate(nov_med_r0 = nov_data$'med_r0',
#                         nov_high_r0 = nov_data$'high_r0',
#                         nov_low_r0 = nov_data$'low_r0',
#                         oct_med_r0 = oct_data$'med_r0',
#                         oct_high_r0 = oct_data$'high_r0',
#                         oct_low_r0 = oct_data$'low_r0') %>%
#   mutate_each(funs(`*`(., high) ), nov_med_r0:oct_low_r0)
#
# oct_plot <- cam_nov_predict  %>%
#             group_by(month) %>%
#             filter(info == min(info))  %>%
#             ggplot(aes(info, oct_med_r0)) + geom_point(size=2) +
#               geom_errorbar(aes(ymax = oct_high_r0, ymin = oct_low_r0)) +
#               coord_cartesian(ylim = c(0,4)) +
#               geom_hline(yintercept=1, lty=2) +
#   labs(x = "Month", y = expression("Predicted October R"[0]))
#
# nov_plot <- cam_nov_predict  %>%
#   group_by(month) %>%
#   filter(info == min(info))  %>%
#   ggplot(aes(info, nov_med_r0)) + geom_point(size=2) +
#   geom_errorbar(aes(ymax = nov_high_r0, ymin = nov_low_r0)) +
#   coord_cartesian(ylim = c(0,4)) +
#   geom_hline(yintercept=1, lty=2) +
#   labs(x = "Month", y = expression("Predicted November R"[0]))
#
# cam_predicted_rnot <- plot_grid(oct_plot, nov_plot, labels = "AUTO")
#
# save_plot("ms_figs/predicted_rnots_cameron.png", cam_predicted_rnot, base_height = 4, base_aspect_ratio = 2.3)
#




########################### Random things for wilke presentation
# load("data_produced/county_r0_distributions.rda")
# county_r0_distributions %>% filter(county=="cameron", month =="Aug") %>% as.numeric() %>%
#   hist(main="", xlab = "R0")
#
# load("data_produced/calculated_tx_county_rnots.rda")
# map_data(map = "county") %>% filter(region=="texas") %>%
#   mutate(subregion = if_else(subregion=="de witt", "dewitt", subregion)) %>%
#   left_join(tx_county_rnots, by=c("subregion" = "county")) %>%
#   filter(month=="Aug") %>% gather(key,value, low_r0:high_r0) %>%
#   mutate(key = factor(key, levels = c("low_r0", "med_r0", "high_r0"))) %>%
#   ggplot(aes(x=long, y=lat, fill = value, group = subregion)) + facet_wrap(~key)+
#   geom_polygon(color = "gray", size=0.2) +
#   scale_fill_gradientn(name = expression("R"[0]), na.value = "white",
#                        colours = c("white", "blue", "yellow", "red"),
#                        # breaks= c(0.1, 1, 10),
#                        values = scales::rescale(c(0, 1, 1.000001, max(tx_county_rnots %>% select(high_r0)))),
#                        guide = guide_colorbar(title=expression("R"[0]), barheight=10)) +
#   theme_nothing() +theme(strip.background = element_blank(),
#                          strip.text.x = element_blank())
#
# hist(rnbinom(10000, mu = 2, size = 1), breaks = 30, right = T, main = "", xlab ="Secondary Cases")
# plot(seq(2,0.25,length.out=1000), dnbinom(x = 0, mu = seq(2,0.25,length.out=1000), size = 1)^10, type = "l", xlim=c(2,0.2), ylab = "Likelihood", xlab = "R0")
# abline(h = 0.05, lty=2)
#
# tx_county_rnots %>% ggplot(aes(month, med_r0, group = county))+geom_line(alpha=0.5) + labs(y="Median R0")

####################################
## First conceptual figure 1 attempt
#####################################
# num_secondary <- 0:5
# like_data <- data_frame(secondary_cases = num_secondary,
#             low = dnbinom(x = num_secondary, mu = 0.5, size = find_overdispersion(0.5)),
#             mid = dnbinom(x = num_secondary, mu = 2, size = find_overdispersion(2)),
#             high = dnbinom(x = num_secondary, mu = 4, size = find_overdispersion(4))) %>%
#           gather(rnot_level, single_intro, low:high) %>%
#           mutate(rnot_level = case_when(.$rnot_level == "low" ~ "0.5",
#                                 .$rnot_level == "mid" ~ "2",
#                                 .$rnot_level == "high" ~ "4"))
#
# secondary_dist_plot <- like_data %>%
#   ggplot(aes(secondary_cases, single_intro, fill=rnot_level)) +
#     geom_bar(stat = "identity", position="dodge", color="black") +
#     geom_hline(yintercept=0.05, lty=2) +
#     scale_x_continuous(breaks = c(0:10)) +
#     scale_y_continuous(expand=c(0,0))+
#     theme(strip.background = element_rect(fill="white"))+
#     scale_fill_brewer(type = "qual", palette = 2) +
#     labs(x = "Number of Secondary Cases", y = "Likelihood")+
#     theme(legend.position = c(0.75,0.75))+
#     guides(fill = guide_legend(title=expression("R"[0]), keywidth = 1.5, keyheight=1.5, direction = "horizontal", label.position = "bottom", label.hjust = 0.5, title.vjust = .75))
# secondary_dist_plot
#
# prob_zero_plot <- like_data %>% filter(secondary_cases==0) %>%
#   mutate(three_intro = single_intro^3, ten_intro = single_intro^10) %>%
#   gather(intro_num, like, three_intro, ten_intro) %>%
#   mutate(mult_case_prob = 1 - like) %>%
#   gather(prob_type, value, like, mult_case_prob) %>%
#   mutate(prob_type = if_else(prob_type=="like", "0", ">0"),
#          intro_num = if_else(intro_num == "ten_intro", "10 Importations", "3 Importations"),
#          prob_type = factor(prob_type, levels = c("0", ">0")),
#          intro_num = factor(intro_num, levels = c("3 Importations", "10 Importations"))) %>%
#   ggplot(aes(prob_type, value, fill = rnot_level)) +
#     facet_wrap(~intro_num)+
#     geom_bar(stat = "identity", position="dodge", color="black") +
#     geom_hline(yintercept=0.05, lty=2) +
#     scale_y_continuous(expand=c(0,0))+
#     theme(strip.background = element_rect(fill="white"))+
#     scale_fill_brewer(type = "qual", palette = 2) +
#     labs(x = "Total Number of Secondary Cases", y = "Likelihood")+
#     theme(legend.position = "none")
# prob_zero_plot
#
# import_example_plot <- plot_grid(secondary_dist_plot, prob_zero_plot, rel_widths = c(1,1), labels = "AUTO")
#
# save_plot("ms_figs/likelihood_ex.png", import_example_plot, base_height = 4, base_aspect_ratio = 2)

#
# ###########################
# # Negative binomial updating R0 figure - supplemental Figure
# ###########################
# load("data_produced/nb_fitod_estimates.rda")
#
# nb_fitod_plot <- nb_fitod_estimates %>% filter(alphas==0.01) %>%
#   mutate(high=ifelse(is.na(high), 10, high)) %>%
#   ggplot(aes(intros, mle)) + geom_point(size=2) +
#   geom_errorbar(aes(ymin=low, ymax=high))+
#   geom_hline(yintercept=1, lty=2) +
#   coord_cartesian(xlim= c(0,20)) +
#   labs(x = "Number of Importations", y = expression("Estimated R"[0]))
#
# save_plot("ms_figs/sf1_nbinom_update.png", nb_fitod_plot, base_height = 4, base_aspect_ratio = 1.8)
#
#
# ############################
# # Individual county scaling alpha plot
# ############################
#
# load("data_produced/county_alphas_rnots.rda")
#
# ## Select counties with highest importations and cameron county where loca transmission occurred
# import_counties <- tx_imports %>% group_by(county) %>%
#   summarize(imports = n()) %>% arrange(-imports) %>%
#   head(8) %>% select(county) %>% add_row(county="cameron")
#
# county_alpha_plot <- county_alphas %>% filter(county %in% import_counties$county) %>%
#   ggplot(aes(date, mle)) + geom_point() + facet_wrap(~county) +
#   geom_errorbar(aes(ymin=low, ymax=high)) +
#   labs(y = "Estimated Scaling Factor", x = "Date") +
#   panel_border()
#
# county_alpha_plot
#
# save_plot("figs/daily_county_alpha.pdf", county_alpha_plot, base_height = 7, base_aspect_ratio = 1.3)
#
#
# ############################
# # County alpha scaling R0 plot
# ############################

################################
## Old inset figure
# August\nImportations\nWithout\nTransmission ## this was earlier legend title
## Expected secondary cases for Importation Month -- inset to previous figure
# exp_secondary_cases %>% filter(month=="Aug") %>%
#   ggplot(aes(secondary_cases, probability)) +
#   geom_bar(stat="identity", position="dodge", fill = "#1b9e77", color = "#1b9e77", alpha=0.6) +
#   coord_cartesian(expand=F) +
#   theme(legend.position = "none", legend.title = element_blank(),plot.background = element_rect(size=1, color="grey", linetype = 1)) +
#   labs(x = expression("2"~degree*" Cases"), y = "Probability") -> secondary_dist_plot
# secondary_dist_plot

# binom_plot <- data_frame(x = 0:15) %>%
#   mutate(dens_low = dbinom(x, size = 37, prob = 0.0574),
#          dens_high = dbinom(x, size = 114, prob = 0.0574),
#          bar_color = if_else(x==2,"yes", "no")) %>%
#   gather(key, value, dens_low:dens_high) %>%
#   mutate(key = factor(if_else(key=="dens_low", "Low", "High"), levels = c("Low", "High"))) %>%
#   ggplot(aes(x,value, fill=bar_color)) + facet_wrap(~key, nrow=2) +
#   geom_histogram(stat="identity") +
#   scale_y_continuous(expand=c(0,0)) +
#   theme(strip.background = element_rect(fill = NA)) +
#   labs(y="Probability", x = "Autochthonous Cases Detected") +
#   scale_fill_manual(values = c("black", "grey80"), guide=FALSE)
# binom_plot
#
# est_ub_r0 <- est_posterior %>% select(-alpha) %>%
#   group_by(county, month) %>%
#   summarise(med_r0 = quantile(x = rnot_samp, probs=0.5),
#             ub_r0 = quantile(x = rnot_samp, probs=0.975))
#
# ub_rnot_plot <- map_data(map = "county") %>% filter(region=="texas") %>%
#   mutate(subregion = if_else(subregion=="de witt", "dewitt", subregion)) %>%
#   left_join(est_ub_r0, by=c("subregion" = "county"))  %>%
#   gather(rnot_level, rnots, med_r0, ub_r0) %>%
#   mutate(month = factor(month, levels=month.abb),
#          rnot_level = factor(if_else(rnot_level == "med_r0", "Median", "97.5%"), levels = c("Median", "97.5%"))) %>%
#   ggplot(aes(x=long, y=lat, fill = rnots, group = subregion)) +
#   facet_wrap(~month+rnot_level, labeller = label_wrap_gen(multi_line=FALSE), ncol = 6)+
#   geom_polygon(color = "gray", size=0.1) +
#   theme_nothing() +
#   scale_fill_gradientn(name = expression("R"[0]), na.value = "white",
#                        colours = c("white", "blue", "yellow", "red"),
#                        values = scales::rescale(c(0, 1, 1.000001, max(est_ub_r0 %>% ungroup() %>% select(ub_r0)))),
#                        guide = guide_colorbar(title=expression("R"[0]), barheight=10))
# ub_rnot_plot
#
# save_plot("ms_figs/sf2_med_and_ub_rnot_posterior.png", plot = ub_rnot_plot, base_height = 7, base_aspect_ratio = 1.3)
#
#
#
#
# ###############################
# ## Prior versus posterior binary whether
# ###############################
#
# bin_post_prior_rnot <- est_posterior %>% select(-alpha) %>%
#   group_by(county, month) %>%
#   summarise(post_high = quantile(x = rnot_samp, probs=0.975)) %>%
#   left_join(tx_county_rnots, by = c("month", "county")) %>%
#   select(county, month, post_high, high_r0) %>%
#   rename(prior_high = high_r0) %>%
#   gather(rnot_level, value, post_high, prior_high) %>%
#   mutate(month = factor(month, levels=month.abb),
#          bin_value = if_else(value >= 1, "0", "1"),
#          rnot_level = factor(if_else(rnot_level == "post_high", "Posterior", "Prior"), levels = c("Prior", "Posterior")))
#
# prior_post_compare_plot <- map_data(map = "county") %>% filter(region=="texas") %>%
#   mutate(subregion = if_else(subregion=="de witt", "dewitt", subregion)) %>%
#   left_join(bin_post_prior_rnot, by=c("subregion" = "county"))  %>%
#   ggplot(aes(x=long, y=lat, fill = bin_value, group = subregion)) +
#   facet_wrap(~month+rnot_level, labeller = label_wrap_gen(multi_line=FALSE), ncol = 6)+
#   geom_polygon(color = "gray", size=0.1) +
#   theme_nothing() +
#   scale_fill_manual(values = c("black", "white"), guide=FALSE)
# prior_post_compare_plot
# save_plot("ms_figs/sf3_prior_vs_post.png", plot = prior_post_compare_plot, base_height = 7, base_aspect_ratio = 1.3)
#
# ############################
# # Statewide daily scaling alpha results - Figure
# ############################
# load("data_produced/post_alpha_ci.rda")
#
#
#
# alpha_plot_fxn <- function(alpha_data, reporting, num_sec_trans){
#   alpha_data %>% filter(reporting_rate == reporting, secondary_trans==num_sec_trans) %>%
#     mutate(date = ymd(date)) %>%
#     ggplot(aes(date, med)) +
#     geom_point(size=2) +
#     geom_errorbar(aes(ymin = low, ymax = high))+
#     labs(y = bquote("Scaling Factor ("~alpha~")"), x = "Date", color = "Reporting\nRate") +
#     theme(legend.position = c(0.2,0.5))+
#     coord_cartesian(ylim=c(0,1)) +
#     scale_x_date(labels = date_format("%b"), breaks=date_breaks("month"))
# }
#
# alpha_plot <- alpha_plot_fxn(post_alpha_ci, 0.0574, 1)
# alpha_plot
#
# fig1 <- plot_grid(prior_rnot_plot,import_plot,  alpha_plot, nrow = 3, align="v", labels = "AUTO")
# save_plot("ms_figs/f1_priorr0_import_alpha.png", fig1, base_height = 8, base_aspect_ratio = 1.2)
#
# daily_alpha <- plot_grid(import_plot, alpha_plot, align="v", nrow=2, labels = "AUTO", rel_heights = c(1,1.3))
#
# save_plot("ms_figs/f3_daily_alpha.png", daily_alpha, base_height = 6, base_aspect_ratio = 1.5)
#
# ########
# ## Supplemental alpha plot figures
# ## Duplicate alpha plots for no secondary transmission, and for more secondary transmission
# alpha_plot_0_trans <- alpha_plot_fxn(post_alpha_ci, 0.0574, 0)
# alpha_plot_5_trans <- alpha_plot_fxn(post_alpha_ci, 0.0574, 5)
#
# alpha_alt_trans_plots <- plot_grid(alpha_plot_0_trans, alpha_plot_5_trans, align = "v", nrow = 2, labels="AUTO")
# alpha_alt_trans_plots
# save_plot("ms_figs/sf4_daily_alpha_alt_trans.png", alpha_alt_trans_plots, base_height = 6, base_aspect_ratio = 1.4)
#
# ## Alpha plot for single secondary transmission that shows effect of reporting rate
# alpha_plot_show_reporting <- post_alpha_ci %>% filter(reporting_rate %in% c("0.0282","0.0574", "0.0866"), secondary_trans==1) %>%
#   mutate(date = ymd(date)) %>%
#   ggplot(aes(date, med)) +
#   geom_point(aes(color = as.factor(reporting_rate)), size=2) +
#   # geom_errorbar(aes(ymin = low, ymax=high))+
#   labs(y = bquote("Scaling Factor ("~alpha~")"), x = "Date", color = "Reporting\nRate") +
#   scale_color_brewer(palette ="YlGnBu", type = "seq") +
#   theme(legend.position = c(0.15,0.8))+
#   coord_cartesian(ylim=c(0,1)) +
#   scale_x_date(labels = date_format("%b"), breaks=date_breaks("month"))+
#   guides(color = guide_legend(override.aes = list(size=5)))
# alpha_plot_show_reporting
# save_plot("ms_figs/sf5_daily_alpha_show_reporting.png", alpha_plot_show_reporting, base_height = 4, base_aspect_ratio = 1.8)

