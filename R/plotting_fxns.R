##############################
## Helper functions for Zika plots
## Spencer Fox October 30, 2017
##############################



## Fxn for adding proper color scaling/outline to maps
tx_outline <- map_data(map="state") %>% filter(region=="texas")
add_map_scale <- function(gmap, max_rnot, texas = tx_outline){
  if(max_rnot > 1){
    gmap <- gmap +
      scale_fill_gradientn(name = expression("R"[0]), na.value = "white",
                           colours = c("white", "#264BAC", "#ffeda0", "#8A1923"),
                           values = scales::rescale(c(0, 1, 1.000001, max_rnot)),
                           guide = guide_colorbar(title=expression("R"[0]), barheight=10))
  } else{
    gmap <- gmap +
      scale_fill_gradient(low="white", high="#264BAC",
                          guide = guide_colorbar(title=expression("R"[0]), barheight=10))
  }
  gmap <- gmap + geom_polygon(data=texas, aes(x=long, y=lat), fill=NA, size=0.1,color="black", inherit.aes = FALSE)
  gmap
}



get_posterior_data <- function(temperature, include_trans, extra_imports, alpha=FALSE){
  ## temperature = "actual" or "historic" - specifies if actualy 2016/17 temps were used or historic ones
  ## include_trans = c(NA, 1, 5) - specifies how many secondary transmission cases should be used in November Cameron county
  ## reporting_rate = c(0.01, 0.0282, 0.0574, 0.0866, 0.1, 0.2) - assumed reporting rate of 0.0574 for most results in paper
  all_post_files <- list.files(path = paste0("data_produced/", data_produced_date,"-posteriors/"), pattern = "*.rda")
  # browser()
  if(alpha){
    path <- paste0("alpha_daily_mcmc_", temperature, "_",
                   ifelse(is.na(include_trans), 0, include_trans),"_",
                   ifelse(extra_imports, "true", "false"), ".rda")
  } else {
    path <- paste0("county_posterior_rnots_", temperature, "_",
                   ifelse(is.na(include_trans), 0, include_trans),"_",
                   ifelse(extra_imports, "true", "false"), ".rda")
  }

  if(!path %in% all_post_files){
    stop("Either your path is incorrect, your parameters are incorrect, or you haven't run the posterior estimation for the parameter combination.")
  }

  load(paste0("data_produced/", data_produced_date,"-posteriors/", path))

  if(alpha){
    est_alphas_df
  } else{
    est_posterior
  }
}

plot_monthly_post_rnots <- function(posterior_rnots, quant = 0.5){
  ## Plots faceted years worth of maps for posterior R0 estimates
  ## quants specifies what quantile to plot 0.5 is median
  est_r0 <- posterior_rnots %>% select(-alpha) %>%
    group_by(county, month, year) %>%
    summarise(r0 = quantile(x = rnot_samp, probs=quant, na.rm=T))

  if(length(unique(est_r0$year)) > 1){
    est_r0 <- est_r0 %>% filter(year==2016)
  }

  rnot_plot <- map_data(map = "county") %>% filter(region=="texas") %>%
    mutate(subregion = if_else(subregion=="de witt", "dewitt", subregion)) %>%
    left_join(est_r0, by=c("subregion" = "county"))  %>%
    mutate(month = factor(month, levels=month.abb)) %>%
    ggplot(aes(x=long, y=lat, fill = r0, group = subregion)) + facet_wrap(~month)+
    geom_polygon(color = "gray", size=0.1) +
    theme_void(base_size = 14)

  rnot_plot <- add_map_scale(rnot_plot, max_rnot = max(est_r0$r0, na.rm=T))
  rnot_plot
}

c_trans <- function(a, b, breaks = b$breaks, format = b$format) {
  a <- as.trans(a)
  b <- as.trans(b)

  name <- paste(a$name, b$name, sep = "-")

  trans <- function(x) a$trans(b$trans(x))
  inv <- function(x) b$inverse(a$inverse(x))

  trans_new(name, trans, inv, breaks, format)
}

alpha_plot_fxn <- function(alpha_data){
  rev_date <- c_trans("reverse", "date")

  tx_temps <- read_csv("data/tx_actual_temps.csv")
  tx_temps <- tx_temps %>% mutate(month = factor(month, levels = month.abb)) %>%
    filter(year==2016) %>%
    group_by(month) %>%
    summarise(avg_temp = mean(avg_temp))

  month_dates <- seq(ymd("2016-01-01"),ymd("2017-03-01"), "month")

  alpha_data %>% gather(date, alpha_samp, 1:ncol(alpha_data)) %>%
    mutate(date = ymd(date),
           month = month(date, label = TRUE)) %>%
    mutate(days_since = as.numeric(date - min(date))) %>%
    group_by(date) %>%
    left_join(tx_temps, by="month") %>%
    ggplot(aes(alpha_samp, date, group = date, fill = avg_temp)) +
    geom_density_ridges(scale=50, rel_min_height = 0.01) +
    scale_y_continuous(trans = rev_date, breaks = month_dates, labels = date_format("%b-%y")) +
    scale_x_continuous(breaks = c(0,0.5,1, 1.5, 2), expand = c(0.01, 0)) +
    theme_ridges() + theme(panel.grid.major = element_line(color="gray50")) +
    scale_fill_gradient(low="#FEEDED", high="#8A1923") +
    guides(fill = guide_colorbar(title = "Average State\nTemperature", barheight = 10,barwidth = 2)) +
    labs(x = bquote("Scaling Factor ("~alpha~")"), y = "Date")
}



plot_med_ub_post_map <- function(temperature){
  ## temperature has to be "historic" or "actual"
  est_posterior = get_posterior_data(temperature = temperature, include_trans = 1, extra_imports = FALSE, alpha = FALSE)

  posterior_med_ub <- est_posterior %>% select(-alpha) %>%
    group_by(county, month) %>%
    summarise(med_r0 = quantile(x = rnot_samp, probs=0.5, na.rm=T),
              ub_r0 = quantile(x = rnot_samp, probs=0.99, na.rm=T))




  post_med_ub <- map_data(map = "county") %>% filter(region=="texas") %>%
    mutate(subregion = if_else(subregion=="de witt", "dewitt", subregion)) %>%
    left_join(posterior_med_ub, by=c("subregion" = "county"))  %>%
    gather(rnot_level, rnots, med_r0, ub_r0) %>%
    mutate(month = factor(month, levels=month.abb),
           rnot_level = factor(if_else(rnot_level == "med_r0", "Median", "99%"), levels = c("Median", "99%"))) %>%
    ggplot(aes(x=long, y=lat, fill = rnots, group = subregion)) +
    facet_wrap(~month+rnot_level, labeller = label_wrap_gen(multi_line=FALSE), ncol = 6) +
    geom_polygon(color = "gray", size=0.1) +
    theme_void(base_size = 14)

  post_med_ub <- add_map_scale(post_med_ub, max_rnot = max(posterior_med_ub$ub_r0, na.rm=T))
  post_med_ub
}
