#################################################
## All data can be downloaded from the repository specified in the README.md
#################################################


####################################################
## Processes the historic temperature data
library(tidyverse)
library(rgdal)
library(raster)
library(ncdf4)
library(lubridate)

extract_county_avg_temp <- function(county, temperature_raster){
  county_map <- map_data(map = "county") %>% filter(region=="texas", subregion == county)

  if(nrow(county_map)==0){
    stop("Problem getting county data, likely spelled a county name incorrectly.")
  }
  p <- Polygon(coords = county_map[,c("long", "lat")])
  ps <- Polygons(list(p), ID = 1)
  sps <- SpatialPolygons(list(ps))

  mean(extract(temperature_raster, sps)[[1]]/10, na.rm=T)
}

compile_temp_df <- function(raster_location){
  ## Takes in the raster location and compiles all of the average county temperatures from that raster

  temperature_raster <- raster(raster_location) # "data/tmean_2-5m_bil/tmean1.bil")

  counties <- map_data(map="county") %>% filter(region=="texas") %>% dplyr::select(subregion) %>% unique()
  county_temps <- purrr::map(counties$subregion, extract_county_avg_temp, temperature_raster)
  counties$avg_temp <- unlist(county_temps)

  print("Done with one of them.")
  
  ## Finds the number of the file being called, and sets to that month abbreviation
  counties$month <- month.abb[as.numeric(strsplit(strsplit(raster_location, split = ".", fixed = T)[[1]], split="n", fixed=T)[[1]][3])]
  counties
}


raster_locations <- list.files(path="data/tmean_2-5m_bil", pattern = "*.bil", full.names = T)

tx_temperatures <- purrr::map(raster_locations, compile_temp_df)

tx_temperatures <- tx_temperatures %>% bind_rows() %>% mutate(month = factor(month, levels = month.abb)) %>%
  arrange(month)

write_csv(tx_temperatures, path = "data/tx_historic_temps.csv")


#########################################################
## Processes the actual temperature data from 2016-17

t_2016 <- brick("data/air.2m.2016.nc")
t_2017 <- brick("data/air.2m.2017.nc")

get_county_temps <- function(current_county, brick_data, current_year){
  ## Extract the specified county spatial information
  county_map <- map_data(map = "county") %>% filter(region=="texas", subregion == current_county)

  ## First turn the county into a spatial object on the correct projection
  p <- Polygon(coords = county_map[,c("long", "lat")])
  ps <- Polygons(list(p), ID = 1)
  sps <- SpatialPolygons(list(ps), proj4string = CRS("+proj=longlat"))
  sps <- spTransform(sps, projection(t_2016))

  ## Then extract the temperature data
  cty_temps <- extract(brick_data, sps)

  ## If list is longer than 1 need to make sure rest properly works
  if(length(cty_temps)!=1){
    browser()
  }

  ## Average all readings within county for each date
  cty_temps <- apply(cty_temps[[1]], MARGIN = 2, FUN = mean)

  ## Convert to data_frame with proper units, and extract and average across month
  data_frame(date = names(cty_temps), temperature = cty_temps-273.15) %>%
    separate(col = date, sep = "\\.", into = c("year", "month", "day", "hour", "minute", "second")) %>%
    mutate(year = substring(year, 2)) %>%
    unite(col = date, year:day, sep="-") %>%
    mutate(date = ymd(date)) %>%
    filter(year(date) == current_year) %>%
    group_by(month = month(date, label=TRUE)) %>%
    summarise(avg_temp = mean(temperature)) %>%
    mutate(subregion = current_county, year = current_year) %>%
    dplyr::select(subregion, year, month, avg_temp)
}

## extract counties
counties <- map_data(map = "county") %>% filter(region=="texas") %>%
  group_by(subregion) %>%
  summarize(county=unique(subregion)) %>%
  ungroup() %>%
  dplyr::select(county)

## Run for all counties of 2016
county_temps_2016 <- counties[[1]] %>% purrr::map(~get_county_temps(current_county =.x, brick_data=t_2016, current_year=2016)) %>% bind_rows()
county_temps_2017 <- counties[[1]] %>% purrr::map(~get_county_temps(current_county =.x, brick_data=t_2017, current_year=2017)) %>% bind_rows()

write_csv(bind_rows(county_temps_2016, county_temps_2017), path = "data/tx_actual_temps.csv")




