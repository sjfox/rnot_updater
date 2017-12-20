base_url <- "zikaEstimatoR"
if(grepl('spencerfox', Sys.info()['login'])) setwd(file.path("~", "projects", base_url))
if(grepl('vagrant', Sys.info()['user'])) setwd( file.path("/vagrant", base_url) )
if(grepl('stampede', Sys.info()['nodename'])) setwd(file.path('/home1/02958/sjf826', base_url))
if(grepl('wrangler', Sys.info()['nodename'])) setwd(file.path('/home/02958/sjf826', base_url))


# rnot_values <- c("low","med", "high")
# reporting_rate <- c(0.01, 0.0282, 0.0574, 0.0866, 0.1, 0.2)
include_trans <- c("NA", "1", "5")
temperature <- c("historic", "actual")
extra_imports <- c(0,1)

sink('launcher/alpha_mcmc_runs.txt')
for(ind in seq_along(extra_imports)){
  for(trans_ind in seq_along(include_trans)){
    for(temp_ind in seq_along(temperature)){
      startCmd <- "R CMD BATCH --no-restore --no-save '--args"
      #fileCmd <- paste0(' single_rnot=', single_rnot, ' rnot_value="', rnot_values[ind], '"')
      fileCmd <- paste0(' temperature="', temperature[temp_ind], '" include_trans="', include_trans[trans_ind], '" extra_imports="', extra_imports[ind], '"')
      endCmd <- "' ../R/calc_alpha_posterior.R"
      full_cmd <- paste0(startCmd, fileCmd, endCmd)
      # print(full_cmd)
      cat(full_cmd)               # add command
      cat('\n')              # add new line
    }
  }
}
sink()
