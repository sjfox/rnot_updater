##############################################################
## R0 functions for calculating distributions for each county
## Spencer Fox
## January 18th, 2017
##############################################################



rnot_calc_dist <- function(mosq_abundance, gdp, temperature,
                           a, b, c.r, mort.fun.list, eip.fun.list, scam.est.list){
  # Function that returns the full distribution of R0s for n county data provided
  require(scam)

  index = seq(1:length(scam.est.list))

  # Mortality function actually gives lifespan, so take 1 / lifespan to get mortality
  # calculates distribution of mortality rates
  g <- 1 / sapply(index, FUN = function(x) {
    mort.fun.list[[x]](temperature)
  })

  # calculates distribution of eip lengths
  e <- sapply(index, FUN = function(x) {
    eip.fun.list[[x]](temperature)
  })

  # calculates distribution of gdp scaling multipliers
  gdp_scaling <- unname(sapply(index, FUN = function(x) {
    predict(scam.est.list[[x]], newdata=data.frame(econ=log(gdp)))
  }))
  # calculates and returns distributions of rnott
  as.data.frame(exp(gdp_scaling) * mosq_abundance * a ^ 2 * b * c.r * exp(-g * e) / g)
}
