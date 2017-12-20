##########################################
## Generate correct scam.est.list for R0 estimation
## Original scam.est.list includes a R0 scaling parameter
## That might distort the estimates
## Perkins Nature Micro 2016 methods
##########################################

reps = 1000

set.seed(808)

# get random samples from confidence interval for EIP-temperature relationship
library(pomp)
beta.0 = c(6,10)
beta.T = c(-.29,-.12)
scalars = sobol(vars = list(beta.0 = beta.0, beta.T = beta.T), 1e5)
eip = function(T,beta.0,beta.T){exp(beta.0 + beta.T * T)}
eip.vec = sapply(1:nrow(scalars),function(ii)eip(30,scalars[ii,1],scalars[ii,2]))
eip.sd = optimize(
  f = function(sd){(qnorm(0.05,6.1,sd) - 3.4) ^ 2 + (qnorm(0.95,6.1,sd) - 9.9) ^ 2},
  interval = c(0,10)
)$minimum
beta.samples = matrix(0,reps,2)
eip.fun = list()
for(rr in 1:reps){
  repeat{
    eip.draw = rnorm(1,6.1,eip.sd)
    if(eip.draw > 0){
      break
    }
  }
  which.range = which(eip.vec > (eip.draw - 0.05) & eip.vec < (eip.draw + 0.05))
  beta.samples[rr,] = colMeans(scalars[which.range,])
  eip.fun[[rr]] =  approxfun(seq(0,50,.1), sapply(seq(0,50,.1), function(T){exp(beta.samples[rr,1] + beta.samples[rr,2] * T)}))
}


# get random samples from confidence interval of mortality-temperature relationship
load('data/vector_suitability/algam_85re.Rdata')
lifespan <- function(Temperature,se.mult){
  dd <- seq(0,120,length.out=(120*24+2))
  nwdd <- data.frame(Days=dd,Temperature=rep(Temperature, (120*24+2)), Study_number=5, Feed_B=2, Feed_S=1)
  nwdd <- cbind(nwdd,logDay=log(nwdd$Days+1), logTemp=log(nwdd$Temperature+1))  ## +1 avoids log(0)
  prediction <- as.vector(unlist(predict(algam,newdata = (nwdd), se.fit = TRUE, type = "response")$fit))
  prediction.se <- as.vector(unlist(predict(algam,newdata = (nwdd), se.fit = TRUE, type = "response")$se.fit))
  prediction <- prediction + se.mult * prediction.se
  prediction <- prediction[-1]
  prediction[1:24] <- prediction[1:24]/prediction[1]
  prediction[which(prediction>1)] <- 1
  prediction[which(prediction<=0.001)] <- 0
  diffDeath <- (prediction[1: (length(prediction)-1)] - prediction[2:length(prediction)])
  diffDeath <- diffDeath/sum(diffDeath)
  return(pmax(0, sum(dd[2:length(prediction)]*diffDeath)))
}
fielddata <- rbind(c(20, 34, 1 - 0.91))
tgrd <- 0.1
lifespan.range <- sapply(seq(fielddata[1],fielddata[2],tgrd), function (tr) {lifespan(tr,0)})
fieldcorxn <- fielddata[3] - 1/mean(lifespan.range)
mort.fun = list()
for(rr in 1 : reps){
  normdev = rnorm(1,0,1)
  mort.fun[[rr]] = approxfun(seq(0,50,tgrd), sapply(seq(0,50,tgrd), function(TT) 1/(1/lifespan(TT,normdev)+fieldcorxn)))
}

if(!dir.exists("data_produced")){
  dir.create("data_produced")
}

if(!dir.exists("data_produced/vector_suitability")){
  dir.create("data_produced/vector_suitability")
} 

save(reps,eip.fun,mort.fun,file='data_produced/vector_suitability/params.RData')


load('data_produced/vector_suitability/params.Rdata')

# constant parameters
b = 0.4
c.r = 3.5
a = 1 / 1.5

# R0 as a function of mosquito-human ratio and temperature
R0 = function(m,T,rep2,rep3){
  g = 1 / mort.fun[[rep2]](T)
  e = eip.fun[[rep3]](T)
  m * a ^ 2 * b * c.r * exp(-g * e) / g
}

# load replicate aegypti occurrence probabilities
aegypti.maps.pts.data = as.matrix(read.csv('data/vector_suitability/aegypti_pt_vals_all.csv')[,-(1:4)])

# seroprevalence data and metadata
x.data = read.csv('data/vector_suitability/seroprev_metadata.csv')
x.data$factor.adj = 1

# attack rate as a function of R0
R0.vec = seq(0,1000,.01)
AR.vec = numeric(length(R0.vec))
for(ii in 1:length(AR.vec)){
  AR.vec[ii] = 1 - optimize(f=function(S){(S-exp(-R0.vec[ii]*(1-S)))^2},interval=c(0,1))$minimum
}
AR.fun.vec = approxfun(R0.vec,AR.vec)
AR.fun = function(R0,h){
  if(R0 < 1){
    return(0)
  }
  if(R0 >= 1){
    return(AR.fun.vec(R0 ^ h))
  }
}

# R0 for a given location
R0.base = function(row,rep1,rep2,rep3){
  vecs = c(0,0)
  vecs[1]=-log(1-pmax(aegypti.maps.pts.data[row,rep1],x.data$albopictus))
  m.fun = max(vecs)
  mean(sort(sapply(x[
    row,c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')],
    function(TT) R0(m.fun,TT,rep2,rep3)),decreasing=TRUE)[1:6])
}

# attack rate for a given location
AR.avg = function(row,h,rep1,rep2,rep3){
  R0.fun = x$factor.adj[row] * R0.base(row,rep1,rep2,rep3)
  ifelse(R0.fun<1000,AR.fun(R0.fun,h),1)
}

# R0 needed to make attack rate consistent with seroprevalence data
R0.needed = function(row,h){
  optim(par=1,fn=function(par){abs(AR.fun(par,h)-x$seroprev[row])})$par
}


# fit a SCAM for the multiplication factor - economic index relationship
# for each of 1,000 random draws of EIP, mortality, and Aedes occurrence
library(scam)
scam.est.list = list()
h.list = list()
h.list[[reps+1]] = 1
rep.master = cbind(
  sample(ncol(aegypti.maps.pts.data),reps,replace=T),
  1:reps,
  1:reps)
for(rr in 1 : reps){
  repeat{
    if(!is.null(h.list[[rr]])){
      break
    }

    # select random combinations of random draws from each of the EIP, occurrence, and mortality distributions
    rep.master[rr,] = c(
      sample(ncol(aegypti.maps.pts.data),1,replace=T),
      sample(reps,1,replace=T),
      sample(reps,1,replace=T))
    rep1 = rep.master[rr,1]
    rep2 = rep.master[rr,2]
    rep3 = rep.master[rr,3]

    # get serological data and Aedes aegypti occurrence probabilities by site
    sites = 1:nrow(x.data)
    x = x.data[sites,]
    aegypti.maps.pts = aegypti.maps.pts.data[sites,rep1]

    # factor change in R0 necessary to make each attack rate == seroprevalence data for each h value
    h.vec = c(1)
    factor.needed = matrix(0,nrow(x),length(h.vec))
    for(hh in 1:length(h.vec)){
      factor.needed[,hh] =
        sapply(1:nrow(x),function(rr)R0.needed(rr,h.vec[hh])) /
        sapply(1:nrow(x),function(rr)R0.base(rr,rep1,rep2,rep3))
    }

    # # fit SCAMs for each h value in the range from 0.01 to 1 and record attack rate residuals
    # res.vec = numeric(length=length(h.vec))
    # for(hh in 1:length(h.vec)){
    #   h = h.vec[hh]
    #   df = data.frame(econ=log(x$gdp_pcppp2005),factor=log(factor.needed[,hh]))
    #   x$factor.adj = 1
    #   if(sum(is.infinite(df$factor))){
    #     next
    #   }
    #   scam.est = scam(factor~s(econ,bs='mdcx'),data=df)
    #   x$factor.adj = exp(scam.est$fitted.values)
    #   res.vec[hh] = sum((sapply(1:nrow(x),function(rr)AR.avg(rr,h,rep1,rep2,rep3)) - x$seroprev)^2)
    # }

    # select h value that minimizes residuals and re-fit associated SCAM
    h = 1
    df = data.frame(econ=log(x$gdp_pcppp2005),factor=log(factor.needed[,1]))
    x$factor.adj = 1
    if(sum(is.infinite(df$factor))){
      next
    }
    scam.est = scam(factor~s(econ,bs='mdcx'),data=df)
    x$factor.adj = exp(scam.est$fitted.values)

    # store best-fit SCAM model and associated h value
    scam.est.list[[rr]] = scam.est
    h.list[[rr]] = h
  }
}

library(scam)
econ = seq(6.5,10.5,0.1)
plot(
  econ,predict(scam.est.list[[1]],newdata=data.frame(econ=econ)),
  type='l',col=rgb(0,0,0,.05),ylim = c(-3,2), xlab='',ylab='')
for(ii in 2:length(scam.est.list)){
  lines(econ,predict(scam.est.list[[ii]],newdata=data.frame(econ=econ)),col=rgb(0,0,0,.1))
}

plot(mort.fun[[1]](1:40), type="l")
for(ii in 2:length(mort.fun)){
  lines(1:40, mort.fun[[ii]](1:40))
}


# save data
save(reps,rep.master,scam.est.list,mort.fun,eip.fun,a,b,c.r,h.list,file='data_produced/vector_suitability/parms_fxns_r0.RData')
