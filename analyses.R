library(R2jags)
library(loo)
library(plyr)

rm(list=ls())
source('paths.R')
source('dataManipulation.R')
source('R/functions-analyses.R')

# run models in lm in order to obtain model matrices of contrasts for jags
lmModel1  <-  lm(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + Temp:lnMass + Temp:lnMass:Species, data=metRates)
lmModel2  <-  lm(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + Temp:lnMass, data=metRates)
lmModel3  <-  lm(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp, data=metRates)
lmModel4  <-  lm(lnRate ~ Species + lnMass + Temp + Species:lnMass + Temp:lnMass, data=metRates)
lmModel5  <-  lm(lnRate ~ Species + lnMass + Temp + Species:Temp + Temp:lnMass, data=metRates)
lmModel6  <-  lm(lnRate ~ Species + lnMass + Temp + Species:lnMass, data=metRates)
lmModel7  <-  lm(lnRate ~ Species + lnMass + Temp + Species:Temp, data=metRates)
lmModel8  <-  lm(lnRate ~ Species + lnMass + Temp + Temp:lnMass, data=metRates)
lmModel9  <-  lm(lnRate ~ Species + lnMass + Temp, data=metRates)

# first use full model to do a model selection for random effects
r            <-  c(TRUE, FALSE)
randomSelec  <-  vector(mode='list', length=length(r))
for(i in seq_along(r)) {
	model             <-  fitJags(lmModel1, run=r[i])
	log_lik           <-  model$BUGSoutput$sims.matrix[, grep('logLik', colnames(model$BUGSoutput$sims.matrix))]
	randomSelec[[i]]  <-  list(jagsFit = model, loo = loo(log_lik))
}; rm(r, i, model, log_lik)

# compare nested models one with another
nested        <-  data.frame(complexModel = 1, nestedModel = 2)
RandEffComps  <-  ddply(nested, .(complexModel, nestedModel), modelComparisonPVals, modelList=randomSelec)

# so, run is significant, now run nested models for assessment of fixed effects
fixedSelec       <-  vector(mode='list', length=9)
fixedSelec[[1]]  <-  randomSelec[[1]]
for(i in 2:9) {
	model            <-  fitJags(get(paste0('lmModel', i)))
	log_lik          <-  model$BUGSoutput$sims.matrix[, grep('logLik', colnames(model$BUGSoutput$sims.matrix))]
	fixedSelec[[i]]  <-  list(jagsFit = model, loo = loo(log_lik))
}; rm(i, model, log_lik)

# compare nested models one with another
# start with full model
nested       <-  data.frame(complexModel = c(rep(1, 8), rep(2, 7), rep(3, 2)), nestedModel = c(2:9, 3:9, 6:7))
FixEffComps  <-  ddply(nested, .(complexModel, nestedModel), modelComparisonPVals, modelList=fixedSelec)

rm(list=ls()[!(ls() %in% c('metRates', 'lmModel1', 'randomSelec', 'fixedSelec', 'RandEffComps', 'FixEffComps'))])
save.image('output/RDatafiles/analyses.RData')
