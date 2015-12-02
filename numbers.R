library(lme4)

rm(list=ls())
source('paths.R')
load('output/RDatafiles/analyses.RData')

# compare q10's between both models
r10  <-  exp(fixef(modelLmer)['(Intercept)'] + fixef(modelLmer)['lnMass']*1)
r25  <-  exp((fixef(modelLmer)['(Intercept)'] + fixef(modelLmer)['Temp25']) + fixef(modelLmer)['lnMass']*1)
q10  <-  (r25 / r10) ^ (10/(25-10))

r10b  <-  exp(fixef(boltzLmer)['(Intercept)'] + fixef(boltzLmer)['lnMass']*1 + (fixef(boltzLmer)['invKT'] / 8.62e-5 * (1 / 293.15 - 1 / 283.15)))
r25b  <-  exp(fixef(boltzLmer)['(Intercept)'] + fixef(boltzLmer)['lnMass']*1 + (fixef(boltzLmer)['invKT'] / 8.62e-5 * (1 / 293.15 - 1 / 298.15)))
q10b  <-  (r25b / r10b) ^ (10/(25-10))
