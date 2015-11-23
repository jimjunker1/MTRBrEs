library(lme4)

rm(list=ls())
source('paths.R')
source('R/functions-analyses.R')

# growth forms
# Bugula and Stolinifera = 3D arborescent
# Bryo and Sponge = 2D and encrusting
# Bugula = Bugula neritina
# Stolinifera = Bugula stolinifera
# Bryo = Hippopodina sp.	
# Sponge = No one knows (Hart & Marshall 2009 - Ecology)

metRates  <-  readFile('data/data.csv')
metRates  <-  metRates[complete.cases(metRates), ]

# correct log values from log 10 to ln
metRates$lnMass  <-  log(10^metRates$logM)
metRates$lnRate  <-  log(10^metRates$logVo2)

# create inverse temperature column
metRates$TempK  <-  metRates$Temp + 273.15
metRates$invKT  <-  1 / 8.62e-5 * (1 / 293.15 - 1 / metRates$TempK)

modelLmer  <-  lmer(lnRate ~ lnMass + invKT + (1 + lnMass + invKT | Species), data=metRates)

# 95% CI for slopes
quantile(fixef(modelLmer)['lnMass'] + ranef(modelLmer)$Species$lnMass, probs=c(.025, 0.975), type=2)
quantile(fixef(modelLmer)['invKT'] + ranef(modelLmer)$Species$invKT, probs=c(.025, 0.975), type=2)
