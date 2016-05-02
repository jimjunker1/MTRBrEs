library(lme4)
library(R2jags)

rm(list=ls())
source('paths.R')
source('dataManipulation.R')
source('R/functions-analyses.R')

# compare random effects among models
lmerModela  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + Temp:lnMass + Temp:lnMass:Species + (1|Run) + (1|Plate), data=metRates, REML=FALSE)
lmerModelb  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + Temp:lnMass + Temp:lnMass:Species + (1|Plate), data=metRates, REML=FALSE)
lmerModelc  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + Temp:lnMass + Temp:lnMass:Species + (1|Run), data=metRates, REML=FALSE)
modelSelectionRandomTab  <-  anova(lmerModela, lmerModelb, lmerModelc)

# lmerModelc most significant. Now run models to asses significance of fixed effects
lmerModel1  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + Temp:lnMass + Temp:lnMass:Species + (1|Run), data=metRates, REML=FALSE)
lmerModel2  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + Temp:lnMass + (1|Run), data=metRates, REML=FALSE)
lmerModel3  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + (1|Run), data=metRates, REML=FALSE)
lmerModel4  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Temp:lnMass + (1|Run), data=metRates, REML=FALSE)
lmerModel5  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:Temp + Temp:lnMass + (1|Run), data=metRates, REML=FALSE)
lmerModel6  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + (1|Run), data=metRates, REML=FALSE)
lmerModel7  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:Temp + (1|Run), data=metRates, REML=FALSE)
lmerModel8  <-  lmer(lnRate ~ Species + lnMass + Temp + Temp:lnMass + (1|Run), data=metRates, REML=FALSE)
lmerModel9  <-  lmer(lnRate ~ Species + lnMass + Temp + (1|Run), data=metRates, REML=FALSE)
modelSelectionFixedTab  <-  anova(lmerModel1,lmerModel2,lmerModel3,lmerModel4,lmerModel5,lmerModel6,lmerModel7,lmerModel8,lmerModel9)

# full model is the most significant one. Refit in jags
model       <-  fitJags(lmerModel1)

rm(list=ls()[!(ls() %in% c('metRates', paste0('lmerModel', 1:9), 'model', 'modelSelectionRandomTab', 'modelSelectionFixedTab'))])
save.image('output/RDatafiles/analyses.RData')
