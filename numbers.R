library(plyr)
library(lme4)

rm(list=ls())
source('paths.R')
source('R/functions-analyses.R')
load('output/RDatafiles/analyses.RData')

# model selection result
deltaDF   <-  anova(modelLmer1, modelLmer2)['Chi Df'][2,]
deltaDIC  <-  tfit2$BUGSoutput$DIC - tfit$BUGSoutput$DIC

# compare deltas between temperatures for each species
species     <-  unique(metRates$Species)
tempDeltas  <-  ldply(species, getTempDeltas)

save.image('output/RDatafiles/numbers.RData')