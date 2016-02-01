library(plyr)
library(lme4)

rm(list=ls())
source('paths.R')
source('R/functions-analyses.R')
load('output/RDatafiles/analyses.RData')

# model selection result
deltaDF   <-  anova(modelLmer1, modelLmer2)['Chi Df'][2,]
deltaDIC  <-  tfit2$BUGSoutput$DIC - tfit$BUGSoutput$DIC

# get average estimates and 95% CI for each parameter
species           <-  unique(metRates$Species)
averageEstimates  <-  ldply(species, getAverageEstimates)
ci95Estimates     <-  ldply(species, get95CIEstimates)

# compare deltas between temperatures for each species
tempDeltas  <-  ddply(averageEstimates, .(species), getTempDeltas)

# q10-equivalents [ q10 == (r25 / r10) ^ (10/(25-10)) ] for Microcionidae and Bugula stolonifera
q10andErs  <-  ddply(averageEstimates[averageEstimates$species %in% c('Sponge', 'Stolonifera'), ], .(species), getq10andEr)

rm(list=ls()[!(ls() %in% c('deltaDF', 'deltaDIC', 'species', 'averageEstimates', 'ci95Estimates', 'tempDeltas', 'q10andErs'))])
save.image('output/RDatafiles/numbers.RData')
