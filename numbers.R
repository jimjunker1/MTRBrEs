library(plyr)
library(lme4)

rm(list=ls())
source('paths.R')
source('R/functions-analyses.R')
load('output/RDatafiles/analyses.RData')

# get average estimates and 95% CI for each parameter
species           <-  unique(metRates$Species)
averageEstimates  <-  ldply(species, getAverageEstimates)
ci95Estimates     <-  ldply(species, get95CIEstimates)

# compare deltas between temperatures for each species
tempDeltas  <-  ddply(averageEstimates, .(species), getTempDeltas)

# metabolic difference (in n-fold) between different temperatures for similar-sized animals
metabolicDiff  <-  ddply(averageEstimates, .(species), getMetabolicDiff, lnMass=log(400))

# q10-equivalents [ q10 == (r25 / r10) ^ (10/(25-10)) ] for Microcionidae and Bugula stolonifera
q10andErs  <-  ddply(averageEstimates[averageEstimates$species %in% c('Sponge', 'Stolonifera'), ], .(species), getq10andEr)

rm(list=ls()[!(ls() %in% c('species', 'averageEstimates', 'ci95Estimates', 'tempDeltas', 'metabolicDiff', 'q10andErs'))])
save.image('output/RDatafiles/numbers.RData')
