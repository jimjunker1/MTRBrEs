library(lme4)
library(extrafont)
library(fontcm)
loadfonts()

rm(list=ls())
source('paths.R')
source('R/functions-figures.R')
source('dataManipulation.R')
load('output/RDatafiles/analyses.RData')

to.pdf(posteriorsBoltzmann(), 'output/figures/posteriorsBoltzmann.pdf', width=11, height=5.5)
to.pdf(modelComparison(), 'output/figures/modelComparison.pdf', width=11, height=5)
embed_fonts('output/figures/modelComparison.pdf')
