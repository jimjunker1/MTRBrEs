library(lme4)
library(MASS)
library(RColorBrewer)
library(extrafont)
library(fontcm)
loadfonts()

rm(list=ls())
source('paths.R')
source('R/functions-figures.R')
load('output/RDatafiles/analyses.RData')

to.pdf(fig1(), 'output/figures/fig1.pdf', width=8, height=7)
embed_fonts('output/figures/fig1.pdf')
