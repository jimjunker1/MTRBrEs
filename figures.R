library(extrafont)
library(fontcm)
loadfonts()
source('R/functions-figures.R')

to.pdf(plotResp(), 'output/figures/metRates.pdf', width=12, height=5.5)
embed_fonts('output/figures/metRates.pdf')
