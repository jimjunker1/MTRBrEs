library(plyr)
rm(list=ls())
source('paths.R')
source('dataManipulation.R')
source('R/functions-figures.R')

table1  <-  ddply(metRates, .(Species, Temp), function(x) {
	data.frame(minMass = rounded(min(exp(x$lnMass))),
			   maxMass = rounded(max(exp(x$lnMass))),
			   n         = nrow(x), stringsAsFactors=FALSE)
})

write.csv(table1, 'output/data/table1.csv', row.names=FALSE)
