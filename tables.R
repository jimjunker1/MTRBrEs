library(lme4)
library(plyr)
rm(list=ls())
source('paths.R')
source('R/functions-analyses.R')
source('R/functions-figures.R')
load('output/RDatafiles/analyses.RData')

# table 1
table1  <-  ddply(metRates, .(Species, Temp), function(x) {
	data.frame(minMass = rounded(min(exp(x$lnMass))),
			   maxMass = rounded(max(exp(x$lnMass))),
			   n         = nrow(x), stringsAsFactors=FALSE)
})
write.csv(table1, 'output/data/table1.csv', row.names=FALSE)

# table 2
table2  <-  data.frame()
for(j in 2:9) {
	refModel  <-  get(paste0('lmerModel', j))
	aovTab    <-  data.frame(anova(lmerModel1, refModel))['lmerModel1',]
	table2    <-  rbind(table2, data.frame('d.f.'=aovTab[,'Chi.Df'], 'Chisq'=aovTab[,'Chisq'], 'P'=aovTab[,'Pr..Chisq.'], stringsAsFactors=FALSE))
}
table2$P  <-  sapply(table2$P, cleanPvals)
write.csv(table2, 'output/data/table2.csv', row.names=FALSE)
