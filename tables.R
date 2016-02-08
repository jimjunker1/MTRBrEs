library(plyr)
rm(list=ls())
source('paths.R')
source('dataManipulation.R')
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
# modelComparisonPVals  <-  function(data, modelList, verbose=TRUE) {
#FixEffComps  <-  ddply(nested, .(complexModel, nestedModel), modelComparisonPVals, modelList=fixedSelec)

table2  <-  ldply(2:9, function(x) {
	modDiff  <-  compare(fixedSelec[[1]]$loo, fixedSelec[[x]]$loo)
	data.frame(elpdDiff  =  rounded(modDiff[['elpd_diff']], 2),
			   se        =  rounded(modDiff[['se']], 2),
			   zPVal     =  cleanPvals(2*pnorm(-abs((modDiff[['elpd_diff']] - 0) / modDiff[['se']]))))
})
write.csv(table2, 'output/data/table2.csv', row.names=FALSE)