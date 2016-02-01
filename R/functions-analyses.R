#########
# GENERAL
#########
readFile  <-  function(path, ...) {
	read.csv(path, header=TRUE, stringsAsFactors=FALSE, ...)
}

rmFcts  <-  function(filename) {
  source(filename, local=TRUE)
  fct_files  <-  ls()
  rm(list=ls(.GlobalEnv)[ls(.GlobalEnv) %in% fct_files], envir=.GlobalEnv)
}

#############################
# DEALING WITH BIB REFERENCES
#############################
getIndividualBibs  <-  function(path='~/bibtex_library/') {
  listRefs     <-  dir(path)[grep('.bib', dir(path))]
  listRefs[listRefs != 'library.bib']
}

readBibRefs  <-  function(individualBibs, path='~/bibtex_library/', verbose=TRUE) {
  lapply(individualBibs, function(x, path, verbose) {
    if(verbose) {
      print(x)
    }
    read.bib(file.path(path, x))
  }, path=path, verbose=verbose)
}

listBibs  <-  function() {
  listRefs        <-  getIndividualBibs()
  theRefs         <-  readBibRefs(listRefs)
  names(theRefs)  <-  gsub('.bib', '', listRefs)
  theRefs
}

exporBibs  <-  function(theRefs=listBibs(), libraryPath='../library.bib', erase=TRUE, verbose=TRUE) {
  if(erase)
    system(paste('rm -r', libraryPath))
  if(!file.exists(libraryPath))
    l_ply(theRefs, function(x, libraryPath, verbose) {
      if(verbose) {
        print(x)
      }
      write.bibtex(x, file=libraryPath, append=TRUE)
    }, libraryPath=libraryPath, verbose=verbose)
}

myCite  <-  function(citationsVec) {
  paste0('[@', paste0(citationsVec, collapse=';@'),']')
}

###############
# PAPER NUMBERS
###############
cleanJagsSummary  <-  function(jagsSummary=tfit$BUGSoutput$summary, modelMatrix=modelLmer1) {
  jagsMat            <-  jagsSummary[paste0('beta[', seq_along(fixef(modelMatrix)), ']'), ]
  rownames(jagsMat)  <-  names(fixef(modelMatrix))
  jagsMat
}

cleanJagsMatrix  <-  function(jagsMatrix=tfit$BUGSoutput$sims.matrix, modelMatrix=modelLmer1) {
  jagsMat            <-  jagsMatrix[, paste0('beta[', seq_along(fixef(modelMatrix)), ']')]
  colnames(jagsMat)  <-  names(fixef(modelMatrix))
  jagsMat
}

get95CIEstimates  <-  function(species) {
  jagsMat  <-  cleanJagsMatrix()
  isReferenceSpecies  <-  species == 'Bryo'
  if(isReferenceSpecies) {
    lnB10  <-  jagsMat[, '(Intercept)']
    s10    <-  jagsMat[, 'lnMass']
    lnB25  <-  rowSums(jagsMat[, c('(Intercept)', 'Temp25')])
    s25    <-  rowSums(jagsMat[, c('lnMass', 'lnMass:Temp25')])
  } else {
    lnB10  <-  rowSums(jagsMat[, c('(Intercept)', paste0('Species', species))])
    s10    <-  rowSums(jagsMat[, c('lnMass', paste0('Species', species, ':', 'lnMass'))])
    lnB25  <-  rowSums(jagsMat[, c('(Intercept)', paste0('Species', species), 'Temp25', paste0('Species', species, ':Temp25'))])
    s25    <-  rowSums(jagsMat[, c('lnMass', paste0('Species', species, ':', 'lnMass'), 'lnMass:Temp25', paste0('Species', species, ':', 'lnMass:Temp25'))])
  }
  data.frame(species=species, 
             'lnB10_2.5'  = quantile(lnB10, probs=0.025, type=2),
             'lnB10_97.5' = quantile(lnB10, probs=0.975, type=2),
             's10_2.5'    = quantile(s10,   probs=0.025, type=2),
             's10_97.5'   = quantile(s10,   probs=0.975, type=2),
             'lnB25_2.5'  = quantile(lnB25, probs=0.025, type=2),
             'lnB25_97.5' = quantile(lnB25, probs=0.975, type=2),
             's25_2.5'    = quantile(s25,   probs=0.025, type=2),
             's25_97.5'   = quantile(s25,   probs=0.975, type=2),
             stringsAsFactors=FALSE)
}

getAverageEstimates  <-  function(species) {
  jagsMat  <-  cleanJagsSummary()
  isReferenceSpecies  <-  species == 'Bryo'
  if(isReferenceSpecies) {
    lnB10  <-  jagsMat['(Intercept)', 'mean']
    s10    <-  jagsMat['lnMass', 'mean']
    lnB25  <-  sum(jagsMat[c('(Intercept)', 'Temp25'), 'mean'])
    s25    <-  sum(jagsMat[c('lnMass', 'lnMass:Temp25'), 'mean'])
  } else {
    lnB10  <-  sum(jagsMat[c('(Intercept)', paste0('Species', species)), 'mean'])
    s10    <-  sum(jagsMat[c('lnMass', paste0('Species', species, ':', 'lnMass')), 'mean'])
    lnB25  <-  sum(jagsMat[c('(Intercept)', paste0('Species', species), 'Temp25', paste0('Species', species, ':Temp25')), 'mean'])
    s25    <-  sum(jagsMat[c('lnMass', paste0('Species', species, ':', 'lnMass'), 'lnMass:Temp25', paste0('Species', species, ':', 'lnMass:Temp25')), 'mean'])
  }
  data.frame(species=species, lnB10=lnB10, s10=s10, lnB25=lnB25, s25=s25, stringsAsFactors=FALSE)
}

getTempDeltas   <-  function(averageEstimates) {
	# deltas between temperatures
	data.frame(species=averageEstimates$species, tempDelta=(exp(mean(averageEstimates$lnB25) + mean(averageEstimates$s25)*2) / exp(mean(averageEstimates$lnB10) + mean(averageEstimates$s10)*2)) / (exp(mean(averageEstimates$lnB25) + mean(averageEstimates$s25)*5) / exp(mean(averageEstimates$lnB10) + mean(averageEstimates$s10)*5)), stringsAsFactors=FALSE)
}

getErEquivalent  <-  function(temp1Kelvin, temp2Kelvin, q10, k=8.62e-5) {
  (k*temp1Kelvin*temp2Kelvin*((temp2Kelvin-temp1Kelvin)/10)*log(q10))/(temp2Kelvin-temp1Kelvin)
}

getq10andEr  <-  function(averageEstimates) {
  # this function, as currently written, cannot be used for Bugula neritina nor Hippopodina sp., because they have interacting slopes between 10 and 25 ˚C
  q10       <-  (exp(averageEstimates$lnB25) / exp(averageEstimates$lnB10)) ^ (10/(25-10))
  er        <-  getErEquivalent(283.15, 298.15, q10)
  q10Boltz  <-  (exp(-er / (8.62e-5 * 298.15)) / exp(-er / (8.62e-5 * 283.15))) ^ (10/(25-10)) # equivalent in Boltzmann terms
  data.frame(species=averageEstimates$species, q10=q10, er=er, stringsAsFactors=FALSE)
}

