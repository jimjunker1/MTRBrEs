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
# DATA ANALYSIS
###############
jagsModel  <-  function(run=TRUE, outputFileName) {
  jmod                                <-  vector(mode='character')
  jmod[(length(jmod) + 1)]  <-  'model {'
  jmod[(length(jmod) + 1)]  <-  '\t# priors for fixed effects in metabolic rate analysis'
  jmod[(length(jmod) + 1)]  <-  '\tfor (i in 1:K) {'
  jmod[(length(jmod) + 1)]  <-  '\t\tbeta[i] ~ dnorm(0, 0.0001)'
  jmod[(length(jmod) + 1)]  <-  '\t}'
  jmod[(length(jmod) + 1)]  <-  '\ttauB  ~   dgamma(1.0E-3,1.0E-3)'
  jmod[(length(jmod) + 1)]  <-  '\tvarB  <-  1/tauB # variance, equivalent to inverse(tauB)'
  jmod[(length(jmod) + 1)]  <-  '\n\t#linear regression for metabolic rate'
  jmod[(length(jmod) + 1)]  <-  '\tfor(i in 1:length(lnRate)) {'
  jmod[(length(jmod) + 1)]  <-  '\t\tlnRate[i]  ~  dnorm(muB[i], tauB)'
  jmod[(length(jmod) + 1)]  <-  '\t\tmuB[i]  <-  inprod(beta[], tjagsModelMatrix[i,])'
  jmod[(length(jmod) + 1)]  <-  '\t\tlogLik[i]  <- - log(tauB) + log(2*pi)  + pow(lnRate[i]-muB[i],2) * tauB'
  jmod[(length(jmod) + 1)]  <-  '\t}'
  
  if(run) {
      jmod                      <-  sub('muB[i]  <-  ', 'muB[i]  <-  r[runJ[i]] + ', jmod, fixed=TRUE)
      jmod[(length(jmod) + 1)]  <-  '\n\t# priors for random effects at the run level'
      jmod[(length(jmod) + 1)]  <-  '\tfor(i in 1:max(runJ)) {'
      jmod[(length(jmod) + 1)]  <-  '\t\tr[i] ~ dnorm(0, tauR)'
      jmod[(length(jmod) + 1)]  <-  '\t}'
      jmod[(length(jmod) + 1)]  <-  '\ttauR  ~   dgamma(1.0E-3,1.0E-3)'
      jmod[(length(jmod) + 1)]  <-  '\tvarR  <-  1/tauR # variance, equivalent to inverse(tauR)'
  }
  jmod[(length(jmod) + 1)]  <-  '\n}'
  write(jmod, outputFileName)
}

fitJags  <-  function(lmModelMatrix, run=TRUE, data=metRates, outputFileName='mod.bug') {
  # fitting an interaction model
  # first extract model matrix from lm model
  tjagsModelMatrix  <-  model.matrix(lmModelMatrix)
  K                 <-  ncol(tjagsModelMatrix)
  set.seed(1)
  parToSave  <-  c('varB', paste0('beta[', seq_len(K), ']'), paste0('logLik[', seq_len(nrow(tjagsModelMatrix)), ']'))
  lnRate     <-  data$lnRate
  tJags      <-  list('lnRate', 'tjagsModelMatrix', 'K', 'pi')
  if(run) {
    parToSave  <-  c(parToSave, c('varR', paste0('r[', sort(unique(as.numeric(as.factor(data$Run)))), ']')))
    runJ       <-  as.numeric(as.factor(data$Run))
    tJags      <-  append(tJags, list('runJ'))
  }
  jagsModel(run, outputFileName)
  tfit        <-  jags.parallel(data=tJags, parameters.to.save=parToSave, model.file=outputFileName, n.chains=5, n.iter=1e6, n.thin=250, envir=environment())
  system(paste0('rm ', outputFileName))
  recompile(tfit)
  tfit        <-  autojags(tfit, n.iter=5e5, n.thin=250, n.update=100)
  tfit
}

modelComparisonPVals  <-  function(data, modelList, verbose=TRUE) {
  modDiff  <-  compare(modelList[[data$complexModel]]$loo, modelList[[data$nestedModel]]$loo)
  if(verbose)
    print(modDiff)
  data.frame(zPVal = 2*pnorm(-abs((modDiff[['elpd_diff']] - 0) / modDiff[['se']])),
         tPVal = 2*pt(abs((modDiff[['elpd_diff']] - 0) / modDiff[['se']]), df=402-1, lower.tail=FALSE))
}

cleanPvals  <-  function(x) {
  if(x <= 0.05 & x > 0.01)
    return('<= 0.05')
  else if(x <= 0.01 & x > 0.001)
    return('<= 0.01')
  else if(x <= 0.001 & x > 0.0001)
    return('<= 0.001')
  else if(x <= 0.0001)
    return('<= 0.0001')
  else
    paste0('= ', round(x, 2))
}

###############
# PAPER NUMBERS
###############
cleanJagsSummary  <-  function(jagsSummary=fixedSelec[[1]]$jagsFit$BUGSoutput$summary, modelMatrix=lmModel1) {
  jagsMat            <-  jagsSummary[paste0('beta[', seq_along(coef(modelMatrix)), ']'), ]
  rownames(jagsMat)  <-  names(coef(modelMatrix))
  jagsMat
}

cleanJagsMatrix  <-  function(jagsMatrix=fixedSelec[[1]]$jagsFit$BUGSoutput$sims.matrix, modelMatrix=lmModel1) {
  jagsMat            <-  jagsMatrix[, paste0('beta[', seq_along(coef(modelMatrix)), ']')]
  colnames(jagsMat)  <-  names(coef(modelMatrix))
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

getMetabolicDiff  <-  function(averageEstimates, lnMass) {
  data.frame(species       =  averageEstimates$species, 
            metabolicDiff  =  exp(averageEstimates$lnB25 + averageEstimates$s25 * lnMass) / exp(averageEstimates$lnB10 + averageEstimates$s10 * lnMass),
            stringsAsFactors=FALSE)
}
