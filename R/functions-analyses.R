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

############
# MODEL FITS
############
perSpeciesLmRun  <-  function(dat) {
	enoughRuns  <-  length(unique(dat$Run)) > 1
	dat$Plate   <-  as.factor(dat$Plate)
	if(enoughRuns) {
		dat$Run    <-  as.factor(dat$Run)		
		lm(lnRate ~ Plate + Run + lnMass * Temp, data=dat, contrasts=list(Run='contr.sum', Plate='contr.sum'))
	} else {
		lm(lnRate ~ Plate + lnMass * Temp, data=dat, contrasts=list(Plate='contr.sum'))
	}
}

modelJags  <-  function(species) {
	bug.file                      <-  vector(mode='character')
	bug.file[length(bug.file)+1]  <-  	'model {'
	bug.file[length(bug.file)+1]  <-  	'\n\t# priors for fixed effects in metabolic rate analysis'
	bug.file[length(bug.file)+1]  <-  	'\tfor(i in 1:K) { '
	bug.file[length(bug.file)+1]  <-  	'\t\tbeta[i] ~ dnorm(0, 0.0001)'
	bug.file[length(bug.file)+1]  <-  	'\t}'
	bug.file[length(bug.file)+1]  <-  	'\ttauB  ~   dgamma(1.0E-3,1.0E-3)'
	bug.file[length(bug.file)+1]  <-  	'\tvarB  <-  1/tauB # variance, equivalent to inverse(tauB)'
	bug.file[length(bug.file)+1]  <-  	'\tfor(i in 1:max(plate)) { '
	bug.file[length(bug.file)+1]  <-  	'\t\tr[i] ~ dnorm(0, tauR)'
	bug.file[length(bug.file)+1]  <-  	'\t}'
	bug.file[length(bug.file)+1]  <-  	'\ttauR  ~   dgamma(1.0E-3,1.0E-3)'
	bug.file[length(bug.file)+1]  <-  	'\tvarR  <-  1/tauR # variance, equivalent to inverse(tauR)'
	bug.file[length(bug.file)+1]  <-  	'\n\t#linear regression for metabolic rate'
	bug.file[length(bug.file)+1]  <-  	'\tfor(i in 1:length(lnRate)) {'
	bug.file[length(bug.file)+1]  <-  	'\t\tlnRate[i]  ~  dnorm(muB[i], tauB)'
	bug.file[length(bug.file)+1]  <-  	'\t\tmuB[i]     <- r[plate[i]] + inprod(beta[], jagsModelMatrix[i,])'
	bug.file[length(bug.file)+1]  <-  	'\t}'
	bug.file[length(bug.file)+1]  <-  	'\n}'
	write(bug.file, paste0('model_two_way_interaction_', species, '.bug'))
}

perSpeciesJagsRun  <-  function(dat) {
	enoughRuns  <-  length(unique(dat$Run)) > 1
	dat$Plate   <-  as.factor(dat$Plate)
	if(enoughRuns) {
		dat$Run      <-  as.factor(dat$Run)
		modelMatrix  <-  model.matrix(~ Run + lnMass * Temp, data=dat, contrasts=list(Run='contr.sum'))
	} else {
		modelMatrix  <-  model.matrix(~ lnMass * Temp, data=dat)
	}

	species   <-  unique(dat$Species)
	bugsFile  <-  paste0('model_two_way_interaction_', species, '.bug')
	
	modelJags(species)
	set.seed(1)
	K      <-  ncol(modelMatrix)
	tJags  <-  list('lnRate'=dat$lnRate, 'jagsModelMatrix'=modelMatrix, 'K'=K, 'plate'=as.numeric(dat$Plate))
	tfit   <-  jags(data=tJags, parameters.to.save=c('varB', paste0('beta[', seq_len(K), ']'), paste0('r[', seq_len(length(unique(dat$Plate))), ']')), model.file=bugsFile, n.chains=3, n.iter=5e5, DIC=TRUE, n.thin=250)
	tfit   <-  autojags(tfit, n.iter=5e5, n.thin=250, n.update=100)
	simsMatrix  <-  tfit$BUGSoutput$sims.matrix
	colNames    <-  grep('beta', colnames(simsMatrix), value=TRUE)
	colnames(simsMatrix)[match(colNames, paste0('beta[', seq_len(K), ']'))]  <-  colnames(modelMatrix)
	jagsOutput  <-  tfit$BUGSoutput$summary
	rowNames    <-  grep('beta', rownames(jagsOutput), value=TRUE)
	rownames(jagsOutput)[match(rowNames, paste0('beta[', seq_len(K), ']'))]  <-  colnames(modelMatrix)
	list(simsMatrix=simsMatrix, jagsOutput=jagsOutput)
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

