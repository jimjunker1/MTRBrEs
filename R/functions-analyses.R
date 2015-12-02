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
	if(enoughRuns) {
		dat$Run  <-  as.factor(dat$Run)
		lm(lnRate ~ Run + lnMass * Temp, data=dat, contrasts=list(Run='contr.sum'))
	} else {
		lm(lnRate ~ lnMass * Temp, data=dat)
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
	bug.file[length(bug.file)+1]  <-  	'\n\t#linear regression for metabolic rate'
	bug.file[length(bug.file)+1]  <-  	'\tfor(i in 1:length(lnRate)) {'
	bug.file[length(bug.file)+1]  <-  	'\t\tlnRate[i]  ~  dnorm(muB[i], tauB)'
	bug.file[length(bug.file)+1]  <-  	'\t\tmuB[i]     <- inprod(beta[], jagsModelMatrix[i,])'
	bug.file[length(bug.file)+1]  <-  	'\t}'
	bug.file[length(bug.file)+1]  <-  	'\n}'
	write(bug.file, paste0('model_two_way_interaction_', species, '.bug'))
}

perSpeciesJagsRun  <-  function(dat) {
	enoughRuns  <-  length(unique(dat$Run)) > 1
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
	tJags  <-  list('lnRate'=dat$lnRate, 'jagsModelMatrix'=modelMatrix, 'K'=K)
	tfit   <-  jags(data=tJags, parameters.to.save=c('varB', paste0('beta[', seq_len(K), ']')), model.file=bugsFile, n.chains=3, n.iter=5e5, DIC=TRUE, n.thin=250)
	tfit   <-  autojags(tfit, n.iter=5e5, n.thin=250, n.update=100)
	simsMatrix  <-  tfit$BUGSoutput$sims.matrix
	colNames    <-  grep('beta', colnames(simsMatrix), value=TRUE)
	colnames(simsMatrix)[match(colNames, paste0('beta[', seq_len(K), ']'))]  <-  colnames(modelMatrix)
	jagsOutput  <-  tfit$BUGSoutput$summary
	rowNames    <-  grep('beta', rownames(jagsOutput), value=TRUE)
	rownames(jagsOutput)[match(rowNames, paste0('beta[', seq_len(K), ']'))]  <-  colnames(modelMatrix)
	list(simsMatrix=simsMatrix, jagsOutput=jagsOutput)
}
