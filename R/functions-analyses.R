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
getTempDeltas   <-  function(species) {
	tmat  <-  tfit$BUGSoutput$sims.matrix
	tmat  <-  tmat[, paste0('beta[', seq_along(fixef(modelLmer1)), ']')]
	colnames(tmat)  <-  names(fixef(modelLmer1))

	isReferenceSpecies  <-  species == 'Bryo'
	if(isReferenceSpecies) {
		i10  <-  tmat[, '(Intercept)']
		s10  <-  tmat[, 'lnMass']
		i25  <-  rowSums(tmat[, c('(Intercept)', 'Temp25')])
		s25  <-  rowSums(tmat[, c('lnMass', 'lnMass:Temp25')])
	} else {
		i10  <-  rowSums(tmat[, c('(Intercept)', paste0('Species', species))])
		s10  <-  rowSums(tmat[, c('lnMass', paste0('Species', species, ':', 'lnMass'))])
		i25  <-  rowSums(tmat[, c('(Intercept)', paste0('Species', species), 'Temp25', paste0('Species', species, ':Temp25'))])
		s25  <-  rowSums(tmat[, c('lnMass', paste0('Species', species, ':', 'lnMass'), 'lnMass:Temp25', paste0('Species', species, ':', 'lnMass:Temp25'))])
	}
	
	# deltas between temperatures
	data.frame(species=species, tempDelta=(exp(mean(i25) + mean(s25)*2) / exp(mean(i10) + mean(s10)*2)) / (exp(mean(i25) + mean(s25)*5) / exp(mean(i10) + mean(s10)*5)), stringsAsFactors=FALSE)
}
