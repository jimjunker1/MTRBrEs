######
# DATA
######
readAndCleanData  <-  function(dataPath) {
    metRates  <-  read.csv(dataPath, header=TRUE, stringsAsFactors=FALSE)
    metRates  <-  metRates[complete.cases(metRates), ]
    
    # correct log values from log 10 to ln
    metRates$lnMass  <-  log(10^metRates$logM)
    metRates$lnRate  <-  log(10^metRates$logVo2)
    
    # create inverse temperature column
    metRates$TempK  <-  metRates$Temp + 273.15
    metRates$invKT  <-  1 / 8.62e-5 * (1 / 293.15 - 1 / metRates$TempK)
    
    # factor temperature column for ancova
    metRates$Temp         <-  as.factor(metRates$Temp)
    metRates
}

########
# TABLES
########
makeTable1  <-  function(data, dest) {
    tab1  <-  ddply(data, .(Species, Temp), function(x) {
        data.frame(minMass = rounded(min(exp(x$lnMass))),
                   maxMass = rounded(max(exp(x$lnMass))),
                   n         = nrow(x), stringsAsFactors=FALSE)
    })
    write.csv(tab1, dest, row.names=FALSE)
}

makeTable2  <-  function(data, dest) {
    tab2  <-  data.frame()
    for(j in 2:9) {
        refModel  <-  data[[paste0('lmerModel', j)]]
        aovTab    <-  data.frame(anova(data$lmerModel1, refModel))['data$lmerModel1',]
        tab2      <-  rbind(tab2, data.frame('d.f.'=aovTab[,'Chi.Df'], 'Chisq'=aovTab[,'Chisq'], 'P'=aovTab[,'Pr..Chisq.'], stringsAsFactors=FALSE))
    }
    tab2$P  <-  sapply(tab2$P, cleanPvals)
    write.csv(tab2, dest, row.names=FALSE)
}

cleanPvals  <-  function(x) {
    if(x < 0.05 & x > 0.01)
        return('< 0.05')
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
# DATA ANALYSIS
###############
produceAnalysesOutputs  <-  function(metRates) {
    # compare random effects among models
    lmerModela  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + Temp:lnMass + Temp:lnMass:Species + (1|Run) + (1|Plate), data=metRates, REML=FALSE)
    lmerModelb  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + Temp:lnMass + Temp:lnMass:Species + (1|Plate), data=metRates, REML=FALSE)
    lmerModelc  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + Temp:lnMass + Temp:lnMass:Species + (1|Run), data=metRates, REML=FALSE)
    modelSelectionRandomTab  <-  anova(lmerModela, lmerModelb, lmerModelc)
    
    # lmerModelc most significant. Now run models to asses significance of fixed effects
    lmerModel1  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + Temp:lnMass + Temp:lnMass:Species + (1|Run), data=metRates, REML=FALSE)
    lmerModel2  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + Temp:lnMass + (1|Run), data=metRates, REML=FALSE)
    lmerModel3  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Species:Temp + (1|Run), data=metRates, REML=FALSE)
    lmerModel4  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + Temp:lnMass + (1|Run), data=metRates, REML=FALSE)
    lmerModel5  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:Temp + Temp:lnMass + (1|Run), data=metRates, REML=FALSE)
    lmerModel6  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:lnMass + (1|Run), data=metRates, REML=FALSE)
    lmerModel7  <-  lmer(lnRate ~ Species + lnMass + Temp + Species:Temp + (1|Run), data=metRates, REML=FALSE)
    lmerModel8  <-  lmer(lnRate ~ Species + lnMass + Temp + Temp:lnMass + (1|Run), data=metRates, REML=FALSE)
    lmerModel9  <-  lmer(lnRate ~ Species + lnMass + Temp + (1|Run), data=metRates, REML=FALSE)
    modelSelectionFixedTab  <-  anova(lmerModel1,lmerModel2,lmerModel3,lmerModel4,lmerModel5,lmerModel6,lmerModel7,lmerModel8,lmerModel9)
    
    # full model is the most significant one. Refit in jags
    model       <-  fitJags(lmerModel1, data=metRates)
    list('lmerModela' = lmerModela, 'lmerModelb' = lmerModelb, 'lmerModelc' = lmerModelc, 'modelSelectionRandomTab' = modelSelectionRandomTab, 'lmerModel1' = lmerModel1, 'lmerModel2' = lmerModel2, 'lmerModel3' = lmerModel3, 'lmerModel4' = lmerModel4, 'lmerModel5' = lmerModel5, 'lmerModel6' = lmerModel6, 'lmerModel7' = lmerModel7, 'lmerModel8' = lmerModel8, 'lmerModel9' = lmerModel9, 'modelSelectionFixedTab' = modelSelectionFixedTab, 'model' = model)
}

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

fitJags  <-  function(lmModelMatrix, data, run=TRUE, outputFileName='mod.bug') {
    # fitting an interaction model
    # first extract model matrix from lm model
    tjagsModelMatrix  <-  model.matrix(lmModelMatrix)
    K                 <-  ncol(tjagsModelMatrix)
    set.seed(1)
    parToSave  <-  c('varB', paste0('beta[', seq_len(K), ']'))
    lnRate     <-  data$lnRate
    tJags      <-  list('lnRate', 'tjagsModelMatrix', 'K')
    if(run) {
        parToSave  <-  c(parToSave, c('varR', paste0('r[', sort(unique(as.numeric(as.factor(data$Run)))), ']')))
        runJ       <-  as.numeric(as.factor(data$Run))
        tJags      <-  append(tJags, list('runJ'))
    }
    jagsModel(run, outputFileName)
    tfit        <-  jags(data=tJags, parameters.to.save=parToSave, model.file=outputFileName, n.chains=3, n.iter=5e5, n.thin=250)
    system(paste0('rm ', outputFileName))
    tfit        <-  autojags(tfit, n.iter=5e5, n.thin=250, n.update=100)
    tfit
}

###############
# PAPER NUMBERS
###############
extractNumbersList  <-  function(metRates, output) {
    # get average estimates and 95% CI for each parameter
    species           <-  unique(metRates$Species)
    averageEstimates  <-  ldply(species, getAverageEstimates, jagsSummary=output$model$BUGSoutput$summary, modelMatrix=output$lmerModel1)
    ci95Estimates     <-  ldply(species, get95CIEstimates, jagsMatrix=output$model$BUGSoutput$sims.matrix, modelMatrix=output$lmerModel1)
    
    list('metRates'         =  metRates,
         'species'          =  species,          
         'averageEstimates' =  averageEstimates, 
         'ci95Estimates'    =  ci95Estimates,
         # compare deltas between temperatures for each species
         'tempDeltas'       =  ddply(averageEstimates, .(species), getTempDeltas),
         # metabolic difference (in n-fold) between different temperatures for similar-sized animals   
         'metabolicDiff'    =  ddply(averageEstimates, .(species), getMetabolicDiff, lnMass=log(400)),
         # q10-equivalents [ q10 == (r25 / r10) ^ (10/(25-10)) ] for Microcionidae and Bugula stolonifera
         'q10andErs'        =  ddply(averageEstimates[averageEstimates$species %in% c('Sponge', 'Stolonifera'), ], .(species), getq10andEr)
    )
}

getAverageEstimates  <-  function(species, ...) {
    jagsMat  <-  cleanJagsSummary(...)
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

cleanJagsSummary  <-  function(jagsSummary, modelMatrix) {
    jagsMat            <-  jagsSummary[paste0('beta[', seq_len(ncol(coef(modelMatrix)$Run)), ']'), ]
    rownames(jagsMat)  <-  names(coef(modelMatrix)$Run)
    jagsMat
}

get95CIEstimates  <-  function(species, ...) {
    jagsMat  <-  cleanJagsMatrix(...)
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

cleanJagsMatrix  <-  function(jagsMatrix, modelMatrix) {
    jagsMat            <-  jagsMatrix[, paste0('beta[', seq_len(ncol(coef(modelMatrix)$Run)), ']')]
    colnames(jagsMat)  <-  names(coef(modelMatrix)$Run)
    jagsMat
}

getTempDeltas   <-  function(averageEstimates) {
    # deltas between temperatures
    data.frame(species=averageEstimates$species, tempDelta=(exp(mean(averageEstimates$lnB25) + mean(averageEstimates$s25)*2) / exp(mean(averageEstimates$lnB10) + mean(averageEstimates$s10)*2)) / (exp(mean(averageEstimates$lnB25) + mean(averageEstimates$s25)*5) / exp(mean(averageEstimates$lnB10) + mean(averageEstimates$s10)*5)), stringsAsFactors=FALSE)
}

getMetabolicDiff  <-  function(averageEstimates, lnMass) {
    data.frame(species       =  averageEstimates$species, 
               metabolicDiff  =  exp(averageEstimates$lnB25 + averageEstimates$s25 * lnMass) / exp(averageEstimates$lnB10 + averageEstimates$s10 * lnMass),
               stringsAsFactors=FALSE)
}

getq10andEr  <-  function(averageEstimates) {
    # this function, as currently written, cannot be used for Bugula neritina nor Hippopodina sp., because they have interacting slopes between 10 and 25 ˚C
    q10       <-  (exp(averageEstimates$lnB25) / exp(averageEstimates$lnB10)) ^ (10/(25-10))
    er        <-  getErEquivalent(283.15, 298.15, q10)
    q10Boltz  <-  (exp(-er / (8.62e-5 * 298.15)) / exp(-er / (8.62e-5 * 283.15))) ^ (10/(25-10)) # equivalent in Boltzmann terms
    data.frame(species=averageEstimates$species, q10=q10, er=er, stringsAsFactors=FALSE)
}

getErEquivalent  <-  function(temp1Kelvin, temp2Kelvin, q10, k=8.62e-5) {
    (k*temp1Kelvin*temp2Kelvin*((temp2Kelvin-temp1Kelvin)/10)*log(q10))/(temp2Kelvin-temp1Kelvin)
}

#############################
# DEALING WITH BIB REFERENCES
#############################
myCite  <-  function(citationsVec) {
    paste0('[@', paste0(citationsVec, collapse=';@'),']')
}

