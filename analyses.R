library(plyr)
library(R2jags)
library(lme4)
library(bbmle)

rm(list=ls())
source('paths.R')
source('dataManipulation.R')
source('R/functions-analyses.R')

# fit model selection - is run a significant random effect?
# http://stackoverflow.com/questions/24019807/how-to-compare-a-model-with-no-random-effects-to-a-model-with-a-random-effect-us
modelLmer1  <-  lmer(lnRate ~ Species * lnMass * Temp + (1 | Run), data=metRates, REML=FALSE)
modelLm1    <-  lm(lnRate ~ Species * lnMass * Temp, data=metRates)
AICtab(modelLmer1, modelLm1) # Delta AIC > 2, run is a significant random effect

# refit modelLmer1 using REML=TRUE
modelLmer  <-  lmer(lnRate ~ Species * lnMass * Temp + (1 | Run), data=metRates, REML=TRUE)
modelLmer2  <-  lmer(lnRate ~ lnMass * Temp + (1 | Run), data=metRates, REML=TRUE)
modelLmer3  <-  lmer(lnRate ~ Temp + (1 | Run), data=metRates, REML=TRUE)
anova(modelLmer, modelLmer2, modelLmer3) # all significant fixed effects

# asses confidence intervals for model parameters
confint.merMod(modelLmer)

# refit modelLmer using Boltzmann function instead
boltzLmer       <-  lmer(lnRate ~ Species + lnMass + invKT + (1 | Run), data=metRates)

#############
# FIT IN JAGS
#############
# fitting the Boltzmann function
set.seed(1)
bJags       <-  list('lnRate'=metRates$lnRate, 'lnMass'=metRates$lnMass, 'species'=as.numeric(as.factor(metRates$Species)), 'run'=as.numeric(as.factor(metRates$Run)), 'invKT'=metRates$invKT)
binits1     <-  list(boTs = -2, A = 0.35, Er = 0.5)
binits2     <-  list(boTs = -1, A = 0.55, Er = 0.6)
binits3     <-  list(boTs =  0, A = 0.75, Er = 0.7)
bfit        <-  jags(data=bJags, inits=list(binits1, binits2, binits3), parameters.to.save=c('boTs', 'A', 'Er', 'varB', 'varS', 'varR', as.character(sapply(1:3, function(x)paste0('s[', sort(unique(as.numeric(as.factor(metRates$Species)))), ',', x,']'))), paste0('r[', sort(unique(as.numeric(as.factor(metRates$Run)))), ']')), model.file='model_boltz.bug', n.chains=3, n.iter=5e5, DIC=TRUE, n.thin=250)
bfit        <-  autojags(bfit, n.iter=5e5, n.thin=250, n.update=100)
bjagsout    <-  bfit$BUGSoutput$summary #jags output

# fitting a three-way interaction function using run as a random effect on intercept
# first create model matrix for three-way interaction - same order as modelLmer
tjagsModelMatrix  <-  model.matrix(modelLmer)
K                 <-  ncol(tjagsModelMatrix)
set.seed(1)
tJags       <-  list('lnRate' = metRates$lnRate,
					 'jagsModelMatrix' = tjagsModelMatrix,
					 'K'  =  K,
					 'run' = as.numeric(as.factor(metRates$Run)))
tfit        <-  jags(data=tJags, parameters.to.save=c('varB', 'varR', paste0('r[', sort(unique(as.numeric(as.factor(metRates$Run)))), ']'), paste0('beta[', seq_len(K), ']')), model.file='model_three_way_interaction.bug', n.chains=3, n.iter=5e5, DIC=TRUE, n.thin=250)
tfit        <-  autojags(tfit, n.iter=5e5, n.thin=250, n.update=100)
tjagsout    <-  tfit$BUGSoutput$summary #jags output

# Compare output between lme4 and jags
plot(tjagsout[paste0('beta[', seq_along(fixef(modelLmer)), ']'), 'mean'], fixef(modelLmer))

# Now fit data for each species separately. Since each species only has 2 runs maximum, treat as a run-specific deviation on intercept (i.e. contr.sum) -- somehow analogously equivalent to random effects above
perSpeciesLmFits    <-  dlply(metRates, .(Species), perSpeciesLmRun)
perSpeciesJagsFits  <-  dlply(metRates, .(Species), perSpeciesJagsRun)

# Compare output between lm and jags
comparedOutputs  <-  data.frame()
for(i in 1:4) {
	species    <-  names(perSpeciesLmFits)[i]
	lme4Coefs  <-  coef(summary(perSpeciesLmFits[[i]]))[,'Estimate']
	jagsCoefs  <-  perSpeciesJagsFits[[i]]$jagsOutput[names(lme4Coefs), 'mean']
	comparedOutputs  <-  rbind(comparedOutputs, data.frame(species=species, lme4Coefs=lme4Coefs, jagsCoefs=jagsCoefs, stringsAsFactors=FALSE))
}
plot(comparedOutputs$lme4Coefs, comparedOutputs$jagsCoefs)

rm(list=ls()[!(ls() %in% c('bfit', 'bjagsout', 'boltzLmer', 'modelLmer', 'tfit', 'tjagsModelMatrix', 'tjagsout', 'perSpeciesLmFits', 'perSpeciesJagsFits'))])
save.image('output/RDatafiles/analyses.RData')
