library(R2jags)
library(lme4)

rm(list=ls())
source('paths.R')
source('dataManipulation.R')
source('R/functions-analyses.R')

# quickly test for significance of fixed effects using lme4
modelLmer1  <-  lmer(lnRate ~ Species * lnMass * Temp + (1 | Run) + (1 | Plate), data=metRates, REML=TRUE)
modelLmer2  <-  lmer(lnRate ~ Species + lnMass * Temp + (1 | Run) + (1 | Plate), data=metRates, REML=TRUE)
modelLmer3  <-  lmer(lnRate ~ lnMass * Temp + (1 | Run) + (1 | Plate), data=metRates, REML=TRUE)
modelLmer4  <-  lmer(lnRate ~ Temp + (1 | Run) + (1 | Plate), data=metRates, REML=TRUE)
modelLmer5  <-  lmer(lnRate ~ Species * Temp + lnMass + (1 | Run) + (1 | Plate), data=metRates, REML=TRUE)
modelLmer6  <-  lmer(lnRate ~ Species * lnMass + Temp + (1 | Run) + (1 | Plate), data=metRates, REML=TRUE)
anova(modelLmer1, modelLmer2, modelLmer3, modelLmer4) # all significant fixed effects
anova(modelLmer1, modelLmer5) # three-way significantly better
anova(modelLmer1, modelLmer6) # three-way significantly better

# asses confidence intervals for model parameters
confint.merMod(modelLmer1)

# fitting a three-way interaction function using run and plate as random effects on intercept
# first create model matrix for three-way interaction - same order as modelLmer1
tjagsModelMatrix  <-  model.matrix(modelLmer1)
K                 <-  ncol(tjagsModelMatrix)
set.seed(1)
tJags       <-  list('lnRate' = metRates$lnRate,
					 'jagsModelMatrix' = tjagsModelMatrix,
					 'K'     =  K,
					 'run'   =  as.numeric(as.factor(metRates$Run)),
					 'plate' =  as.numeric(as.factor(metRates$Plate)))
tfit        <-  jags(data=tJags, parameters.to.save=c('varB', 'varR', 'varP', paste0('r[', sort(unique(as.numeric(as.factor(metRates$Run)))), ']'), paste0('p[', sort(unique(as.numeric(as.factor(metRates$Plate)))), ']'), paste0('beta[', seq_len(K), ']')), model.file='model_interaction.bug', n.chains=3, n.iter=5e5, DIC=TRUE, n.thin=250)
tfit        <-  autojags(tfit, n.iter=5e5, n.thin=250, n.update=100)
tjagsout    <-  tfit$BUGSoutput$summary #jags output

# Compare output between lme4 and jags
plot(tjagsout[paste0('beta[', seq_along(fixef(modelLmer1)), ']'), 'mean'], fixef(modelLmer1))

# fitting a two-way interaction function using run and plate as random effects on intercept
# first create model matrix for two-way interaction - same order as modelLmer2
tjagsModelMatrix2  <-  model.matrix(modelLmer2)
K                  <-  ncol(tjagsModelMatrix2)
set.seed(1)
tJags2      <-  list('lnRate' = metRates$lnRate,
					 'jagsModelMatrix' = tjagsModelMatrix2,
					 'K'     =  K,
					 'run'   =  as.numeric(as.factor(metRates$Run)),
					 'plate' =  as.numeric(as.factor(metRates$Plate)))
tfit2       <-  jags(data=tJags2, parameters.to.save=c('varB', 'varR', 'varP', paste0('r[', sort(unique(as.numeric(as.factor(metRates$Run)))), ']'), paste0('p[', sort(unique(as.numeric(as.factor(metRates$Plate)))), ']'), paste0('beta[', seq_len(K), ']')), model.file='model_interaction.bug', n.chains=3, n.iter=5e5, DIC=TRUE, n.thin=250)
tfit2       <-  autojags(tfit2, n.iter=5e5, n.thin=250, n.update=100)
tjagsout2   <-  tfit2$BUGSoutput$summary #jags output

# Compare output between lme4 and jags
plot(tjagsout2[paste0('beta[', seq_along(fixef(modelLmer2)), ']'), 'mean'], fixef(modelLmer2))

rm(list=ls()[!(ls() %in% c('metRates', 'modelLmer1', 'modelLmer2', 'tfit', 'tjagsModelMatrix', 'tjagsout', 'tfit2', 'tjagsModelMatrix2', 'tjagsout2'))])
save.image('output/RDatafiles/analyses.RData')
