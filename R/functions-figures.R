######################
# AUXILLIARY FUNCTIONS
######################
toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, family='CM Roman', ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

rounded  <-  function(value, precision=1, change=FALSE) {
  if(change) {
    value  <-  value * -1
  }
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}

###############
# PAPER FIGURES
###############
figureProp  <-  function(x, side='x') {
  if(side == 'x') {
    n  <-  1:2
  } else {
    n  <-  3:4
  }
  usr  <-  par('usr')
  (x-usr[n[1]])/(usr[n[2]]-usr[n[1]])
}

getProb  <-  function(densObjt) {
  densObjt$y  <-  (densObjt$x[2:length(densObjt$x)] - densObjt$x[1:(length(densObjt$x)-1)]) * rowMeans(cbind(densObjt$y[2:length(densObjt$y)], densObjt$y[1:(length(densObjt$y)-1)]))
  densObjt$x  <-  rowMeans(cbind(densObjt$x[2:length(densObjt$x)], densObjt$x[1:(length(densObjt$x)-1)]))
  densObjt
}

plotTriad  <-  function(species=species, spNames=spNames, tmat, x1=TRUE, y1=TRUE, x2=TRUE, y2=TRUE, x3=TRUE, y3=TRUE, x4=TRUE, y4=TRUE) {
  col10     <-  colorRampPalette(rev(brewer.pal(9, 'Blues')))
  col10     <-  col10(32)
  col10[1]  <-  '#FFFFFF'
  col25     <-  colorRampPalette(rev(brewer.pal(9, 'Reds')))
  col25     <-  col25(32)
  col25[1]  <-  '#FFFFFF'

  # all species together - 3-way interaction
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
 
  xDens10    <-  getProb(density(s10))
  xDens25    <-  getProb(density(s25))
  
  if(x1) {
    x1l  <-  substitute('Mass scaling exponent, '*italic(alpha))
  } else {
    x1l  <-  ''
  }

  if(y1) {
    y1l  <-  substitute('Posterior prob. '%*%' 10'^{-3})
  } else {
    y1l  <-  ''
  }

  if(x2) {
    x2l  <-  substitute('ln Mass, '*italic('M'['i'])*' (mg)')
  } else {
    x2l  <-  ''
  }

  if(y2) {
    y2l  <-  substitute('ln Rates, '*italic('B'['i'])*' ('*mu*'l O'[2]*' h'^{-1}*')')
  } else {
    y2l  <-  ''
  }

  if(x3) {
    x3l  <-  substitute('Mass scaling exponent, '*italic(alpha))
  } else {
    x3l  <-  ''
  }

  if(y3) {
    y3l  <-  substitute('ln '*italic('B'['o']))
  } else {
    y3l  <-  ''
  }

  if(x4) {
    x4l  <-  substitute('Posterior prob. '%*%' 10'^{-3})
  } else {
    x4l  <-  ''
  }

  if(y4) {
    y4l  <-  substitute('ln '*italic('B'['o']))
  } else {
    y4l  <-  ''
  }

  plot(NA, xlab='', ylab='', xlim=c(-1, 2), ylim=c(0,0.01), xpd=NA, type='l', axes=FALSE)
  proportionalLabel(0.5, 1.4, x1l, adj=c(0.5, 1), xpd=NA)
  proportionalLabel(-0.4, 0.5, y1l, adj=c(0.5, 0.5), xpd=NA, srt=90)
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  if(x1) axis(3, at=seq(-1, 2, 1))
  axis(2, at=seq(0, 0.008, 0.004), labels=seq(0, 8, 4), las=1)
  box()
  polygon(c(xDens10$x, min(xDens10$x)), c(xDens10$y, min(xDens10$y)), border='dodgerblue2', col=transparentColor('dodgerblue2', 0.3))
  polygon(c(xDens25$x, min(xDens25$x)), c(xDens25$y, min(xDens25$y)), border='tomato', col=transparentColor('tomato', 0.3))
  quants  <-  quantile(s10, probs=c(0.025, 0.975), type=2)
  proportionalLabel(figureProp(quants), rep(0.95, 2), text=FALSE, type='l', col='dodgerblue2', lwd=1.2)
  proportionalLabel(0.97, 0.95, rounded(mean(s10), 2), adj=c(1, 0.5), col='dodgerblue2', cex=0.8)
  quants  <-  quantile(s25, probs=c(0.025, 0.975), type=2)
  proportionalLabel(figureProp(quants), rep(0.86, 2), text=FALSE, type='l', col='tomato', lwd=1.2)
  proportionalLabel(0.97, 0.86, rounded(mean(s25), 2), adj=c(1, 0.5), col='tomato', cex=0.8)

  dat         <-  metRates[metRates$Species == species, ]
  dat$colors  <-  c('dodgerblue2', 'tomato')[match(dat$TempK, c(283.15, 298.15))]
  plot(NA, xlab='', ylab='', xlim=c(-1, 7), ylim=c(-6, 6), axes=FALSE, type='n')
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()  
  points(dat$lnRate ~ dat$lnMass, pch=16, col=transparentColor(dat$colors, 0.6))
  proportionalLabel(0.5, 1.4, x2l, adj=c(0.5, 1), xpd=NA)
  proportionalLabel(1.4, 0.5, y2l, adj=c(0.5, 1), xpd=NA, srt=270)
  if(x2) axis(3, at=seq(0, 6, 2))
  axis(4, at=seq(-5, 5, 5), las=1)
  box()
  xpts  <-  range(dat$lnMass[dat$TempK == 283.15])
  lines(xpts, mean(i10) + mean(s10)*xpts, lwd=2.3)
  lines(xpts, mean(i10) + mean(s10)*xpts, col='dodgerblue2', lwd=2)
  xpts  <-  range(dat$lnMass[dat$TempK == 298.15])
  lines(xpts, mean(i25) + mean(s25)*xpts, lwd=2.3)
  lines(xpts, mean(i25) + mean(s25)*xpts, col='tomato', lwd=2)
  proportionalLabel(0.03, 0.9, spNames, adj=c(0, 0.5), xpd=NA, cex=1)

  den3d10  <-  kde2d(s10, i10, n=500)
  den3d25  <-  kde2d(s25, i25, n=500)
  plot(NA, xlim=c(-1, 2), ylim=c(-10, 7), las=1, ylab='', xlab='', axes=FALSE)
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  proportionalLabel(0.5, -0.4, x3l, adj=c(0.5, 0.5), xpd=NA)
  proportionalLabel(-0.4, 0.5, y3l, adj=c(0.5, 0.5), xpd=NA, srt=90)
  points(s10, i10, pch=16, col=transparentColor('dodgerblue2', 0.1))
  points(s25, i25, pch=16, col=transparentColor('tomato', 0.1))
  image(den3d10, col=transparentColor(col10, c(0, rep(0.1, 31))), add=TRUE)
  image(den3d25, col=transparentColor(col25, c(0, rep(0.1, 31))), add=TRUE)
  if(x3) axis(1, at=seq(-1, 2, 1))
  if(y3) axis(2, at=seq(-10, 5, 5), las=1)
  box()

  xDens10    <-  getProb(density(i10))
  xDens25    <-  getProb(density(i25))
  plot(NA, xlab='', ylab='', xlim=c(0,0.03), ylim=c(-10, 7), xpd=NA, type='l', axes=FALSE)
  proportionalLabel(0.5, -0.4, x4l, adj=c(0.5, 0.5), xpd=NA)
  proportionalLabel(1.4, 0.5, y4l, adj=c(0.5, 1), xpd=NA, srt=270)
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  if(x4) axis(1, seq(0, 0.02, 0.01), labels=seq(0, 20, 10))
  if(y4) axis(4, at=seq(-10, 5, 5), las=1)
  box()
  polygon(c(xDens10$y, min(xDens10$y)), c(xDens10$x, min(xDens10$x)), border='dodgerblue2', col=transparentColor('dodgerblue2', 0.3))
  polygon(c(xDens25$y, min(xDens25$y)), c(xDens25$x, min(xDens25$x)), border='tomato', col=transparentColor('tomato', 0.3))
  quants  <-  quantile(i10, probs=c(0.025, 0.975), type=2)
  proportionalLabel(rep(0.95, 2), figureProp(quants, 'y'), text=FALSE, type='l', col='dodgerblue2', lwd=1.2)
  proportionalLabel(0.95, 0.15, rounded(mean(i10), 2), adj=c(1, 0.5), col='dodgerblue2', cex=0.8)
  quants  <-  quantile(i25, probs=c(0.025, 0.975), type=2)
  proportionalLabel(rep(0.86, 2), figureProp(quants, 'y'), text=FALSE, type='l', col='tomato', lwd=1.2)
  proportionalLabel(0.86, 0.06, rounded(mean(i25), 2), adj=c(1, 0.5), col='tomato', cex=0.8)
}

fig1  <-  function() {
  tmat  <-  model$BUGSoutput$sims.matrix
  tmat  <-  tmat[,paste0('beta[', seq_len(ncol(coef(lmerModel1)$Run)), ']')]
  colnames(tmat)  <-  names(coef(lmerModel1)$Run)
  species  <-  unique(metRates$Species)
  # species real names
  spNames  <-  list(substitute(italic('Hippopodina iririkiensis')),
                    'Microcionidae',
                    substitute(italic('Bugula neritina')),
                    substitute(italic('Bugula stolonifera')))

  figMat  <-  matrix(
                      c(
                      rep(c(rep(1, 4), rep(2, 4), rep(5, 2), rep(6, 4), rep(7, 4)), 4),
                      rep(c(rep(3, 4), rep(4, 4), rep(5, 2), rep(8, 4), rep(9, 4)), 4),
                      rep(10, 36),
                      rep(c(rep(11, 4), rep(12, 4), rep(15, 2), rep(16, 4), rep(17, 4)), 4),
                      rep(c(rep(13, 4), rep(14, 4), rep(15, 2), rep(18, 4), rep(19, 4)), 4)
                      ),
                      nrow=18, ncol=18, byrow=TRUE)

  par(omi=rep(1, 4), mai=rep(0, 4), cex=1)
  layout(figMat)
  plotTriad(species[1], spNames[[1]], tmat, x1=TRUE, y1=TRUE, x2=TRUE, y2=FALSE, x3=FALSE, y3=TRUE, x4=FALSE, y4=FALSE)
  plot(0, 0, axes=FALSE, type='n', xlab='', ylab='')
  plotTriad(species[2], spNames[[2]], tmat, x1=TRUE, y1=FALSE, x2=TRUE, y2=TRUE, x3=FALSE, y3=FALSE, x4=FALSE, y4=TRUE)
  plot(0, 0, axes=FALSE, type='n', xlab='', ylab='')
  plotTriad(species[3], spNames[[3]], tmat, x1=FALSE, y1=TRUE, x2=FALSE, y2=FALSE, x3=TRUE, y3=TRUE, x4=TRUE, y4=FALSE)
  plot(0, 0, axes=FALSE, type='n', xlab='', ylab='')
  plotTriad(species[4], spNames[[4]], tmat, x1=FALSE, y1=FALSE, x2=FALSE, y2=TRUE, x3=TRUE, y3=FALSE, x4=TRUE, y4=TRUE)
}
