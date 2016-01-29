######################
# AUXILLIARY FUNCTIONS
######################
make.transparent <- function(col, opacity=0.5) {
  if (length(opacity) > 1 && any(is.na(opacity))) {
    n <- max(length(col), length(opacity))
    opacity <- rep(opacity, length.out=n)
    col <- rep(col, length.out=n)
    ok <- !is.na(opacity)
    ret <- rep(NA, length(col))
    ret[ok] <- Recall(col[ok], opacity[ok])
    ret
  } else {
    tmp <- col2rgb(col)/255
    rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
  }
}

label <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
  usr  <-  par('usr')
  x.p  <-  usr[1] + px*(usr[2] - usr[1])
  y.p  <-  usr[3] + py*(usr[4] - usr[3])
  if(log=='x'){x.p<-10^(x.p)}
  if(log=='y'){y.p<-10^(y.p)}
  if(log=='xy'){x.p<-10^(x.p);y.p<-10^(y.p)}
  if(text){
    text(x.p, y.p, lab, adj=adj, ...)
  } else {
    points(x.p, y.p, ...)
  }
}

to.pdf <- function(expr, filename, ...) {
  to.dev(expr, pdf, filename, ...)
}

fig.path  <-  function(name) {
  file.path('output/figures', name)
}

to.dev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, family='CM Roman', ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

to.pdf <- function(expr, filename, ...) {
  to.dev(expr, pdf, filename, ...)
}

rounded  <-  function(value, precision=1, change=FALSE) {
  if(change) {
    value  <-  value * -1
  }
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}

whiteGrid  <-  function() {
  label(rep(0.2, 2), c(0,1), text=FALSE, type='l', col='white', lwd=0.5)
  label(rep(0.4, 2), c(0,1), text=FALSE, type='l', col='white', lwd=0.5)
  label(rep(0.6, 2), c(0,1), text=FALSE, type='l', col='white', lwd=0.5)
  label(rep(0.8, 2), c(0,1), text=FALSE, type='l', col='white', lwd=0.5)
  label(c(0,1), rep(0.2, 2), text=FALSE, type='l', col='white', lwd=0.5)
  label(c(0,1), rep(0.4, 2), text=FALSE, type='l', col='white', lwd=0.5)
  label(c(0,1), rep(0.6, 2), text=FALSE, type='l', col='white', lwd=0.5)
  label(c(0,1), rep(0.8, 2), text=FALSE, type='l', col='white', lwd=0.5)
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

plotTriad  <-  function(species=species, spNames=spNames, tmat, y1=TRUE, x1=TRUE, y2=TRUE, x2=TRUE, y3=TRUE, x3=TRUE) {
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
  
  if(y1) {
    y1  <-  substitute('Posterior probability '%*%' 10'^{-3})
  } else {
    y1  <-  ''
  }

  if(x1) {
    x1  <-  substitute('Mass scaling exponent, '*italic(alpha))
  } else {
    x1  <-  ''
  }

  if(y2) {
    y2  <-  substitute('Metabolic normalisation, '*italic('B'['o']))
  } else {
    y2  <-  ''
  }

  if(x2) {
    x2  <-  substitute('Posterior probability '%*%' 10'^{-3})
  } else {
    x2  <-  ''
  }

  if(y3) {
    y3  <-  substitute('Metabolic rates, '*italic('B'['i']))
  } else {
    y3  <-  ''
  }

  if(x3) {
    x3  <-  'ln Mass (mg)'
  } else {
    x3  <-  ''
  }

  plot(NA, xlab='', ylab='', xlim=c(-1, 2), ylim=c(0,0.01), xpd=NA, type='l', axes=FALSE)
  label(-0.4, 0.5, y1, adj=c(0.5, 0.5), xpd=NA, srt=90)
  label(0.5, 1.4, x1, adj=c(0.5, 1), xpd=NA)
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  axis(2, at=seq(0, 0.008, 0.004), labels=seq(0, 8, 4), las=1)
  axis(3, at=seq(-1, 2, 1))
  box()
  polygon(c(xDens10$x, min(xDens10$x)), c(xDens10$y, min(xDens10$y)), border='dodgerblue2', col=make.transparent('dodgerblue2', 0.3))
  polygon(c(xDens25$x, min(xDens25$x)), c(xDens25$y, min(xDens25$y)), border='tomato', col=make.transparent('tomato', 0.3))
  quants  <-  quantile(s10, probs=c(0.025, 0.975), type=2)
  label(figureProp(quants), rep(0.95, 2), text=FALSE, type='l', col='dodgerblue2', lwd=1.2)
  label(0.97, 0.95, rounded(mean(s10), 2), adj=c(1, 0.5), col='dodgerblue2', cex=0.8)
  quants  <-  quantile(s25, probs=c(0.025, 0.975), type=2)
  label(figureProp(quants), rep(0.86, 2), text=FALSE, type='l', col='tomato', lwd=1.2)
  label(0.97, 0.86, rounded(mean(s25), 2), adj=c(1, 0.5), col='tomato', cex=0.8)

  dat         <-  metRates[metRates$Species == species, ]
  dat$colors  <-  c('dodgerblue2', 'tomato')[match(dat$TempK, c(283.15, 298.15))]
  plot(NA, xlab='', ylab='', xlim=c(-1, 7), ylim=c(-6, 6), axes=FALSE, type='n')
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()  
  points(dat$lnRate ~ dat$lnMass, pch=16, col=make.transparent(dat$colors, 0.6))
  label(1.4, 0.5, y3, adj=c(0.5, 1), xpd=NA, srt=270)
  label(0.5, 1.4, x3, adj=c(0.5, 1), xpd=NA)
  axis(4, at=seq(-5, 5, 5), las=1)
  axis(3, at=seq(0, 6, 2))
  box()
  xpts  <-  range(dat$lnMass[dat$TempK == 283.15])
  lines(xpts, mean(i10) + mean(s10)*xpts, lwd=2.3)
  lines(xpts, mean(i10) + mean(s10)*xpts, col='dodgerblue2', lwd=2)
  xpts  <-  range(dat$lnMass[dat$TempK == 298.15])
  lines(xpts, mean(i25) + mean(s25)*xpts, lwd=2.3)
  lines(xpts, mean(i25) + mean(s25)*xpts, col='tomato', lwd=2)
  label(0.03, 0.9, spNames, adj=c(0, 0.5), xpd=NA, cex=1)

  den3d10  <-  kde2d(s10, i10, n=500)
  den3d25  <-  kde2d(s25, i25, n=500)
  plot(NA, xlim=c(-1, 2), ylim=c(-10, 5), las=1, ylab='', xlab='', axes=FALSE)
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  points(s10, i10, pch=16, col=make.transparent('dodgerblue2', 0.1))
  points(s25, i25, pch=16, col=make.transparent('tomato', 0.1))
  image(den3d10, col=make.transparent(col10, c(0, rep(0.1, 31))), add=TRUE)
  image(den3d25, col=make.transparent(col25, c(0, rep(0.1, 31))), add=TRUE)
  box(bty='l', lty=2)
  box(bty='7')

  xDens10    <-  getProb(density(i10))
  xDens25    <-  getProb(density(i25))
  plot(NA, xlab='', ylab='', xlim=c(0,0.025), ylim=c(-10, 5), xpd=NA, type='l', axes=FALSE)
  label(0.5, -0.4, x2, adj=c(0.5, 0.5), xpd=NA)
  label(1.4, 0.5, y2, adj=c(0.5, 1), xpd=NA, srt=270)
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  axis(4, at=seq(-10, 5, 5), las=1)
  axis(1, seq(0, 0.02, 0.01), labels=seq(0, 20, 10))
  box()
  polygon(c(xDens10$y, min(xDens10$y)), c(xDens10$x, min(xDens10$x)), border='dodgerblue2', col=make.transparent('dodgerblue2', 0.3))
  polygon(c(xDens25$y, min(xDens25$y)), c(xDens25$x, min(xDens25$x)), border='tomato', col=make.transparent('tomato', 0.3))
  quants  <-  quantile(i10, probs=c(0.025, 0.975), type=2)
  label(rep(0.95, 2), figureProp(quants, 'y'), text=FALSE, type='l', col='dodgerblue2', lwd=1.2)
  label(0.95, 0.15, rounded(mean(i10), 2), adj=c(1, 0.5), col='dodgerblue2', cex=0.8)
  quants  <-  quantile(i25, probs=c(0.025, 0.975), type=2)
  label(rep(0.86, 2), figureProp(quants, 'y'), text=FALSE, type='l', col='tomato', lwd=1.2)
  label(0.86, 0.06, rounded(mean(i25), 2), adj=c(1, 0.5), col='tomato', cex=0.8)
}

fig2  <-  function() {
  tmat  <-  tfit$BUGSoutput$sims.matrix
  tmat  <-  tmat[,paste0('beta[', seq_along(fixef(modelLmer1)), ']')]
  colnames(tmat)  <-  names(fixef(modelLmer1))
  species  <-  unique(metRates$Species)
  # species real names
  spNames  <-  list(substitute(italic('Hippopodina')*' sp.'),
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
  plotTriad(species[1], spNames[[1]], tmat, x1=TRUE, y1=TRUE, x2=FALSE, y2=FALSE, x3=TRUE, y3=FALSE)
  plot(0, 0, axes=FALSE, type='n', xlab='', ylab='')
  plotTriad(species[2], spNames[[2]], tmat, x1=TRUE, y1=FALSE, x2=FALSE, y2=TRUE, x3=TRUE, y3=TRUE)
  plot(0, 0, axes=FALSE, type='n', xlab='', ylab='')
  plotTriad(species[3], spNames[[3]], tmat, x1=FALSE, y1=TRUE, x2=TRUE, y2=FALSE, x3=FALSE, y3=FALSE)
  plot(0, 0, axes=FALSE, type='n', xlab='', ylab='')
  plotTriad(species[4], spNames[[4]], tmat, x1=FALSE, y1=FALSE, x2=TRUE, y2=TRUE, x3=FALSE, y3=TRUE)
}
