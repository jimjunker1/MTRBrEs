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

linearRescale <- function(x, r.out) {
  p <- (x - min(x)) / (max(x) - min(x))
  r.out[[1]] + p * (r.out[[2]] - r.out[[1]])
}

rounded  <-  function(value, precision=1, change=FALSE) {
  if(change) {
    value  <-  value * -1
  }
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}

###########################
# SPECIES AS RANFOM EFFECTS
###########################
posteriorsBoltzmann  <-  function() {
  par(mfrow=c(1, 2), mar=c(5.1,4.1,4.1,3), omi=c(0.5,1,0.5,1), cex=1)
  theCol  <-  c('Hippopodina sp.'='dodgerblue2', 'Bugula neritina'='darkolivegreen3', 'Sponge'='tomato', 'Bugula stolinifera'='darkgoldenrod')

  # mass
  plot(NA, xlim=c(-1.5,2.5), ylim=c(0,1.5), xlab='Mass scaling', ylab='Posterior density', las=1)
  for(i in 1:4) {
    xDens  <-  density(rowSums(bfit$BUGSoutput$sims.matrix[,c('A', paste0('s[', i, ',1]'))]))
    polygon(c(xDens$x, min(xDens$x)), c(xDens$y, min(xDens$y)), border=theCol[i], col=make.transparent(theCol[i], 0.3))
    label(c(0.02, 0.07), rep(1-i*0.1, 2), text=FALSE, type='l', lwd=2, col=theCol[i])
      label(0.08, 1-i*0.1, round(xDens$x[which.max(xDens$y)[1]], 2), adj=c(0, 0.5), cex=0.8)
  }

  # temperature
  par(mar=c(5.1,1.1,4.1,6))
  plot(NA, xlim=c(0,1.1), ylim=c(0,9), xlab='Temperature dependence', ylab='', las=1)
  for(i in 1:4) {
    xDens  <-  density(rowSums(bfit$BUGSoutput$sims.matrix[,c('Er', paste0('s[', i, ',3]'))]))
    polygon(c(xDens$x, min(xDens$x)), c(xDens$y, min(xDens$y)), border=theCol[i], col=make.transparent(theCol[i], 0.3))
    label(c(0.02, 0.07), rep(1-i*0.1, 2), text=FALSE, type='l', lwd=2, col=theCol[i])
      label(0.08, 1-i*0.1, round(xDens$x[which.max(xDens$y)[1]], 2), adj=c(0, 0.5), cex=0.8)
  }

  # species color code
  for(k in seq_along(theCol)) {
    label(1.1, 1-k*0.1, text=FALSE, pch=21, col=theCol[k], bg=make.transparent(theCol[k], .5), adj=c(0.5, 0.5), xpd=NA)
    label(1.15, 1-k*0.1, names(theCol)[k], adj=c(0, 0.5), xpd=NA, font=3)
  }
}

modelComparison  <-  function() {
  tmat  <-  tfit$BUGSoutput$sims.matrix
  tmat  <-  tmat[,paste0('beta[', seq_along(fixef(modelLmer)), ']')]
  colnames(tmat)  <-  names(fixef(modelLmer))
  species  <-  unique(metRates$Species)
  # species real names
  spNames  <-  c('Hippopodina sp.'='Bryo', 'Bugula neritina'='Bugula', 'Microcionidae'='Sponge', 'Bugula stolonifera'='Stolonifera')

  par(mfcol=c(2,4), omi=c(0.5, 0.5, 0, 0), mai=rep(0.3, 4))
  for(j in seq_along(species)) {
    # all species together - 3-way interaction
    isReferenceSpecies  <-  species[j] == 'Bryo'
    if(isReferenceSpecies) {
      ylb  <-  'Rescaled [0,1] posterior density'
      t10  <-  tmat[, 'lnMass']
      t25  <-  rowSums(tmat[, c('lnMass', 'lnMass:Temp25')])
    } else {
      ylb  <-  ''
      t10  <-  rowSums(tmat[, c('lnMass', paste0('Species', species[j], ':', 'lnMass'))])
      t25  <-  rowSums(tmat[, c('lnMass', 'lnMass:Temp25', paste0('Species', species[j], ':', 'lnMass:Temp25'))])
    }
    
    xDens10    <-  density(t10)
    xDens10$y  <-  xDens10$y / max(xDens10$y)
    xDens25    <-  density(t25)
    xDens25$y  <-  xDens25$y / max(xDens25$y)

    # species separate
    tmatb  <-  perSpeciesJagsFits[[species[j]]]$simsMatrix
    t10b   <-  tmatb[, 'lnMass']
    t25b   <-  rowSums(tmatb[, c('lnMass', 'lnMass:Temp25')])
    
    xDens10b    <-  density(t10b)
    xDens10b$y  <-  xDens10b$y / max(xDens10b$y)
    xDens25b    <-  density(t25b)
    xDens25b$y  <-  xDens25b$y / max(xDens25b$y)

    # fix xlim
    xli        <-  range(c(xDens10$x, xDens10b$x, xDens25$x, xDens25b$x))
    xli        <-  c(xli[1] - 0.05*(xli[2]-xli[1]), xli[2] + 0.05*(xli[2]-xli[1]))
    plot(NA, xlab='', ylab=ylb, xlim=xli, ylim=c(0,1.2), las=1, xpd=NA, type='l')
    polygon(c(xDens10$x, min(xDens10$x)), c(xDens10$y, min(xDens10$y)), border='dodgerblue2', col=make.transparent('dodgerblue2', 0.3))
    polygon(c(xDens25$x, min(xDens25$x)), c(xDens25$y, min(xDens25$y)), border='tomato', col=make.transparent('tomato', 0.3))
    label(c(0.03, 0.1), rep(0.93, 2), text=FALSE, col='tomato', lwd=1.2, adj=c(0, 0.5), type='l')
    label(0.12, 0.93, round(xDens25$x[which.max(xDens25$y)[1]], 2), adj=c(0, 0.5), cex=0.7)
    label(c(0.03, 0.1), rep(0.83, 2), text=FALSE, col='dodgerblue2', lwd=1.2, adj=c(0, 0.5), type='l')
    label(0.12, 0.83, round(xDens10$x[which.max(xDens10$y)[1]], 2), adj=c(0, 0.5), cex=0.7)
    label(0.95, 0.93, names(spNames)[grep(species[j], spNames)], font=3, adj=c(1, 0.5), cex=0.9)

    plot(NA, xlab='Mass scaling', ylab=ylb, xlim=xli, ylim=c(0,1.2), las=1, xpd=NA, type='l')
    polygon(c(xDens10b$x, min(xDens10b$x)), c(xDens10b$y, min(xDens10b$y)), border='dodgerblue2', col=make.transparent('dodgerblue2', 0.3))
    polygon(c(xDens25b$x, min(xDens25b$x)), c(xDens25b$y, min(xDens25b$y)), border='tomato', col=make.transparent('tomato', 0.3))
    label(c(0.03, 0.1), rep(0.93, 2), text=FALSE, col='tomato', lwd=1.2, adj=c(0, 0.5), type='l')
    label(0.12, 0.93, round(xDens25b$x[which.max(xDens25b$y)[1]], 2), adj=c(0, 0.5), cex=0.7)
    label(c(0.03, 0.1), rep(0.83, 2), text=FALSE, col='dodgerblue2', lwd=1.2, adj=c(0, 0.5), type='l')
    label(0.12, 0.83, round(xDens10b$x[which.max(xDens10b$y)[1]], 2), adj=c(0, 0.5), cex=0.7)
  }
}
