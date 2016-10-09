######################
# AUXILLIARY FUNCTIONS
######################

extrafont::loadfonts(quiet=TRUE)

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

circle  <-  function(x, y, radius, start=0, end=pi, nsteps=100,...){  
    pinRatio  <-  par('pin')[2]/par('pin')[1]
    rs        <-  seq(start, end, len=nsteps) 
    xc        <-  x + radius * cos(rs)
    yc        <-  (y + radius * sin(rs)) / pinRatio
    lines(xc,yc,...) 
}

proportionalPng  <-  function(logo, px, py, log = FALSE, ...) {
    if(!is.numeric(px) | !is.numeric(py) | length(px) != 2) {
        stop('wrong position coordinates [0,1]')
    }
    usr  <-  par('usr')
    pin  <-  par('pin')
    # first get proportions of coordinates right
    pxRg    <-  px[2] - px[1]
    xProp   <-  usr[1] + px * (usr[2] - usr[1]) # x range from relative px 
    yProp   <-  usr[3] + py * (usr[4] - usr[3]) # minimum y from relative py
    # now get aspect ratio to calculate maximum y
    pinRatio  <-  pin[2] / pin[1] # aspect ratio of actual plot region, depends on device and plot size
    dims      <-  dim(logo)[1:2] # number of x-y pixels for the logo (aspect ratio)
    AR        <-  dims[1] / dims[2]
    yProp     <-  c(yProp, usr[3] + (py + pxRg * AR / pinRatio) * (usr[4] - usr[3])) # maximum y from relative py correcting for x and plot ratios
    if (log == 'x') xProp <- 10^(xProp)
    if (log == 'y') yProp <- 10^(yProp)
    if (log == 'xy') {xProp <- 10^(xProp); yProp <- 10^(yProp)}
    rasterImage(logo, xProp[1], yProp[1], xProp[2], yProp[2], interpolate=TRUE, ...)
}

changePngColour  <-  function(pngObject, col, ...) {
    rgbVals  <-  col2rgb(col, ...) / 255
    for(i in 1:3)
        pngObject[,,i]  <-  rgbVals[i]
    pngObject
}

whiteGrid  <-  function (...) {
    proportionalLabel(rep(0.2, 2), c(0, 1), text = FALSE, type = 'l', 
        col = 'grey70', lwd = 0.3, ...)
    proportionalLabel(rep(0.4, 2), c(0, 1), text = FALSE, type = 'l', 
        col = 'grey70', lwd = 0.3, ...)
    proportionalLabel(rep(0.6, 2), c(0, 1), text = FALSE, type = 'l', 
        col = 'grey70', lwd = 0.3, ...)
    proportionalLabel(rep(0.8, 2), c(0, 1), text = FALSE, type = 'l', 
        col = 'grey70', lwd = 0.3, ...)
    proportionalLabel(c(0, 1), rep(0.2, 2), text = FALSE, type = 'l', 
        col = 'grey70', lwd = 0.3, ...)
    proportionalLabel(c(0, 1), rep(0.4, 2), text = FALSE, type = 'l', 
        col = 'grey70', lwd = 0.3, ...)
    proportionalLabel(c(0, 1), rep(0.6, 2), text = FALSE, type = 'l', 
        col = 'grey70', lwd = 0.3, ...)
    proportionalLabel(c(0, 1), rep(0.8, 2), text = FALSE, type = 'l', 
        col = 'grey70', lwd = 0.3, ...)
}

makeFigure1  <-  function() {
    toPdf(fig1(), 'output/figures/fig1.pdf', width=8, height=5)
    embed_fonts('output/figures/fig1.pdf')
}

makeFigure2  <-  function(...) {
    toPdf(fig2(...), 'output/figures/fig2.pdf', width=8, height=7)
    embed_fonts('output/figures/fig2.pdf')
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

plotTriad  <-  function(species=species, spNames=spNames, tmat, dat, x1=TRUE, y1=TRUE, x2=TRUE, y2=TRUE, x3=TRUE, y3=TRUE, x4=TRUE, y4=TRUE) {
    col10     <-  colorRampPalette(rev(brewer.pal(9, 'Greys')))
    col10     <-  col10(32)
    col10[1]  <-  '#FFFFFF'
    col25     <-  colorRampPalette(rev(brewer.pal(9, 'Greys')))
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
    rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
    whiteGrid()
    if(x1) axis(3, at=seq(-1, 2, 1))
    axis(2, at=seq(0, 0.008, 0.004), labels=seq(0, 8, 4), las=1)
    box()
    polygon(c(xDens10$x, min(xDens10$x)), c(xDens10$y, min(xDens10$y)), border='grey50', col=transparentColor('grey50', 0.3))
    polygon(c(xDens25$x, min(xDens25$x)), c(xDens25$y, min(xDens25$y)), border='black', col=transparentColor('black', 0.6))
    quants  <-  quantile(s10, probs=c(0.025, 0.975), type=2)
    proportionalLabel(figureProp(quants), rep(0.95, 2), text=FALSE, type='l', col='grey50', lwd=1.2)
    proportionalLabel(0.97, 0.95, rounded(mean(s10), 2), adj=c(1, 0.5), col='grey50', cex=0.8)
    quants  <-  quantile(s25, probs=c(0.025, 0.975), type=2)
    proportionalLabel(figureProp(quants), rep(0.86, 2), text=FALSE, type='l', col='black', lwd=1.2)
    proportionalLabel(0.97, 0.86, rounded(mean(s25), 2), adj=c(1, 0.5), col='black', cex=0.8)
    
    dat$colors  <-  c('grey50', 'black')[match(dat$TempK, c(283.15, 298.15))]
    plot(NA, xlab='', ylab='', xlim=c(-1, 7), ylim=c(-6, 6), axes=FALSE, type='n')
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
    whiteGrid()  
    points(dat$lnRate ~ dat$lnMass, pch=16, col=transparentColor(dat$colors, 0.6))
    proportionalLabel(0.5, 1.4, x2l, adj=c(0.5, 1), xpd=NA)
    proportionalLabel(1.4, 0.5, y2l, adj=c(0.5, 1), xpd=NA, srt=270)
    if(x2) axis(3, at=seq(0, 6, 2))
    axis(4, at=seq(-5, 5, 5), las=1)
    box()
    xpts  <-  range(dat$lnMass[dat$TempK == 283.15])
    lines(xpts, mean(i10) + mean(s10)*xpts, lwd=2.3)
    lines(xpts, mean(i10) + mean(s10)*xpts, col='grey50', lwd=2)
    xpts  <-  range(dat$lnMass[dat$TempK == 298.15])
    lines(xpts, mean(i25) + mean(s25)*xpts, lwd=2.3)
    lines(xpts, mean(i25) + mean(s25)*xpts, col='black', lwd=2)
    proportionalLabel(0.03, 0.9, spNames, adj=c(0, 0.5), xpd=NA, cex=1)
    
    den3d10  <-  kde2d(s10, i10, n=500)
    den3d25  <-  kde2d(s25, i25, n=500)
    plot(NA, xlim=c(-1, 2), ylim=c(-10, 7), las=1, ylab='', xlab='', axes=FALSE)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
    whiteGrid()
    proportionalLabel(0.5, -0.4, x3l, adj=c(0.5, 0.5), xpd=NA)
    proportionalLabel(-0.4, 0.5, y3l, adj=c(0.5, 0.5), xpd=NA, srt=90)
    points(s10, i10, pch=16, col=transparentColor('grey50', 0.1))
    points(s25, i25, pch=16, col=transparentColor('black', 0.1))
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
    rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
    whiteGrid()
    if(x4) axis(1, seq(0, 0.02, 0.01), labels=seq(0, 20, 10))
    if(y4) axis(4, at=seq(-10, 5, 5), las=1)
    box()
    polygon(c(xDens10$y, min(xDens10$y)), c(xDens10$x, min(xDens10$x)), border='grey50', col=transparentColor('grey50', 0.3))
    polygon(c(xDens25$y, min(xDens25$y)), c(xDens25$x, min(xDens25$x)), border='black', col=transparentColor('black', 0.6))
    quants  <-  quantile(i10, probs=c(0.025, 0.975), type=2)
    proportionalLabel(rep(0.95, 2), figureProp(quants, 'y'), text=FALSE, type='l', col='grey50', lwd=1.2)
    proportionalLabel(0.95, 0.15, rounded(mean(i10), 2), adj=c(1, 0.5), col='grey50', cex=0.8)
    quants  <-  quantile(i25, probs=c(0.025, 0.975), type=2)
    proportionalLabel(rep(0.86, 2), figureProp(quants, 'y'), text=FALSE, type='l', col='black', lwd=1.2)
    proportionalLabel(0.86, 0.06, rounded(mean(i25), 2), adj=c(1, 0.5), col='black', cex=0.8)
}

fig1  <-  function() {
    par(mai=rep(0,4))
    plot(0, 0, type='n')
    pinRatio  <-  par('pin')[2]/par('pin')[1]
    lines(c(-0.8, 0.8), c(0,0), lty=2)
    polygon(c(-0.8, 0.8, 0.8, -0.8, -0.8), c(-0.9, -0.9, 0.9, 0.9, -0.9), lty=2)
    
    # Erect species
    bugulaWhole  <-  readPNG('data/bugula_whole.png')
    bugulaCut1   <-  readPNG('data/bugula_cut1.png')
    bugulaCut2   <-  readPNG('data/bugula_cut2.png')
    proportionalPng(bugulaWhole, c(0.14,0.39), c(0.51))
    proportionalPng(changePngColour(bugulaCut1, 'grey30'), c(0.44,0.6), c(0.51))
    proportionalPng(changePngColour(bugulaCut2, 'grey60'), c(0.64,0.85), c(0.51))
    lines(c(-0.5115759, -0.4814338), c(0.2616007, 0.1159565), lwd=2.5, col='grey30', lend=2)
    lines(c(-0.4292721, -0.4514349), c(0.3390424, 0.1743852), lwd=2.5, col='grey60', lend=2)
    
    # Flat species
    ## full circle
    circle(-0.35, -0.3, 0.2, end=pi*2, lwd=1.5)
    text(-0.35, -0.3/pinRatio, 'A = 2D Area', adj=c(0.5, 0.5), cex=1.3, font=3)
    text(-0.35, -0.1, 'G = Growing edge', adj=c(0.5, 0.5), cex=1.3, font=3)
    ## 1/4 circle
    circle(0.35, -0.3, 0.2, end=pi*2, lwd=1.5)
    circle(0.35, -0.3, 0.2, end=pi/2, col='grey60', lwd=1.5)
    lines(c(0.35,0.35),c(-0.3,-0.1)/pinRatio, lwd=1.5, lty=2, col='grey60')
    lines(c(0.35,0.55),c(-0.3,-0.3)/pinRatio, lwd=1.5, lty=2, col='grey60')
    text(0.365, -0.23/pinRatio, 'A / 4', adj=c(0, 0.5), col='grey60', cex=1.3, font=3)
    text(0.5, -0.14/pinRatio, 'G / 4', adj=c(0, 0.5), col='grey60', cex=1.3, font=3)
}

fig2  <-  function(metRates, output) {
    tmat  <-  output$model$BUGSoutput$sims.matrix
    tmat  <-  tmat[,paste0('beta[', seq_len(ncol(coef(output$lmerModel1)$Run)), ']')]
    colnames(tmat)  <-  names(coef(output$lmerModel1)$Run)
    species         <-  unique(metRates$Species)
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
    plotTriad(species[1], spNames[[1]], tmat, metRates[metRates$Species == species[1], ], x1=TRUE, y1=TRUE, x2=TRUE, y2=FALSE, x3=FALSE, y3=TRUE, x4=FALSE, y4=FALSE)
    plot(0, 0, axes=FALSE, type='n', xlab='', ylab='')
    plotTriad(species[2], spNames[[2]], tmat, metRates[metRates$Species == species[2], ], x1=TRUE, y1=FALSE, x2=TRUE, y2=TRUE, x3=FALSE, y3=FALSE, x4=FALSE, y4=TRUE)
    plot(0, 0, axes=FALSE, type='n', xlab='', ylab='')
    plotTriad(species[3], spNames[[3]], tmat, metRates[metRates$Species == species[3], ], x1=FALSE, y1=TRUE, x2=FALSE, y2=FALSE, x3=TRUE, y3=TRUE, x4=TRUE, y4=FALSE)
    plot(0, 0, axes=FALSE, type='n', xlab='', ylab='')
    plotTriad(species[4], spNames[[4]], tmat, metRates[metRates$Species == species[4], ], x1=FALSE, y1=FALSE, x2=FALSE, y2=TRUE, x3=TRUE, y3=FALSE, x4=TRUE, y4=TRUE)
}
