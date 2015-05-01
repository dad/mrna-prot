
histMiss <- function(state, layout=NULL, type="b",
                     key=TRUE, names=TRUE, pdffile=NULL, truth=NULL,
                     width=20, height=16, logit=TRUE, ...) {
  
  invlogit <- function(x, eta0, eta1) {
    1/(1+exp(-(eta0+eta1*x)))
  }

  if (!is.null(pdffile)) {
    pdf(pdffile, width, height, ...)
  }

  if (is.null(layout)) {
    layout <- matrix(1:56, nrow=8, byrow=TRUE)
  }
  layout(layout)

  doplot <- function(I, imputed, eta0, eta1, name, truth, teta0, teta1) {
    h0 <- hist(c(imputed,range(truth)), 30, plot=FALSE)
    h1 <- hist(imputed, breaks=h0$breaks, plot=FALSE)
    h2 <- hist(imputed[I==0], breaks=h0$breaks, plot=FALSE)
    h3 <- hist(imputed[I!=0], breaks=h0$breaks, plot=FALSE)
    par(mar=c(3,3,1,1)+.1)
    plot(h1$mids[h1$count!=0], h1$count[h1$count!=0], type=type,
         ylim=c(0, max(h1$count)), pch=1,
         xlab="data", ylab="count")
    lines(h2$mids[h2$count!=0], h2$count[h2$count!=0],
          type=type, col=2, pch=2)
    lines(h3$mids[h3$count!=0], h3$count[h3$count!=0],
          type=type, col=3, pch=2)
    if (!is.null(truth)) {
      h4 <- hist(truth, breaks=h0$breaks, plot=FALSE)      
      lines(h4$mids[h4$count!=0], h4$count[h4$count!=0],
            type=type, col=5, pch=3)
    }

    if (logit) {
      cens <- sapply(h1$mids, invlogit, eta0=eta0, eta1=eta1)
      lines(h1$mids, cens*max(h1$count)/max(cens), col=4)
      at=seq(0, 1, length=11)
      axis(4, at=at*max(h1$count)/max(cens), label=at, col=4)
      if (!is.null(truth)) {
        tcens <- sapply(h1$mids, invlogit, eta0=teta0, eta1=teta1)
        lines(h1$mids, tcens*max(h1$count)/max(tcens), col=5)
      }
    }

    if (names) { title(main=name) }
    
    if (key) {
      if (logit) {
        legend("topright", bty="n", pch=c(1,2,2,NA), col=1:4, lty=1,
               c("total", "missing", "non-missing", "censoring"))
      } else {
        legend("topright", bty="n", pch=c(1,2,2), col=1:3, lty=1,
               c("total", "missing", "non-missing"))
      }
    }
  }  

  for (i in seq_len(ncol(state$X))) {
    doplot(state$I[,i], state$X[,i], state$eta0[state$kj[i]],
           state$eta1[state$kj[i]], name=colnames(state$X)[i],
           truth=truth$X[,i], truth$eta0[truth$kj[i]],
           truth$eta1[truth$kj[i]])
  }

  if (!is.null(pdffile)) {
    dev.off()
  }

  invisible(NULL)
}
  
getHist <- function(samples, par=c("cor", "G", "Xi", "Theta","Tau", "RelErr","eta0", "eta1", "nu","Lvar","Lmean","Tvar","Tmean"), truth=NULL) {

  par <- match.arg(par)
  print(par)
  if (par=="cor") {
    ev <- function(psi) {
        psi[upper.tri(psi, diag=FALSE)]
      }
    top <- sapply(samples, function(x) ev(x$psi))
    tru <- NULL ; if (!is.null(truth)) { tru <- ev(truth$psi) }
    
  } else if (par=="psi") {
    idx <- upper.tri(samples[[1]]$psi, diag=TRUE)
    top <- sapply(samples, function(x) x$psi[idx])
    tru <- if (is.null(truth)) NULL else truth$psi[idx]
  } 
    else if (par=="G") {
    top <- sapply(samples, function(x) x$G)
    tru <- truth$G
  } else if (par=="Xi") {
    top <- sapply(samples, function(x) diag(x$Xi))
    tru <- diag(truth$Xi)
  } else if (par=="Theta") {
    top <- sapply(samples, function(x) diag(x$Theta))
    tru <- diag(truth$Theta)
  } else if(par=="RelErr"){
    top <- sapply(samples, function(x) getRelErr(x))
    tru <- getRelErr(truth)
  } else if (par=="eta0") {
    top <- sapply(samples, function(x) x$eta0)
    tru <- truth$eta0
  } else if (par=="eta1") {
    top <- sapply(samples, function(x) x$eta1)
    tru <- truth$eta1
  } else if (par=="nu") {
    top <- sapply(samples, function(x) x$nu)
    tru <- truth$nu
  } else if (par=="Tau"){
    top <- sapply(samples, function(x) diag(x$Tau))
    tru <- diag(truth$Tau)
  } else if(par=="Lvar"){
    top <- sapply(samples, function(x) apply(x$L,2,var))
    tru <- rep(1,length(unique(truth$lj)))
  } else if(par=="Tvar"){
    top <- sapply(samples, function(x) apply(x$T,2,var))
    tru <- diag(truth$Tau)
  } else if(par=="Lmean"){
    top <- sapply(samples, function(x) apply(x$L,2,mean))
    tru <- rep(0,length(unique(truth$lj)))
  } else if(par=="Tmean"){
    top <- sapply(samples, function(x) apply(x$T,2,mean))
    tru <- rep(0,length(unique(truth$tj)))
  } 
  
  if (is.list(top)) { top <- unlist(top) } # workaround 
  structure(t(rbind(top)), truth=tru)
}

addRightTri <- function(pos, col=1:6, border=NA, xpd=TRUE, cex=1, ...) {
  col <- rep(col, length=length(pos))
  usr <- par("usr")
  cex <- 1 * cex
  xfac <- (usr[2]-usr[1])/20
  yfac <- (usr[4]-usr[3])/20
  if (all(pos == pos[1])) {
    xc <- usr[2]
    xo <- cex*sqrt(3)/2*xfac
    polygon(c(xc, xc+xo, xc+xo),
            c(pos[1], pos[1]+cex/2*yfac, pos[1]-cex/2*yfac),
            xpd=xpd, col=NA, border=NULL, ...)
  } else {
    for (i in seq_along(pos)) {
      xc <- usr[2]
      xo <- cex*sqrt(3)/2*xfac
      if (i > 1) { xc <- xc + xo/2 *
                     sum(abs(pos[1:(i-1)] - pos[i]) < cex*yfac/4) }
      polygon(c(xc, xc+xo, xc+xo),
              c(pos[i], pos[i]+cex/2*yfac, pos[i]-cex/2*yfac),
              xpd=xpd, col=col[i], border=border, ...)
    }
  }
}

tracePlot <- function(samples, par, truth=NULL, truthCex=1,
                      pdffile=NULL, mar=c(5,4,2,3)+.1, ...) {

  pars <- eval(formals(getHist)$par)  
  par <- match.arg(par, pars)  
  
  if (!is.null(pdffile)) {
    pdf(pdffile, ...)
  }

  top <- getHist(samples, par, truth)

  ylim <- if (is.null(attr(top, "truth"))) {
    r <- range(top)
    r + c(-.05,.05) * abs(r)
  } else {
    range(c(top, attr(top, "truth")))
  }
  par(mar=mar)
  if (length(top) != 0) {
    matplot(top, type="l", xlab="iterations", ylab=par, ylim=ylim)
    if (!is.null(attr(top, "truth"))) {
      addRightTri(attr(top, "truth"), col=1:6, cex=truthCex)
    }
  }
  
  if (!is.null(pdffile)) {
    dev.off()
  }

  invisible(NULL)
}

## Does it make sense to plot 'par' for 'sample'?

shouldPlot <- function(sample, par) {
  if (par %in% c("cor", "G", "Xi", "Theta", "Tau", "RelErr","nu",
                 "Lvar", "Lmean", "Tvar", "Tmean")) {
    TRUE
  } else if (par %in% c("eta0", "eta1")) {
    is.null(sample$I) || any(sample$I == 0)
  }
}

allTracePlots <- function(samples, truth=NULL, pdffile=NULL, layout=NULL,
                          pdfargs=NULL, width=20, height=16,
                          pars=NULL, ...) {

  if (is.null(pars)) {
    pars <- eval(formals(getHist)$par)
    pars <- pars[sapply(pars, shouldPlot, sample=samples[[1]])]
  }
  
  if (!is.null(pdffile)) {
    do.call(pdf, c(list(pdffile, width=width, height=height), pdfargs))
  }
  if (!is.null(layout)) {
    layout(layout)
  }

  sapply(pars, function(p)
         tracePlot(samples=samples, par=p, truth=truth, ...))
  
  if (!is.null(pdffile)) {
    dev.off()
  }
  invisible(NULL)
}

acfPlot <- function(samples, par, chain=1, lag.max=NULL, pdffile=NULL, ...) {

  pars <- eval(formals(getHist)$par)  
  par <- match.arg(par, pars)  
  
  if (!is.null(pdffile)) {
    pdf(pdffile, ...)
  }

  top <- getHist(samples, chain, par)
  n <- ncol(top)

  if (!is.null(lag.max)) { lag.max <- rep(lag.max, length=n) }  
  
  nr <- ceiling(sqrt(n))
  nc <- ceiling(n/nr)
  layout(matrix(1:(nr*nc), nr=nr, nc=nc, byrow=TRUE))

  for (i in 1:n) {
    acf(top[,i], lag.max=lag.max[i], main=paste(par, i))
  }
  
  if (!is.null(pdffile)) {
    dev.off()
  }

  invisible(NULL)
}


histRes <- function(sample, preps, cens1, cens2, layout=NULL, type="b",
                    key=TRUE, names=TRUE, pdffile=NULL,
                     width=20, height=16, ...) {
  
  if (!is.null(pdffile)) {
    pdf(pdffile, width, height, ...)
  }

  if (is.null(layout)) {
    layout <- matrix(1:56, nrow=8, byrow=TRUE)
  }
  layout(layout)

  for (i in seq_len(ncol(sample$ResM1))) {
    cols <- getPrep(colnames(cens1))==colnames(sample$ResM1)[i]
    obsGenes <- !apply(cens1[,cols],1,function(x) any(is.na(x)))

    hist(sample$ResM1[,i],xlim=c(-6,6),main=colnames(sample$ResM1)[i])
    hist(sample$ResM1[obsGenes,i],xlim=c(-6,6),main=colnames(sample$ResM1)[i],col="red",add=TRUE)
  }
  for (i in seq_len(ncol(cens2))) {
    obsGenes <- !is.na(cens2[,i])
    hist(sample$ResM2[,i],xlim=c(-6,6),main=colnames(sample$ResM2)[i])
    hist(sample$ResM2[obsGenes,i],xlim=c(-6,6),col="red",add=TRUE)
  }
  

  if (!is.null(pdffile)) {
    dev.off()
  }

  invisible(NULL)
}

grImage <- function(sample, extrapolate=.1, resolution=200, CC=NULL,
                    markFit=TRUE, markZero=TRUE, verbose=TRUE, ...) {
  
  require(Matrix)
  attach(sample)
  
  nz <- function(x) x[x!=0]
  extrapol <- function(r, frac) {
    if (r[1] >= r[2]) { stop("Invalid range") }
    if (frac < 0) { stop("Invalid fraction") }
    ex <- (r[2]-r[1]) * frac
    r + c(-ex, ex)
  }

  allga <- nz(c(G1["S.abund",], G2["S.abund",]))
  allge <- nz(c(G1["S.expr",], G2["S.expr",]))
  arange <- extrapol(range(allga), extrapolate)
  erange <- extrapol(range(allge), extrapolate)

  aseq <- seq(arange[1], arange[2], length=resolution)
  eseq <- seq(erange[1], erange[2], length=resolution)

  if (is.null(CC)) {
    CC <- matrix(0, length(aseq), length(eseq))
    for (i in seq_along(aseq)) {
      if (verbose) { cat(".") }
      for (j in seq_along(eseq)) {
        ga <- aseq[i] ; ge <- eseq[j]
        CC[i,j] <- CC[j,i] <- cor(L[,1] + S[,1] * ga, L[,2] + S[,2] * ge)
      }
    }
    if (verbose) { cat("\n") }
    rownames(CC) <- aseq
    colnames(CC) <- eseq
  }

  detach(sample)

  rescale <- function(x, xmin, xmax, omin=min(x), omax=max(x)) {
    (x-omin) / (omax-omin) * (xmax-xmin) + xmin
  }
  
  alabel <- pretty(aseq)
  elabel <- pretty(eseq)
  aat <- rescale(alabel, omin=min(aseq), omax=max(aseq),xmin=0.5, xmax=resolution+0.5)
  eat <- rescale(elabel, omin=min(eseq), omax=max(eseq),xmin=0.5, xmax=resolution+0.5)

  anms <- sapply(names(allga),function(x) paste("\n",strsplit(x,".",fixed=TRUE)[[1]][2]))
  enms <- sapply(names(allge),function(x) strsplit(x,".",fixed=TRUE)[[1]][2])

  alabs <- mapply(function(x,y) paste(x,y,sep="\n"),anms,as.character(1:length(anms)))
 elabs <-  mapply(function(x,y) paste(x,y,sep="\n"),enms,as.character(1:length(enms)))

  eat <- c(eat,rescale(allge, omin=min(eseq), omax=max(eseq),
                     xmin=0.5, xmax=resolution+0.5))
  aat <- c(aat,rescale(allga, omin=min(aseq), omax=max(aseq),
                     xmin=0.5, xmax=resolution+0.5))
  alabs <- c(alabel,anms)
  elabs <- c(elabel,enms)
  
  myplot <- image(Matrix(CC), xlab="abundance", ylab="expression", lwd=0,
                  col.regions=heat.colors(30), colorkey=TRUE,
                  scales=list(x=list(at=aat, labels=alabs, rot=45),
                    y=list(at=eat, labels=elabs), rot=45), sub="", ...)

  if (markFit || markZero) {
    opan <- myplot$panel
    myplot$panel <- function(...) {
      opan(...)
      if (markFit) {
        panel.abline(v=rescale(allga, omin=min(aseq), omax=max(aseq),
                       xmin=0.5, xmax=resolution+0.5))
        panel.abline(h=rescale(allge, omin=min(eseq), omax=max(eseq),
                       xmin=0.5, xmax=resolution+0.5))
      }
      if (markZero) {
        panel.abline(v=rescale(0, omin=min(aseq), omax=max(aseq),
                       xmin=0.5, xmax=resolution+0.5))
        panel.abline(h=rescale(0, omin=min(eseq), omax=max(eseq),
                       xmin=0.5, xmax=resolution+0.5))
      }

##       panel.text(x=rescale(allga,omin=min(aseq),omax=max(aseq),xmin=0.5,xmax=resolution+0.5),
##       y=0, labels=names(allga),xpd=FALSE)
##       panel.text(x=0,y=rescale(allge,omin=min(eseq),omax=max(eseq),xmin=0.5,xmax=resolution+0.5),
##                        labels=names(allga))
      
    }
  }
  
  print(myplot)
  
  invisible(CC)
}

postDist <- function(fits, par, truth=NULL, main="") {
  pars <- eval(formals(getHist)$par)  
  par <- match.arg(par, pars)  
  toplot <- getHist(fits, par=par, truth=truth)

  for (i in 1:ncol(toplot)) {
    toplot.i <- toplot[,i]
    tru <- attr(toplot, "truth")[i]
    if (!is.null(tru) && is.na(tru)) tru <- NULL
    h <- hist(toplot[,i], plot=FALSE)
    xlim <- range(c(h$breaks, tru))
    plot(h, xlim=xlim, main=main, xlab=par)
    if (!is.null(tru)) { text(tru, 0, "*", col="red", cex=3) }
  }
  invisible(NULL)
}

allPostDist <- function(fits, pars=NULL, pdffile=NULL, layout=NULL,
                        pdfargs=NULL, width=20, height=15, ...) {

  if (is.null(pars)) {
    pars <- eval(formals(getHist)$par)
    pars <- pars[sapply(pars, shouldPlot, sample=fits[[1]][[1]])]
  }

  if (!is.null(pdffile)) {
    do.call(pdf, c(list(pdffile, width=width, height=height), pdfargs))
  }
  if (!is.null(layout)) {
    layout(layout)
  }

  sapply(pars, function(p) postDist(fits=fits, par=p, ...))
  
  if (!is.null(pdffile)) {
    dev.off()
  }
  invisible(NULL)
}
  
plotAll <- function(dir, include=NULL, drop=1:50,
                    truth=sprintf("%s/truth.Rdata", dir), truthname="data",
                    patternPrefix="thinned", patternSuffix=".RData",
                    sampleName="nstart",
                    pdffile=sprintf("%s/plots.pdf", dir),
                    layout=matrix(1:24, nrow=4, byrow=TRUE),
                    pdfwidth=20, pdfheight=15,
                    tracePars=NULL,
                    skipTracePars=character(),
                    postPars=NULL,
                    skipPostPars=character(),
                    pdfargs=NULL) {

  ## load all data first
  if (!is.null(include)) {
    files <- sprintf("%s-%i.Rdata", patternPrefix, include)
  } else {
    files <- list.files(dir, pattern=sprintf("%s-.*\\%s",patternPrefix,patternSuffix))
    files <- setdiff(files, sprintf("%s-%i%s", patternPrefix,patternSuffix,drop))
  }
  files <- files[order(as.numeric(gsub("[^0-9]","", files)))]
  fits <- lapply(files, function(f) {
    print(f) ; load(sprintf("%s/%s", dir, f)) ;
    get(sampleName, inherits=FALSE)
  })

  ## open output file, set layout
  if (!is.null(pdffile)) {
    do.call(pdf, c(list(pdffile, width=pdfwidth, height=pdfheight), pdfargs))
  }
  if (!is.null(layout)) { layout(layout) }

  ## Load truth, if given
  if (!is.null(truth)) { load(truth); }
  
  ## trace plots first
  allTracePlots(samples=fits, pars=tracePars,
                truth=if (!is.null(truth)) get(truthname) else NULL)

  ## posterior distributions next
  #allPostDist(fits=fits, pars=postPars,
  #            truth=if (!is.null(truth)) get(truthname) else NULL)
  
  ## close file
  if (!is.null(pdffile)) { dev.off() }
}
