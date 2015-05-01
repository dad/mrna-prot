
obsCor <- function(dataMat, method=c("pearson", "kendall", "spearman")) {
  nc <- ncol(dataMat)
  res <- matrix(0, nc, nc)
  for (i in 1:nc) {
    for (j in i:nc) {
      ind <- which(!is.na(dataMat[,i]) & !is.na(dataMat[,j]))
      res[i,j] <- res[j,i] <- cor(dataMat[ind,i], dataMat[ind,j],
                                  method=method)
    }
  }
  res
}

## Correlation of mean expression and mean abundance
meanCor <- function(dataMatRep, dataMatObs, cur){

  means <- cur$nu
  means <- means[sort(names(means))]

  dataVar <- diag(V.X(cur))
  dataVar <- dataVar[sort(names(dataVar))]

  scaledMat <- t((t(dataMatObs$X)-means)/dataVar)
  abundMeanObs <- apply(scaledMat[,grep("abund", colnames(scaledMat))], 1,
                        function(x) mean(x,na.rm=TRUE))
  exprMeanObs <-  apply(scaledMat[,grep("expr", colnames(scaledMat))], 1,
                        function(x) mean(x,na.rm=TRUE))
  ind1 <- which(!is.na(exprMeanObs) & !is.na(abundMeanObs))
  
  scaledMat <- t((t(dataMatRep)-means)/dataVar)
  abundMeanRep <- apply(scaledMat[,grep("abund", colnames(scaledMat))], 1,
                        function(x) mean(x,na.rm=TRUE))
  exprMeanRep <-  apply(scaledMat[,grep("expr", colnames(scaledMat))], 1,
                        function(x) mean(x,na.rm=TRUE))
  ind2 <- which(!is.na(exprMeanRep) & !is.na(abundMeanRep))
  
  c(cor(abundMeanObs[ind1],exprMeanObs[ind1]),
    cor(abundMeanRep[ind2],exprMeanRep[ind2]))

}

## Mean of all pairwise correlations
meanCor2 <- function(dataMat){

  dataMat <- dataMat[,sort(colnames(dataMat))]
  exprCens <- dataMat[,grep("expr",colnames(dataMat))]
  abundCens <- dataMat[,grep("abund",colnames(dataMat))]

  cors <- c()
  for(i in 1:ncol(exprCens)){
    for(j in 1:ncol(abundCens)){
      exprCol <- exprCens[,i]
      abundCol <- abundCens[,j]
      ind <- which(!is.na(exprCol) & !is.na(abundCol))
      cors <- c(cors,cor(exprCol[ind],abundCol[ind]))
    }
  }

  mean(cors)
}

ppc <- function(samplefiles, par=c("obsCov", "obsCor", "meanCor",
                               "meanCor2","var", "iqr", "obsMean",
                               "kendall", "spearman"),
                verbose=TRUE, samplename="samp") {

  par <- match.arg(par)

  getval <- function(data, dataObs, cur, par){
    if (par=="obsCov") {
      cov(data, use="pairwise.complete.obs")
    } else if (par=="obsCor") {
      cor(data, use="pairwise.complete.obs")
    } else if (par=="meanCor") {
      meanCor(data, dataObs, cur)
    } else if (par=="meanCor2"){
      meanCor2(data)
    } else if (par=="var") {
      apply(data, 2, var, na.rm=TRUE)
    } else if (par=="iqr") {
      apply(data, 2, IQR, na.rm=TRUE)
    } else if (par=="obsMean") {
      apply(data, 2, mean, na.rm=TRUE)
    } else if (par=="kendall") {
      obsCor(data, method="kendall")
    } else if (par=="spearman") {
      obsCor(data, method="spearman")
    }
  }

  ## Observed data
  myenv <- new.env()
  load(samplefiles[[1]], envir=myenv)
  data <- myenv[[samplename]]
  data$X[ data$I == 0 ] <- NA

  ## Generate artificial data sets
  repl <- lapply(samplefiles, function(sf) {
    if (verbose) { print(sf) }
    myenv <- new.env()
    load(sf, envir=myenv)
    cur <- myenv[[samplename]]

    getnoise <- function(covmat, len) {
      sapply(diag(covmat), function(v) rnorm(len, 0, sqrt(v)))
    }

    T <- getnoise(cur$Tau, cur$N)
    E <- getnoise(cur$Xi, cur$N)
    R <- getnoise(cur$Theta, cur$N)

    X <- (cur$L[,cur$lj] %*% diag(cur$G[cur$kj]) +
          T[,cur$tj] + E[,cur$kj] + R +
          rep(cur$nu, each=cur$N))
    colnames(X) <- cur$Rnames
    
    cens.prob <- 1/(1+exp(-rep(cur$eta0[cur$kj], each=cur$N) -
                          rep(cur$eta1[cur$kj], each=cur$N) * X))
    I <- 1 - (runif(length(X)) > cens.prob)
    X[ I == 0 ] <- NA
    
    getval(X, data, cur, par)
  })  

  pp <- if (par=="obsCov") {
    dataval <- getval(data$X, par=par)
    covReps <- t(sapply(repl, function(x) x[upper.tri(x, diag=TRUE)]))
    covObsVec <- dataval[upper.tri(dataval, diag=TRUE)]
    sapply(1:length(covObsVec), function(i) {
      sum(covReps[,i] > covObsVec[i]) / nrow(covReps)
    })
  } else if (par %in% c("obsCor", "kendall", "spearman")) {
    dataval <- getval(data$X, par=par)
    corReps <- t(sapply(repl, function(x) x[upper.tri(x, diag=FALSE)]))
    corObsVec <- dataval[upper.tri(dataval, diag=FALSE)]
    sapply(1:length(corObsVec), function(i) {
      sum(corReps[,i] > corObsVec[i]) / nrow(corReps)
    })
  } else if (par=="meanCor") {
    dataval <- NULL
    mean(sapply(repl, function(x) x[1]<x[2]))
  }  else if (par=="meanCor2") {
    dataval <- getval(data$X, par=par)
    sum(dataval > repl)/length(repl)
  } else if (par=="var") {
    dataval <- getval(data$X, par=par)
    varReps <- do.call(rbind, repl)
    sapply(1:length(dataval), function(i) {
      sum(varReps[,i] > dataval[i]) / nrow(varReps)
    })
  } else if (par=="iqr") {
    dataval <- getval(data$X, par=par)
    iqrReps <- do.call(rbind, repl)
    sapply(1:length(dataval), function(i) {
      sum(iqrReps[,i] > dataval[i]) / nrow(iqrReps)
    })
  } else if (par=="obsMean") {
    dataval <- getval(data$X, par=par)
    meanReps <- do.call(rbind, repl)
    sapply(1:length(dataval), function(i) {
      sum(meanReps[,i] > dataval[i]) / nrow(meanReps)
    })
  }

  list(repl=repl, dataval=dataval, pp=pp, par=par)
}

ppcPlot <- function(PP, data) {
  if (! PP$par %in% c("obsCov", "obsCor", "meanCor", "var", "iqr", "obsMean",
                      "kendall", "spearman")) {
    stop("Unknown `par'")
  }

  if (PP$par %in% c("var", "iqr", "obsMean")) {
    layout(matrix(1:63, nr=7, nc=9, byrow=TRUE))
    for (i in 1:length(nstart$preps)) {
      top <- do.call(cbind, PP$repl)[i,]
      xlim <- range(c(top, PP$dataval[i]))
      hist(top, xlim=xlim, main=data$preps[i], xlab=PP$par)
      text(PP$dataval[i], 0, "*", cex=4, col=2, xpd=NA)
    }
    
  } else if (PP$par %in% c("obsCov", "obsCor", "kendall", "spearman")) {
    ## TODO
    stop("This `par' is not implemented yet")
  } else if (PP$par %in% c("meanCor")) {
    top <- do.call(rbind, PP$repl)
    lim <- range(top)
    plot(top, xlim=lim, ylim=lim, xlab="Data mean cor", ylab="Model mean cor")
    lines(lim, lim)
    text(lim[1], lim[1], adj=c(0,0), label=sprintf("p=%g", round(PP$pp, 3)),
         cex=2)
  }
  invisible(NULL)
}
