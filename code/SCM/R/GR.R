
grfit <- function(sample, data=NULL, what=c("X", "exp", "L",
                                       "Xobs", "Xobs.sc"),
                  scale=TRUE, dropAbund=TRUE) {

  what <- match.arg(what)
  
  if (is.null(sample$X) && is.null(data)) {
    stop("No `X' in `sample' and `data' is not given")
  }

  if (!is.null(sample$X) && !is.null(data)) {
    stop("`X' in `sample' and `data' is also given")
  }

  if (!is.null(data)) {
    sample$X <- data$X
    sample$X[data$I==0] <- sample$X.imp
    sample$I <- data$I
  }
  
  GRfile <- system.file("table_s1_gene_parameters_08-06-07.csv",
                        package="SCM")
  GRtab <- read.csv(GRfile)

  GRparamsfile <- system.file("GRparams.RData", package="SCM")
  load(GRparamsfile)

  genes <- intersect(GRtab$ORF, rownames(sample$X))
  calibration_genes <- intersect(lsCalibration, genes)

  GRtab <- GRtab[ match(genes, GRtab$ORF), c("ORF", "Growth.Rate.Slope")]
  gidx <- match(genes, rownames(sample$X))
  gidx2 <- match(calibration_genes, rownames(sample$X))

  if (what == "X") {
    full.dat <- sample$X
    exprData <- full.dat[gidx, ]
    if (dropAbund) { exprData <- full.dat[, sample$lj == "expr"] }
    if (scale) { exprData <- myscale(exprData) }
    grStats <- calculateGrowthStats(exprData, frmeGRParameters,
                                    calibration_genes, NULL)
    grStats$vdRates
  } else if (what == "exp") {
    prep <- sample$L %*% t(sample$G.mat) +
      sample$T %*% t(sample$B.mat) + sample$E
    rownames(prep) <- rownames(sample$X)
    exprPrep <- prep[gidx,]
    if (dropAbund) {
      exprPrep <- exprPrep[, unique(sample$kj[sample$lj == "expr"])]
    }
    if (scale) { exprPrep <- scale(exprPrep) }
    grStats <- calculateGrowthStats(exprPrep, frmeGRParameters,
                                    calibration_genes, NULL)
    grStats$vdRates
  } else if (what == "L") {
    L <- sample$L
    rownames(L) <- rownames(sample$X)
    if (dropAbund) { L <- L[,"expr",drop=FALSE] }
    if (scale) L <- scale(L)
    grStats <- calculateGrowthStats(L, frmeGRParameters,
                                    calibration_genes, NULL)
    grStats$vdRates
  } else if (what == "Xobs") {
    X <- sample$X
    X[sample$I==0] <- NA
    X <- X[gidx,]
    if (dropAbund) { X <- X[, sample$lj == "expr"] }
    if (scale) X <- scale(X)
    grStats <- calculateGrowthStats(X, frmeGRParameters,
                                    calibration_genes, NULL)
    grStats$vdRates
  } else if (what == "Xobs.sc") {
    X <- sample$X
    if (scale) X <- scale(X)
    X[sample$I==0] <- NA
    if (dropAbund) { X <- X[, sample$lj == "expr"] }
    X <- X[gidx,]
    grStats <- calculateGrowthStats(X, frmeGRParameters,
                                    calibration_genes, NULL)
    grStats$vdRates
  }
}

grcor <- function(sample, data=NULL, err=c("noise/signal", "relerr")) {

  err <-  match.arg(err)
  
  if (is.null(sample$X) && is.null(data)) {
    stop("No `X' in `sample' and `data' is not given")
  }

  if (!is.null(sample$X) && !is.null(data)) {
    stop("`X' in `sample' and `data' is also given")
  }

  if (!is.null(data)) {
    sample$X <- data$X
    sample$X[data$I==0] <- sample$X.imp
  }

  prepErr <- function(samp) {
    err <- (diag(samp$Tau)[samp$tj] +
            diag(samp$Xi)[samp$kj])/samp$G[samp$kj]^2
    err <- unique(err)
    names(err) <- unique(samp$kj)
    err
  }
  
  if (err == "noise/signal") {
    relerr <- getRelErr(sample)
    perr <- prepErr(sample)
  } else if (err == "relerr") {
    relerr <- getRelErr(sample)
    relerr <- relerr / (relerr+1)
    perr <- prepErr(sample)
    perr <- perr / (perr+1)
  }
  
  res <- list()  
  
  res$G.cor.rep <- cor.test(grfit(sample),
         sample$G[sample$kj[(sample$lj == "expr")]])
  res$G.cor.exp <- cor.test(grfit(sample, what="exp"),
         sample$G[unique(sample$kj[sample$lj == "expr"])])
  res$TE.cor.rep <- cor.test(grfit(sample),
         relerr[sample$Rnames[sample$lj == "expr"]])
  res$TE.cor.exp <- cor.test(grfit(sample, what="exp"),
         perr[unique(sample$kj[sample$lj == "expr"])])
  
  res
}
