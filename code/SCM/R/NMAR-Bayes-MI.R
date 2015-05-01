
# This is from MCMCpack and GPL-3

riwish <- function (v, S)  {
  solve(rwish(v, solve(S)))
}

# This is from MCMCpack and GPL-3

rwish <- function (v, S) {
  if (!is.matrix(S)) 
    S <- matrix(S)
  if (nrow(S) != ncol(S)) {
    stop(message = "S not square in rwish().\n")
  }
  if (v < nrow(S)) {
    stop(message = "v is less than the dimension of S in rwish().\n")
  }
  p <- nrow(S)
  CC <- chol(S)
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
  if (p > 1) {
    pseq <- 1:(p - 1)
    Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <-
      rnorm(p * (p - 1)/2)
  }
  crossprod(Z %*% CC)
}

startRandom <- function(data) {
  start <- list()
  start$data <- data
  start$data[ is.na(data) ] <- rnorm(sum(is.na(data)),
                                     mean(as.vector(data), na.rm=TRUE),
                                     sd=sd(as.vector(data), na.rm=TRUE))
  start$obsParam <- list(phi=NA, rho=NA) # recalculated anyway
  start$dataMean <- colMeans(start$data)
  start$dataCov <- cov(start$data)
  start
}

startEM <- function(data) {
  start <- list()
  start$data <- amelia(data, m=1)$imputations[[1]]
  start$obsParam <- list(phi=NA, rho=NA) # recalculated anyway
  start$dataMean <- colMeans(start$data)
  start$dataCov <- cov(start$data)
  start
}

updatePrecalcMI <- function(cur) {
  cur$covInv <- list()
  for (j in 1:ncol(cur$data)) {
    cur$covInv[[j]] <- cur$dataCov[j,-j,drop=FALSE] %*%
      solve(cur$dataCov[-j,-j,drop=FALSE])
  }
  cur
}

## Draw the phi and rho observation model parameters
## We use a logit now
drawObsParamMI <- function(data, current, preps) {
  lev <- levels(preps)
  current$obsParam <- list(phi=rep(NA, length(lev)),
                           rho=rep(NA, length(lev)))
  for (l in seq_along(lev)) {
    y <- as.vector(is.na(data[, lev[l] == preps]))
    x <- as.vector(current$data[, lev[l] == preps])
    lfit <- speedglm(y ~ x, data=data.frame(x=x, y=y),
                     family=binomial(link="logit"))
    lfitCoef <- as.matrix(coef(summary(lfit)))
    ## fix broken speedglm
    mode(lfitCoef) <- "numeric"
    current$obsParam$phi[l] <- rnorm(1, lfitCoef[1,1], lfitCoef[1,2])
    current$obsParam$rho[l] <- rnorm(1, lfitCoef[2,1], lfitCoef[2,2])
  }
  current
}

## Fill in the missing values, based on the observation model and
## the mean and covariance of the data. This is a
## Metropolis-Hastings step.
drawMissingMI <- function(data, current, preps, missIdx, propWidth) {
  imputed <- .Call("R_SCM_drawMissing", current, data,
                   as.integer(preps), length(levels(preps)), missIdx,
                   as.numeric(propWidth),
                   PACKAGE="SCM")
  current$data <- imputed$sample
  list(sample=current, param=imputed$param)
}

## Draw the mean of the data, conjugate (=Normal) prior,
## with large variance
## Draw the covariance matrix of the data, Jeffrys prior
drawMeanVar <- function(data, current) {
  n <- nrow(data)
  popMean <- colMeans(current$data)
  popCov <- cov(current$data)
  current$dataMean <- rmvnorm(1, mean=popMean, sigma=popCov/n)
  current$dataCov <- riwish(n-1, (n-1)*popCov)
  current
}

indepProp <- function(current, data, idx) {
  .Call("R_SCM_indepProp", current, data, as.integer(idx), PACKAGE="SCM")
}

## Gelman-Rubin convergence criteria
## sampledPars must be a matrix, each row corresponding to
## an MCMC chain, and each column corresponding to a sample
## The 'R' convergence statistics is returned

checkConv.GR <- function(sampledPars) {
  n <- ncol(sampledPars)
  m <- nrow(sampledPars)
  psi.j <- rowMeans(sampledPars)
  psi.. <- mean(psi.j)
  B <- n/(m-1) * sum( (psi.j - psi..)^2 )
  s.j.sq <- sapply(seq_len(m), function(j)
                   sum( (sampledPars[j,] - psi.j[[j]])^2 ) / (n-1))
  W <- mean(s.j.sq)
  var.psi <- (n-1)/n * W + B/n
  R <- sqrt(var.psi / W)
  R
}

checkConvObs <- function(samples) {
  idx <- as.integer(length(samples)/2+1):length(samples)
  nopreps <- length(samples[idx][[1]][[1]]$obsParam$phi)
  phi.R <- sapply(1:nopreps, function(x) {
    pp <- sapply(samples[idx], sapply, function(y) y$obsParam$phi[x])
    checkConv.GR(pp)
  })
  rho.R <- sapply(1:nopreps, function(x) {
    pp <- sapply(samples[idx], sapply, function(y) y$obsParam$rho[x])
    checkConv.GR(pp)
  })

  print(paste(phi.R, rho.R))
  if ( all(phi.R <= 1.1) && all(rho.R <= 1.1)) {
    samples <- unlist(samples[idx], recursive=FALSE)
    attr(samples, "break") <- TRUE
  } else {
    attr(samples, "break") <- FALSE
    attr(samples, "checkAt") <- attr(samples, "checkAt") * 2
  }

  samples
}

## data: must be a matrix, with missing values as NAs
## start: a function that takes data and gives back a
##        named list with slots:
##          'data': data with the missing values filled in
##          'obsParam': observation model parameters
##          'dataMean': mean values for each columns in 'data'
##          'dataCov': the covariance matrix of 'data'
##        Alternatively, if it is not a function, then it must be
##        the list of starting values, one for each chain.
## noChains: the number of MCMC chains to run. This must be at least two
##        if we want to access the convergence of the chains.
## stopMode: how to stop the simulation:
##          'fixed': take a fixed number of steps
##          'convObs': check the convergence of the observation
##             parameters. At least two chains are needed for this.
## propWidth: a multiplication factor for the normal Metropolis proposal.
## noSteps: number of steps to take, if 'stopMode' is 'fixed', if
##          'stopMode' is 'convObs', then it gives the minimum number
##          of steps to make, half of which will be ignored. If there
##          is no convergence, then the number of steps is doubled,
##          again dropping half of the sequence.

NMARsampler <- function(data, preps=rep(1, ncol(data)), start, noChains=2,
                        stopMode=c('convObs', 'fixed'), propWidth=1.0,
                        noSteps, verbose=FALSE) {

  ## Check arguments (a bit)
  if (!is.matrix(data) || !is.numeric(data)) {
    stop("`data' must be a numeric matrix")
  }
  if (sum(is.na(data)) == 0) {
    stop("No missing values in `data'")
  }
  noChains <- as.integer(noChains)
  if (noChains <= 1) {
    stop("Illegal `noChains' value")
  }
  stopMode <- match.arg(stopMode)
  if (stopMode == 'convObs' && noChains == 1) {
    stop("More than 1 chain is needed to check convergence")
  }
  if (!is.function(start) && length(start) != noChains) {
    stop("Length of 'start' must match 'noChains'")
  }
  
  ## Calculate starting values
  if (is.function(start)) {
    start <- lapply(1:noChains, function(x) start(data))
  }

  ## Which rows have missing values
  missIdx <- which(is.na(rowSums(data)))
  
  ## Main loop
  preps <- as.factor(preps)
  samples <- list()
  attr(samples, "checkAt") <- noSteps

  param <- list(noAccepted=0, noRejected=0)
  current <- lapply(start, updatePrecalcMI)
  stepsTaken <- 0L
  
  while (TRUE) {
    stepsTaken <- stepsTaken + 1L
    ## Draw the phi & rho parameters of the missingness logit model
    current <- lapply(current, drawObsParamMI, data=data, preps=preps)
    ## Draw the missing data
    current <- lapply(current, drawMissingMI, data=data, preps=preps,
                      missIdx=missIdx, propWidth=propWidth)
    newParam <- lapply(current, "[[", "param")
    for (p in newParam) {
      param$noAccepted <- param$noAccepted + p$noAccepted
      param$noRejected <- param$noRejected + p$noRejected
    }
    current <- lapply(current, "[[", "sample")
    ## Draw the mean and variance 
    current <- lapply(current, drawMeanVar, data=data)
    ## Update some temporary variables
    current <- lapply(current, updatePrecalcMI)
    ## Store sample
    samples[[ length(samples)+1 ]] <- current
    if (verbose) { cat(".") }
    ## Check stopping criterion
    if (stopMode == "fixed" && stepsTaken == noSteps) {
      break;
    } else if (stopMode == "convObs") {
      ## All this is from Bayesian Data Analysis, 2nd ed., page 296
      if (stepsTaken == attr(samples, "checkAt")) {
        samples <- checkConvObs(samples)
        if (attr(samples, "break")) { break; } 
      }
    }
  }

  if (verbose) { cat("\n") }
  list(samples=samples, param=param)
}

