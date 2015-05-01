
loadExprAbund <- function(file="data/scer-raw-extended-Apr.20.2012.Rdata",
                          transformation=c("asinh", "log"),
                          preps=TRUE, ...) {

  transformation <- match.arg(transformation)
  load(file)
  dat <- extended.data$trans
  rm(extended.data)

  rownames(dat) <- dat[,1]
  dat <- dat[, -(1:3)] # Remove first three columns

  if (preps) {
    ## Create three measurements from ingolia
    colnames(dat)[colnames(dat) == "expr.ingolia.1"] <- "expr.ingolia1.1"
    colnames(dat)[colnames(dat) == "expr.ingolia.2"] <- "expr.ingolia1.2"
    colnames(dat)[colnames(dat) == "expr.ingolia.rq"] <- "expr.ingolia2.1"
    colnames(dat)[colnames(dat) == "expr.ingolia.rq.late"] <-
      "expr.ingolia2.2"
    colnames(dat)[colnames(dat) == "expr.ingolia.ca"] <- "expr.ingolia3.1"
    colnames(dat)[colnames(dat) == "expr.ingolia.ca.late"] <-
      "expr.ingolia3.2"
    
    ## Create two measurements for Lipson
    colnames(dat)[colnames(dat) == "expr.lipson.ma"] <- "expr.lipson2.ma"
    
    newNames <- matrix(nc=2, byrow=TRUE,
                       c("abund.gygi", "abund.gygi.1",
                         "abund.futcher", "abund.futcher.1",
                         "abund.washburn", "abund.washburn.1",
                         "abund.peng", "abund.peng.1",
                         "abund.ghaem", "abund.ghaem.1",
                         "abund.newman", "abund.newman.1",
                         "abund.lu", "abund.lu.1",
                         "abund.degodoy", "abund.degodoy.1",
                         "expr.velc", "expr.velc.1",
                         "expr.holstege", "expr.holstege.1",
                         "expr.causton.acid.a", "expr.causton.acida",
                         "expr.causton.acid.b", "expr.causton.acidb",
                         "expr.causton.peroxide.a", "expr.causton.peroxidea",
                         "expr.causton.peroxide.b", "expr.causton.peroxideb",
                         "expr.causton.nacl", "expr.causton.nacl",
                         "expr.nagalakshmi", "expr.nagalakshmi.1", 
                         "expr.yassour.ypd0.1", "expr.yassour.ypd01",
                         "expr.yassour.ypd0.2", "expr.yassour.ypd02",
                         "expr.yassour.ypd15.1", "expr.yassour.ypd151",
                         "expr.yassour.ypd15.2", "expr.yassour.ypd152",
                         "expr.yassour", "expr.yassour.1",
                         "expr.mackay", "expr.mackay.1", 
                         "expr.garcia", "expr.garcia.1",
                         "expr.pelechano", "expr.pelechano.1"))
    for (i in 1:nrow(newNames)) {
      colnames(dat)[colnames(dat)==newNames[i,1]] <- newNames[i,2]
    }

    limits <- c(abund.lee.1=2, abund.lee.2=2, abund.lee.3=2,
              expr.causton.acida=3.8, expr.causton.acidb=3.8,
              expr.causton.nacl=3.8, expr.causton.peroxidea=3.8,
              expr.causton.peroxideb=3.8, expr.nagalakshmi.1=1.0)
     
  } else {
    ## no preps, replace all but the first dots with _
    cn <- strsplit(colnames(dat), ".", fixed=TRUE)
    colnames(dat) <- sapply(cn, function(x) {
      paste(sep="", x[1], ".", paste(x[-1], collapse="_"), ".1")})

    limits <- c(abund.lee_1.1=2, abund.lee_2.1=2, abund.lee_3.1=2,
                expr.causton_acid_a.1=3.8, expr.causton_acid_b.1=3.8,
                expr.causton_nacl.1=3.8, expr.causton_peroxide_a.1=3.8,
                expr.causton_peroxide_b.1=3.8, expr.nagalakshmi.1=1.0)
  }
  
  ## get measurement type for each column
  col.types <- sapply(colnames(dat), function(cname)
                      strsplit(cname, split=".", fixed=T)[[1]][1])
  unique.types <- unique(col.types)
  
  ## pred.types <- c("int", "disp", "abund", "expr", "cai", "gc", "mw")
  pred.types <- c("abund", "expr")
  
  pred <- dat[, col.types %in% pred.types]
  predmat.types <- sapply(strsplit(colnames(pred), "\\."), "[", 1)

  for (l in seq_along(limits)) {
    idx <- which(colnames(pred)==names(limits)[l])
    pred[,idx][ pred[,idx] < limits[l] ] <- NA
  }
  
  ## transformation
  if (transformation=="asinh") {
    pred <- asinh(exp(pred))
    for (l in seq_along(limits)) {
      idx <- which(colnames(pred)==names(limits)[l])
      pred[,idx][ pred[,idx] < limits[l] ] <- NA
    }
  } else if (transformation=="log") {
    for (l in seq_along(limits)) {
      idx <- which(colnames(pred)==names(limits)[l])
      pred[,idx][ pred[,idx] < log(limits[l]) ] <- NA
    }
  }
  
  ## Remove the completely missing genes
  pred <- pred[!is.na(rowMeans(pred, na.rm=TRUE)),]

  ## Make it a matrix instead of a data frame
  pred <- as.matrix(pred)
  
  ## sort the columns
  preps <- sort(colnames(pred))
  chkpreps(preps)
  pred <- pred[,preps,drop=FALSE]

  ## separate it based on number of levels
  multiMeasNames <- multiMeasNames(preps)
  uniMeasNames <- uniMeasNames(preps)
  list(preps=preps, data1=pred[,multiMeasNames, drop=FALSE],
       data2=pred[,uniMeasNames, drop=FALSE])
}

## Load the data into three latent variables, only keep
## microarray and RNAseq measurements for mRNA.

loadData3 <- function(file="data/scer-raw-extended-Apr.20.2012.Rdata",
                      transformation=c("asinh", "log"), preps=TRUE, ...) {

  data <- loadExprAbund(file=file, transformation=transformation,
                        preps=preps, ...)

  type <- c(velc="SAGE", roth="MA", holstege="MA", causton="MA",
            dudley="cMA", miura="cPCR", nagalakshmi="RS",
            yassour="RS", lipson="RS", lipson2="MA", ingolia1="RS",
            ingolia2="RS", ingolia3="RS", mackay="cMA", garcia="cMA",
            pelechano="cMA")

  ma.pn <- paste("expr", sep=".", names(type)[type %in% c("MA", "cMA")])
  ma.rn <- paste("expr", sep=".", names(type)[type %in% c("RS")])
  ma.preps <- grepl(paste(ma.pn, collapse="|"), data$preps)
  rn.preps <- grepl(paste(ma.rn, collapse="|"), data$preps)
  ab.preps <- grepl("^abund", data$preps)

  alldata <- cbind(data$data1, data$data2)[, data$preps]

  ma.names <- data$preps[ma.preps]
  rn.names <- data$preps[rn.preps]

  oldpreps <- preps <- data$preps
  preps[ma.preps] <- sub("^expr", "ma", preps[ma.preps])
  preps[rn.preps] <- sub("^expr", "rnaseq", preps[rn.preps])
  colnames(alldata) <- preps

  preps <- preps[ ma.preps | rn.preps | ab.preps ]
  oldpreps <- oldpreps[ ma.preps | rn.preps | ab.preps ]

  oldpreps <- oldpreps[order(preps)]
  preps <- sort(preps)
  alldata <- alldata[, preps]
  
  list(preps=preps, oldpreps=oldpreps,
       data1=alldata[, multiMeasNames(preps), drop=FALSE],
       data2=alldata[, uniMeasNames(preps), drop=FALSE])
}

splitLatent <- function(state, newpreps) {

  ## Check newpreps

  ## All prep names must stay equal after the first dot, or NA
  valid <- !is.na(newpreps)
  stopifnot(all(sub("[^\\.]*\\.", "", state$preps[valid]) ==
                sub("[^\\.]*\\.", "", newpreps[valid])))

  ## TODO: other checks
  
  ## Mapping of the latent variables
  latentNames <- getLatent(state$preps)
  newLatentNames <- sub("([^\\.]*).*$", "\\1", newpreps)
  latentNames <- latentNames[!duplicated(newLatentNames) &
                             !is.na(newLatentNames)]
  newLatentNames <- newLatentNames[!duplicated(newLatentNames) &
                                   !is.na(newLatentNames)]
  latentNames <- latentNames[order(newLatentNames)]
  newLatentNames <- sort(newLatentNames)
  
  ## Mapping of the preps
  prepNames <- getPrep(state$preps)
  newPrepNames <- sub("\\.[^\\.]*$", "", newpreps)
  newPrepNames <- newPrepNames[!duplicated(prepNames)]
  prepNames <- prepNames[!duplicated(prepNames)]

  cens <- cbind(state$cens1, state$cens2)[, state$preps]
  data <- cbind(state$data1, state$data2)[, state$preps]
  rho <- c(state$rho1, state$rho2)[state$preps]
  phi <- c(state$phi1, state$phi2)[state$preps]
  nu <- c(state$nu1, state$nu2)[state$preps]  
  
  newcens <- cens[, !is.na(newpreps)]
  newdata <- data[, !is.na(newpreps)]
  colnames(newcens) <- colnames(newdata) <- na.omit(newpreps)

  newrho <- rho[!is.na(newpreps)]
  newphi <- phi[!is.na(newpreps)]
  newnu <- nu[!is.na(newpreps)]
  names(newrho) <- names(newphi) <- names(newnu) <- na.omit(newpreps)  
  
  fpreps <- sort(na.omit(newpreps))

  newstate <- state[c("numGenes")]
  
  newstate$preps <- fpreps

  newstate$cens1 <- newcens[, multiMeasNames(fpreps)]
  newstate$cens2 <- newcens[, uniMeasNames(fpreps)]
  newstate$data1 <- newdata[, multiMeasNames(fpreps)]
  newstate$data2 <- newdata[, uniMeasNames(fpreps)]

  newstate$rho1 <- newrho[multiMeasNames(fpreps)]
  newstate$rho2 <- newrho[uniMeasNames(fpreps)]
  newstate$phi1 <- newphi[multiMeasNames(fpreps)]
  newstate$phi2 <- newphi[uniMeasNames(fpreps)]

  ## We just duplicate the latent variable(s) that were split
  newstate$L <- state$L[, paste("L.", sep="", latentNames)]
  colnames(newstate$L) <- paste("L.", sep="", newLatentNames)
  newstate$L <- newstate$L[, order(colnames(newstate$L))]
  
  newstate$S <- state$S[, paste("S.", sep="", latentNames)]
  colnames(newstate$S) <- paste("S.", sep="", newLatentNames)
  newstate$S <- newstate$S[, order(colnames(newstate$S))]
  
  newstate$P <- state$P[, prepNames[match(multiPrepNames(fpreps),
                                          newPrepNames)]]
  colnames(newstate$P) <- multiPrepNames(fpreps)

  ## We need to duplicate some values here, plus jitter the 1.0 values
  ## a bit, otherwise the psi matrix will be singular
  rnz <- function(mat) {
    mat[mat != 0] <- mat[mat != 0] - runif(sum(mat != 0)) / 10
    diag(mat)[diag(mat)!=0] <- 1.0
    (mat + t(mat)) / 2
  }
  
  newstate$psi <- matrix(NA, nrow=2*noLatent(fpreps), ncol=2*noLatent(fpreps))
  rownames(newstate$psi) <- colnames(newstate$psi) <- psiNames(fpreps)
  newPsiIdx <- c(paste("L.", sep="", newLatentNames),
                 paste("S.", sep="", newLatentNames))
  oldPsiIdx <- c(paste("L.", sep="", latentNames),
                 paste("S.", sep="", latentNames))
  newstate$psi[newPsiIdx, newPsiIdx] <- state$psi[oldPsiIdx, oldPsiIdx]
  newstate$psi <- rnz(newstate$psi)
  
  updateScaleMatrix <- function(new, old) {
    nz <- function(x) x[x!=0]
    for (c in seq_len(ncol(new))) {
      nv <- nz(old[, prepNames[match(colnames(new)[c], newPrepNames)]])
      new[which(new[,c] != 0), c] <- nv
    }
    new
  }
  updateScaleMatrix2 <- function(new, old) {
    nz <- function(x) x[x!=0]
    for (c in seq_len(ncol(new))) {
      nv <- nz(old[, state$preps[match(colnames(new)[c], newpreps)]])
      new[which(new[,c] != 0), c] <- nv
    }
    new
  }
  
  newstate$B <- updateScaleMatrix(createUnitB(fpreps), state$B)
  newstate$G1 <- updateScaleMatrix(createUnitG1(fpreps), state$G1)
  newstate$G2 <- updateScaleMatrix2(createUnitG2(fpreps), state$G2)
  newstate$Lambda1 <- updateScaleMatrix2(createUnitLambda1(fpreps),
                                         state$Lambda1)
  newstate$Lambda2 <- updateScaleMatrix2(createUnitLambda2(fpreps),
                                         state$Lambda2)
  
  newstate$Eps1 <- diag(noMeasMulti(fpreps))
  colnames(newstate$Eps1) <- rownames(newstate$Eps1) <- multiMeasNames(fpreps)
  diag(newstate$Eps1) <-
    diag(state$Eps1)[ state$preps[match(colnames(newstate$Eps1), newpreps)] ]

  newstate$Eps2 <- diag(noMeasUni(fpreps))
  colnames(newstate$Eps2) <- rownames(newstate$Eps2) <- uniMeasNames(fpreps)
  diag(newstate$Eps2) <-
    diag(state$Eps2)[ state$preps[match(colnames(newstate$Eps2), newpreps)] ]

  newstate$Chi <- diag(noPrepsMulti(fpreps))
  colnames(newstate$Chi) <- rownames(newstate$Chi) <- multiPrepNames(fpreps)
  diag(newstate$Chi) <-
    diag(state$Chi)[ prepNames[match(colnames(newstate$Chi), newPrepNames)] ]

  newstate$nu1 <- newnu[multiMeasNames(fpreps)]
  newstate$nu2 <- newnu[uniMeasNames(fpreps)]

  newstate$B.mask <- getBmask(fpreps, scaling=FALSE)
  newstate$Lambda1.mask <- getLambda1mask(fpreps, scaling=FALSE)
  newstate$Lambda2.mask <- getLambda2mask(fpreps, scaling=FALSE)
  newstate$G1.mask <- getmask(createUnitG1(fpreps), scaling=TRUE,
                              fixfirst=FALSE, fixsecond=FALSE)
  newstate$G2.mask <- getmask(createUnitG2(fpreps), scaling=TRUE,
                              fixfirst=FALSE, fixsecond=FALSE)

  newstate
}

mergeIngolia <- function(state) {

  newstate <- state[c("numGenes")]
  
  newstate$preps <- state$preps
  newstate$preps[grep("ingolia[12]", newstate$preps)] <-
    sprintf("expr.ingolia.%s", c("1", "2", "rq", "rq_late"))
  ord <- order(newstate$preps)
  newstate$preps <- newstate$preps[ord]
  
  map <- structure(newstate$preps, names=state$preps)
  
  updateColNames <- function(mat) {
    colnames(mat) <- unname(map[colnames(mat)])
    mat
  }

  ord1 <- order(map[colnames(state$cens1)])
  newstate$cens1 <- updateColNames(state$cens1)[, ord1]
  newstate$cens2 <- state$cens2
  newstate$data1 <- updateColNames(state$data1)[, ord1]
  newstate$data2 <- state$data2

  updateNames <- function(vec) {
    names(vec) <- unname(map[names(vec)])
    vec
  }

  newstate$rho1 <- updateNames(state$rho1)[ord1]
  newstate$rho2 <- state$rho2
  newstate$phi1 <- updateNames(state$phi1)[ord1]
  newstate$phi2 <- state$phi2

  newstate$L <- state$L
  newstate$S <- state$S

  meanColMat <- function(mat) {
    ingP <- grepl("ingolia", colnames(mat))
    newmat <- mat[, !ingP]
    newmat <- cbind(newmat, "expr.ingolia"=rowMeans(mat[, ingP]))
    newmat[, order(colnames(newmat))]
  }

  newstate$P <- meanColMat(state$P)
  
  newstate$psi <- state$psi

  newstate$B <- meanColMat(state$B)
  newstate$G1 <- meanColMat(state$G1)
  
  newstate$G2 <- state$G2

  colMax <- function(x) apply(x, 2, max)
  
  maxRowMat <- function(mat) {
    ingP <- grepl("ingolia", rownames(mat))
    newmat <- mat[!ingP, ]
    newmat <- rbind(newmat, "expr.ingolia"=colMax(mat[ingP, ]))
    newmat[order(rownames(newmat)), ]
  }    
  
  newstate$Lambda1 <- updateColNames(maxRowMat(state$Lambda1))[, ord1]

  newstate$Lambda2 <- state$Lambda2
  
  newstate$Eps1 <- t(updateColNames(t(updateColNames(state$Eps1))))[ord1,
                                                                    ord1]
  
  newstate$Eps2 <- state$Eps2

  meanVec <- function(vec) {
    ingP <- grepl("ingolia", names(vec))
    newvec <- vec[!ingP]
    newvec <- c(newvec, "expr.ingolia"=mean(vec[ingP]))
    newvec[order(names(newvec))]
  }
  
  newstate$Chi <- diag(meanVec(diag(state$Chi)))
  rownames(newstate$Chi) <- colnames(newstate$Chi) <- colnames(newstate$B)

  newstate$nu1 <- updateNames(state$nu1)[ord1]
  newstate$nu2 <- state$nu2

  newstate$B.mask <- getBmask(newstate$preps, scaling=FALSE)
  newstate$Lambda1.mask <- getLambda1mask(newstate$preps, scaling=FALSE)
  newstate$Lambda2.mask <- getLambda2mask(newstate$preps, scaling=FALSE)
  newstate$G1.mask <- getmask(createUnitG1(newstate$preps), scaling=TRUE,
                              fixfirst=FALSE, fixsecond=FALSE)
  newstate$G2.mask <- getmask(createUnitG2(newstate$preps), scaling=TRUE,
                              fixfirst=FALSE, fixsecond=FALSE)

  newstate
}

#' Get thresholds for droppong a given number of genes
#'
#' @param data As loaded with \code{loadExprAbund}.
#' @param no_to_drop Single integer, number of points to drop
#' @return Named vector of the thresholds.
#'
#' @export

get_drop_thresholds <- function(data, no_to_drop) {

  stopifnot(is.matrix(data))
  stopifnot(length(no_to_drop) == 1 ||
              length(no_to_drop) == ncol(data),
            all(no_to_drop >= 0), all(no_to_drop <= nrow(data)))

  no_to_drop <- rep(no_to_drop, length.out = ncol(data))

  get_drop_threshold <- function(x, no) {
    if (no == 0) {
      -Inf
    } else {
      thr <- rev(sort(x, na.last = TRUE, decreasing = TRUE))[no]
      if (is.na(thr)) thr <- min(x, na.rm = TRUE)
      thr
    }
  }

  thr <- rep(NA_real_, ncol(data))
  for (i in seq_len(ncol(data))) {
    thr[i] <- get_drop_threshold(data[,i], no = no_to_drop[i])
  }
  names(thr) <-colnames(data)
  thr
}

threshold_obs <- function(data, thresholds) {

  stopifnot(is.matrix(data))
  stopifnot(all(colnames(data) == names(thresholds)))

  thr_col <- function(x, thr) ifelse(x < thr | is.na(x), 0, 1)

  thr_cols <- function(df, thrs) {
    for (i in seq_len(ncol(df))) {
      thr <- thrs[ colnames(df)[i] ]
      df[,i] <- thr_col(df[,i], thr)
    }
    df
  }

  thr_cols(data, thresholds)
}

get_start_state <- function(file = "data/nonpooled-tech.RData",
                            name = "samp") {
  e <- new.env()
  load(file, envir = e)
  e[[name]]
}

threshold_state <- function(state, no_to_drop) {
  data <- state$X
  data[state$I == 0] <- NA
  state$thresholds <- get_drop_thresholds(data, no_to_drop)
  state$I <- threshold_obs(data, state$thresholds)
  state
}
