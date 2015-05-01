
"~" <- function(...) UseMethod("~")
"~.default" <- .Primitive("~")
"~.character" <- function(...) {
  if (any(! sapply(list(...), mode) == "character")) {
    stop("Only characters can be concatenated with '~'")
  }
  paste(..., sep="")
}

chr <- as.character

defaultLNames <- function(NL) { "L" ~ chr(1:NL) }
defaultTNames <- function(NT) { "T" ~ chr(1:NT) }
defaultENames <- function(NE) { "E" ~ chr(1:NE) }

defaultRNames <- function(NL, NT, NE, NR) {
  (chr(defaultLJ(NL, NT, NE, NR)) ~ "."  ~
   chr(defaultTJ(NL, NT, NE, NR)) ~ "."  ~
   chr(defaultKJ(NL, NT, NE, NR)) ~ ".R" ~ chr(1:NR))
}

defaultLJ <- function(NL, NT, NE, NR) { rep(defaultLNames(NL), each=NR / NL) }
defaultTJ <- function(NL, NT, NE, NR) { rep(defaultTNames(NT), each=NR / NT) }
defaultKJ <- function(NL, NT, NE, NR) { rep(defaultENames(NE), each=NR / NE) }

defaultNu <- function(names) {
  structure(rep(0, length(names)), names=names)
}

defaultEta0 <- function(cEta0, NL, NT, NE, NR) {
  structure(rep(cEta0, NE),
            names=as.character(unique(defaultKJ(NL, NT, NE, NR))))
}

defaultEta1 <- function(cEta1, NL, NT, NE, NR) {
  structure(rep(cEta1, NE),
            names=as.character(unique(defaultKJ(NL, NT, NE, NR))))
}

defaultG <- function(NL, NT, NE, NR) {
  structure(rep(1, NE), names=as.character(unique(defaultKJ(NL, NT, NE, NR))))
}

defaultPsi <- function(NL) {
  psi <- matrix(0.9, NL, NL)
  diag(psi) <- 1
  colnames(psi) <- rownames(psi) <- defaultLNames(NL)
  psi
}

setNames <- function(state){
  newstate <- with(state, {
    if(is.null(names(G)))  names(G) <- unique(kj)
    if(is.null(dimnames(Theta)))
      rownames(Theta) <- colnames(Theta) <- Rnames
    if(is.null(dimnames(Xi)))
      rownames(Xi) <- colnames(Xi) <- unique(kj)
    if(is.null(dimnames(Tau))) 
      rownames(Tau) <- colnames(Tau) <- unique(tj)
    if(is.null(dimnames(TechWeights)))
      rownames(TechWeights) <- colnames(TechWeights) <- unique(tj)
    if(is.null(dimnames(psi)))
      rownames(psi) <- colnames(psi) <- unique(lj)
    if(is.null(names(nu)))  names(nu) <- Rnames
    if (is.null(names(eta0))) names(eta0) <- unique(kj)
    if (is.null(names(eta1))) names(eta1) <- unique(kj)
    state$G <- G
    state$Theta <- Theta
    state$Xi <- Xi
    state$Tau <- Tau
    state$TechWeights <- TechWeights
    state$psi <- psi
    state$nu <- nu
    state$eta0 <- eta0
    state$eta1 <- eta1
    state
  })
  newstate
}

ndiag <- function(size, names, diag=1) {
  structure(diag(rep(diag, size)), dimnames=list(names, names))
}

scaleCols <- function(M, scale) {
  t(t(M) * scale)
}

getX <- function(state) {
  with(state, {
    X <- scaleCols(L[, lj], G[kj]) + T[, tj] + E[, kj] + R + rep(nu, each=N)
    colnames(X) <- colnames(R)
    X
  })
} 

posrnorm <- function(...) {
  x <- rnorm(...)
  ifelse(x > 0, x, 0)
}

homoTheta <- normTheta <- function(mean, sd=0) {
  function(lj, tj, kj, Rnames) {
    diag(posrnorm(length(Rnames), mean, sd))
  }
}

homoXi <- normXi <- function(mean, sd=0) {
  function(lj, tj, kj, Rnames) {
    diag(posrnorm(length(unique(kj)), mean, sd))
  }
}

homoTau <- normTau <- function(mean, sd=0) {
  function(lj, tj, kj, Rnames) {
    diag(posrnorm(length(unique(tj)), mean, sd))
  }
}
homoNu <- normNu <- function(mean, sd=0) {
  function(lj, tj, kj, Rnames) {
    rnorm(length(unique(Rnames)), mean, sd)
  }
}

homoEta <- normEta <- function(mean, sd=0) {
  function(lj, tj, kj, Rnames) {
    posrnorm(length(unique(kj)), mean, sd)
  }
}

E2R <- function(state,enames){
  rnames <- with(state,{
    Rnames[which(kj%in%enames)]
  })
  rnames
}

T2E <- function(state,tnames){
  enames <- with(state,{
    unique(kj[which(tj%in%tnames)])
  })
  enames
}
  
T2R <- function(state,tnames){
  E2R(state,T2E(state,tnames))
}

removeReps <- function(state,rNames){
    newstate <- with(state, {
        newstate <- state

        toRemove <- match(rNames,Rnames)
        kj <- kj[-toRemove]
        tj <- tj[-toRemove]
        lj <- lj[-toRemove]
        Rnames <- Rnames[-toRemove]

        newstate$kj <- kj
        newstate$tj <- tj
        newstate$lj <- lj
        newstate$Rnames <- Rnames
        
        newstate$NR <- length(Rnames)
        newstate$NE <- length(unique(kj))
        newstate$NT <- length(unique(tj))
        newstate$NL <- length(unique(lj))
        newstate$nu <- nu[Rnames]
        newstate$Theta <- Theta[Rnames,Rnames]
        newstate$X <- X[,Rnames]
        newstate$I <- I[,Rnames]
        newstate$R <- R[,Rnames]

        newstate$G <- G[unique(kj)]
        newstate$eta0 <- eta0[unique(kj)]
        newstate$eta1 <- eta1[unique(kj)]
        newstate$Xi <- Xi[unique(kj),unique(kj)]
        newstate$E <- E[,unique(kj)]
        
        newstate$TechWeights <- TechWeights[unique(tj),unique(tj)]
        newstate$Tau <- Tau[unique(tj),unique(tj)]
        newstate$T <- T[,unique(tj)]
        
        newstate$psi <- psi[unique(lj),unique(lj)]
        newstate$L <- L[,unique(lj)]

        newstate$G.mat <- G.mat(newstate)
        newstate$B.mat <- B.mat(newstate)
        newstate$Lambda.mat <- Lambda.mat(newstate)
        
        newstate
    })

    newstate
}

getRelErr <- function(state){
  err <- (diag(state$Tau)[state$tj]+diag(state$Xi)[state$kj]+diag(state$Theta))/state$G[state$kj]^2
  names(err) <- state$Rnames
  err
}

genDataSpec <- function(N=1000, NL=2, tList=c(3,2), eList=c(2,2,2,1,4),
                        rList=rep(3,sum(eList)),
                        Lnames=defaultLNames(NL),
                        Tnames=defaultTNames(NT),
                        Enames=defaultENames(NE),
                        psi=defaultPsi(NL),
                        Theta=homoTheta(.1), Xi=homoXi(.3),
                        TechWeights=homoTau(1), cTau=homoTau(0),
                        nu=homoNu(0), eta0=homoEta(Inf), eta1=homoEta(0)) {

  NT <- sum(tList)
  NE <- sum(eList)
  NR <- sum(rList)
  kj <- rep(Enames, rList)
  tj <- rep(rep(Tnames, eList), rList)
  lj <- rep(rep(rep(Lnames, tList), eList), rList)
  Rnames <- lj ~ "." ~ tj ~ "." ~ kj ~ ".R" ~ chr(seq_along(kj))

  if (is.function(Theta)) { Theta <- Theta(lj, tj, kj, Rnames) }
  if (is.function(Xi))    { Xi    <- Xi   (lj, tj, kj, Rnames) }
  if (is.function(eta0))  { eta0  <- eta0 (lj, tj, kj, Rnames) }
  if (is.function(eta1))  { eta1  <- eta1 (lj, tj, kj, Rnames) }

  G <- structure(rep(1, NE), names=Enames)
  Xi=ndiag(NE, Enames, diag=.3)
  Theta=ndiag(NR, Rnames, diag=.1)

  ## Don't draw theta if not identifiable
  Theta.draw <- unname(!(table(kj)==1)[kj])
  names(Theta.draw) <- Rnames
  diag(Theta)[!Theta.draw] <- 0

  ## Don't draw xi if not identifiable
  E2T <- tj[match(unique(kj),kj)]
  Xi.draw <- unname(!(table(E2T)==1)[E2T])
  names(Xi.draw) <- unique(kj)
  diag(Xi)[!Xi.draw] <- 0

  data <- genData(N=N, NL=NL, NT=NT, NE=NE, NR=NR, Rnames=Rnames,
                  lj=lj, tj=tj, kj=kj, psi=psi, nu=nu, G=G, TechWeights=TechWeights,
                  eta0=eta0, eta1=eta1, Xi=Xi, Theta=Theta, cTau=cTau)


  data
}
  


genData <- function(N=1000,  NL=2, NT=NL*2, NE=NT*3, NR=NE*3,
                    Rnames=defaultRNames(NL, NT, NE, NR),
                    psi=defaultPsi(NL),
                    lj=defaultLJ(NL, NT, NE, NR),
                    tj=defaultTJ(NL, NT, NE, NR),
                    kj=defaultKJ(NL, NT, NE, NR),
                    nu=defaultNu(Rnames),
                    G=defaultG(NL, NT, NE, NR),
                    eta0=homoEta(Inf), eta1=homoEta(0),
                    Xi=homoXi(.3), Theta=homoTheta(.1),
                    TechWeights=diag(NT),
                    cTau=homoTau(0)) {

  result <- list()
  result$N <- N
  result$NL <- NL
  result$NT <- NT
  result$NE <- NE
  result$NR <- NR
  result$Rnames <- Rnames
  result$lj <- lj
  result$tj <- tj
  result$kj <- kj
  result$psi <- psi
  result$G  <- G

  if (is.function(Theta)) { Theta <- Theta(lj, tj, kj, Rnames) }
  if (is.function(Xi))    { Xi    <- Xi   (lj, tj, kj, Rnames) }
  if (is.function(cTau))  { cTau  <- cTau (lj, tj, kj, Rnames) }
  if (is.function(eta0))  { eta0  <- eta0 (lj, tj, kj, Rnames) }
  if (is.function(eta1))  { eta1  <- eta1 (lj, tj, kj, Rnames) }
  if (is.function(nu))    { nu    <- nu   (lj, tj, kj, Rnames) }
  if (is.function(TechWeights)) {
    TechWeights <- TechWeights(lj, tj, kj, Rnames)
  }

  result$eta0 <- eta0
  result$eta1 <- eta1
  result$Xi <- Xi
  result$Theta <- Theta
  result$nu <- nu

  result$noAccImp <- result$noRejImp <- 0L
  result$noAccCor <- result$noRejCor <- 0L
  
  techsPerLatent <- table(lj[match(unique(tj),tj)])

  Tau <- cTau %*% solve(TechWeights)
  result$Tau <- Tau
  result$TechWeights <- TechWeights

  result <- setNames(result)
  
  ## Latent variables
  result$L <- rmvnorm(N, mean=rep(0, NL), sigma=psi)
  colnames(result$L) <- unique(as.character(lj))
  
  result$T <- rmvnorm(N, mean=rep(0, NT), sigma=Tau)
  colnames(result$T) <- unique(as.character(tj))
  
  result$E <- rmvnorm(N, mean=rep(0, NE), sigma=Xi)
  colnames(result$E) <- unique(as.character(kj))
  
  result$R <- rmvnorm(N, mean=rep(0, NR), sigma=Theta)
  colnames(result$R) <- Rnames

  result$Theta.draw <- unname(!(table(result$kj)==1)[result$kj])
  names(result$Theta.draw) <- result$Rnames
  
  E2T <- result$tj[match(unique(result$kj),result$kj)]
  result$Xi.draw <- unname(!(table(E2T)==1)[E2T])
  names(result$Xi.draw) <- unique(result$kj)

  ## Censoring

  X <- getX(result)
  cens.prob <- 1/(1+exp(-rep(result$eta0[kj], each=N) -
                        rep(result$eta1[kj], each=N) * X))
  result$I <- 1 - (runif(length(X)) > cens.prob)
  result$X <- X

  ## Mapping matrices 
  result$G.mat <- G.mat(result)
  result$B.mat <- B.mat(result)
  result$Lambda.mat <- Lambda.mat(result)
  
  ## Done
  
  result
}

chkState <- function(state) {

  with(state, {
  
    ## N
    stopifnot(is.numeric(N), length(N) == 1)
    
    ## NL
    stopifnot(is.numeric(NL), length(NL) == 1, NL <= NT)

    ## NT
    stopifnot(is.numeric(NT), length(NT) == 1, NT >= NL, NT <= NE)
              
    ## NE
    stopifnot(is.numeric(NE), length(NE) == 1, NE >= NT, NR <= NR)
    
    ## NR
    stopifnot(is.numeric(NR), length(NR) == 1, NR >= NE)
      
    ## lj
    stopifnot(is.character(lj), length(lj) == NR, length(unique(lj)) == NL)
    stopifnot(all(sapply(tapply(lj, tj, unique, simplify=FALSE), length) == 1))
    
    ## tj
    stopifnot(is.character(tj), length(tj) == NR, length(unique(tj)) == NT)
    stopifnot(all(sapply(tapply(tj, kj, unique, simplify=FALSE), length) == 1))
    
    ## kj
    stopifnot(is.character(kj), length(kj) == NR, length(unique(kj)) == NE)
              
    ## Rnames
    stopifnot(is.character(Rnames), length(Rnames) == NR,
              all(sub("^([^\\.]+)\\..*$", "\\1", Rnames) == lj),
              all(sub("^[^\\.]+\\.([^\\.]+)\\..*$", "\\1", Rnames) == tj),
              all(sub("^[^\\.]+\\.[^\\.]+\\.([^\\.]+).*$", "\\1", Rnames)
                  == kj))
    
    ## G
    stopifnot(is.numeric(G), length(G) == NE, all(G > 0),
              all(names(G) == unique(kj)))
              
    ## G.mat
    stopifnot(is.numeric(G.mat), is.matrix(G.mat), nrow(G.mat) == NE,
              ncol(G.mat) == NL, all(rownames(G.mat) == unique(kj)),
              all(colnames(G.mat) == unique(lj)),
              all(rowSums(G.mat) == G),
              all(lj[match(unique(kj), kj)] ==
                  colnames(G.mat)[apply(G.mat,1, function(x) which(x!=0))]))
    
    ## psi
    stopifnot(is.numeric(psi), is.matrix(psi), nrow(psi) == NL,
              ncol(psi) == NL, all(diag(psi) == 1), isSymmetric(psi),
              all(psi[lower.tri(psi)] < 1),
              all(rownames(psi) == unique(lj)),
              all(colnames(psi) == unique(lj)))
    
    ## Tau
    stopifnot(is.numeric(Tau), is.matrix(Tau), nrow(Tau) == NT,
              ncol(Tau) == NT, all(rownames(Tau) == unique(tj)),
              all(colnames(Tau) == unique(tj)), all(diag(Tau) > 0),
              isSymmetric(Tau), all(Tau[lower.tri(Tau)] == 0))
      
    ## Xi
    stopifnot(is.numeric(Xi), is.matrix(Xi), nrow(Xi) == NE, ncol(Xi) == NE,
              isSymmetric(Xi), all(rownames(Xi) == unique(kj)),
              all(colnames(Xi) == unique(kj)),
              all(diag(Xi) >= 0), all(Xi[lower.tri(Xi)] == 0),
              all((diag(Xi) != 0) == Xi.draw))
    
    ## Theta
    stopifnot(is.numeric(Theta), is.matrix(Theta), nrow(Theta) == NR,
              ncol(Theta) == NR, all(rownames(Theta) == Rnames),
              all(colnames(Theta) == Rnames), isSymmetric(Theta),
              all(diag(Theta) >= 0), all(Theta[lower.tri(Theta)] == 0),
              all((diag(Theta) != 0) == Theta.draw))
    
    ## eta0
    stopifnot(is.numeric(eta0), length(eta0) == NE,
              all(names(eta0) == unique(kj)))
    
    ## eta1
    stopifnot(is.numeric(eta1), length(eta1) == NE, all(eta1 >= 0),
              all(names(eta1) == unique(kj)))
    
    ## L
    stopifnot(is.numeric(L), is.matrix(L), nrow(L) == N, ncol(L) == NL,
              all(colnames(L) == unique(lj)))
              
    ## T
    stopifnot(is.numeric(T), is.matrix(T), nrow(T) == N, ncol(T) == NT,
              all(colnames(T) == unique(tj)))
    
    ## E
    stopifnot(is.numeric(E), is.matrix(X), nrow(E) == N, ncol(E) == NE,
              all(colnames(E) == unique(kj)),
              all(E[,!Xi.draw] == 0))
    
    ## R
    stopifnot(is.numeric(R), is.matrix(R), nrow(R) == N, ncol(R) == NR,
              all(colnames(R) == Rnames),
              all(abs(R[,!Theta.draw]) < 1e-10))
    
    ## X
    stopifnot(is.numeric(X), is.matrix(X), nrow(X) == N, ncol(X) == NR,
              all(colnames(X) == Rnames))
    
    ## I
    stopifnot(is.numeric(I), is.matrix(I), nrow(I) == N, ncol(I) == NR,
              all(colnames(I) == Rnames), all(I %in% c(0,1)))
    
    ## nu
    stopifnot(is.numeric(nu), length(nu) == NR, names(nu) == Rnames)
    
    ## TechWeights
    stopifnot(is.numeric(TechWeights), is.matrix(TechWeights),
              isSymmetric(TechWeights),
              all(colnames(TechWeights) == unique(tj)),
              all(rownames(TechWeights) == unique(tj)),
              all(TechWeights[lower.tri(TechWeights)] == 0),
              all(diag(TechWeights) > 0))
    
    ## Theta.draw
    stopifnot(is.logical(Theta.draw), length(Theta.draw) == NR,
              names(Theta.draw) == Rnames,
              all(Theta.draw == unname(!(table(kj)==1)[kj])))
                  
    ## Xi.draw
    E2T <- tj[match(unique(kj), kj)]
    stopifnot(is.logical(Xi.draw), length(Xi.draw) == NE,
              names(Xi.draw) == unique(kj),
              all(Xi.draw == unname(!(table(E2T)==1)[E2T])))
    
    ## B.mat
    stopifnot(is.numeric(B.mat), is.matrix(B.mat), nrow(B.mat) == NE,
              ncol(B.mat) == NT, all(B.mat %in% c(0,1)),
              all(rownames(B.mat) == unique(kj)),
              all(colnames(B.mat) == unique(tj)),
              all(rowSums(B.mat) == 1),
              all(colnames(B.mat)[apply(B.mat, 1,
                                        function(x) which(x!=0))] ==
                  tj[match(unique(kj), kj)]))

    ## Lambda.mat
    stopifnot(is.numeric(Lambda.mat), is.matrix(Lambda.mat),
              nrow(Lambda.mat) == NR, ncol(Lambda.mat) == NE,
              all(rownames(Lambda.mat) == Rnames),
              all(colnames(Lambda.mat) == unique(kj)),
              all(Lambda.mat %in% c(0,1)),
              all(rowSums(Lambda.mat) == 1),
              all(colnames(Lambda.mat)[apply(Lambda.mat, 1,
                             function(x) which(x!=0))] == kj))

    ## Check that the residuals sum up
    stopifnot(max(abs(getX(state) - X)) < 1e-10)
  })
}
     
chkSampState <- function(start, samp) {
  samp$X <- start$X
  samp$I <- start$I
  samp$X[samp$I == 0] <- samp$X.imp
  samp$R <- samp$X - (samp$L[,samp$lj] %*% diag(samp$G[samp$kj]) +
                      samp$E[,samp$kj] + samp$T[,samp$tj] +
                      rep(samp$nu, each=samp$N))
  chkState(samp)
}

datamatrix <- function(data) {
  res <- data[c("N", "NL", "NT", "NE", "NR", "Rnames", "lj", "kj", "tj")]
  res$X <- getX(data)
  res
}
