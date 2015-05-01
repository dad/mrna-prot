
'
The model is

X(i,j) = L(i,l[j]) G(k[j]) + T(i,t[j]) + E(i,k[j]) + R(i,j) + nu(j)

And we have the following notation:

L: N x NL
T: N x NT
E: N x NE
R: N x NR
X: N x NR

G: vector of length NE
nu: vector of length NR
eta0, eta1: vectors of length NR

tj: vector of length NR, factor with tech names as levels
kj: vector of length NR, factor with experiment names as levels
lj: vector of length NR, factor with latent names

Psi = Cov(L), NL x NL matrix
Xi: covariance matrix for E, diagonal, NE x NE
Theta: covariance matrix for Theta, diagonal, NR x NR
Tau: covariance matrix for T, diagonal, NT x NT 
'

rmvnorm <- function (n, mean = rep(0, nrow(sigma)),
                     sigma = diag(length(mean)),
                     method = c("eigen", "svd", "chol"),
                     pre0.9_9994 = FALSE, checks=FALSE) {
  if (checks) {
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE)) {
      stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) {
      stop("mean and sigma have non-conforming size")
    }
    sigma1 <- sigma
    dimnames(sigma1) <- NULL
    if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
      warning("sigma is numerically not symmetric")
    }
  }
  method <- match.arg(method)
  if (method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive definite")
    }
    ev$values <- zapsmall(ev$values)
    retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*%
      t(ev$vectors)
  } else if (method == "svd") {
    sigsvd <- svd(sigma)
    if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  } else if (method == "chol") {
    retval <- chol(sigma, pivot = TRUE)
    o <- order(attr(retval, "pivot"))
    retval <- retval[, o]
  }
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n,
                   byrow = !pre0.9_9994) %*% retval
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}

myrmvnorm <- function(n=nrow(mean),
                      mean=matrix(0, nrow=n, ncol=nrow(sigma)),
                      sigma=diag(length(mean)),
                      method=c("eigen", "svd", "chol"),
                      checks=TRUE) {
  method <- match.arg(method)
  if (checks) {
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                     check.attributes = FALSE)) {
      stop("sigma must be a symmetric matrix")
    }
    if (nrow(mean) != n || ncol(mean) != nrow(sigma)) {
      stop("mean and n and sigma have non-conforming size")
    }
    sigma1 <- sigma
    dimnames(sigma1) <- NULL
    if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
      warning("sigma is numerically not symmetric")
    }
  }
  if (method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive definite")
    }
    ev$values <- zapsmall(ev$values)
    retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*%
      t(ev$vectors)
  } else if (method == "svd") {
    sigsvd <- svd(sigma)
    if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  } else if (method == "chol") {
    retval <- chol(sigma, pivot = TRUE)
    o <- order(attr(retval, "pivot"))
    retval <- retval[, o]
  }
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
  retval <- retval + mean
  colnames(retval) <- colnames(mean)
  retval
}

## map from NL to NE, scale by G
G.mat <- function(state){

  with(state, {
    G.mat <- matrix(0,nrow=length(G),ncol=state$NL)
    colnames(G.mat) <- unique(lj)
    rownames(G.mat) <- unique(kj)
    for(l in unique(lj)){
      kidx <- as.character(kj[lj==l])
      G.mat[kidx,l] <- G[kidx]
    }
    name(G.mat,state)
  })
}

## map from NT to NE
B.mat <- function(state){
  with(state, {
    B.mat <- matrix(0,nrow=NE,ncol=NT)
    for(k in 1:NE){
      B.mat[k,tj[match(unique(kj)[k],kj)]==unique(tj)] <- 1
    }
    name(B.mat,state)
  })
}

## map from NE to NR
Lambda.mat <- function(state){
  with(state, {
    Lambda.mat <- matrix(0,nrow=NR,ncol=NE)
    for(j in 1:NR){
      Lambda.mat[j,kj[j]==unique(kj)] <- 1
    } 
    name(Lambda.mat,state)
  })
}


V.LTE <- function(state){
  with(state, {
    d <- nrow(psi)+nrow(Tau)+nrow(Xi)
    v.lte <- matrix(0,nrow=d,ncol=d)
    v.lte[1:NL,1:NL] <- psi
    v.lte[(NL+1):(NL+NT),(NL+1):(NL+NT)] <- Tau
    v.lte[(NL+NT+1):d,(NL+NT+1):d] <- Xi
    v.lte
  })
}

V.X <- function(state){
  with(state, {
    T2R <- Lambda.mat%*%B.mat
    v.x <- Lambda.mat%*%G.mat%*%psi%*%t(G.mat)%*%t(Lambda.mat)+T2R%*%Tau%*%t(T2R)+
      Lambda.mat%*%Xi%*%t(Lambda.mat)+Theta
    rownames(v.x) <- Rnames
    v.x })
}
V.X.inv <- function(state){
  solve(V.X(state))
}

V.LTEX <- function(state){
   with(state, {
       v.ltex <- t(cbind(Lambda.mat%*%G.mat%*%psi,Lambda.mat%*%B.mat%*%Tau,Lambda.mat%*%Xi))
       colnames(v.ltex) <- colnames(R)
       rownames(v.ltex) <- c(as.character(unique(lj)),as.character(unique(tj)),as.character(unique(kj)))
       v.ltex})
}

createLavaanModel <- function(datamatrix,
                              tau=c("response", "full", "diag")) {

  tau <- match.arg(tau)

  ## Unit variance for response variables
  lvar <- datamatrix$lj ~ " ~~ 1 * " ~ datamatrix$lj
  ##
  l1 <- datamatrix$lj ~ " =~ 1 * " ~ datamatrix$tj
  ##
  l2 <- datamatrix$tj ~ " =~ 1 * " ~ datamatrix$kj
  ##
  l3 <- datamatrix$kj ~ " =~ 1 * " ~ datamatrix$Rnames
  if (tau == "full") {
    tjvar <- combn(unique(D$tj), 2)
    l22 <- tjvar[1,] ~ " ~~ " ~ tjvar[2,]
  } else if (tau == "response") {
    tjvar <- tapply(datamatrix$tj, datamatrix$lj,
                    function(x) combn(unique(x), 2))
    tjvar <- do.call(cbind, tjvar)
    l22 <- tjvar[1,] ~ " ~~ " ~ tjvar[2,]
  } else {
    l22 <- ""
  }
  ## Intercept terms
  lnu <- datamatrix$Rnames ~ " ~ 1"
  
  paste(unique(c(lvar, l1, l2, l3, l22, lnu)), collapse="\n")
}

createStart <- function(datamatrix, method="EM",
                        tau=c("response", "full", "diag"), ...) {
  method <- match.arg(method)
  if (method == "EM") {
    X <- as.data.frame(datamatrix$X)
    model1 <- createLavaanModel(datamatrix, tau=tau)
    fit1 <- cfa(model1, data=X, std.lv=FALSE, meanstructure=FALSE)
    fit1
  }
}

name <- function(what, base) {
  for (n in names(what)) {
    names(what[[n]]) <- names(base[[n]])
    dimnames(what[[n]]) <- dimnames(base[[n]])
  }
  what
}

drawSingleCor <- function(newpsi, state, idx, V) {
  with(state, {
    f <- function(newval) {
      cm2 <- newpsi
      cm2[idx[1], idx[2]] <- cm2[idx[2], idx[1]] <- newval
      det(cm2)
    }

    ## Likelihood
    like <- function(mycor) {
      invpsi <- solve(mycor)
      t0 <- - N * NL * log(2*pi)
      t1 <- - N * log(det(mycor))
      t2 <- - sum(invpsi * V)
      (t0 + t1 + t2) /2
    }
    like.idx <- function(scor) {
      cm2 <- newpsi
      cm2[ idx[1], idx[2] ] <- cm2[ idx[2], idx[1] ] <- scor
      like(cm2)
    }

    ## Calculate boundaries
    f1 <- f(1) ; f0 <- f(0); fm1 <- f(-1)
    a <- (f1 + fm1 - 2*f0) / 2
    b <- (f1 - fm1) / 2
    c <- f0
    if (abs(a) > 1e-14) {
      ## Quadratic
      di <- b*b - 4*a*c
      if (di < 0) {
        i1 <- -1 ; i2 <- 1
      } else {
        i1 <- min(max((-b - sqrt(di)) / (2*a), -1), 1)
        i2 <- min(max((-b + sqrt(di)) / (2*a), -1), 1)
        if (i1 > i2) { tmp <- i1 ; i1 <- i2 ; i2 <- tmp }
      }
    } else {
      ## Linear
      if (abs(b) < 1e-14) {
        i1 <- -1 ; i2 <- 1;
      } else if (b < 0) {
        i1 <- -1 ; i2 <- min(-c/b, 1)
      } else if (b > 0) {
        i1 <- max(-c/b, -1) ; i2 <- 1
      }
    }

    cur <- newpsi[idx[1], idx[2]]
    propsd <- sqrt(1/N)/10

    prop <- rtruncnorm(1, mean=cur, sd=propsd, a=i1, b=i2)
    cur.like <- like(newpsi)
    prop.like <- like.idx(prop)

    trans.like <- dtruncnorm(prop, a=i1, b=i2, mean=cur, sd=propsd)
    rev.like <- dtruncnorm(cur, a=i1, b=i2, mean=prop, sd=propsd)

    if (prop.like > cur.like ||
        runif(1) < exp(prop.like - cur.like) * rev.like / trans.like) {
      newpsi[idx[1], idx[2]] <- newpsi[idx[2], idx[1]] <- prop
    }
    newpsi
  })
}

drawPsi <- function(state, priorScale=diag(0, state$NL),
                    priorDF=0) {
  newstate <- with(state, {
    idx <- expand.grid(1:NL, 1:NL)
    idx <- idx[idx[,1] < idx[,2],]
    idx <- lapply(1:nrow(idx), function(i) unlist(idx[i,]))
    V <- rowSums(sapply(1:N, function(i) {
      as.vector(cbind(L[i,]) %*% rbind(L[i,]))
    }))
    newpsi <- psi
    for (i in idx) {
      newpsi <- drawSingleCor(newpsi, state, i, V)
    }
    name(list(psi=newpsi),state)
  })
  state[names(newstate)] <- newstate
  state
}

drawThetaNu <- function(state, nu0=3, sigma02=1/5) {
  newstate <- with (state, {
    ## Dependent variable to fit intercept only
    NR.draw <- sum(Theta.draw)
    DV <- X - L[, lj] %*% diag(G[kj]) - T[, tj] - E[, kj]
    DV <- DV[,Theta.draw]
    nu.hat <- apply(DV, 2, mean)
    df <- nu0 + (N - 1)
    
    scale <- (N * apply(DV, 2, var) + nu0 * sigma02) / (nu0 + N)
    diag(Theta)[Theta.draw]  <- as.numeric(df * scale / rchisq(NR.draw, df))

    nu[Theta.draw] <- as.vector(rnorm(NR.draw, mean=nu.hat, diag(Theta)[Theta.draw] / N))
    R[,Theta.draw] <- t(t(DV) - nu[Theta.draw])
    name(list(Theta=Theta, nu=nu, R=R), state)
  })
  
  state[names(newstate)] <- newstate
  
  state
}

## Draw global G
drawG <- function(state, nu0=3,sigma02=1/5){
    newstate <- with(state, {
        topresR <-  Theta.draw
        topresE <- !Theta.draw &  Xi.draw[kj]
        topresT <- !Theta.draw & !Xi.draw[kj]

        ## TODO: Check this
        DV <- (X - rep(!topresE, each=N) * E[,kj] -
               rep(!topresT, each=N) * T[,tj]) - rep(nu, each=N)

        weights <- diag(Theta)+diag(Xi)[kj]*(!topresT&!topresR)+diag(Tau)[tj]*topresT

        for(l in unique(lj)){

            lweights <- rep(weights[lj==l],each=N)
            Y <- as.vector(DV[,lj==l])/sqrt(lweights)
            W <- as.vector(L[, lj[lj==l]])/sqrt(lweights)

            V.beta <- solve(t(W)%*%W)
            Beta.hat <- V.beta%*%t(W)%*%Y

            Beta.samp <- rnorm(1,mean=Beta.hat,sd=sqrt(V.beta))
            G[unique(kj[lj==l])] <- Beta.samp
        }

        name(list(G=G), state)
    })

    state[names(newstate)] <- newstate
    state$G.mat <- G.mat(state)
    state

}

draw_stupidG <- function(state, nu0=3,sigma02=1/5){
    newstate <- with(state, {
        topresR <-  Theta.draw
        topresE <- !Theta.draw &  Xi.draw[kj]
        topresT <- !Theta.draw & !Xi.draw[kj]

        ## TODO: Check this
        DV <- (X - rep(!topresE, each=N) * E[,kj] -
               rep(!topresT, each=N) * T[,tj]) - rep(nu, each=N)

        weights <- diag(Theta)+diag(Xi)[kj]*(!topresT&!topresR)+diag(Tau)[tj]*topresT

        lweights <- rep(weights,each=N)
        Y <- as.vector(DV[,lj==lj])/sqrt(lweights)
        W <- as.vector(L[, lj])/sqrt(lweights)
        
        V.beta <- solve(t(W)%*%W)
        Beta.hat <- V.beta%*%t(W)%*%Y
        
        Beta.samp <- rnorm(1,mean=Beta.hat,sd=sqrt(V.beta))
        G[unique(kj)] <- Beta.samp
        
        name(list(G=G), state)
    })

    state[names(newstate)] <- newstate
    state$G.mat <- G.mat(state)
    state

}

drawXi <- function(state, nu0=3,sigma02=1/5){
    newstate <- with (state, {

    ## Fit slope only
    NE.draw <- sum(Xi.draw)
    E.idx <- match(unique(kj), kj)
    res <- X[, E.idx] - L[, lj[E.idx]]*rep(G[kj[E.idx]],each=N) - T[, tj[E.idx]] - R[, E.idx] - rep(nu[E.idx], each=N)
    res <- res[,Xi.draw]
    E[,Xi.draw] <- res

    df <- nu0 + N
    scale <- (diag(t(res) %*% res) + nu0 * sigma02) / (nu0 + N)
    Xi.tmp <- diag(as.numeric(df * scale / rchisq(NE.draw, df)))

    Xi[Xi.draw,Xi.draw] <- Xi.tmp
    name(list(E=E, Xi=Xi), state)
  })

    state[names(newstate)] <- newstate
    state
}

drawXiNu <- function(state, nu0=3,sigma02=1/5){
    newstate <- with (state, {
      
    ## Find single replicate, single experiment case 
    toDraw <- kj[!Theta.draw&Xi.draw[kj]]
    NE.draw <- length(toDraw)
    E.idx <- match(unique(kj), kj)

    ## R should be zero but subtract anyway
    DV <- X[, E.idx,drop=FALSE] - L[, lj[E.idx]]*rep(G[kj[E.idx]],each=N) -
        T[, tj[E.idx],drop=FALSE] - R[, E.idx,drop=FALSE]
    colnames(DV) <- unique(kj)
    DV <- DV[,toDraw,drop=FALSE]

    for(idx in toDraw){
      ## Fit slope / intercept
      W <- rep(1,N)
      V.beta <- 1/N
      Beta.hat <- mean(DV[,idx])
      
      res <- (DV[,idx] - rep(Beta.hat,N))

      df <- nu0 + N
      scale <- (diag(t(res) %*% res) + nu0 * sigma02) / (nu0 + N)

      Xi.samp <- as.numeric(df * scale / rchisq(1, df))
      
      Beta.samp <- rnorm(1, mean=Beta.hat, sd=sqrt(V.beta*Xi.samp))
      nu[kj==idx] <- Beta.samp
      E[,idx] <- DV[,idx]-Beta.samp
      diag(Xi)[idx] <- Xi.samp
  }
    name(list(E=E, Xi=Xi,nu=nu), state)
})
  state[names(newstate)] <- newstate
  state$G.mat <- G.mat(state)
  state    
}

drawTauNu <- function(state, nu0=3,sigma02=1/5,poolTau=FALSE){
      newstate <- with (state, {
    ## Find single replicate, single experiment cases
    toDraw <- unique(kj[!Theta.draw&!Xi.draw[kj]])
    NE.draw <- length(toDraw)
    E.idx <- match(toDraw, kj)
    t.names <- tj[E.idx]
    
    ## Fit slope only
    ## E and R should be zero but subtract anyway
    DV <- X[, E.idx,drop=FALSE] - L[, lj[E.idx]]*rep(G[kj[E.idx]],each=N) -
        E[,kj[E.idx],drop=FALSE] - R[, E.idx,drop=FALSE]

    for(idx in toDraw){
      W <- cbind(rep(1,N))
      V.beta <- 1/N
      DVcol <- DV[,E2R(state,idx),drop=FALSE]
      Beta.hat <- mean(DVcol)
      res <- (DVcol - rep(Beta.hat,1))
      l <- lj[kj==idx]
      t <- tj[kj==idx]

      if(poolTau){
          ltechs <- unique(tj[lj==l])
          T.tmp <- T[,ltechs,drop=FALSE]
          T.tmp[,t] <- res
          df <- nu0 + (N-1)
          NT.l <- length(ltechs)
          NN <- N * NT.l
          scale <- (sum(as.vector(T.tmp%*%sqrt(TechWeights[ltechs,ltechs,drop=FALSE]))^2) +
                    nu0*sigma02) / (nu0+NN)
          cTau.samp <- as.numeric(df * scale / rchisq(1, df))
          Tau[ltechs,ltechs] <- cTau.samp*solve(TechWeights[ltechs,ltechs])
      }else{
          df <- nu0 + (N-1)
          scale <- (sum(res^2) + nu0*sigma02) / (nu0+N)
          Tau[t,t] <- as.numeric(df * scale / rchisq(1, df))
      }
      
      Beta.samp <- rnorm(1, mean=Beta.hat, sd=sqrt(Tau[t,t]/N))
      nu[E2R(state,idx)] <- Beta.samp
      T[,t] <- DVcol - nu[E2R(state,idx)]

    }
    
    name(list(T=T, nu=nu, Tau=Tau), state)
    
  })
  state[names(newstate)] <- newstate
  state$G.mat <- G.mat(state)
  state
}

## Single replicate, single experiments
drawGXiNu <- function(state, nu0=3, sigma02=1/5){
  newstate <- with (state, {
      
    ## Find single replicate, single experiment case 
    toDraw <- kj[!Theta.draw&Xi.draw[kj]]
    NE.draw <- length(toDraw)
    E.idx <- match(unique(kj), kj)

    ## R should be zero but subtract anyway
    DV <- X[, E.idx,drop=FALSE] - T[, tj[E.idx],drop=FALSE] - R[, E.idx,drop=FALSE]
    colnames(DV) <- unique(kj)
    DV <- DV[,toDraw,drop=FALSE]
    LMat <- L[, lj[E.idx]]
    colnames(LMat) <- unique(kj)
    LMat <- LMat[,toDraw,drop=FALSE]

    for(idx in toDraw){
      ## Fit slope / intercept
      W <- cbind(rep(1,N),LMat[,idx])
      V.beta <- solve(t(W)%*%W)
      Beta.hat <- V.beta%*%t(W)%*%DV[,idx]
      
      res <- (DV[,idx] - W %*% Beta.hat)

      df <- nu0 + N
      scale <- (diag(t(res) %*% res) + nu0 * sigma02) / (nu0 + N)
      Xi.samp <- as.numeric(df * scale / rchisq(1, df))

      Beta.samp <- rmvnorm(1, mean=Beta.hat, sigma=V.beta*Xi.samp)
      nu[kj==idx] <- Beta.samp[1]
      G[idx] <- Beta.samp[2]
      E[,idx] <- DV[,idx] - W%*%t(Beta.samp)
      diag(Xi)[idx] <- Xi.samp
  }
    name(list(G=G, E=E, Xi=Xi,nu=nu), state)
})
  state[names(newstate)] <- newstate
  state$G.mat <- G.mat(state)
  state

  
}

drawGTau <- function(state, poolTau, nu0=3, sigma02=1/5) {
  newstate <- with (state, {

    ## Find multiple replicate, single experiment cases
    toDraw <- unique(kj[Theta.draw&!Xi.draw[kj]])
    NE.draw <- length(toDraw)
    E.idx <- match(toDraw, kj)
    t.names <- unique(tj[E.idx])
    
    ## Fit slope only
    ## E should be zero but subtract anyway
    DV <- X[, E.idx,drop=FALSE] - E[,kj[E.idx],drop=FALSE] - R[, E.idx,drop=FALSE] - rep(nu[E.idx], each=N)
    LMat <- L[, lj[E.idx],drop=FALSE]

    G.hat <- diag(t(LMat) %*% DV) / diag(t(LMat) %*% LMat)
    names(G.hat) <- tj[E.idx]
    
    res <- (DV - LMat %*% diag(x=G.hat,nrow=length(G.hat)))
    colnames(res) <- tj[E.idx]

    if(poolTau){
        for(l in unique(lj[E.idx])){
            ltechs <- unique(tj[lj==l])
            T.tmp <- T[,ltechs,drop=FALSE]
            t.draw <- intersect(t.names,ltechs)
            T.tmp[,t.draw] <- res[,t.draw]
            df <- nu0 + (N-1)
            NT.l <- length(ltechs)
            NN <- N * NT.l
            scale <- (sum(as.vector(T.tmp%*%sqrt(TechWeights[ltechs,ltechs,drop=FALSE]))^2) +
                      nu0*sigma02) / (nu0+NN)
            cTau.samp <- as.numeric(df * scale / rchisq(1, df))
            Tau[ltechs,ltechs] <- cTau.samp*solve(TechWeights[ltechs,ltechs])

            G.var <- cTau.samp/diag(TechWeights)[t.draw] / as.numeric(t(LMat) %*% LMat)
            G.samp <- rnorm(length(t.draw), mean=G.hat[t.draw], sd=sqrt(G.var))
            G[kj[match(t.draw,tj)]] <- G.samp

            T[,t.draw] <- DV - (LMat %*% diag(x=G.samp,nrow=length(G.samp)))
            
        }} else{

            T.tmp <- res[,t.names,drop=FALSE]
            df <- nu0 + (N-1)
            scale <- apply(res,2,function(T.tmp) (sum(T.tmp^2)+nu0*sigma02)/(nu0+N))
            diag(Tau)[t.names] <- as.numeric(df * scale / rchisq(length(t.names), df))

            G.var <- diag(Tau)[t.names] / as.numeric(t(LMat) %*% LMat)
            G.samp <- rnorm(length(t.names), mean=G.hat, sd=sqrt(G.var))
            G[kj[match(t.names,tj)]] <- G.samp

            T[,t.names] <- DV - (LMat %*% diag(x=G.samp,nrow=length(G.samp)))
        }
    
    name(list(G=G, T=T, Tau=Tau), state)
    
  })
  state[names(newstate)] <- newstate
  state$G.mat <- G.mat(state)
  state
}


drawGTauNu <- function(state, poolTau, nu0=3, sigma02=1/5) {
  newstate <- with (state, {
    ## Find single replicate, single experiment cases
    toDraw <- unique(kj[!Theta.draw&!Xi.draw[kj]])
    NE.draw <- length(toDraw)
    E.idx <- match(toDraw, kj)
    t.names <- tj[E.idx]
    
    ## Fit slope only
    ## E and R should be zero but subtract anyway
    DV <- X[, E.idx,drop=FALSE] - E[,kj[E.idx],drop=FALSE] - R[, E.idx,drop=FALSE]
    LMat <- L[, lj[E.idx],drop=FALSE]
    colnames(LMat) <- kj[E.idx]

    for(idx in toDraw){
      W <- cbind(rep(1,N),LMat[,idx])
      V.beta <- solve(t(W)%*%W)
      DVcol <- DV[,E2R(state,idx),drop=FALSE]
      Beta.hat <- V.beta%*%t(W)%*%DVcol
      res <- (DVcol - W %*% Beta.hat)
      l <- lj[kj==idx]
      t <- tj[kj==idx]

      if(poolTau){
          ltechs <- unique(tj[lj==l])
          T.tmp <- T[,ltechs,drop=FALSE]
          T.tmp[,t] <- res
          df <- nu0 + (N-1)
          NT.l <- length(ltechs)
          NN <- N * NT.l
          scale <- (sum(as.vector(T.tmp%*%sqrt(TechWeights[ltechs,ltechs,drop=FALSE]))^2) +
                    nu0*sigma02) / (nu0+NN)
          cTau.samp <- as.numeric(df * scale / rchisq(1, df))
          Tau[ltechs,ltechs] <- cTau.samp*solve(TechWeights[ltechs,ltechs])
      }else{
          df <- nu0 + (N-1)
          scale <- (sum(res^2) + nu0*sigma02) / (nu0+N)
          Tau[t,t] <- as.numeric(df * scale / rchisq(1, df))
      }
      
      Beta.samp <- rmvnorm(1, mean=Beta.hat, sigma=V.beta*Tau[t,t])
      nu[E2R(state,idx)] <- Beta.samp[1]
      G[idx] <- Beta.samp[2]
      T[,t] <- DVcol - (LMat[,idx] * Beta.samp[2]) - nu[E2R(state,idx)]

    }
    
    name(list(G=G, T=T, nu=nu, Tau=Tau), state)
    
  })
  state[names(newstate)] <- newstate
  state$G.mat <- G.mat(state)
  state
}

drawGXi <- function(state, nu0=3, sigma02=1/5) {
  newstate <- with (state, {
  
    ## Fit slope only
    NE.draw <- sum(Xi.draw)
    E.idx <- match(unique(kj), kj)
    DV <- X[, E.idx] - T[, tj[E.idx]] - R[, E.idx] - rep(nu[E.idx], each=N)
    DV <- DV[,Xi.draw]
    LMat <- L[, lj[E.idx]]
    LMat <- LMat[,Xi.draw]


    G.hat <- diag(t(LMat) %*% DV) / diag(t(LMat) %*% LMat)

    res <- (DV - LMat %*% diag(G.hat))

    df <- nu0 + N
    scale <- (diag(t(res) %*% res) + nu0 * sigma02) / (nu0 + N)
    Xi.tmp <- diag(as.numeric(df * scale / rchisq(NE.draw, df)))
    G.tmp <- rnorm(NE.draw, mean=G.hat, sd=sqrt(diag(Xi.tmp) / diag(t(LMat) %*% LMat)) )
    G[Xi.draw] <- G.tmp
    E[,Xi.draw] <- DV - (LMat %*% diag(G.tmp))
    Xi[Xi.draw,Xi.draw] <- Xi.tmp
    name(list(G=G, E=E, Xi=Xi), state)
  })
  state[names(newstate)] <- newstate
  state$G.mat <- G.mat(state)
  state
}


drawLTE <- function(state){

  newstate <- with(state, {

    L2Draw <- unique(lj)
    ## Only draw E & R are nonzero (mult reps, mult exps)
    E2Draw <- unique(kj[Theta.draw&Xi.draw[kj]])
    ## Draw T if either E or R is non-zero
    T2Draw <- unique(tj[Theta.draw|Xi.draw[kj]])

    drawIdx <- c(L2Draw,T2Draw,E2Draw)
    
    mu.LTE <- t(V.LTEX(state) %*% V.X.inv(state) %*% ( t(X)-nu ))[,drawIdx]
    V.LTE.X <- V.LTE(state) - V.LTEX(state) %*% V.X.inv(state) %*% t(V.LTEX(state))
    V.LTE.X <- V.LTE.X[drawIdx,drawIdx]

    LTE <- myrmvnorm(n=nrow(mu.LTE), mean=mu.LTE, sigma=V.LTE.X)
    colnames(LTE) <- drawIdx
    L[,L2Draw] <- LTE[,L2Draw]

    T[,T2Draw] <- LTE[,T2Draw]
    ## For the complete collapsed cases, T = X-LG-nu (E and R should be zero)
    ## Ie E and R are zero (non-identifiable) for this technology
    Tsingle <- setdiff(unique(tj),T2Draw)
    T[,Tsingle] <- t(t(X[,T2R(state,Tsingle),drop=FALSE] -
      (L%*%t(G.mat))[,T2E(state,Tsingle),drop=FALSE])-nu[T2R(state,Tsingle)])
    
    E[,E2Draw] <- LTE[,E2Draw]
    ## Single rep experiments, multiple experiments in a tech
    Esingle <- unique(kj[!Theta.draw&Xi.draw[kj]])
    ## Either when there is only a single experiment for a technology
    Ezero <- unique(kj[!Xi.draw[kj]])
    E[,Ezero] <- 0
    ## For single rep experiments, E = T-LG-T-nu (R = zero)
    E[,Esingle] <- t(t(X[,E2R(state,Esingle)]-(L%*%t(G.mat))[,Esingle]-(T%*%t(B.mat))[,match(Esingle,unique(kj))])-nu[E2R(state,Esingle)])

    R <- t(t(X-(L%*%t(G.mat))%*%t(Lambda.mat)- T%*%t(B.mat)%*%t(Lambda.mat) - E%*%t(Lambda.mat))-nu)
    colnames(T) <- unique(tj)
    colnames(E) <- unique(kj)
    name(list(L=L,T=T, E=E,R=R),state)
  })
  state[names(newstate)] <- newstate
  state

}

drawTau <- function(state, FitTechWeights=state$TechWeights, poolTau,
                    nu0=3, sigma02=1/5){
  newstate <- with(state, {
      if(poolTau){
          for(l in unique(lj)){
              ltechs <- unique(tj[lj==l])
              df <- nu0 + (N-1)
              NT.l <- length(ltechs)
              NN <- N * NT.l
              scale <- (sum(as.vector(T[,ltechs,drop=FALSE]%*%sqrt(TechWeights[ltechs,ltechs,drop=FALSE]))^2)
                        + nu0 * sigma02) / (nu0 + NN)
              cTau <- as.numeric(df * scale / rchisq(1, df))
              Tau[ltechs,ltechs] <- cTau*solve(TechWeights[ltechs,ltechs])
          }
      } else{
          for(t in unique(tj)){
              df <- nu0 + (N-1)
              scale <- (sum(T[,t,drop=FALSE]^2)+ nu0 * sigma02) / (nu0 + N)
              Tau[t,t] <- as.numeric(df * scale / rchisq(1, df))
          }
      }

      name(list(Tau=Tau * diag(NT)), state)
  })

  state[names(newstate)] <- newstate
  state
}


drawObs.myfit <- function(x, y, method=c("bayesglm", "speedglm", "glm"),
                          obs_scale, obs_scale_intercept, ...) {
  method <- match.arg(method)
  yx <- as.formula("y ~ x")
  if (all(y==0) || all(y==1)) {
    c(eta1=NA, eta0=NA)
  } else {
    if (method=="bayesglm") {
      require(arm)
      lfit <- bayesglm(yx, data=data.frame(x=x, y=y),
                       family=binomial(link="logit"), Warning=FALSE,
                       prior.scale = obs_scale,
                       prior.scale.for.intercept = obs_scale_intercept, ...)
      lfitCoef <- as.matrix(coef(summary(lfit)))
      c(eta0=rnorm(1, lfitCoef[1,1], lfitCoef[1,2]),
        eta1=rnorm(1, lfitCoef[2,1], lfitCoef[2,2]))
    } else if (method=="glm") {
      lfit <- glm(yx, data=data.frame(x=x, y=y),
                  family=binomial(link="logit"), ...)
      lfitCoef <- as.matrix(coef(summary(lfit)))
      c(eta0=rnorm(1, lfitCoef[1,1], lfitCoef[1,2]),
        eta1=rnorm(1, lfitCoef[2,1], lfitCoef[2,2]))
    } else if (method=="speedglm") {
      require(speedglm)
      lfit <- speedglm(yx, data=data.frame(x=x, y=y),
                       family=binomial(link="logit"), ...)
      lfitCoef <- as.matrix(coef(summary(lfit)))
      mode(lfitCoef) <- "numeric"
      c(eta0=rnorm(1, lfitCoef[1,1], lfitCoef[1,2]),
        eta1=rnorm(1, lfitCoef[2,1], lfitCoef[2,2]))
    }
  }
}

## untested 
drawObsParam <- function(state, method=c("bayesglm", "speedglm", "glm"),
                         obs_scale, obs_scale_intercept, ...) {
  method <- match.arg(method)
  newstate <- with(state, {
    for (k in unique(kj)) {

      k.indices <- kj==k
      y <- c(as.vector(I[,k.indices]))
      x <- c(as.vector(X[,k.indices]))
      rp <- drawObs.myfit(x, y, method=method, obs_scale = obs_scale,
                          obs_scale_intercept = obs_scale_intercept, ...)

      eta0[k] <- rp["eta0"]
      eta1[k] <- rp["eta1"]
    }
    name(list(eta0=eta0,eta1=eta1), state)
  })

  state[names(newstate)] <- newstate
  state
}

drawMissing <- function(state) {
  newstate <- with(state, {
    missIdx <- which(I == 0, arr.ind=TRUE)
    topresR <-  Theta.draw
    topresE <- !Theta.draw &  Xi.draw[kj]
    topresT <- !Theta.draw & !Xi.draw[kj]

    dataMean <- (X - R - rep(topresE, each=N) * E[,kj] -
                 rep(topresT, each=N) * T[,tj])
    dataVar  <- diag(Theta) + topresE * diag(Xi)[kj] + topresT * diag(Tau)[tj]

    .Call("R_SCM_drawMissing", state, missIdx, as.double(eta0[kj]),
          as.double(eta1[kj]), as.double(thresholds), dataMean,
          dataVar, topresR, topresE, topresT, noAccImp, noRejImp,
          match(kj, colnames(E)), match(tj, colnames(T)), PACKAGE="SCM")
  })
  
  state[names(newstate)] <- newstate
  state
} 

BayesCFA <- function(start, noSteps=50, TechWeights=start$TechWeights,
                     draw=c(), poolTau=FALSE, singleG=FALSE,
                     stupidG=FALSE, obs_scale = NULL,
                     obs_scale_intercept = NULL, verbose=TRUE) {

  samples <- list()

  stepsTaken <- 0L

  current <- start

  current$G.mat <- G.mat(current)
  current$B.mat <- B.mat(current)
  current$Lambda.mat <- Lambda.mat(current)

  current$noAccImp <- current$noRejImp <- 0L
  current$noAccCor <- current$noRejCor <- 0L

  ## Don't draw theta if not identifiable
  current$Theta.draw <- unname(!(table(current$kj)==1)[current$kj])
  names(current$Theta.draw) <- current$Rnames
  #diag(!current$Theta)[current$Theta.draw] <- 0
  
  ## Don't draw xi if not identifiable
  E2T <- current$tj[match(unique(current$kj),current$kj)]
  current$Xi.draw <- unname(!(table(E2T)==1)[E2T])
  names(current$Xi.draw) <- unique(current$kj)
  #diag(current$Xi)[!current$Xi.draw] <- 0

  ## If thresholds are not given, then we assume -Inf
  if (is.null(current$thresholds)) {
    current$thresholds <- rep(-Inf, length(current$Rnames))
  }

  draw <- c(draw, LTE=TRUE, Psi=TRUE, ThetaNu=TRUE, GXi=TRUE,
            GXiNu=TRUE, GTau=TRUE, GTauNu=TRUE, Tau=TRUE,
            ObsParam=TRUE, Missing=TRUE, G=TRUE, Xi=TRUE,
            XiNu=TRUE, TauNu=TRUE)
  draw <- draw[unique(names(draw))]

  while (TRUE) {
    stepsTaken <- stepsTaken + 1

    if (draw["LTE"])      current <- drawLTE(current)
    if (draw["Psi"])      current <- drawPsi(current)
    if (draw["ThetaNu"])  current <- drawThetaNu(current)
    if (!singleG) {
        if (draw["GXi"])      current <- drawGXi(current)
        if (draw["GXiNu"])    current <- drawGXiNu(current)
        if (draw["GTau"])     current <- drawGTau(current, poolTau=poolTau)
        if (draw["GTauNu"])   current <- drawGTauNu(current, poolTau=poolTau)

    } else {
        if (draw["G"]) {
          if(stupidG)
            current <- draw_stupidG(current)
          else
            current <- drawG(current)
        }
        if (draw["Xi"])    current <- drawXi(current)
        if (draw["XiNu"])    current <- drawXiNu(current)
        if (draw["TauNu"])   current <- drawTauNu(current, poolTau=poolTau)
    }
    
    if (draw["Tau"])      current <- drawTau(current, TechWeights,
                                             poolTau=poolTau)
    if (draw["ObsParam"]) current <- drawObsParam(current,
                            obs_scale = obs_scale,
                            obs_scale_intercept = obs_scale_intercept)
    if (draw["Missing"])  current <- drawMissing(current)

    samples[[ length(samples)+1 ]] <- current
    if (verbose) { cat(".") }

    if (stepsTaken == noSteps) { break }
  }

  if (verbose) { cat("\n") }
  samples
}

## TODO: handle arguments

odyCFA <- function(outdir, startdir, samplefiles="sample.*.RData",
                   startfile="start.RData", noRuns=1000,
                   bsub=list(q="airoldi", r="", J="scm[1]",
                     oo=file.path(outdir, "%I.out"),
                     eo=file.path(outdir, "%I.err"))
                   ) {

  ## Create outdir
  dir.create(outdir)

  ## Get the last sample from startdir
  sf <- list.files(path=startdir, pattern=samplefiles)
  sf1 <- sf[ which.max(as.numeric(gsub("^[^0-9]*([0-9]+).*$", "\\1", sf))) ]
  ee <- new.env()
  load(file.path(startdir, sf1), envir=ee)
  samp <- ee[[ ls(ee)[1] ]]

  ## Also get the start sample to have all the fields
  ee <- new.env()
  load(file.path(startdir, startfile), envir=ee)
  start <- ee[[ ls(ee)[1] ]]

  ## Put them together
  samp$I <- start$I
  samp$X <- start$X
  samp$X[ samp$I == 0 ] <- samp$X.imp
  samp$R <- with(samp, { X - L[, lj] %*% diag(G[kj]) - T[, tj] - E[, kj] - rep(nu, each=N) })                 

  ## Put the start file in the output directory
  start <- samp
  save(start, file=file.path(outdir, "start.RData"))

  ## Create an R script that will be run
  run <- function() {

    library(SCM)
    load(file.path(outdir, "start.RData"))
    samp <- start
    for (i in 1:noRuns) {
      next50 <- BayesCFA(samp, noSteps=50)
      samp <- next50[[50]]
      save(samp, file=file.path(outdir, sprintf("sample-%i.RData", i)))
    }
  }

  cat("outdir <- ", paste("'", sep="", outdir, "'"),
      "noRuns <- ", noRuns,
      "run <- ", deparse(run), "run()", sep="\n",
      file=file.path(outdir, "script.R"))

  ## Create a bsub file
  bsub.txt <- paste(sep="\n", paste(sep="", "# BSUB -", names(bsub), " ", bsub, collapse="\n"),
                    paste("Rscript", file.path(outdir, "script.R")))
  cat(bsub.txt, file=file.path(outdir, "run.bsub"))

  ## Submit
  system(paste("bsub < ", file.path(outdir, "run.bsub")))

  invisible(NULL)
}
