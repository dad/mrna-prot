
## Generate fake data and check the V.* covariance matrices

library(SCM)

set.seed(42)

preps <- getPrepsDef(2, 4, 4)

B <- createUnitB(preps)
Lambda1 <- createUnitLambda1(preps)
Lambda2 <- createUnitLambda2(preps)
G1 <- randScaling(createUnitG1(preps), 1.5)
G2 <- randScaling(createUnitG2(preps), 1.5)

data <- genData(numGenes=5000,
                B=B, Lambda1=Lambda1, Lambda2=Lambda2, G1=G1, G2=G2,
                cChi=.01, cEps1=0.001, cEps2=.011)

V.L(data) %plus% - cov(data$L)
V.S(data) %plus% - cov(data$S)
V.L.S(data) %plus% - cov(data$L, data$S)
V.L.P(data) %plus% - cov(data$L, data$P)
max(sapply(V.P(data), function(x) max(abs(x %plus% - cov(data$P)))))
V.S.P(data) %plus% - cov(data$S, data$P)
V.L.X1(data) %plus% - cov(data$L, data$data1)
V.S.X1(data) %plus% - cov(data$S, data$data1)
max(sapply(V.P.X1(data), function(x)
           max(abs(x %plus% - cov(data$P, data$data1)))))
max(sapply(V.X1(data), function(x)
           max(abs(x %plus% - cov(data$data1, data$data1)))))
V.L.X2(data) %plus% - cov(data$L, data$data2)
V.S.X2(data) %plus% - cov(data$S, data$data2)
V.P.X2(data) %plus% - cov(data$P, data$data2)
V.X1.X2(data) %plus% - cov(data$data1, data$data2)
V.X2(data) %plus% - cov(data$data2)

dim(dataMean1(data))

dataCov1(data)

## A data set with some measurements in data2
## and non-symmetric psi matrix

library(SCM)

set.seed(42)

preps <- c("L1.P1.M1", "L1.P1.M2", "L1.P1.M3",
           "L1.P2.M1",
           "L1.P3.M1", "L1.P3.M2",
           "L2.P1.M1",
           "L2.P2.M1", "L2.P2.M2", "L2.P2.M3",
           "L2.P3.M1", "L2.P3.M2", "L2.P3.M3",
           "L2.P4.M1")

psi <- matrix(NA, 4, 4)
psi[1,] <- psi[,1] <- c(2  , 2,   1.6, 1.6)
psi[2,] <- psi[,2] <- c(2  , 3  , 2.6, 2.5)
psi[3,] <- psi[,3] <- c(1.6, 2.6, 2.8, 2.7)
psi[4,] <- psi[,4] <- c(1.6, 2.5, 2.7, 2.8)
colnames(psi) <- rownames(psi) <- psiNames(preps)

B <- randScaling(createUnitB(preps), 1.5)
Lambda1 <- randScaling(createUnitLambda1(preps), 1.1)
Lambda2 <- randScaling(createUnitLambda2(preps), 1.6)
G1 <- randScaling(createUnitG1(preps), 1.5)
G2 <- randScaling(createUnitG2(preps), 1.5)

data <- genData(numGenes=5000, preps=preps, psi=psi,
                B=B, Lambda1=Lambda1, Lambda2=Lambda2, G1=G1, G2=G2,
                cChi=.5, cEps1=.1, cEps2=.6)

V.L(data) %plus% - cov(data$L)
V.S(data) %plus% - cov(data$S)
V.L.S(data) %plus% - cov(data$L, data$S)
V.L.P(data) %plus% - cov(data$L, data$P)
V.P(data) %plus% - cov(data$P)
V.S.P(data) %plus% - cov(data$S, data$P)
V.L.X1(data) %plus% - cov(data$L, data$data1)
V.S.X1(data) %plus% - cov(data$S, data$data1)
V.P.X1(data) %plus% - cov(data$P, data$data1)
V.X1(data) %plus% - cov(data$data1, data$data1)
V.L.X2(data) %plus% - cov(data$L, data$data2)
V.S.X2(data) %plus% - cov(data$S, data$data2)
V.P.X2(data) %plus% - cov(data$P, data$data2)
V.X1.X2(data) %plus% - cov(data$data1, data$data2)
V.X2(data) %plus% - cov(data$data2)

dim(dataMean1(data))
dim(dataMean2(data))

dataCov1(data)
dataCov2(data)

#############################################################
## Test the individual draws

library(SCM)

set.seed(42)

preps <- c("L1.P1.M1", "L1.P1.M2", "L1.P1.M3",
           "L1.P2.M1",
           "L1.P3.M1", "L1.P3.M2",
           "L2.P1.M1",
           "L2.P2.M1", "L2.P2.M2", "L2.P2.M3",
           "L2.P3.M1", "L2.P3.M2", "L2.P3.M3",
           "L2.P4.M1")

psi <- matrix(NA, 4, 4)
psi[1,] <- psi[,1] <- c(2, 2,   0, 0)
psi[2,] <- psi[,2] <- c(2, 3, 0, 0)
psi[3,] <- psi[,3] <- c(0, 0, 2.8, 2.7)
psi[4,] <- psi[,4] <- c(0, 0, 2.7, 2.8)
colnames(psi) <- rownames(psi) <- psiNames(preps)

B <- randScaling(createUnitB(preps), 1.5)
Lambda1 <- randScaling(createUnitLambda1(preps), 1.1)
Lambda2 <- randScaling(createUnitLambda2(preps), 1.6)
G1 <- randScaling(createUnitG1(preps), 1.5)
G2 <- randScaling(createUnitG2(preps), 1.5)

data <- genData(numGenes=5000, preps=preps, psi=psi,
                B=B, Lambda1=Lambda1, Lambda2=Lambda2, G1=G1, G2=G2,
                cChi=.5, cEps1=.1, cEps2=.6)

## Check latent draws

dL <- lapply(1:100, function(x) drawLatent(data)$L)
max(abs(Reduce("+", dL) / length(dL) - data$L))
max(abs(sapply(dL, colMeans) - 0))
max(abs(sapply(dL, apply, 2, var) - diag(psi)[1:2]))

layout(rbind(1:2))
hist(sapply(dL, apply, 2, var)[1,])
abline(v=var(data$L[,1]), col=2)
abline(v=mean(sapply(dL, function(x) var(x[,1]))), col=3)
hist(sapply(dL, apply, 2, var)[2,])
abline(v=var(data$L[,2]), col=2)
abline(v=mean(sapply(dL, function(x) var(x[,2]))), col=3)

## Check S draws

dS <- lapply(1:100, function(x) drawS(data)$S)
max(abs(Reduce("+", dS) / length(dS) - data$S))
max(abs(sapply(dS, colMeans) - 0))
max(abs(sapply(dS, apply, 2, var) - diag(psi)[3:4]))

layout(rbind(1:2))
hist(sapply(dS, apply, 2, var)[1,])
abline(v=var(data$S[,1]), col=2)
abline(v=mean(sapply(dS, function(x) var(x[,1]))), col=3)
hist(sapply(dS, apply, 2, var)[2,])
abline(v=var(data$S[,2]), col=2)
abline(v=mean(sapply(dS, function(x) var(x[,2]))), col=3)

## Check P draws

dP <- lapply(1:100, function(x) drawPrepLatent(data)$P)
max(abs(Reduce("+", dP) / length(dP) - data$P))
max(abs(sapply(dP, colMeans) - 0))
max(abs(sapply(dP, apply, 2, var) - apply(data$P, 2, var)))

layout(rbind(1:4))
hist(sapply(dP, apply, 2, var)[1,])
abline(v=var(data$P[,1]), col=2)
abline(v=mean(sapply(dP, function(x) var(x[,1]))), col=3)
hist(sapply(dP, apply, 2, var)[2,])
abline(v=var(data$P[,2]), col=2)
abline(v=mean(sapply(dP, function(x) var(x[,2]))), col=3)

hist(sapply(dP, apply, 2, var)[3,])
abline(v=var(data$P[,3]), col=2)
abline(v=mean(sapply(dP, function(x) var(x[,3]))), col=3)
hist(sapply(dP, apply, 2, var)[4,])
abline(v=var(data$P[,4]), col=2)
abline(v=mean(sapply(dP, function(x) var(x[,4]))), col=3)

data <- addMasks(data)


new.data <- lapply(1:100, function(x) drawLambda1Eps1Nu1(data, 0, 0))
# Lambda1 test
idx <- which(data$Lambda1.mask==0)
nc <- ceiling(sqrt(length(idx))); nr <- ceiling(length(idx)/nc)
layout(matrix(1:(nr*nc),nrow=nr,ncol=nc))
Lambda1vals <- sapply(new.data,function(x) x$Lambda1[idx])
for(i in 1:length(idx)){
  hist(Lambda1vals[i,],breaks=10)
}
# Eps1
nr <- ceiling(sqrt(nrow(data$Eps1))); nc <- ceiling(nrow(data$Eps1)/nr)
layout(matrix(1:(nr*nc),nrow=nr,ncol=nc))
Eps1.vals <- sapply(new.data,function(x) diag(x$Eps1))
for(i in 1:nrow(data$Eps1)){
  hist(Eps1.vals[i,],breaks=10)
}
# nu1
nr <- ceiling(sqrt(length(data$nu1))); nc <- ceiling(length(data$nu1)/nr)
layout(matrix(1:(nr*nc),nrow=nr,ncol=nc))
nu1.vals <- sapply(new.data,function(x) x$nu1)
for(i in 1:length(data$nu1)){
  hist(nu1.vals[i,],breaks=10)
}




## drawLambda2G2Eps2Nu2 ##
new.data <- lapply(1:100, function(x) drawLambda2G2Eps2Nu2(data, 0, 0))
# Lambda2 test
idx <- which(data$Lambda2.mask==0)
nc <- ceiling(sqrt(length(idx))); nr <- ceiling(length(idx)/nr)
layout(matrix(1:(nr*nc),nrow=nr,ncol=nc))
Lambda2vals <- sapply(new.data,function(x) x$Lambda2[idx])
for(i in 1:length(idx)){
  hist(Lambda2vals[i,],breaks=10)
}
# G2 test
idx <- which(data$G2.mask==0)
nc <- ceiling(sqrt(length(idx))); nr <- ceiling(length(idx)/nc)
layout(matrix(1:(nr*nc),nrow=nr,ncol=nc))
G2vals <- sapply(new.data,function(x) x$G2[idx])
for(i in 1:length(idx)){
  hist(G2vals[i,],breaks=10)
}
# Eps2
nr <- ceiling(sqrt(nrow(data$Eps2))); nc <- ceiling(nrow(data$Eps2)/nr)
layout(matrix(1:(nr*nc),nrow=nr,ncol=nc))
Eps2.vals <- sapply(new.data,function(x) diag(x$Eps2))
for(i in 1:nrow(data$Eps2)){
  hist(Eps2.vals[i,],breaks=10)
}
# Nu2
nr <- ceiling(sqrt(length(data$nu2))); nc <- ceiling(length(data$nu2)/nr)
layout(matrix(1:(nr*nc),nrow=nr,ncol=nc))
nu2.vals <- sapply(new.data,function(x) x$nu2)
for(i in 1:length(data$nu2)){
  hist(nu2.vals[i,],breaks=10)
}

## drawBG1Chi ##
new.data <- lapply(1:100, function(x) drawBG1Chi(data, 0, 0))
# B test
idx <- which(data$B.mask==0)
nc <- ceiling(sqrt(length(idx))); nr <- ceiling(length(idx)/nr)
layout(matrix(1:length(idx),nrow=nr,ncol=nc))
Bvals <- sapply(new.data,function(x) x$B[idx])
for(i in 1:length(idx)){
  hist(Bvals[i,],breaks=10)
}
# G1 test
idx <- which(data$G1.mask==0)
nc <- ceiling(sqrt(length(idx))); nr <- ceiling(length(idx)/nr)
layout(matrix(1:length(idx),nrow=nr,ncol=nc))
G1vals <- sapply(new.data,function(x) x$G1[idx])
for(i in 1:length(idx)){
  hist(G1vals[i,],breaks=10)
}
# Chi
nr <- ceiling(sqrt(nrow(data$Chi))); nc <- ceiling(nrow(data$Chi)/nr)
layout(matrix(1:nrow(data$Chi),nrow=nr,ncol=nc))
Chi.vals <- sapply(new.data,function(x) diag(x$Chi))
for(i in 1:nrow(data$Chi)){
  hist(Chi.vals[i,],breaks=10)
}

##############################
# Test the starting point generator based on the old model

library(SCM)

load("data/thinned-100.Rdata")
oldstate <- nstart ; rm(nstart)
newstate <- BayesCFAStart(oldstate, S.mask=c(FALSE, FALSE))

##############################
# Test the full model on a small data set,
# starting from the truth

library(SCM)

set.seed(42)

preps <- c("L1.P1.M1", "L1.P1.M2", "L1.P1.M3",
           "L1.P2.M1",
           "L1.P3.M1", "L1.P3.M2",
           "L2.P1.M1",
           "L2.P2.M1", "L2.P2.M2", "L2.P2.M3",
           "L2.P3.M1", "L2.P3.M2", "L2.P3.M3",
           "L2.P4.M1")

psi <- matrix(NA, 4, 4)
psi[1,] <- psi[,1] <- c(2  , 2,   1.6, 1.6)
psi[2,] <- psi[,2] <- c(2  , 3  , 2.6, 2.5)
psi[3,] <- psi[,3] <- c(1.6, 2.6, 2.8, 2.7)
psi[4,] <- psi[,4] <- c(1.6, 2.5, 2.7, 2.8)
colnames(psi) <- rownames(psi) <- psiNames(preps)

B <- createUnitB(preps)
Lambda1 <- createUnitLambda1(preps)
Lambda2 <- createUnitLambda2(preps)
G1 <- randScaling(createUnitG1(preps), 1.5)
G2 <- randScaling(createUnitG2(preps), 1.5)

data <- genData(numGenes=5000, preps=preps, psi=psi,
                B=B, Lambda1=Lambda1, Lambda2=Lambda2, G1=G1, G2=G2,
                cChi=.5, cEps1=.1, cEps2=.6)

data <- addMasks(data)
data <- addResTerms(data)

fitted <- BayesCFA(data$cens1, data$cens2, preps, start=list(data),
                   scaling=FALSE, noSteps=50, verbose=TRUE, draw=c(L=TRUE))

for (i in 1:10000) {
  allTracePlots(fitted$samples, truth=data,
                pars=c("cor", "psi", "Chi",
                  "Eps1", "Eps2", "Lmean", "Smean",
                  "Pmean", "Lvar", "Svar", "Pvar",
                  "nu1", "nu2", "Ldot", "Sdot", "Lcor", "Scor"))

  fitted <- BayesCFA(data$cens1, data$cens2, preps,
                start=list(fitted$samples[[length(fitted$samples)]][[1]]),
                scaling=FALSE, noSteps=50, verbose=TRUE, draw=c(L=TRUE))
}

#################################
## Start from a random starting point
## For this we essentially need to generate another data set,
## with different parameters

library(SCM)


set.seed(42)

preps <- c("L1.P1.M1", "L1.P1.M2", "L1.P1.M3",
           "L1.P2.M1",
           "L1.P3.M1", "L1.P3.M2",
           "L2.P1.M1", "L2.P1.M2", "L2.P1.M3",
           "L2.P2.M1",
           "L2.P3.M1", "L2.P3.M2")

psi <- matrix(NA, 4, 4)
psi[1,] <- c(2, 1.8, 0  , 0)
psi[2,] <- c(1.8, 2, 0  , 0)
psi[3,] <- c(0, 0, 2.8, 2.7)
psi[4,] <- c(0, 0, 2.7, 2.8)
colnames(psi) <- rownames(psi) <- psiNames(preps)

B <- createUnitB(preps)
Lambda1 <- createUnitLambda1(preps)
Lambda2 <- createUnitLambda2(preps)
G1 <- randScaling(createUnitG1(preps), 1.5)
G2 <- randScaling(createUnitG2(preps), 1.5)

data <- genData(numGenes=5000, preps=preps, psi=psi,
                B=B, Lambda1=Lambda1, Lambda2=Lambda2, G1=G1, G2=G2,
                cChi=.5, cEps1=.1, cEps2=.6)

data <- addMasks(data)

psi2 <- matrix(NA, 4, 4)
psi2[1,] <- c(2  , 1.2, 0  , 0)
psi2[2,] <- c(1.2, 2  , 0  , 0)
psi2[3,] <- c(0  , 0  , 3  , 2.4)
psi2[4,] <- c(0  , 0  , 2.4, 3)
colnames(psi2) <- rownames(psi2) <- psiNames(preps)

G12 <- randScaling(createUnitG1(preps), 1.5)
G22 <- randScaling(createUnitG2(preps), 1.5)

Chi <- diag(noPrepsMulti(preps))
diag(Chi) <- runif(nrow(Chi))
rownames(Chi) <- colnames(Chi) <- multiPrepNames(preps)

Eps1 <- diag(noMeasMulti(preps))
diag(Eps1) <- runif(nrow(Eps1))
rownames(Eps1) <- colnames(Eps1) <- multiMeasNames(preps)

Eps2 <- diag(noMeasUni(preps))
diag(Eps2) <- runif(nrow(Eps2))
rownames(Eps2) <- colnames(Eps2) <- uniMeasNames(preps)

data2 <- genData(numGenes=5000, preps=preps, psi=psi2,
                 B=B, Lambda1=Lambda1, Lambda2=Lambda2, G1=G1, G2=G2,
                 Chi=Chi, Eps1=Eps1, Eps2=Eps2)

data2 <- addMasks(data2)
data$G1.mask[] <- 1 ; data2$G1.mask[] <- 1
data$G2.mask[] <- 1 ; data2$G2.mask[] <- 1

## data2$S <- data$S
## data2$L <- data$L
## data2$P <- data$P
## data2$Eps1 <- data$Eps1
## data2$nu1 <- data$nu1

data2$data1 <- data$data1
data2$data2 <- data$data2

fitted <- BayesCFA(data$cens1, data$cens2, preps, start=list(data2),
                   scaling=FALSE, noSteps=50, verbose=TRUE)

for (i in 1:10000) {
  allTracePlots(fitted$samples, truth=data,
                pars=c("cor", "psi", "Chi",
                  "Eps1", "Eps2", "Lmean", "Smean",
                  "Pmean", "Lvar", "Svar", "Pvar",
                  "nu1", "nu2", "Ldot", "Sdot", "Lcor", "Scor"))

  fitted <- BayesCFA(data$cens1, data$cens2, preps,
                start=list(fitted$samples[[length(fitted$samples)]][[1]]),
                scaling=FALSE, noSteps=50, verbose=TRUE)

}

#################################################

## Bigger data set, start from truth

library(SCM)

set.seed(42)

preps <- c("L1.P1.M1", "L1.P1.M2", "L1.P1.M3", "L1.P1.M4", "L1.P1.M5",
           "L1.P2.M1", "L1.P2.M2", "L1.P2.M3", "L1.P2.M4",
           "L1.P3.M1", "L1.P3.M2", "L1.P3.M3",
           "L1.P4.M1", "L1.P4.M2", "L1.P4.M3", "L1.P4.M4",
           "L1.P5.M1",
           "L1.P6.M1",
           "L1.P7.M1",
           "L1.P8.M1",
           "L2.P1.M1", "L2.P1.M2", "L2.P1.M3", "L2.P1.M4",
           "L2.P2.M1", "L2.P2.M2", "L2.P2.M3",
           "L2.P3.M1", "L2.P3.M2", "L2.P3.M3",
           "L2.P4.M1", "L2.P4.M2", "L2.P4.M3", "L2.P4.M4",
           "L2.P5.M1",
           "L2.P6.M1",
           "L2.P7.M1")

psi <- matrix(NA, 4, 4)
psi[1,] <- psi[,1] <- c(2, 2, 0,   0)
psi[2,] <- psi[,2] <- c(2, 3, 0,   0)
psi[3,] <- psi[,3] <- c(0, 0, 2.8, 2.7)
psi[4,] <- psi[,4] <- c(0, 0, 2.7, 2.8)
colnames(psi) <- rownames(psi) <- psiNames(preps)

B <- createUnitB(preps)
Lambda1 <- createUnitLambda1(preps)
Lambda2 <- createUnitLambda2(preps)
G1 <- randScaling(createUnitG1(preps), 1.5)
G2 <- randScaling(createUnitG2(preps), 1.5)

data <- genData(numGenes=5000, preps=preps, psi=psi,
                B=B, Lambda1=Lambda1, Lambda2=Lambda2, G1=G1, G2=G2,
                cChi=.5, cEps1=.1, cEps2=.6)

data <- addMasks(data)

fitted <- BayesCFA(data$cens1, data$cens2, preps, start=list(data),
                   scaling=FALSE, noSteps=50, verbose=TRUE)

for (i in 1:10000) {
#  samp <- fitted$samples[[length(fitted$samples)]]
#  save(samp, file=sprintf("bigger-test/thinned-%i.Rdata", i))
#  pdf("bigger-test/plot.pdf", width=20, height=15)
  allTracePlots(fitted$samples, truth=data,
                pars=c("cor", "psi", "Chi", "G1", "G2",
                  "Eps1", "Eps2", "Lmean", "Smean",
                  "Pmean", "Lvar", "Svar", "Pvar",
                  "nu1", "nu2", "Ldot", "Sdot", "Lcor", "Scor", "Pcor"))
#  dev.off()

  fitted <- BayesCFA(data$cens1, data$cens2, preps,
                start=list(fitted$samples[[length(fitted$samples)]][[1]]),
                scaling=FALSE, noSteps=50, verbose=TRUE)
}

## Start from a random starting point

psi2 <- matrix(NA, 4, 4)
psi2[1,] <- c(2  , 1.2, 0  , 0)
psi2[2,] <- c(1.2, 2  , 0  , 0)
psi2[3,] <- c(0  , 0  , 3  , 2.4)
psi2[4,] <- c(0  , 0  , 2.4, 3)
colnames(psi2) <- rownames(psi2) <- psiNames(preps)

G12 <- randScaling(createUnitG1(preps), 1.5)
G22 <- randScaling(createUnitG2(preps), 1.5)
G12[as.logical(data$G1.mask)] <- data$G1[as.logical(data$G1.mask)]
G22[as.logical(data$G2.mask)] <- data$G2[as.logical(data$G2.mask)]

Chi <- diag(noPrepsMulti(preps))
diag(Chi) <- runif(nrow(Chi))
rownames(Chi) <- colnames(Chi) <- multiPrepNames(preps)

Eps1 <- diag(noMeasMulti(preps))
diag(Eps1) <- runif(nrow(Eps1))
rownames(Eps1) <- colnames(Eps1) <- multiMeasNames(preps)

Eps2 <- diag(noMeasUni(preps))
diag(Eps2) <- runif(nrow(Eps2))
rownames(Eps2) <- colnames(Eps2) <- uniMeasNames(preps)

data2 <- genData(numGenes=5000, preps=preps, psi=psi2,
                 B=B, Lambda1=Lambda1, Lambda2=Lambda2, G1=G12, G2=G22,
                 Chi=Chi, Eps1=Eps1, Eps2=Eps2)

data2 <- addMasks(data2)

## data$G1.mask[] <- 1 ; data2$G1.mask[] <- 1
## data$G2.mask[] <- 1 ; data2$G2.mask[] <- 1

data2$data1 <- data$data1
data2$data2 <- data$data2

fitted <- BayesCFA(data$cens1, data$cens2, preps, start=list(data2),
                   scaling=FALSE, noSteps=50, verbose=TRUE)

for (i in 1:10000) {
  samp <- fitted$samples[[length(fitted$samples)]]
  save(samp, file=sprintf("bigger-test2/thinned-%i.Rdata", i))
  pdf("bigger-test2/plot.pdf", width=20, height=15)
  allTracePlots(fitted$samples, truth=data,
                pars=c("cor", "psi", "Chi", "G1", "G2",
                  "Eps1", "Eps2", "Lmean", "Smean",
                  "Pmean", "Lvar", "Svar", "Pvar",
                  "nu1", "nu2", "Ldot", "Sdot", "Lcor", "Scor"))
  dev.off()
  
  fitted <- BayesCFA(data$cens1, data$cens2, preps,
                start=list(fitted$samples[[length(fitted$samples)]][[1]]),
                scaling=FALSE, noSteps=50, verbose=TRUE)

}

#######################################
## A data set with gene-specific noise

library(SCM)

set.seed(42)

preps <- c("L1.P1.M1", "L1.P1.M2", "L1.P1.M3",
           "L1.P2.M1",
           "L1.P3.M1", "L1.P3.M2",
           "L2.P1.M1",
           "L2.P2.M1", "L2.P2.M2", "L2.P2.M3",
           "L2.P3.M1", "L2.P3.M2", "L2.P3.M3",
           "L2.P4.M1")

psi <- matrix(NA, 4, 4)
psi[1,] <- psi[,1] <- c(2  , 2,   1.6, 1.6)
psi[2,] <- psi[,2] <- c(2  , 3  , 2.6, 2.5)
psi[3,] <- psi[,3] <- c(1.6, 2.6, 2.8, 2.7)
psi[4,] <- psi[,4] <- c(1.6, 2.5, 2.7, 2.8)
colnames(psi) <- rownames(psi) <- psiNames(preps)

B <- createUnitB(preps)
Lambda1 <- createUnitLambda1(preps)
Lambda2 <- createUnitLambda2(preps)
G1 <- randScaling(createUnitG1(preps), 1.5)
G2 <- randScaling(createUnitG2(preps), 1.5)

Omega <- runif(5000)
data <- genData(numGenes=5000, preps=preps, psi=psi,
                B=B, Lambda1=Lambda1, Lambda2=Lambda2, G1=G1, G2=G2,
                cChi=.5, cEps1=.1, cEps2=.6, Omega=Omega)

data <- addMasks(data)
data <- addResTerms(data)

fitted <- BayesCFA(data$cens1, data$cens2, preps, start=list(data),
                   scaling=FALSE, noSteps=10, verbose=TRUE,
                   draw=c(Omega=TRUE))
                   
for (i in 1:10000) {
  allTracePlots(fitted$samples, truth=data,
                pars=c("cor", "psi", "Chi",
                  "Eps1", "Eps2", "Lmean", "Smean",
                  "Pmean", "Lvar", "Svar", "Pvar",
                  "nu1", "nu2", "Ldot", "Sdot", "Lcor", "Scor"))

  fitted <- BayesCFA(data$cens1, data$cens2, preps,
                start=list(fitted$samples[[length(fitted$samples)]][[1]]),
                scaling=FALSE, noSteps=50, verbose=TRUE, draw=c(L=TRUE))
}
