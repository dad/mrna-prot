rm(list=ls())
unloadNamespace("SCM")
library(SCM)

## Test 1 (homogeneous variances)
data <- genData(NT=4, cTau=homoTau(.2))

start <- data
samps <- BayesCFA(start)

allTracePlots(samps,truth=start,layout=matrix(1:16,4)) 

## Test 2 (inhomogeneous variances)
rm(list=ls())
unloadNamespace("SCM")
library(SCM)
set.seed(42)
NT <- 4
NE <- NT*4
NR <- NE*4
Theta <- diag(abs(rnorm(NR,mean=.2,sd=.2)))
Xi <- diag(abs(rnorm(NE,mean=1,sd=.5)))
G <- rnorm(NE,mean=1,sd=.5)
data <- genData(N=1000,NT=NT,NE=NE,NR=NR, Theta=Theta,Xi=Xi,G=G, cTau=0.5)
samps <- BayesCFA(data,noStep=100)
allTracePlots(samps,truth=data,layout=matrix(1:12,4))

## Test 3 (inhomogeneous variances, random start)
start <- data
start$psi <- matrix(c(1,.8,.8,1),2)
start$Tau <- .2*solve(start$TechWeights)
start$Xi=ndiag(start$NE, unique(as.character(start$kj)), diag=.3)
start$Theta=ndiag(start$NR, start$Rnames, diag=.1)

start$G <- rep(1,start$NE)
names(start$G) <- unique(start$kj)

start$L <- rmvnorm(start$N,sigma=start$psi)
colnames(start$L) <- unique(start$lj)

start$T <- matrix(0,nrow=start$N,ncol=start$NT)
colnames(start$T) <- unique(start$tj)

start$E <- matrix(0,nrow=start$N,ncol=start$NE)
colnames(start$E) <- unique(start$kj)
start$R <- start$X-(start$L[,start$lj[match(unique(start$kj),start$kj)]]%*%diag(start$G))%*%t(start$Lambda.mat)

colnames(start$R) <- colnames(start$X)

samps <- BayesCFA(start,noStep=200)

allTracePlots(samps,truth=data,layout=matrix(1:12,4))


#### Test 4 (some single replicate experiments)
rm(list=ls())
unloadNamespace("SCM")
library(SCM)
set.seed(42)
data <- genDataSpec(eList=c(2, 2, 2, 2, 2), rList=c(1,rep(3, 9) ))
start <- data
samps <- BayesCFA(noSteps=100,start)

allTracePlots(samps,truth=start,layout=matrix(1:16,4)) 


#### Test 5 (a multi rep, single experiment, tech)
rm(list=ls())
unloadNamespace("SCM")
library(SCM)

set.seed(41)

data <- genDataSpec()
start <- data
samps <- BayesCFA(noSteps=100,start)
allTracePlots(samps,truth=start,layout=matrix(1:16,4)) 

#### Test 6 (complete collpase-- single rep, single tech)
rm(list=ls())
unloadNamespace("SCM")
library(SCM)
set.seed(41)
tList <- c(3,2)
eList <- c(1, 3, 2, 2, 4)
rList <- c(1,rep(3, sum(eList)-1))
data <- genDataSpec(tList=tList,eList=eList,rList=rList,cTau=homoTau(.5),TechWeights=diag(1:5), nu=rnorm(sum(rList),0,1))
start <- data
samps <- BayesCFA(noSteps=100,start)
allTracePlots(samps,truth=start,layout=matrix(1:16,4)) 

#### Test 7 nu
rm(list=ls())
unloadNamespace("SCM")
library(SCM)
set.seed(41)

tList <- c(2,2)
eList <- c(2, 2, 2, 2)
rList <- rep(3, sum(eList))
data <- genDataSpec(tList=tList,eList=eList,rList=rList,cTau=homoTau(.5),TechWeights=diag(1:4), nu=rnorm(sum(rList),0,1))
start <- data
samps <- BayesCFA(noSteps=100,start)
allTracePlots(samps,truth=start,layout=matrix(1:16,4)) 


#### Test 8 all collapse types and non-zero nu
rm(list=ls())
unloadNamespace("SCM")
library(SCM)
set.seed(41)
tList <- c(3,2)
eList <- c(1, 3, 2, 1, 4)
rList <- c(1,rep(3, sum(eList)-2),1)
data <- genDataSpec(tList=tList,eList=eList,rList=rList,cTau=homoTau(.5),TechWeights=diag(1:5), nu=rnorm(sum(rList),0,1))
start <- data
samps <- BayesCFA(noSteps=100,start)
allTracePlots(samps,truth=start,layout=matrix(1:16,4)) 


#### Test 9, missing data
rm(list=ls())
unloadNamespace("SCM")
library(SCM)
set.seed(42 * 42)
data <- genData(eta0=homoEta(2), eta1=homoEta(1))
start <- data
samps <- BayesCFA(noSteps=100,start)
allTracePlots(samps,truth=start,layout=matrix(1:16,4)) 

#### Test 10, missing data + special cases
rm(list=ls())
unloadNamespace("SCM")
library(SCM)
set.seed(42)
data <- genDataSpec(eta0=homoEta(2), eta1=homoEta(1),
                    tList=c(3,2), eList=c(1,1,2,2,2),
                    rList=c(1,2,1,2,2,2,2,2),
                    cTau=homoTau(.2))
start <- data
samps <- BayesCFA(noSteps=100, start)
allTracePlots(samps,truth=start,layout=matrix(1:16,4)) 

#### Test 11, missing data + special cases + no pooling Tau
rm(list=ls())
unloadNamespace("SCM")
library(SCM)
set.seed(42)
data <- genDataSpec(eta0=homoEta(2), eta1=homoEta(1),
                    tList=c(3,2), eList=c(1,1,2,2,2),
                    rList=c(1,2,1,2,2,2,2,2),
                    cTau=normTau(.2, .1))
start <- data
samps <- BayesCFA(noSteps=100, start, poolTau=FALSE)
allTracePlots(samps,truth=start,layout=matrix(1:16,4))
