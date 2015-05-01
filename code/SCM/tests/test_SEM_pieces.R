
library(testthat)
library(SCM)

context("SEM draws")

test_that("Drawing Theta and Nu works", {
  set.seed(42)
  options(digits=10)
  data <- genData(NT=2, NE=8, NR=32)
  
  newdata <- drawThetaNu(data)
  expect_that(newdata, is_a("list"))
  expect_that(attributes(newdata$Theta),
              is_identical_to(attributes(data$Theta)))
  expect_that(attributes(newdata$nu),
              is_identical_to(attributes(data$nu)))
  expect_that(getX(newdata), equals(getX(data)))

  ThetaNu <- lapply(1:100, function(x) { drawThetaNu(data)[c("Theta", "nu")] })
  sumTheta <- summary(as.vector(sapply(lapply(ThetaNu, "[[", "Theta"), diag)))
  sumNu <- summary(as.vector(sapply(ThetaNu, "[[", "nu")))

  expect_that(unclass(sumTheta),
              is_equivalent_to(c(0.08020894, 0.09771378, 0.10211410,
                                 0.10204210, 0.10639740, 0.12280670)))
  expect_that(unclass(sumNu),
              is_equivalent_to(c(-0.017235830, -0.005374416, 0.001176128,
                                 0.001077403, 0.008241299, 0.025906420)))
})

######

test_that("Drawing G and Xi works", {
  set.seed(42)
  data <- genData(NT=2, NE=8, NR=32)
  
  newdata <- drawGXi(data)
  expect_that(newdata, is_a("list"))
  expect_that(attributes(newdata$G),
              is_identical_to(attributes(data$G)))
  expect_that(attributes(newdata$Xi),
              is_identical_to(attributes(data$Xi)))
  expect_that(getX(newdata), equals(getX(data)))

  GXi <- lapply(1:100, function(x) { drawGXi(data)[c("G", "Xi")] })

  sumXi <- summary(as.vector(sapply(lapply(GXi, "[[", "Xi"), diag)))
  sumG <- summary(as.vector(sapply(GXi, "[[", "G")))
  
  expect_that(unclass(sumG),
              is_equivalent_to(c(0.9018559, 0.9873730, 1.0054300,
                                 1.0022570, 1.0198940, 1.0692080)))
  expect_that(unclass(sumXi),
              is_equivalent_to(c(0.2558644, 0.2936982, 0.3035309,
                                 0.3037948, 0.3142314, 0.3542335)))
})

#######

test_that("Drawing Tau works", {
  set.seed(42)
  data <- genData(NT=4, NE=16, NR=64, cTau=homoTau(.2))

  newdata <- drawTau(data, poolTau=TRUE)
  expect_that(newdata, is_a("list"))
  expect_that(attributes(newdata$Tau),
              is_identical_to(attributes(data$Tau)))
  expect_that(getX(newdata), equals(getX(data)))
  
  Tau <- lapply(1:100, function(x) { drawTau(data, poolTau=TRUE)["Tau"] })

  sumTau <- summary(as.vector(sapply(lapply(Tau, "[[", "Tau"), diag)))
  
  expect_that(unclass(sumTau),
              is_equivalent_to(c(0.1798716, 0.2038385, 0.2091183,
                                 0.2093361, 0.2151560, 0.2307163)))
})

test_that("Drawing unpooled Tau works", {
  set.seed(42)
  options(digits=7)
  data <- genData(NT=4, NE=16, NR=64, cTau=homoTau(.2))

  newdata <- drawTau(data, poolTau=FALSE)
  expect_that(newdata, is_a("list"))
  expect_that(attributes(newdata$Tau),
              is_identical_to(attributes(data$Tau)))
  expect_that(getX(newdata), equals(getX(data)))

  Tau <- lapply(1:100, function(x) { drawTau(data, poolTau=FALSE)["Tau"] })

  sumTau <- summary(as.vector(sapply(lapply(Tau, "[[", "Tau"), diag)))

  expect_that(unclass(sumTau),
              is_equivalent_to(c(0.1779, 0.2025, 0.2086,
                                 0.2089, 0.2154, 0.2374)))
})

test_that("Drawing Tau runs", {
  data <- genDataSpec(eta0=homoEta(2), eta1=homoEta(1),
                      tList=c(3,2), eList=c(1,1,2,2,2),
                      rList=c(1,2,1,2,2,2,2,2),
                      cTau=homoTau(.2))
  tmp <- drawTau(data, poolTau=TRUE)
  tmp <- drawTau(data, poolTau=FALSE)
  tmp <- drawLTE(data)
  tmp <- drawGTau(data, poolTau=TRUE)
  tmp <- drawGTau(data, poolTau=FALSE)
  tmp <- drawGTauNu(data, poolTau=TRUE)
  tmp <- drawGTauNu(data, poolTau=FALSE)
  tmp <- drawGXi(data)
  tmp <- drawGXiNu(data)
  tmp <- drawThetaNu(data)
})

test_that("zero Tau works (well, runs)", {
  data <- genData(NT=4)
  tmp <- drawTau(data, poolTau=TRUE)
  tmp <- drawTau(data, poolTau=FALSE)
  tmp <- drawLTE(data)
  tmp <- drawGTau(data, poolTau=TRUE)
  tmp <- drawGTau(data, poolTau=FALSE)
  tmp <- drawGTauNu(data, poolTau=TRUE)
  tmp <- drawGTauNu(data, poolTau=FALSE)
  tmp <- drawGXi(data)
  tmp <- drawGXiNu(data)
  tmp <- drawThetaNu(data)
})

########

test_that("Drawing psi works", {
  set.seed(42)
  options(digits=10)
  data <- genData(NT=2, NE=8, NR=32)

  newdata <- drawPsi(data)
  expect_that(newdata, is_a("list"))
  expect_that(attributes(newdata$psi),
              is_identical_to(attributes(data$psi)))
  expect_that(getX(newdata), equals(getX(data)))

  psi <- lapply(1:100, function(x) { drawPsi(data)["psi"] })

  sumpsi <- summary(sapply(psi, function(x) x$psi[1,2]))

  expect_that(unclass(sumpsi),
              is_equivalent_to(c(0.8917299, 0.8980496, 0.9000000,
                                 0.8995953, 0.9008159, 0.9060382)))
})

#########

test_that("Drawing LTER works", {
  set.seed(42)
  options(digits=7)
  data <- genData(NT=4, NE=16, NR=64, cTau=homoTau(.2))

  newdata <- drawLTE(data)
  expect_that(newdata, is_a("list"))
  expect_that(attributes(newdata$L),
              is_identical_to(attributes(data$L)))
  expect_that(attributes(newdata$T),
              is_identical_to(attributes(data$T)))
  expect_that(attributes(newdata$E),
              is_identical_to(attributes(data$E)))
  expect_that(attributes(newdata$R),
              is_identical_to(attributes(data$R)))
  expect_that(getX(newdata), equals(getX(data)))

  LTER <- lapply(1:100, function(x) { drawLTE(data)[c("L", "T", "E", "R")] })

  sumL1 <- summary(as.vector(sapply(LTER, function(x) diag(var(x$L)))))
  sumL2 <- summary(sapply(LTER, function(x) var(x$L)[1,2]))
  sumT  <- summary(as.vector(sapply(LTER, function(x) apply(x$T, 2, var))))
  sumE  <- summary(as.vector(sapply(LTER, function(x) apply(x$E, 2, var))))
  sumR  <- summary(as.vector(sapply(LTER, function(x) apply(x$R, 2, var))))

  expect_that(unclass(sumL1), is_equivalent_to(c(0.8918, 0.9224, 0.9342,
                                                 0.9341, 0.9454, 0.9878)))
  expect_that(unclass(sumL2), is_equivalent_to(c(0.8068, 0.8258, 0.8376,
                                                 0.8370, 0.8453, 0.8837)))
  expect_that(unclass(sumT),  is_equivalent_to(c(0.1778, 0.1954, 0.2024,
                                                 0.2026, 0.2082, 0.2359)))
  expect_that(unclass(sumE),  is_equivalent_to(c(0.2669, 0.2957, 0.3039,
                                                 0.3037, 0.3114, 0.3442)))
  expect_that(unclass(sumR),  is_equivalent_to(c(0.08511, 0.09747, 0.10040,
                                                 0.10040, 0.10330, 0.11630)))
})

#########

test_that("Observation parameters are drawn properly", {
  set.seed(42)
  data <- genData(NT=2, NE=8, NR=32, eta0=homoEta(2), eta1=homoEta(1))

  newdata <- drawObsParam(data)
  expect_that(newdata, is_a("list"))
  expect_that(getX(data), equals(getX(newdata)))

  eta <- lapply(1:10, function(x) { drawObsParam(data)[c("eta0", "eta1")] })

  sumEta0 <- summary(as.vector(sapply(eta, "[[", "eta0")))
  sumEta1 <- summary(as.vector(sapply(eta, "[[", "eta1")))

  expect_that(unclass(sumEta0), is_equivalent_to(c(1.824,  1.947,  1.995,
                                                   1.992,  2.050,  2.186)))
  expect_that(unclass(sumEta1), is_equivalent_to(c(0.8368, 0.9606, 1.0110,
                                                   1.0100, 1.0550, 1.1580)))
})

#########

test_that("Imputation works", {
  set.seed(42)
  options(digits=7)
  data <- genData(NT=2, NE=8, NR=32, eta0=homoEta(2), eta1=homoEta(1))
  newdata <- drawMissing(data)
  expect_that(newdata, is_a("list"))
  expect_that(newdata$X, equals(getX(newdata)))
  expect_that(newdata$noAccImp, equals(5467))
  expect_that(newdata$noRejImp, equals(3))
  expect_that(newdata$X[newdata$I==1], equals(data$X[data$I==1]))

  sumX <- summary(newdata$X[newdata$I==0])
  expect_that(unclass(sumX), is_equivalent_to(c(-4.2970, -1.7470, -1.0190,
                                                -0.9985, -0.2794, 3.5500)))
})

#########

test_that("Imputation works for special cases", {
  set.seed(42)
  options(digits=7)
  data <- genDataSpec(eta0=homoEta(2), eta1=homoEta(1),
                      tList=c(3,2), eList=c(1,1,2,2,2),
                      rList=c(1,2,1,2,2,2,2,2),
                      cTau=homoTau(.2))
  newdata <- drawMissing(data)
  expect_that(newdata, is_a("list"))
  expect_that(newdata$X, equals(getX(newdata)))
  expect_that(newdata$noAccImp, equals(2446))
  expect_that(newdata$noRejImp, equals(4))
  expect_that(newdata$X[newdata$I==1], equals(data$X[data$I==1]))
  
  sumX <- summary(newdata$X[newdata$I==0])
  expect_that(unclass(sumX), is_equivalent_to(c(-4.8250, -1.8150, -1.0570,
                                                -1.0500, -0.3081,  3.3020)))
})

test_that("Drawing G and Tau works (pooled)", {
  set.seed(42)
  options(digits=7)
  data <- genDataSpec(eta0=homoEta(2), eta1=homoEta(1),
                      tList=c(3,2), eList=c(1,1,2,2,2),
                      rList=c(1,2,1,2,2,2,2,2),
                      cTau=homoTau(.2))
  newdata <- drawGTau(data, poolTau=TRUE)
  expect_that(unique(newdata$G[-2]), equals(1.0))
  expect_that(unique(diag(newdata$Tau)[1:3]), equals(0.2147358303839194682))
  expect_that(unique(diag(newdata$Tau)[4:5]), equals(0.2))

  GTau <- lapply(1:100, function(x) { drawGTau(data, poolTau=TRUE)[c("G", "Tau")] })
  sumTau <- summary(as.vector(sapply(lapply(GTau, "[[", "Tau"), diag)))
  sumG <- summary(as.vector(sapply(GTau, "[[", "G")))
  expect_that(unclass(sumTau), is_equivalent_to(c(0.1832, 0.1996, 0.2000,
                                                  0.2015, 0.2034, 0.2270)))
  expect_that(unclass(sumG), is_equivalent_to(c(0.9852, 1.0000, 1.0000,
                                                1.0020, 1.0000, 1.0460)))
})

test_that("Drawing G and Tau works (unpooled)", {
  set.seed(42)
  options(digits=7)
  data <- genDataSpec(eta0=homoEta(2), eta1=homoEta(1),
                      tList=c(3,2), eList=c(1,1,2,2,2),
                      rList=c(1,2,1,2,2,2,2,2),
                      cTau=homoTau(.2))
  newdata <- drawGTau(data, poolTau=FALSE)
  expect_that(newdata$G[2], equals(c(E2=0.9896782806781062591384)))
  expect_that(unique(newdata$G[-2]), equals(1.0))
  expect_that(unique(diag(newdata$Tau)[2]), equals(0.228727879664697))
  expect_that(unique(diag(newdata$Tau)[-2]), equals(0.2))

  GTau <- lapply(1:100, function(x) { drawGTau(data, poolTau=FALSE)[c("G", "Tau")] })
  sumTau <- summary(as.vector(sapply(lapply(GTau, "[[", "Tau"), diag)))
  sumG <- summary(as.vector(sapply(GTau, "[[", "G")))
  expect_that(unclass(sumTau), is_equivalent_to(c(0.1952, 0.2000, 0.2000,
                                                  0.2031, 0.2000, 0.2418)))
  expect_that(unclass(sumG), is_equivalent_to(c(0.9844, 1.0000, 1.0000,
                                                1.0020, 1.0000, 1.0470)))
})

test_that("Drawing G, Tau and Nu works", {
  set.seed(42)
  options(digits=7)
  data <- genDataSpec(eta0=homoEta(2), eta1=homoEta(1),
                      tList=c(3,2), eList=c(1,1,2,2,2),
                      rList=c(1,2,1,2,2,2,2,2),
                      cTau=homoTau(.2))
  newdata <- drawGTauNu(data, poolTau=TRUE)
  expect_that(newdata$G[1], equals(c(E1=0.980993907336250026141)))
  expect_that(unique(newdata$G[-1]), equals(1))
  expect_that(unique(diag(newdata$Tau)[1:3]), equals(0.2147404466740395201363))
  expect_that(unique(diag(newdata$Tau)[4:5]), equals(0.2))
  expect_that(newdata$nu[1], equals(c(L1.T1.E1.R1=-0.02263175080508259343071)))
  expect_that(unique(newdata$nu[-1]), equals(0))

  newdata <- drawGTauNu(data, poolTau=FALSE)
  expect_that(newdata$G[1], equals(c(E1=0.9884619429827676828637)))
  expect_that(unique(newdata$G[-1]), equals(1))
  expect_that(diag(newdata$Tau)[1], equals(c(T1=0.1974965133207789991232)))
  expect_that(unique(diag(newdata$Tau)[2:5]), equals(0.2))
  expect_that(newdata$nu[1], equals(c(L1.T1.E1.R1=-0.006901617501362249974817)))
  expect_that(unique(newdata$nu[-1]), equals(0))

  GTauNu <- lapply(1:100, function(x) { drawGTauNu(data, poolTau=FALSE)[c("G", "Tau", "nu")] })
  sumTau <- summary(as.vector(sapply(lapply(GTauNu, "[[", "Tau"), diag)))
  sumG <- summary(as.vector(sapply(GTauNu, "[[", "G")))
  sumNu <- summary(as.vector(sapply(GTauNu, "[[", "nu")))
  expect_that(unclass(sumTau), is_equivalent_to(c(0.1856, 0.2000, 0.2000,
                                                  0.2018, 0.2000, 0.2310)))
  expect_that(unclass(sumG), is_equivalent_to(c(0.9589, 1.0000, 1.0000,
                                                0.9986, 1.0000, 1.0240)))
  expect_that(unclass(sumNu), is_equivalent_to(c(-3.149e-02, 0.000e+00, 0.000e+00,
                                                 -1.127e-05, 0.000e+00, 3.236e-02)))

  GTauNu <- lapply(1:100, function(x) { drawGTauNu(data, poolTau=TRUE)[c("G", "Tau", "nu")] })
  sumTau <- summary(as.vector(sapply(lapply(GTauNu, "[[", "Tau"), diag)))
  sumG <- summary(as.vector(sapply(GTauNu, "[[", "G")))
  sumNu <- summary(as.vector(sapply(GTauNu, "[[", "nu")))
  expect_that(unclass(sumTau), is_equivalent_to(c(0.1788, 0.2000, 0.2000,
                                                  0.2016, 0.2037, 0.2245)))
  expect_that(unclass(sumG), is_equivalent_to(c(0.9585, 1.0000, 1.0000,
                                                0.9986, 1.0000, 1.0160)))
  expect_that(unclass(sumNu), is_equivalent_to(c(-2.969e-02, 0.000e+00, 0.000e+00,
                                                 -4.169e-05, 0.000e+00, 3.269e-02)))
})
