
context("genData and genDataSpec")

test_that("homogenous variances work", {
  library(SCM)

  D <- genData(Theta=homoTheta(.5), Xi=homoXi(.2), Tau=homoTau(.1))
  expect_that(D$Theta, is_equivalent_to(.5 * diag(D$NR)))
  expect_that(D$Xi, is_equivalent_to(.2 * diag(D$NE)))
  expect_that(D$Tau, is_equivalent_to(.1 * diag(D$NT)))

  D2 <- genData(Theta=homoTheta(0), Xi=homoXi(0), Tau=homoTau(0))
  expect_that(D2$Theta, is_equivalent_to(0 * diag(D2$NR)))
  expect_that(D2$Xi,    is_equivalent_to(0 * diag(D2$NE)))
  expect_that(D2$Tau,   is_equivalent_to(0 * diag(D2$NT)))
})

test_that("normally distributed variances", {
  library(SCM)
  set.seed(42)

  D <- genData(Theta=normTheta(.5,.1), Xi=normXi(.3,.1),
               Tau=normTau(.1,.1))
  expect_that(unclass(summary(as.vector(D$Theta))),
              is_equivalent_to(c(0.00000, 0.00000, 0.00000,
                                 0.01594, 0.00000, 0.72870)))
  expect_that(unclass(summary(as.vector(D$Xi))),
              is_equivalent_to(c(0.00000, 0.00000, 0.00000,
                                 0.0300, 0.0000, 0.4035)))
  expect_that(unclass(summary(as.vector(D$Tau))),
              is_equivalent_to(c(0.00000, 0.00000, 0.03195,
                                 0.04612, 0.07807, 0.12060)))
  
})

test_that("genDataSpec works", {
  library(SCM)
  set.seed(42)

  D <- genDataSpec()
  tList <- eval(formals(genDataSpec)$tList)
  eList <- eval(formals(genDataSpec)$eList)
  rList <- eval(formals(genDataSpec)$rList)
  
  expect_that(D$NT, equals(sum(tList)))
  expect_that(D$NE, equals(sum(eList)))
  expect_that(D$NR, equals(sum(rList)))

  expect_that(unique(dim(D$Theta)), equals(D$NR))
  expect_that(unique(dim(D$Xi)), equals(D$NE))
  expect_that(unique(dim(D$Tau)), equals(D$NT))
})
