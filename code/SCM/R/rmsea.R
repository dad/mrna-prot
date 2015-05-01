
fittedvar <- function(sample) {
  cbind(rbind(V.X1(sample), V.X2.X1(sample)),
        rbind(t(V.X2.X1(sample)), V.X2(sample)))
}

impvar <- function(sample) {
  cbind(rbind(cov(sample$data1), cov(sample$data2, sample$data1)),
        rbind(cov(sample$data1, sample$data2), cov(sample$data2)))
}

chisq <- function(state) {
  cov1 <- impvar(state)
  cov2 <- fittedvar(state)
  N <- nrow(nstart$data1)-1
  (N-1) * (sum(diag(cov1 %*% solve(cov2))) -
           nrow(cov1) + log(det(cov2)) - log(det(cov1)))
}

rmsea <- function(state) {
  cov1 <- impvar(state)
  cov2 <- fittedvar(state)
  N <- nrow(nstart$data1)-1
  X2 <- (N-1) * (sum(diag(cov1 %*% solve(cov2))) -
                 nrow(cov1) + log(det(cov2)) - log(det(cov1)))
  df <- ncol(cov1)
  G <- 1
  GG <- 0
  sqrt(max(c((X2/N)/df - 1/(N - GG), 0))) * sqrt(G)
}
