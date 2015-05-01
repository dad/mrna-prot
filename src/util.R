
## Averages over sets of data arising from the same study (manuscript, ms), as indicated
## by the ms.flds vector.
ms.mean <- function(x, ms.flds, flds=NULL, na.rm=TRUE, mean.fxn=mean) {
  if (is.null(flds)) {
    flds <- names(x)
  }
  stopifnot(length(ms.flds) == length(flds))
  if (is.null(ncol(x[flds]))) {
    res <- sapply(unique(ms.flds), function(ms) {
      mean.fxn(x[ms.flds==ms], na.rm=na.rm)
    })
  }
  else {
    y <- as.data.frame(x[,flds])
    msf <- unique(ms.flds)
    res <- as.data.frame(sapply(msf, function(ms) {
      f <- ms.flds==ms
      xf <- y[,f]
      if (is.null(ncol(xf))) {
        # No data to average -- return raw data
        resy <- xf
      } else {
        resy <- as.data.frame(apply(xf, 1, mean.fxn, na.rm=na.rm))
      }
      resy
    }))
    names(res) <- msf
  }
  res[is.na(res)] <- NA # Turns NaN into NA
  res
}

## Standard error on the mean
std.err <- function(x,na.rm=T) {
  y <- x
  if (na.rm) {
    y <- na.omit(x)
  }
  sd(y)/sqrt(length(y))
}

##############
## Test cases
##############
test.ms.mean1 <- function() {
  ms <- as.character(c(2,2,2,4,4,4))
  x <- c(1,2,3,4,5,6)
  names(x) <- ms
  xm <- ms.mean(x, ms)
  res <- xm['2']==2 & xm['4'] == 5
}

test.ms.mean2 <- function() {
  xa <- data.frame(x1=1:6, x2=1:6, y1=7:12, y2=7:12)
  xam <- ms.mean(xa, c('x','x','y','y'))
  res <- all(xam[,'x'] == 1:6 & xam[,'y']==7:12)
}

test.ms.mean3 <- function() {
  xa <- data.frame(x1=c(1:5,NA), x2=c(1:5,NA), y1=7:12, y2=7:12)
  xam <- ms.mean(xa, c('x','x','y','y'), na.rm=T)
  res <- all(xam[1:5,'x'] == 1:5) & is.na(xam$x[6]) & all(xam[,'y']==7:12)
}

test.ms.mean4 <- function() {
  # Single observation for one variable
  xa <- data.frame(x1=c(1:5,NA), x2=c(1:5,NA), y=7:12)
  xam <- ms.mean(xa, c('x','x','y'), na.rm=T)
  res <- all(xam[1:5,'x'] == 1:5) & is.na(xam$x[6]) & all(xam[,'y']==7:12)
}

util.tests <- c(test.ms.mean1,
                test.ms.mean2,
                test.ms.mean3,
                test.ms.mean4)

util.runtest = F
if (util.runtest) {
  test.res <- sapply(util.tests, function(x) {x()})
  if (any(!test.res)) {
    cat("Tests failed\n")
  } else { 
    cat("Tests passed\n")
  }
}
