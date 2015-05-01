
## Create in-silico prep names
getPrepsDef <- function(numLatent, numPrepsPer, numMeasPer) {
  lat <- rep(paste("L", sep="", 1:numLatent), each=numPrepsPer*numMeasPer)
  prep <- rep(paste("P", sep="", 1:numPrepsPer), each=numMeasPer, numLatent)
  meas <- rep(paste("M", sep="", 1:numMeasPer), numLatent*numPrepsPer)
  paste(sep=".", lat, prep, meas)
}

## Create preps to techs map
getTechsDef <- function(preps, noTechsPerLatent){

  prep <- getPrep(preps)
  names(prep) <- preps
  uniquePrep <- unique(prep)
  
  ## There should be at least 3 observations for each technology
  stopifnot(all(3*noTechsPerLatent <= summary(as.factor((getLatentFromPrep(uniquePrep))))))
  
  techsDef <- rep("",length(preps))
  names(techsDef) <- preps
  
  numTechs <- noTechsPerLatent*noLatent(preps)
  techs <- toupper(letters[1:numTechs])
  latents <- getLatent(preps)
  for(i in seq_along(unique(latents))){
    l <- unique(latents)[i]
    ltechs <- techs[(noTechsPerLatent*(i-1)+1) : (noTechsPerLatent*i) ]

    prepL <- unique(prep[latents==l])
    
    techVec <- rep(ltechs,length(prepL))
    count <- 1
    for(p in sample(prepL)){
      techsDef[names(prep[prep==p])] <- paste(l,".",techVec[count],sep="")
      count <- count+1
    }
  }
  techsDef
}

## Check preps, we check that
chkpreps <- function(preps) {
  ## 1. they are ordered
  stopifnot(!is.unsorted(preps))
  ## 2. unique
  stopifnot(!any(duplicated(preps)))
  ## 3. have the right pattern: latent dot prep dot meas
  stopifnot(all(grepl("^[^\\.]+\\.[^\\.]+\\.[^\\.]+$", preps)))
}

## techs is a named list, sorted
## names is "preps"
## values are the corresponding technology
chktechs <- function(techs){
  ## 1. they are ordered by prep
  stopifnot(!is.unsorted(names(techs)))
  ## 2. tech names (the preps) are not duplicated
  stopifnot(!any(duplicated(names(techs))))
  ## 3. have the right pattern: latent dot tech
  ## and names are latent.prep.meas
  stopifnot(all(grepl("^[^\\.]+\\.[^\\.]+$", techs)))
  stopifnot(all(grepl("^[^\\.]+\\.[^\\.]+\\.[^\\.]+$", names(techs))))
}

chkprepnames <- function(prepnames) {
  ## 1. they are ordered
  stopifnot(!is.unsorted(prepnames))
  ## 2. unique
  stopifnot(!any(duplicated(prepnames)))
  ## 3. have the right pattern: latent dot prep 
  stopifnot(all(grepl("^[^\\.]+\\.[^\\.]+$", prepnames)))
}

## Table of contents of functions to query preps
## 
## 1 number of things, a single non-negative integer:
## --------------------------------------------------
##   noMeas - number of measurements
##   noMeasUni - number of measurements from non-identifyable preps
##   noMeasMulti - number of measurements from identifyable preps
##   noLatent - number of latent variables
##   noPreps - number of _all_ prep variables, including the ones
##             that cannot be identified
##   noPrepsUni - number of prep variables that cannot be
##                identified
##   noPrepsMulti - number of prep variables that can be identified
##    
## 2 number of things, grouped by other things
## -------------------------------------------
##   noPrepsPer - number of _all_ preps per latent variables
##   noMeasPerPer - number of measurements per preps per latent
##                  variables. This is a list of named integer
##                  vectors. Both the list and the vectors are
##                  ordered alphabetically. It contains _all_ preps.
##   noMeasPer - number of measurements per preps. It contains _all_
##               preps, in a named integer vector, ordered alphabetically.
## 
## 3 query properties of the measuments
## ------------------------------------
##   getPrep - which preps the different measurements belong to.
##             a character vector, contains _all_ measurements.
##   getLatent - which latent variables the measurements belong to.
##             a character vector, contains _all_ measurements.
##   getPrepUni - like getPrep, but only for measurements from
##                non-identifyable preps
##   getPrepMulti - like getPrep, but only for measurements from
##                  identifyable preps
##   getLatentUni - like getLatent, but only for measurements from
##                  non-identifyable preps
##   getLatentMulti - like getLatent, but only for measurements from
##                  identifyable preps
##   getLLatentFromPrep - which latent variables the given preps belong to
##   getLSatentFromPrep - which latent variables the given preps belong to
## 
## 4 query names of preps and latent variables
## -------------------------------------------
##   latentNames - names of latent variables, sorted alphabetically
##   prepNames - _all_ prep names, sorted alphabetically.
##   uniPrepNames - names of the non-identifiable preps
##   multiPrepNames - names of the identifyable preps
##   uniMeasNames - names of the measurements that belong to
##                  non-identifyable preps
##   multiMeasNames - names of the measurement that belong to
##                  identifiable preps

## Total number of measurements
noMeas <- function(preps) {
  chkpreps(preps)
  length(preps)
}

## Number of measurements that belong to non-identifyable preps
noMeasUni <- function(preps) {
  chkpreps(preps)
  nos <- noMeasPer(preps)
  sel <- names(nos)[nos==1]
  pr <- getPrep(preps)
  sum(pr %in% sel)
}

## Number of measuements that belong to identifyable preps
noMeasMulti <- function(preps) {
  chkpreps(preps)
  nos <- noMeasPer(preps)
  sel <- names(nos)[nos!=1]
  pr <- getPrep(preps)
  sum(pr %in% sel)  
}

## Number of latent variables
noLatent <- function(preps) {
  chkpreps(preps)
  length(unique(sapply(strsplit(preps, ".", fixed=TRUE), "[", 1)))
}

## Number of latent variables
noLatentFromTechs <- function(techs) {
  chkpreps(techs)
  length(unique(getLatentFromTechs(techs)))
}


## Number of technologies
noTechs <- function(techs){
  chktechs(techs)
  length(unique(techs))
}

## Number of _all_ preps
noPreps <- function(preps) {
  chkpreps(preps)
  pr <- sub("\\.[^\\.]+$", "", preps)
  length(unique(pr))
}

## Number or real preps, i.e. preps with more than one measurement
noPrepsUni <- function(preps) {
  chkpreps(preps)
  nom <- noMeasPer(preps)
  sum(nom==1)
}

## Number or real preps, i.e. preps with more than one measurement
noPrepsMulti <- function(preps) {
  chkpreps(preps)
  nom <- noMeasPer(preps)
  sum(nom!=1)
}

## Number of preps per latent variable
noPrepsPer <- function(preps) {
  chkpreps(preps)
  lp <- strsplit(sub("\\.[^\\.]*$", "", preps), ".", fixed=TRUE)
  tapply(sapply(lp, "[", 2), sapply(lp, "[", 1),
         function(x) length(unique(x)))
}

## Number of measurements per preps per latent variables
noMeasPerPer <- function(preps) {
  chkpreps(preps)
  sp <- strsplit(preps, "\\.")
  tapply(lapply(sp, "[", 2:3), sapply(sp, "[", 1), function(lp) {
    tapply(sapply(lp, "[", 2), sapply(lp, "[", 1), function(x)
           length(unique(x)))
  }, simplify=FALSE)
}

## Number of measurements per preps
noMeasPer <- function(preps) {
  chkpreps(preps)
  pr <- sub("\\.[^\\.]+$", "", preps)
  tab <- table(pr)
  structure(as.vector(tab), names=dimnames(tab)[[1]])
}

## Prep ids of the measurements
getPrep <- function(preps) {
  chkpreps(preps)
  sub("\\.[^\\.]*$", "", preps)
}

## Latent variables of the measurements
getLatent <- function(preps) {
  chkpreps(preps)
  sub("\\..*$", "", preps)
}

## Prep ids of the measuments from non-identifyable preps
getPrepUni <- function(preps) {
  chkpreps(preps)
  up <- uniPrepNames(preps)
  pr <- getPrep(preps)
  pr[ pr %in% up ]
}

## Prep ids of the measurements from identifyable preps
getPrepMulti <- function(preps) {
  chkpreps(preps)
  mp <- multiPrepNames(preps)
  pr <- getPrep(preps)
  pr[ pr %in% mp ]
}

## Latent variables for measurements from non-identifyable preps
getLatentUni <- function(preps) {
  chkpreps(preps)
  up <- uniPrepNames(preps)
  pr <- getPrep(preps)
  la <- getLatent(preps)
  la[ pr %in% up ]
}

## Latent variables for mesurements from identifyable preps
getLatentMulti <- function(preps) {
  chkpreps(preps)
  mp <- multiPrepNames(preps)
  pr <- getPrep(preps)
  la <- getLatent(preps)
  la[ pr %in% mp ]
}

## All prep names
prepNames <- function(preps) {
  chkpreps(preps)
  sort(unique(getPrep(preps)))
}

## Names of non-identifyable preps
uniPrepNames <- function(preps) {
  chkpreps(preps)
  no <- noMeasPer(preps)
  names(no)[no==1]
}

## Names of identifyable preps
multiPrepNames <- function(preps) {
  chkpreps(preps)
  no <- noMeasPer(preps)
  names(no)[no!=1]
}

## Which latent variable for the given prep names
getLatentFromPrep <- function(prepnames) {
  chkprepnames(prepnames)
  sub("\\..*$", "", prepnames)
}

getLatentFromTechs <- function(techs) {
  chktechs(techs)
  sub("\\..*$", "", techs)
}


## Names of latent variables
latentNames <- function(preps) {
  chkpreps(preps)
  unique(getLatent(preps))
}

## Names of the measuments from non-identifyable preps
uniMeasNames <- function(preps) {
  chkpreps(preps)
  up <- uniPrepNames(preps)
  preps[ getPrep(preps) %in% up ]
} 

## Names of the measuments from identifyable preps
multiMeasNames <- function(preps) {
  chkpreps(preps)
  mp <- multiPrepNames(preps)
  preps[ getPrep(preps) %in% mp ]
} 

## Names of the rows and columns of the psi matrix
psiNames <- function(preps) {
  latentNames <- latentNames(preps)
  c(paste("L.", sep="", latentNames), paste("S.", sep="", latentNames))
}
