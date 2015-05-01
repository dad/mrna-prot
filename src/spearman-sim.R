
gen.data <- function(n, target.corr, target.rel, sub.n=n, censoring=c("random", "bottom", "logit"), error.corr.x=0.0, error.corr.y=0.0) {
	# Generate n rows of data having two measurments of each variable x and y, with target true correlation between x and y, and target reliability of each variable.
	censoring <- match.arg(censoring)
	sdx <- 1.0
	x <- rnorm(n, mean=0, sd=sdx)
	# What is the standard deviation of zero-mean noise necessary to make the true correlation between x and y equal to target.corr?
	corr.sd.target <- sdx*sqrt(1/target.corr^2 - 1)
	y <- x + rnorm(n,sd=corr.sd.target)
	sdy <- sqrt(sdx^2 + corr.sd.target^2)
	# What is the standard deviation of zero-mean noise necessary to make observations of x have target.rel reliability?
	rel.sd.target.x <- sdx*sqrt(1/target.rel - 1)
	rel.sd.target.y <- sdy*sqrt(1/target.rel - 1)
	# What is the noise 
	# Must have total noise variance in X = rel.sd.target.x^2 = err.corr.sd.target.x^2 + err.uncorr.sd.target.x^2
	corr.noise.x <- 0
	err.uncorr.sd.target.x <- rel.sd.target.x
	corr.noise.y <- 0
	err.uncorr.sd.target.y <- rel.sd.target.y
	if (error.corr.x > 0.0) {
		err.corr.sd.target.x <- rel.sd.target.x*sqrt(1/error.corr.x - 1)
		err.uncorr.sd.target.x <- sqrt(rel.sd.target.x^2 - err.corr.sd.target.x^2)
		corr.noise.x <- rnorm(n, sd=err.corr.sd.target.x)
	}
	if (error.corr.y > 0.0) {
		err.corr.sd.target.y <- rel.sd.target.y*sqrt(1/error.corr.y - 1)
		err.uncorr.sd.target.y <- sqrt(rel.sd.target.y^2 - err.corr.sd.target.y^2)
		corr.noise.y <- rnorm(n, sd=err.corr.sd.target.y)
	}

	#cat(corr.noise.x, err.uncorr.sd.target.x, '\n')
	# Generate the data
	d <- data.frame(x=x, y=y, 
		x1 = x + corr.noise.x + rnorm(n,sd=err.uncorr.sd.target.x), 
		x2 = x + corr.noise.x + rnorm(n,sd=err.uncorr.sd.target.x), 
		y1 = y + corr.noise.y + rnorm(n,sd=err.uncorr.sd.target.y), 
		y2 = y + corr.noise.y + rnorm(n,sd=err.uncorr.sd.target.y))

	#cat('about to optimize','n',err.uncorr.sd.target.x, err.uncorr.sd.target.y,'\n')
		
	# Censor the data
	if (sub.n != n) {
		#cat(sub.n, n, '\n')
		if (censoring == "bottom") {
			d[order(d$x)[(sub.n+1):n], c("x", "x1", "x2")] <- NA
		} else if (censoring == "random") {
			#print('about to optimize')
			d[(sub.n+1):n, c("x", "x1", "x2")] <- NA
		} else if (censoring == "logit") {
			noc <- n - sub.n
			invlogit <- function(x, phi, rho) {
					1/(1+exp(-(phi+rho*x)))
			}
			f <- function(par) {
					phi <- par[1]
					rho <- par[2]
					(sum(invlogit(d[,"x"], phi, rho)) - noc) ^ 2
			}
			res <- optim(c(-2,-2), f, method="L-BFGS-B",
						 lower=c(-20,-5), upper=c(20,0),
						 control=list(factr=1e15))
			cp <- invlogit(d[,"x"], res$par[1], res$par[2])
			d[runif(n) < cp, c("x", "x1", "x2")] <- NA
		}
	}
	d
}

gen.corr <- function(x, method='pearson') {
	# Based on data in x, assumed to have two measurements each of x and y, generate observed correlation, reliabilities, and corrected correlation
	#rsp <- rcorr.sp(x[,c('x1','x2')], x[,c('y1','y2')], method=method)
	rsp <- cor.sp(x=x[,c('x1','x2')], y=x[,c('y1','y2')], method=method)
	c(true=cor(x$x,x$y,method=method,use="complete.obs"), obs=cor(x$x1,x$y1,method=method,use="complete.obs"),
          est=rsp$estimate, rel=sqrt(rsp$rxx*rsp$ryy), relx=rsp$rxx, rely=rsp$ryy)
#          est=rsp$r, rel=rsp$rel, relx=rsp$relx, rely=rsp$rely)
}
gen.reps <- function(m, n, target.corr, target.rel, method='pearson', sub.n=n, censoring="random") {
	# Utility to generate many replicates of obs/est correlations, discarding the underlying data
	x <- as.data.frame(rbind(t(replicate(m, gen.corr(gen.data(n, target.corr, target.rel, sub.n, censoring), method=method)))))
	x
}

n <- 5000
m <- 50
meth='pearson'

repnums <- 10
bias.xlim <- c(100,n)
bias.sub.n <- seq(log2(bias.xlim[1]), log2(bias.xlim[2]), length.out=repnums)
bias.corr.target <- 0.9
# The true reliability target
rel.target <- 0.7

# Do the latent model fits if they are not available
#if (!file.exists("data/spearman-sim-model.Rdata")) {
#        stop("spearman simulation model fits not available")
#}
load("../data/spearman-sim-model2.Rdata")

target <- seq(0.1,1,length.out=repnums)
cols <- rep("darkorange",repnums) #rainbow(repnums)
corr.cols <- rep("blue",repnums) #rainbow(repnums)
dosq <- 1
tp <- 0.6
cex.pts <- 0.75
set.seed(111)
#if (F){
varycor <- list()
for (reps in 1:repnums) {
	varycor[[as.character(target[reps])]] <- gen.reps(m, n, target.corr=target[reps], target.rel=rel.target, method=meth)
}
cat("done varycor\n")
corr.target <- 0.9
varyrel <- list()
for (reps in 1:repnums) {
	varyrel[[as.character(target[reps])]] <- gen.reps(m, n, corr.target, target.rel=target[reps], method=meth)
}
cat("done varyrel\n")
corr.target <- 0.9
mar <- list()
for (l.n in seq(log2(bias.xlim[1]), log2(bias.xlim[2]), length.out=repnums)) {
	sub.n <- 2^l.n
	mar[[as.character(sub.n)]] <- gen.reps(m, n, target.corr=corr.target, target.rel=rel.target, method=meth, sub.n=sub.n, censoring="random")
}
cat("done mar\n")
#}
if (T) {
	bias.corr.target <- 0.9
	nmar <- list()
	for (l.n in bias.sub.n) {
		sub.n <- 2^l.n
		nmar[[as.character(sub.n)]] <- gen.reps(m, n, bias.corr.target, rel.target, method=meth, sub.n=sub.n, censoring="logit")
	}
}
cat("done nmar\n")
corEst.all <- list(varycor=varycor, varyrel=varyrel, mar=mar, nmar=nmar, nmar.inf=corEst)
fname <- "~/research/spearman/data/spearman-model-results.Rdata"
save(corEst.all, file=fname)
cat("# Wrote simulation results to", fname, '\n')
