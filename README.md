## Accounting for experimental noise reveals that mRNA levels, amplified by post-transcriptional processes, largely determine steady-state protein levels in yeast

This repository contains the manuscript mentioned in the title, and associated code and data sets.
Should you need help running our code, please contact us.

### Citation

Csárdi G, Franks AM, Choi DS, Airoldi EM, Drummond DA, "Accounting for experimental noise reveals that mRNA levels, amplified by post-transcriptional processes, largely determine steady-state protein levels in yeast," *PLoS Genetics* (2015).

### Abstract

Cells respond to their environment by modulating protein levels
through mRNA transcription and post-transcriptional control. Modest observed
correlations between global steady-state mRNA and protein measurements
have been interpreted as evidence that mRNA levels determine
roughly 40% of the variation in protein levels, indicating dominant
post-transcriptional effects. However, the techniques underlying these
conclusions, such as correlation and regression, yield biased results
when data are noisy, missing systematically, and collinear---properties
of mRNA and protein measurements---which motivated us to revisit this
subject. Noise-robust analyses of 24 studies of budding yeast reveal
that mRNA levels explain more than 85% of the variation in steady-state
protein levels. Protein levels are not proportional to mRNA levels, but rise much more rapidly. Regulation of translation suffices to explain this nonlinear effect, revealing post-transcriptional amplification of, rather than competition with, transcriptional signals. These results substantially revise widely credited models of protein-level regulation, and introduce multiple noise-aware approaches essential for proper analysis of many biological phenomena.

### Code

Source for the SCM is in `code`, and source used to generate merged datasets and figures is in `src`. The latter relies on code in [`dad:base.git`](http://github.com/dad/base)

Here, for the impatient, is an implementation of Spearman's correction in [R](http://www.r-project.org).

```
# Log-transform x, treating values <= 0 or infinite as NA
log.nozero <- function(x, log.fxn=base::log, ...) {
	x[x<=0 | x==Inf] <- NA
	log.fxn(x, ...)
}


# Spearman correction for correlations between x and y, with each variable being a matrix or data.frame
# of replicates.
cor.sp.matrix <- function(x, y, method='pearson', use='pairwise.complete.obs', log=FALSE, na.rm=FALSE) {
	d <- data.frame(x,y)
	fnstr.beg <- ""
	fnstr.end <- ""
	if (log) {
		fnstr.beg <- "log.nz("
		fnstr.end <- ")"
	}
	if (log) {
		d <- log.nozero(d)
	}
	if (na.rm) {
		d <- na.omit(d)
	}
	if (nrow(na.omit(d))<3) {
		# Correlations with fewer than 3 points throw errors (and are probably garbage anyway)
		warning("Insufficient data to compute correlations")
	}

	# Dimensions: do the right thing if there's only one measurement of x or y.
	nx <- 1
	if (!is.null(dx <- dim(x))) {
		nx <- dx[2]
	}
	ny <- 1
	if (!is.null(dy <- dim(y))) {
		ny <- dy[2]
	}
	# Calculate the full correlation matrix
	r <- cor(d, method=method, use=use)
	# Extract reliabilities and correlations
	if (nx>1) {
		mrelx <- geom.mean(lt(r[1:nx,1:nx]))
	} else {
		mrelx <- r[1,1] # Just one X value
	}

	if (ny>1) {
		mrely <- geom.mean(lt(r[(nx+1):ncol(d),(nx+1):ncol(d),drop=F]))
	} else {
		mrely <- r[ncol(d),ncol(d)] # Just one Y value
	}
	rr <- r[1:nx,(nx+1):ncol(d)]
	mr <- geom.mean(rr)
	res <- mr/sqrt(mrelx*mrely)
	# Correlation estimate is the Spearman-corrected correlation
	# No confidence interval
	list(r=res, n=nrow(d), r.unc=mr, cor=r, nx=nx, ny=ny, relx=mrelx, rely=mrely, estimate=res, conf.int=c(NA,NA))
}
```

### Data

Published data can be downloaded from [Dryad](http://datadryad.org/resource/doi:10.5061/dryad.d644f). Many datafiles, including raw data, are also available here for convenience.

### License

MIT
