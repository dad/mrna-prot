source("util.R")
raw.data.dirname <- "../data/raw"
bg.data.dirname <- "../data/bg"
project.data.dirname <- "../data"
master.fname <- paste(project.data.dirname,"scer-raw-extended.Rdata",sep='/')

assemble.background.data = T
assemble.fit.data = T
compute.spearman.corrections = T
normalize.to.molecules.per.cell = T
assemble.translational.efficiency.data = T
save.data = T
package.data = T

prot.files <- read.delim(paste(project.data.dirname,"/","scer-prot-files.txt",sep=''), comment.char='#', stringsAsFactors=FALSE)
mrna.files <- read.delim(paste(project.data.dirname,"/","scer-mrna-files.txt",sep=''), comment.char='#', stringsAsFactors=FALSE)

exp.flds <- c("prot","mrna")
exp.vars <- list(prot=prot.files, mrna=mrna.files)

load.var <- function(x, master.key, data.cache=NULL) {
	#print(x)
	x <- as.list(x)
	dat <- data.cache[[x$fname]]
	if (is.null(dat)) {
		cat("Loading",x$fname,'\n')
		#print(x$log)
		# Read the data
		dat <- read.table(paste(raw.data.dirname,x$fname,sep='/'), header=TRUE, sep='\t', comment.char="#", quote="\"", fill=FALSE, stringsAsFactors=FALSE)
		data.cache[[x$fname]] <- dat
	} else {
		cat("Loaded",x$fname,'from cache\n')
	}
	# Match the key
	z <- match(master.key, dat[,x$orfname])
	
	# Determine if any transformation must be done to the data
	fn <- noop
	if (!is.na(x$log)) {
		if (x$log == 'log') {
			fn <- exp
		}
		else if (x$log == 'log2') {
			fn <- function(x) {2^x}
		}
	}
	# Select the data and transform it
	d <- fn(dat[z,x$valname])
	#cat(x$valname, x$fname, length(d), '\n')
	return(d)
}

assemble.data <- function(file.var.list, master.key, label) {
	data.cache <- list()
	fs <- apply(file.var.list, 1, load.var, master.key=master.key, data.cache=data.cache)
	colnames(fs) <- p.0(label,file.var.list$alias)
	return(fs)
}

if (assemble.background.data) {
	cat("# Loading background S.cer data...\n")
	options("stringsAsFactors"=FALSE)
	# Use this for updated SGD R64 genome release
	#y.feat <- read.delim(paste(bg.data.dirname,"/scerevisiae-coding-features.txt",sep=''), comment.char='#', quote=NULL)
	# Use this for the older SGD release
	y.feat <- read.delim(paste(bg.data.dirname,"/scer-features.txt",sep=''), comment.char='#', quote=NULL)
	# Paralogs
	para <- read.delim("../data/bg/ygob-fungal-pillars-v7.txt", comment.char='#')
	z1 <- match(y.feat$orf, para$scer1)
	z2 <- match(y.feat$orf, para$scer2)
	y.feat$pair.orf <- para[z1,'scer2']
	y.feat$pair.orf[!is.na(para[z2,'scer1'])] <- para[z2,'scer1'][!is.na(para[z2,'scer1'])]
	# Gene names
	sgd <- read.delim(paste(bg.data.dirname,"/SGD_features.tab",sep=''), comment.char='#', quote=NULL, header=F)
	z <- match(y.feat$orf, sgd$V4)
	y.feat$gene <- sgd[z,'V5']  # common gene name
	y.feat$sgd.type <- sgd[z,'V3'] # dubious, uncharacterized, verified

	bg <- y.feat

	# Order names for easy viewing
    ordered.names <- unique(c("orf","gene","pair.orf",names(bg)))
    bg <- bg[,ordered.names]

	
	cat("# Loaded", nrow(bg), "rows of background S.cer data\n")
	
	raw <- NULL
	for (variable in exp.flds) {
		d <- assemble.data(exp.vars[[variable]], bg$orf, label=variable)
		if (!is.null(raw)) {
			raw <- data.frame(d, raw)
		} else {
			raw <- data.frame(d)
		}
	}
	
	var.flds <- sapply(exp.flds, function(m){paste(m, exp.vars[[m]]$alias, sep='.')})
	var.fld.descs <- list(orf="Open reading frame ID", gene="Gene name")
	ms.descs <- list()
	for (x in exp.flds) {
		xv <- exp.vars[[x]]
		descs <- as.list(xv$desc)
		names(descs) <- p.0(x, xv$alias)
		var.fld.descs <- c(var.fld.descs, descs)
		# Isolate unique manuscripts and store ms descriptions
		xvu <- xv[match(unique(xv$ms),xv$ms),]
		nd <- xvu$ms.desc
		names(nd) <- p.0(x,xvu$ms)
		ms.descs <- c(ms.descs, nd)
	}
}
if (T) {
	yres <- list(bg=bg, descs=descs, ms.descs=ms.descs)
}

if (assemble.fit.data) {
	cat("# Loading fit results...\n")
	# Load "corr"
	#load("../data/nonpooled-tech.RData")
	#summary.byorf <- read.table("../data/nonpooled-tech-sample-summary.txt", header=T, sep='\t')
	#summary.data <- get(load("../data/nonpooled-tech-sample-corrs.RData"))
	load("../data/g2-nonpooled-tech-sample.RData")
	# Loads L1, L2
	load("../data/g2-nonpooled-tech-sample-L.Rdata")

	summary.byorf <- read.table("../data/g2-nonpooled-tech-sample-summary.txt", header=T, sep=' ')
	summary.data <- get(load("../data/g2-nonpooled-tech-sample-corrs.RData"))
	# Load column mapping
	cn <- read.delim("../data/column-mapping.txt")

	model <- samp[c('X','L','N','NL','NT','NE','NR','lj','kj','tj','G','G.mat','nu','psi','Tau','Xi','Theta','eta0','eta1','T','E','I','TechWeights','Theta.draw','Xi.draw','B.mat','Lambda.mat')]

	# Raw data
	res <- list()
	# Imputed data
	res$X.raw <- data.frame(model$X)[,cn$fit] # imputed data
	colnames(res$X.raw) <- cn$orig
	
	
	# Is data imputed?
	res$I <- data.frame(model$I)[,cn$fit]==0
	colnames(res$I) <- cn$orig
	# Normalize the data using inferred means
	res$X <- scale.mat(res$X.raw, scale=F, center=model$nu[cn$fit])
	# Normalized data, omitting imputed data (censored)
	res$Xc <- res$X
	res$Xc[res$I] <- NA
	# Best estimate
	res$L <- data.frame(samp$L[,c('abund','expr')])
	colnames(res$L) <- c('prot','mrna')
	#colnames(raw.norm) <- cn$orig
	
	# Fields and identifiers
	fields <- list()
	fields$all <- cn$orig
	fields$prot <- as.character(subset(cn$orig, substring(cn$orig,1,4)=='prot'))
	fields$mrna <- as.character(subset(cn$orig, substring(cn$orig,1,4)=='mrna'))
	# Reorder experimental variables to put in same order as fields.
	exp.vars.ordered <- lapply(names(exp.vars), function(x) {
		v <- exp.vars[[x]]
		v[match(fields[[x]], paste(x,v$alias,sep='.')),]
		})
	names(exp.vars.ordered) <- names(exp.vars)
	# Per-manuscript fields
	ms.fields=list(mrna=paste('mrna',unique(exp.vars.ordered$mrna$ms),sep='.'), prot=paste('prot',unique(exp.vars.ordered$prot$ms),sep='.'))

	# Standardize all the data to have the same ORF order
	zbg <- match(yres$bg$orf, rownames(samp$X))

	target.sample <- 6
	L.sample <- data.frame(mrna=L1[,target.sample], prot=L2[,target.sample])

	# Error/reliability information
	# err is the noise variance/signal variance 
	e <- err[cn$fit]
	names(e) <- cn$orig
	# Convert to proportion of variance due to noise -- err = 1 - reliability.
	res$err <- e/(1+e)
	e.mean.ms <- c(ms.mean(e, p.0('mrna',exp.vars.ordered$mrna$ms), fields$mrna, na.rm=T, mean.fxn=geo.mean),
		ms.mean(e, p.0('prot',exp.vars.ordered$prot$ms), fields$prot, na.rm=T, mean.fxn=geo.mean))
	res$err.ms <- e.mean.ms/(1+e.mean.ms)

	raw.norm <- log.nozero(raw)
	raw.norm[,fields$mrna] <- scale(raw.norm[,fields$mrna], center=T, scale=T)
	raw.norm[,fields$prot] <- scale(raw.norm[,fields$prot], center=T, scale=T)
  
	# Averages
	raw.ms.avg <- data.frame(
		mrna=ms.mean(raw.norm[fields$mrna], exp.vars.ordered$mrna$ms), 
		prot=ms.mean(raw.norm[fields$prot], exp.vars.ordered$prot$ms))

	raw.avg <- data.frame(
		mrna=rowMeans(raw.ms.avg[ms.fields$mrna], na.rm=T),
		prot=rowMeans(raw.ms.avg[ms.fields$prot], na.rm=T),
		sd.mrna=apply(raw.ms.avg[ms.fields$mrna], 1, std.err, na.rm=T), 
		sd.prot=apply(raw.ms.avg[ms.fields$prot], 1, std.err, na.rm=T),
		n.mrna=apply(raw.ms.avg[ms.fields$mrna], 1, na.len),
		n.prot=apply(raw.ms.avg[ms.fields$prot], 1, na.len))
	
	X.avg <- data.frame(
		mrna=rowMeans(mav <- ms.mean(res$X[zbg,], exp.vars.ordered$mrna$ms, fields$mrna, na.rm=T), na.rm=T),
		prot=rowMeans(pav <- ms.mean(res$X[zbg,], exp.vars.ordered$prot$ms, fields$prot, na.rm=T), na.rm=T),
		sd.mrna=apply(mav, 1, std.err, na.rm=T), 
		sd.prot=apply(pav, 1, std.err, na.rm=T))
	
	X.ms.avg <- data.frame(mrna=mav, prot=pav)
	
	Xc.avg <- data.frame(
		mrna=rowMeans(mav <- ms.mean(res$Xc[zbg,], exp.vars.ordered$mrna$ms, fields$mrna, na.rm=T), na.rm=T),
		prot=rowMeans(pav <- ms.mean(res$Xc[zbg,], exp.vars.ordered$prot$ms, fields$prot, na.rm=T), na.rm=T),
		sd.mrna=apply(mav, 1, std.err, na.rm=T), 
		sd.prot=apply(pav, 1, std.err, na.rm=T),
		n.mrna=apply(mav, 1, na.len),
		n.prot=apply(pav, 1, na.len))

	Xc.ms.avg <- data.frame(mrna=mav, prot=pav)

    # Expanding imputation variable (I) means introducing NA's; set these to TRUE (imputed)
    res$I.expanded <- res$I[zbg,]
    res$I.expanded[is.na(res$I.expanded)] <- TRUE
    rownames(res$I.expanded) <- NULL
	
	head <- yres$bg[,c('gene','orf','pair.orf')]
	
	sums <- summary.byorf[zbg,c('mean.mrna','mean.prot','sd.mrna','sd.prot')]
	colnames(sums) <- c('mrna','prot','sd.mrna','sd.prot')
	est <- cbind(head, sums, Xc.avg[,c('n.mrna','n.prot')])
	est.sample <- cbind(head, L.sample[zbg,c('mrna','prot')], Xc.avg[,c('n.mrna','n.prot')])
}

if (T) {
	# Compile master object from these results
	yres.dat <- list(
		fields=fields,
		ms.fields=ms.fields,
		exp.vars=exp.vars.ordered,
		n.genes=samp$N,
		corrs=summary.data[[1]],
		raw=cbind(head,raw),
		raw.avg=cbind(head,raw.avg),
		raw.ms.avg=cbind(head,raw.ms.avg),
		X.raw=cbind(head,res$X.raw[zbg,]),
		X=cbind(head,res$X[zbg,]),
		Xc=cbind(head,res$Xc[zbg,]),
        I=cbind(head,res$I.expanded),
		L=cbind(head,res$L[zbg,]),
		est.sample.unscaled=est.sample,
		est.unscaled=est,
		X.avg=cbind(head,X.avg),
		Xc.avg=cbind(head,Xc.avg),
		X.ms.avg=cbind(head,X.ms.avg),
		Xc.ms.avg=cbind(head,Xc.ms.avg),
		exp.vars=exp.vars.ordered,
		err=res$err[c(fields$mrna,fields$prot)],
		err.ms=res$err.ms[c(ms.fields$mrna,ms.fields$prot)],
		model=model
		)
	yres <- c(yres, yres.dat)
	
	# Special-case some data
	ribosome.orf <- subset(yres$bg, substring(gene,1,3) %in% c('RPS','RPL','RPP') & !gene=='RPP1')$orf
	histone.gene <- c("HTA1","HTA2","HTB1","HTB2","HHT1","HHT2","HHF1","HHF2")
	histone.orf <- yres$bg[match(histone.gene,yres$bg$gene),'orf']
	# Glycolytic enzymes (core): 
	#  Hexokinase (HXK1,HXK2)
	#  Phosphoglucose isomerase (PGI1)
	#  Phosphofructokinase (PFK1,PFK2)
	#  Fructose bisphosphate aldolase (FBA1)
	#  Triosephosphate isomerase (TPI1)
	#  Glyceraldehyde phosphate dehydrogenase (TDH3,TDH2,TDH1)
	#  Phosphoglycerate kinase (PGK1)
	#  Phosphoglycerate mutase (GPM1 [GPM2,GPM3 cannot complement gpm1 and are believed non-functional])
	#  Enolase (ENO1,ENO2)
	#  Pyruvate kinase (CDC19,PYK2)
	glycolytic.gene <- c('HXK1','HXK2','PGI1','PFK1','PFK2','FBA1','TPI1','TDH3','TDH2','TDH1','PGK1','GPM1','ENO1','ENO2','CDC19','PYK2')
	glycolytic.orf <- yres$bg[match(glycolytic.gene,yres$bg$gene),'orf']

	poortag.orf <- c(ribosome.orf, histone.orf)
	not.poortag.orf <- yres$bg$orf[is.na(match(yres$bg$orf, poortag.orf))]
	orfs.list <- list(ribosome=ribosome.orf, histone=histone.orf, glycolytic=glycolytic.orf, poortag=poortag.orf, other=not.poortag.orf)
	subsets <- lapply(orfs.list, function(m){match(m, yres$bg$orf)})
	yres$subset <- subsets
}

rsp.comb <- function(x, vx, vy, method, min.n) {
	# Compute pairs of Spearman-corrected correlations between variables in v1 and v2.
	# Formula: r_corrected = [( r(x1,y1)*r(x2,y2)*r(x1,y2)*r(x2,y1) )^0.25] / sqrt(r(x1,x2)*r(y1,y2))
	# Complete cases only; enforce n >= min.n
	retvals <- c('r','n','lower','upper','runc','rxx','ryy')
	corrs <- apply(combn(vx,2), 2, function(xi){
		res <- apply(combn(vy,2), 2, function(yi) {
			d <- na.omit(data.frame(x[,c(xi,yi)]))
			colnames(d) <- c('x1','x2','y1','y2')
			r <- NA
			mr <- NA
			rxx <- NA
			ryy <- NA
			lower <- NA
			upper <- NA
			
			n <- nrow(d)
			if (n>=min.n) {
              #print(paste(n,xi,yi))
				rsp <- cor.sp(d$x1, d$x2, d$y1, d$y2, method=method)
				r <- rsp$estimate
				mr <- rsp$r.uncorrected
				rxx <- rsp$raw['x1','x2']
				ryy <- rsp$raw['y1','y2']
				lower <- rsp$range[1]
				upper <- rsp$range[2]
			}
			c(r=r, n=n, lower=lower, upper=upper, runc=mr, rxx=rxx, ryy=ryy)
		})
		res
	})
	full.res <- matrix(corrs, nrow=length(retvals))
	full.res <- na.omit(as.data.frame(t(full.res)))
	colnames(full.res) <- retvals
	full.res
}

if (compute.spearman.corrections) {
	# Compute Spearman correction
	cat("# Computing Spearman-corrected correlations...\n")
	# pairwise Spearman-corrected values
	# Exclude mrna.causton and prot.nagaraj, as they are reliability outliers due to overlap with mrna.holstege and prot.degodoy
	rsp.pairs.p <- rsp.comb(raw.ms.avg, ms.fields$mrna[-c(1)], ms.fields$prot[-c(11)], method='p', min.n=100)
	rsp.pairs.s <- rsp.comb(raw.ms.avg, ms.fields$mrna[-c(1)], ms.fields$prot[-c(11)], method='s', min.n=100)
	cat("# Done computing Spearman-corrected correlations.\n")
}
if (T) {
	# Add in the spearman-corrected results
	yres$rsp.pearson <- rsp.pairs.p
	yres$rsp.spearman <- rsp.pairs.s
}

if (normalize.to.molecules.per.cell) {
	cat("# Normalizing to molecules per cell...\n")
	# Scaling factors from model
	prot.names <- names(yres$model$G)[1:11]
	mrna.names <- names(yres$model$G)[12:25]
	G.prot <- yres$model$G[prot.names]
	G.mrna <- yres$model$G[mrna.names]

    # Normalization functions
    # Normalize protein to mass per cell
    norm.prot <- function(x.in, total.protein.mass, mw) {
    	e.prot <- exp(x.in)
		# Molecular weight is in grams/mole; we want molecules/cell
		# molec./cell = (grams/cell) * (molec./mole) / (grams/mole)
		scaling.const <- total.protein.mass * 6.022e23 / sum(e.prot * mw,na.rm=T)
		as.vector(e.prot * scaling.const)
    }
    # Normalize mRNA to molecules per cell
    norm.mrna <- function(x.in, total.mrna) {
		e.mrna <- exp(x.in)
		scaling.const <- total.mrna / sum(e.mrna, na.rm=T)
		as.vector(e.mrna * scaling.const)   	
    }

	# Normalize
	# 36,139 mRNAs/cell, Miura et al. BMC Genomics 9(574), 2008; http://www.biomedcentral.com/1471-2164/9/574
	number.of.mrnas <- 36139
	# 12.5 ug of protein in 4.8 x 10^6 cells during outgrowth.
	# 4 ug of protein in 1.5 x 10^6 cells, Johnston, G. C., F. R. Pringle, and L. H. Hartwell. 1977. Coordination of 
	# growth with cell division in the yeast S. cerevisiae. Exp. Cell Res. 105:79–98. 
	protein.mass <- (4/1.5e6)/1e6

	# Set up variables
	est <- yres$est.unscaled
	rescaled.mrna <- est$mrna*median(G.mrna)
	rescaled.prot <- est$prot*median(G.prot)
	est$mrna <- norm.mrna(rescaled.mrna, number.of.mrnas)
	est$prot <- norm.prot(rescaled.prot, protein.mass, yres$bg$mw)
	
	# Bounds: scale the upper and lower bounds just like the means
	scale.factor <- est$mrna/exp(yres$est.unscaled$mrna)
	est$mrna.upper.sd <- exp(yres$est.unscaled$mrna + yres$est.unscaled$sd.mrna)*scale.factor
	est$mrna.upper.95 <- exp(yres$est.unscaled$mrna + 1.96*yres$est.unscaled$sd.mrna)*scale.factor
	est$mrna.lower.sd <- exp(yres$est.unscaled$mrna - yres$est.unscaled$sd.mrna)*scale.factor
	est$mrna.lower.95 <- exp(yres$est.unscaled$mrna - 1.96*yres$est.unscaled$sd.mrna)*scale.factor
	
	scale.factor <- est$prot/exp(yres$est.unscaled$prot)
	est$prot.upper.sd <- exp(yres$est.unscaled$prot + yres$est.unscaled$sd.prot)*scale.factor
	est$prot.upper.95 <- exp(yres$est.unscaled$prot + 1.96*yres$est.unscaled$sd.prot)*scale.factor
	est$prot.lower.sd <- exp(yres$est.unscaled$prot - yres$est.unscaled$sd.prot)*scale.factor
	est$prot.lower.95 <- exp(yres$est.unscaled$prot - 1.96*yres$est.unscaled$sd.prot)*scale.factor

	# Sample
	samp <- yres$est.sample.unscaled
	rescaled.mrna.sample <- samp$mrna*median(G.mrna)
	rescaled.prot.sample <- samp$prot*median(G.prot)
	samp$mrna <- norm.mrna(rescaled.mrna.sample, number.of.mrnas)
	samp$prot <- norm.prot(rescaled.prot.sample, protein.mass, yres$bg$mw)
}
if (T) {
	# Add in the normalized results
	yres$est <- est
	yres$est.sample <- samp
}

if (assemble.translational.efficiency.data) {
	cat("# Loading translational efficiency/ribosome footprinting data\n")
	ing1.raw <- read.delim("../data/raw/ingolia-footprinting/ingolia09-fp-rich1.txt", comment.char='#')
	z <- match(yres$bg$orf, ing1.raw$orf)
	ing1 <- ing1.raw[z,]
	ing2.raw <- read.delim("../data/raw/ingolia-footprinting/ingolia09-fp-rich2.txt", comment.char='#')
	z <- match(yres$bg$orf, ing2.raw$orf)
	ing2 <- ing2.raw[z,]
	ing1m.raw <- read.delim("../data/raw/ingolia-footprinting/ingolia09-mrna-rich1.txt", comment.char='#')
	z <- match(yres$bg$orf, ing1m.raw$orf)
	ing1m <- ing1m.raw[z,]
	ing2m.raw <- read.delim("../data/raw/ingolia-footprinting/ingolia09-mrna-rich2.txt", comment.char='#')
	z <- match(yres$bg$orf, ing2m.raw$orf)
	ing2m <- ing2m.raw[z,]
	ger12.raw <- read.delim("../data/raw/gerashchenko12-riboseq/sd02.txt", comment.char='#')
	z <- match(yres$bg$orf, ger12.raw$orf)
	ger12 <- ger12.raw[z,]
	ger.raw <- read.delim("../data/raw/gerashchenko14-riboseq.txt", comment.char='#')
	z <- match(yres$bg$orf, ger.raw$orf)
	ger <- ger.raw[z,]
	sub.raw <- read.delim("../data/raw/subtelny14-riboseq.txt", comment.char='#')
	z <- match(yres$bg$orf, sub.raw$orf)
	sub <- sub.raw[z,]
	mcm.raw <- read.delim("../data/raw/mcmanus14-riboseq.txt", comment.char='#')
	mcm.raw$orf <- mcm.raw$Gene
	z <- match(yres$bg$orf, mcm.raw$orf)
	mcm <- mcm.raw[z,]
	rds <- data.frame(
		rd.ing1=(ing1$norm+ing2$norm)/2,
		rd.ger1=(ger12$rdens.init.rpkm1+ger12$rdens.init.rpkm2)/2,
		rd.ger2=(ger$unstressed_noCHX+ger$unstressed_1x_CHX)/2,
		rd.mcm1=(mcm$Scer_Mix_RPF_1+mcm$Scer_Mix_RPF_2)/2,
		rd.sub1=sub$rpf.rpkm
	)
	rd.desc <- list(
		'rd.ing1'="\\cite{ingolia09}", 
		'rd.ger1'="\\cite{gerashchenko12}", 
		'rd.ger2'="\\cite{Gerashchenko2014}", 
		'rd.mcm1'="\\cite{McManus2014}", 
		'rd.sub1'="\\cite{Subtelny2014}")

	# Mean
	rds.sc <- apply(rds, 2, function(x) {exp(scale(log.nz(x), center=T, scale=F))})
	rd.median <- apply(rds.sc, 1, median, na.rm=T)
	rd.mean <- apply(rds.sc, 1, mean, na.rm=T)
	rd.n <- apply(rds.sc, 1, na.len)

	mrna.desc <- list(
		ing1="\\cite{ingolia09}", 
		ger1="\\cite{gerashchenko12}", 
		mcm1="\\cite{McManus2014}", 
		sub1="\\cite{Subtelny2014}",
		est="SCM")
	mrnas <- data.frame(
		ing1=(ing1m$norm+ing2m$norm)/2, 
		ger1=ger12$mrna.init.rpkm1, 
		mcm1=(mcm$Scer_Mix_mRNA_1 + mcm$Scer_Mix_mRNA_2)/2,
		sub1=sub$mrna.rpkm,
		est=yres$est$mrna)
	names(mrnas) <- names(mrna.desc)
	# Goal is to make a fully independent mRNA dataset using only recent RNA-seq measurements
	mrna.sc <- apply(mrnas[1:(ncol(mrnas)-1)], 2, function(x) {exp(scale(log.nz(x), center=T, scale=F))})
	# Exclude Ingolia 2009 and SCM estimate from the averages
	mrna.sc.new <- mrna.sc[,2:ncol(mrna.sc)]
	mrna.median <- apply(mrna.sc.new, 1, median, na.rm=T)
	mrna.mean <- apply(mrna.sc.new, 1, mean, na.rm=T)
	mrna.n <- apply(mrna.sc.new, 1, na.len)

	te <- list(data=data.frame(mrna.sc, rds.sc, yres$est, mrna.median, mrna.mean, mrna.n, rd.median, rd.mean, rd.n, te=rd.median/mrna.median),
		te.mrna.names=colnames(mrna.sc.new),
		te.rd.names=colnames(rds.sc))

	yres$te <- te
}

if (save.data) {
	save(yres, file=master.fname)
	cat("# Saved data to ",master.fname,'\n')
}


if (package.data) {
	# Load header descriptions
	headers <- read.table("../data/column-descriptions.txt", header=T, sep='\t')

	write.datafile <- function(fname, description, dat) {
		write(paste("#",description,"\n#"), file=fname)
		write(paste("# Generated",date(),"\n#"), file=fname, append=TRUE)
		write(paste("# Citation:\n#   \"Accounting for experimental noise reveals that transcript levels,\n",
			"#   amplified by post-transcriptional regulation, largely determine steady-state protein levels in yeast,\"\n",
			"#   Gábor Csárdi, Alexander Franks, David S. Choi, Eduardo M. Airoldi, and D. Allan Drummond, PLoS Genetics, 2015","\n#",sep=''), 
			file=fname, append=TRUE)
		#flds <- colnames(dat)
		#headers[match(flds,headers$column.name),'description'
		header.list <- apply(headers[match(colnames(dat),headers$column.name),], 1, paste, collapse=' = ')
		write(paste('#',header.list, collapse='\n'), file=fname, append=TRUE)
		write.table(format(dat, trim=TRUE, digits=3, scientific=F), file=fname, row.names=F, sep='\t', quote=F, append=TRUE)
		cat("# Wrote data to", fname, '\n')
	}
	# Suppress warnings about appending data
	old.warning.level <- getOption("warn")
	options(warn=-1)
	# Raw data
	write.datafile(fname="../data/package/scer-mrna-protein-raw.txt", 
		description="Raw mRNA and protein measurements for haploid budding yeast (S.cerevisiae) growing exponentially in rich medium with glucose",
		dat=yres$raw[,c('orf','gene',yres$fields$mrna, yres$fields$prot)])
	# Normalized values
	write.datafile(fname="../data/package/scer-mrna-protein-normalized.txt", 
		description="SCM-normalized mRNA and protein measurements for haploid budding yeast (S.cerevisiae) growing exponentially in rich medium with glucose",
		dat=yres$Xc[,c('orf','gene',yres$fields$mrna, yres$fields$prot)])
	# Normalized and imputed values
	write.datafile(fname="../data/package/scer-mrna-protein-normalized-imputed.txt", 
		description="SCM-normalized mRNA and protein measurements, with imputed values, for haploid budding yeast (S.cerevisiae) growing exponentially in rich medium with glucose",
		dat=yres$X[,c('orf','gene',yres$fields$mrna, yres$fields$prot)])
	# Unscaled estimates
	write.datafile(fname="../data/package/scer-mrna-protein-unscaled-estimate.txt", 
		description="SCM estimates of mRNA and protein levels, log-transformed and unscaled, for haploid budding yeast (S.cerevisiae) growing exponentially in rich medium with glucose",
		dat=yres$est.unscaled[,c('orf','gene','mrna','prot','sd.mrna','sd.prot','n.mrna','n.prot')])
	# Scaled estimates
	write.datafile(fname="../data/package/scer-mrna-protein-absolute-estimate.txt", 
		description="SCM estimates of absolute mRNA and protein levels, scaled to (average) molecules per haploid cell, for haploid budding yeast (S.cerevisiae) growing exponentially in rich medium with glucose",
		dat=yres$est[,c('orf','gene','mrna','prot','sd.mrna','sd.prot','n.mrna','n.prot')])
	# Scaled estimates -- sample
	write.datafile(fname="../data/package/scer-mrna-protein-absolute-estimate-sample.txt", 
		description="mRNA and protein levels, scaled to (average) molecules per haploid cell, sampled according to mean and variance of SCM estimates, for haploid budding yeast (S.cerevisiae) growing exponentially in rich medium with glucose",
		dat=yres$est.sample[,c('orf','gene','mrna','prot','n.mrna','n.prot')])
	# Translational efficiency
	write.datafile(fname="../data/package/scer-translational-efficiency.txt", 
		description="Ribosome density, mRNA level, and translational efficiency measurements for haploid budding yeast (S.cerevisiae) growing exponentially in rich medium with glucose",
		dat=yres$te$data[,c('orf','gene',c('ing1',yres$te$te.mrna.names), yres$te$te.rd.names, 'mrna.median','mrna.n','rd.median','rd.n','te')])
	options(warn=old.warning.level)
}
