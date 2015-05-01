source("util.R")
library(Hmisc)
library(lmodel2)

output.type = 'svg'
file.out = F

load.data = F
# Precalculated data
calc.exponent = F

# Figure 1
fig.reliability = F

# Figure 2
fig.correlations.between.scatter = F

# Figure 3
fig.spearman.sim = F

# Figure 4
fig.averaging.increases.corr = F

# Figure 5
fig.mrna.prot.scatter = F
fig.scm.abs.standard.comparison = F
fig.undetected.prot = F
fig.paxdb.vs.scm = F
fig.mrna.prot.overview.distribution = F

# Figure 6
fig.exponent = F
fig.relative.dynamic.range = F
fig.rdens.mrna.slope.vs.cutoff = F
fig.sma.ma.rma.slopes = F

# Figure 7
fig.toy.model = F

fig.translational.efficiency = F

raw.data.dirname <- "../data/raw"
project.data.dirname <- "../data"
master.fname <- paste(project.data.dirname,"scer-raw-extended.Rdata",sep='/')

if (load.data) {
	yres <- get(load(master.fname))
    cat("# Loaded", nrow(yres$raw.avg), "rows of S.cerevisiae data\n")
}

if (fig.spearman.sim) {
	corEst <- get(load(paste(project.data.dirname,"sim/spearman-sim-scm-results.Rdata",sep='/')))

	# Number of genes
	n <- 5000
	# Number of replicates per level
	m <- 50
	# Number of levels
	repnums <- 10
	meth='pearson'

	bias.xlim <- c(100,n)
	bias.sub.n <- seq(log2(bias.xlim[1]), log2(bias.xlim[2]), length.out=repnums)
	bias.corr.target <- 0.8

	dosq <- 1
	sqstr = ""
	if (dosq == 2) {
		sqstr = "squared "
	}
	tp <- 0.6
	cex.pts <- 0.75
	pch <- 16
	cols <- myrainbow(3)
	obscol <- tcol(cols[1],0.4)
	estcol <- tcol(cols[2],0.4)
	fitcol <- tcol(cols[3],0.4)

	if (file.out) dev.out('fig.spearman.sim', width=6, height=6, output.type=output.type)
	split.screen(c(2,2))
	screen(1)
	# Panel A: Vary true correlation at fixed reliability
	par(mar=c(5,5,1,0))
	plot(c(0,1), c(0,1), type='n', xlim=c(0,1), ylim=c(0,1.1), las=1, xlab=paste('True ',sqstr,'correlation',sep=''), ylab=paste('Observed/corrected\n',sqstr,' correlation',sep=''))
	abline(0,1,lty='dotted')
	#abline(h=1,lty='dashed')
	#abline(h=0,lty='dashed')
	points(rep(as.numeric(names(corEst$varycor)), sapply(corEst$varycor, nrow)),
		sapply(corEst$varycor, function(x){x$obs^dosq}), pch=pch, col=obscol, cex=cex.pts)
	points(rep(as.numeric(names(corEst$varycor)), sapply(corEst$varycor, nrow)),
		sapply(corEst$varycor, function(x){x$est^dosq}), pch=pch, col=estcol, cex=cex.pts)

	screen(2)
	# Panel B: Vary reliability at fixed true correlation
	par(mar=c(5,0,1,5))
	plot(c(0,1), c(0,1), type='n', xlim=c(0,1), ylim=c(0,1.1), las=1, yaxt='n', xlab='Reliability', ylab='')
	# The true correlation target
	corr.target <- 0.9
	abline(h=corr.target^dosq,lty='dotted')
	#abline(h=1,lty='dashed')
	points(rep(as.numeric(names(corEst$varyrel)), sapply(corEst$varyrel, nrow)),
		sapply(corEst$varyrel, function(x){x$obs^dosq}), pch=pch, col=obscol, cex=cex.pts)
	points(rep(as.numeric(names(corEst$varyrel)), sapply(corEst$varyrel, nrow)),
		sapply(corEst$varyrel, function(x){x$est^dosq}), pch=pch, col=estcol, cex=cex.pts)

	screen(3)
	# Panel C: At fixed reliability and true correlation, subsample randomly
	par(mar=c(5,5,1,0))
	xlim <- c(100,n)
	plot(c(0,1), c(0,1), type='n', xlim=xlim, ylim=c(0,1.1), las=1, xlab='Number of detected genes', ylab=paste('Observed/corrected\n',sqstr,' correlation',sep=''), log='x')
	corr.target <- 0.9
	abline(h=corr.target^dosq,lty='dotted')
	#abline(h=1,lty='dashed')
	points(rep(as.numeric(names(corEst$mar)), sapply(corEst$mar, nrow)),
		sapply(corEst$mar, function(x){x$obs^dosq}), pch=pch, col=obscol, cex=cex.pts)
	points(rep(as.numeric(names(corEst$mar)), sapply(corEst$mar, nrow)),
		sapply(corEst$mar, function(x){x$est^dosq}), pch=pch, col=estcol, cex=cex.pts)

	screen(4)
	# Panel D: At fixed reliability and true correlation, subsample non-randomly (taking top M genes ranked by first variable)
	par(mar=c(5,0,1,5))
	plot(c(0,1), c(0,1), type='n', xlim=bias.xlim, ylim=c(0,1.1), las=1, yaxt='n', xlab='Number of detected genes', ylab='', log='x')
	bias.corr.target <- 0.9
	abline(h=bias.corr.target^dosq,lty='dotted')
	#abline(h=1,lty='dashed')
	points(rep(as.numeric(names(corEst$nmar)), sapply(corEst$nmar, nrow)),
		sapply(corEst$nmar, function(x){x$obs^dosq}), pch=pch, col=obscol, cex=cex.pts)
	points(rep(as.numeric(names(corEst$nmar)), sapply(corEst$nmar, nrow)),
		sapply(corEst$nmar, function(x){x$est^dosq}), pch=pch, col=estcol, cex=cex.pts)
	points(rep(as.numeric(names(corEst$nmar.inf)), sapply(corEst$nmar.inf, length)),
		   unlist(corEst$nmar.inf)^dosq, pch=pch, col=fitcol, cex=cex.pts)
	close.screen(all=TRUE)
	if (file.out) dev.off()
	#source("~/research/spearman/analysis/spearman-sim.R")
}

if (fig.reliability) {
	meth='p'
	pch=16
	within.pch = 16
	within.col = 'gray70'
	between.col = tcol('black',0.7)
	rc.mrna.within <- sapply(unique(yres$exp.vars$mrna$ms), function(m) {
		d <- subset(yres$exp.vars$mrna, ms==m)
		f <- p.0('mrna',d$alias)
		if (length(f)>1) {
			rc <- rcormat(log.nozero(yres$raw[,f]), log=F, meth=meth, na.rm=F)
			cbind(lt(rc$r), lt(rc$n))
		} else {
			c(NA,NA)
		}
		})
	rc.prot.within <- sapply(unique(yres$exp.vars$prot$ms), function(m) {
		d <- subset(yres$exp.vars$prot, ms==m)
		f <- p.0('prot',d$alias)
		if (length(f)>1) {
			rc <- rcormat(log.nozero(yres$raw[,f]), log=F, meth=meth, na.rm=F)
			cbind(lt(rc$r), lt(rc$n))
		} else {
			c(NA,NA)
		}
		})
	mrna.within.r = unlist(sapply(rc.mrna.within, function(x){if (!is.null(nrow(x))) {x[,1]} else {x[1]}}))
	mrna.within.n = unlist(sapply(rc.mrna.within, function(x){if (!is.null(nrow(x))) {x[,2]} else {x[2]}}))
	med.mrna.within = median(mrna.within.r, na.rm=T)
	prot.within.r = unlist(sapply(rc.prot.within, function(x){if (!is.null(nrow(x))) {x[,1]} else {x[1]}}))
	prot.within.n = unlist(sapply(rc.prot.within, function(x){if (!is.null(nrow(x))) {x[,2]} else {x[2]}}))
	med.prot.within = median(prot.within.r, na.rm=T)
	
	rc.mrna <- rcormat(yres$raw.ms.avg[yres$ms.fields$mrna], log=F, meth=meth)
	rc.prot <- rcormat(yres$raw.ms.avg[yres$ms.fields$prot], log=F, meth=meth)

	med.mrna.between = median(lt(rc.mrna$r), na.rm=T)
	med.prot.between = median(lt(rc.prot$r), na.rm=T)
	
	if (file.out) dev.out('fig.reliability', width=10, height=3.5, output.type=output.type)
	split.screen(c(1,3))
	screen(1)
	par(mar=c(4,4,1,0))
	plot(1,1, type='n', pch=pch, las=1, xlim=c(0,6000), ylim=c(0,1), xaxs='r', yaxs='i', xlab='Number of mRNAs detected\nin both studies', ylab='Correlation between studies')
	#points(mrna.within.n, mrna.within.r, col=within.col, pch=pch)
	points(lt(rc.mrna$n), lt(rc.mrna$r), col=between.col, pch=pch)
	#points(rc.mrna$n['mrna.ingolia.2009','mrna.ingolia.2010'], rc.mrna$r['mrna.ingolia.2009','mrna.ingolia.2010'], col='blue', cex=2, pch=pch)
	points(rc.mrna$n['mrna.causton','mrna.holstege'], rc.mrna$r['mrna.causton','mrna.holstege'], col='blue', cex=2, pch=pch)
	abline(h=med.mrna.within, lty='dotted')
	abline(h=med.mrna.between, lty='dashed')
	screen(2)
	par(mar=c(4,0,1,4))
	plot(1, 1, type='n', pch=pch, las=1, xlim=c(0,6000), ylim=c(0,1), xaxs='r', yaxt='n', yaxs='i', xlab='Number of proteins detected\nin both studies', ylab='Correlation between studies')
	#points(prot.within.n, prot.within.r, col=within.col, pch=pch)
	points(lt(rc.prot$n), lt(rc.prot$r), col=between.col, pch=pch)
	points(rc.prot$n['prot.nagaraj','prot.degodoy'], rc.prot$r['prot.nagaraj','prot.degodoy'], col='blue', cex=2, pch=pch)
	abline(h=med.prot.within, lty='dotted')
	abline(h=med.prot.between, lty='dashed')
	screen(3)
	par(mar=c(4,1,1,1))
	lg.raw.set <- na.omit(data.frame(yres$raw.ms.avg[,c('mrna.yassour','mrna.lipson','prot.degodoy','prot.ghaem','orf')]))
	#lg.imp.set <- na.omit(data.frame(yres$imp.ms.avg[,c('mrna.yassour','mrna.lipson','prot.degodoy','prot.ghaem','orf')]))
	d <- lg.raw.set
	gatt <- function(x,attrname) {x[[attrname]]}
	r.est <- function(x) {x$estimate}
	r.low <- function(x) {x$conf.int[1]}
	r.up <- function(x) {x$conf.int[2]}
	g <- list(c11 <- cortest(d$mrna.yassour, d$prot.degodoy, method=meth),
		c12 <- cortest(d$mrna.yassour, d$prot.ghaem, method=meth),
		c21 <- cortest(d$mrna.lipson, d$prot.degodoy, method=meth),
		c22 <- cortest(d$mrna.lipson, d$prot.ghaem, method=meth),
		cortest(d$mrna.lipson+d$mrna.yassour, d$prot.degodoy+d$prot.ghaem, method=meth))
	rsp <- cor.sp(lg.raw.set[,1:4], meth='p', na.rm=T) #geo.mean(c(c11$estimate, c12$estimate, c21$estimate, c22$estimate))/sqrt(cor(d$mrna.yassour,d$mrna.lipson,meth=meth)*cor(d$prot.degodoy,d$prot.ghaem,meth=meth))
	barplot.err(c(sapply(g,r.est),rsp$estimate), lower=c(sapply(g,r.low),rsp$conf.int[1]), upper=c(sapply(g,r.up),rsp$conf.int[2]), space=c(0.1,0.1,0.1,0.1, 0.8, 0.8), ylim=c(0,1.0), las=1)
	close.screen(all=TRUE)
    if (file.out) dev.off()
    
    
}

if (fig.correlations.between.scatter) {
	meth = 's'
	y <- yres$est
	y$mrna[is.na(yres$Xc.avg$mrna)] <- NA
	y$prot[is.na(yres$Xc.avg$prot)] <- NA
	L.noimp.ct <- cortest(y$mrna, y$prot, log=T, meth=meth)
	# Censored estimate: correlations for genes for which we observe at least one protein and mRNA level
	# Spearman-corrected estimate
	n.cutoff <- 3840
	# for the set of datasets which yield n>cutoff, what is the reliability?
	
	
	# DAD: issue is that we do not want to include reliability estimates from the same lab doing different studies.
	# So no Causton/Holstege, and no de Godoy/Thakur/Nagaraj.
	nprot <- count.pairwise(yres$raw.ms.avg[yres$ms.fields$prot], yres$raw.ms.avg[yres$ms.fields$prot])
	rc.prot <- cor(yres$raw.ms.avg[yres$ms.fields$prot], method=meth, use='pairwise.complete.obs')
	# Delete the offending same-lab/different-ms values
	rc.prot['prot.nagaraj','prot.thakur'] <- NA
	rc.prot['prot.nagaraj','prot.degodoy'] <- NA
	rc.prot['prot.degodoy','prot.thakur'] <- NA
	rc.prot['prot.thakur','prot.nagaraj'] <- NA
	rc.prot['prot.degodoy','prot.nagaraj'] <- NA
	rc.prot['prot.thakur','prot.degodoy'] <- NA
	diag(nprot) <- NA
	diag(rc.prot) <- NA
	nprot['prot.nagaraj','prot.thakur'] <- NA
	nprot['prot.nagaraj','prot.degodoy'] <- NA
	nprot['prot.degodoy','prot.thakur'] <- NA
	nprot['prot.thakur','prot.nagaraj'] <- NA
	nprot['prot.degodoy','prot.nagaraj'] <- NA
	nprot['prot.thakur','prot.degodoy'] <- NA
	nmrna <- count.pairwise(yres$raw.ms.avg[yres$ms.fields$mrna], yres$raw.ms.avg[yres$ms.fields$mrna])
	rc.mrna <- cor(yres$raw.ms.avg[yres$ms.fields$mrna], method=meth, use='pairwise.complete.obs')
	# Delete the offending same-lab/different-ms values
	rc.mrna['mrna.causton','mrna.holstege'] <- NA
	rc.mrna['mrna.holstege','mrna.causton'] <- NA
	#rc.mrna['mrna.ingolia.2009','mrna.ingolia.2010'] <- NA
	#rc.mrna['mrna.ingolia.2010','mrna.ingolia.2009'] <- NA
	diag(nmrna) <- NA
	diag(rc.mrna) <- NA
	nmrna['mrna.causton','mrna.holstege'] <- NA
	nmrna['mrna.holstege','mrna.causton'] <- NA

	avg.corrs <- cor(yres$raw.ms.avg[yres$ms.fields$mrna], yres$raw.ms.avg[yres$ms.fields$prot], method=meth, use='pairwise.complete.obs')
	avg.corrs.s <- cor(yres$raw.ms.avg[yres$ms.fields$mrna], yres$raw.ms.avg[yres$ms.fields$prot], method='s', use='pairwise.complete.obs')
	avg.corrs.p <- cor(yres$raw.ms.avg[yres$ms.fields$mrna], yres$raw.ms.avg[yres$ms.fields$prot], method='p', use='pairwise.complete.obs')
	nboth <- count.pairwise(yres$raw.ms.avg[yres$ms.fields$mrna], yres$raw.ms.avg[yres$ms.fields$prot])
	
	# Largest dataset
	lg.raw.set <- na.omit(data.frame(yres$raw.ms.avg[,c('mrna.lipson','mrna.yassour','prot.degodoy','prot.ghaem','orf')]))
	#lg.raw.set <- na.omit(data.frame(yres$raw.ms.avg[,c('mrna.lipson','mrna.yassour','prot.degodoy','prot.ghaem','orf')]))
	colnames(lg.raw.set) <- c('x1','x2','y1','y2','orf')
	lg.rsp <- cor.sp(lg.raw.set[,c('x1','x2','y1','y2')], method=meth)
	lg.n <- nrow(lg.raw.set)
	#lg.imp.avg.ct <- cortest(rowMeans(lg.imp.set[,1:2]),rowMeans(lg.imp.set[,3:4]), meth='p')
	lg.raw.avg.ct <- cortest(rowMeans(lg.raw.set[,1:2]),rowMeans(lg.raw.set[,3:4]), meth='p')
	
	bin.size <- 500
	if (meth=='p') {
		rsp <- yres$rsp.pearson
	} else {
		rsp <- yres$rsp.spearman
	}
	rsp.bins <- data.frame(t(sapply(seq(bin.size,3820,bin.size), function(m) {
		y <- subset(rsp, n > (m-bin.size) & n < m)
		c(n=m-bin.size/2, mean=mean(y$r),sd=sd(y$r), n.points=nrow(y), se.95=1.96*sd(y$r)/sqrt(nrow(y)))
	})))

	pch <- 16
	pch.tri <- 18
    pch.rsp <- 15
    cex <- 0.8
    big.cex <- 1.1
    cex.rsp <- 1.1
	    
	if (file.out) dev.out(p.0('fig.correlations.scatter',meth), width=5, height=4, output.type=output.type)
	plot(1, 1, type='n', xlim=c(0,6000), ylim=c(0,1), pch=pch, cex=cex, las=1, yaxs='i', xlab='Number of genes with quantified mRNA and protein levels', ylab='mRNA-protein correlation')
	points(nboth, avg.corrs, pch=pch, col=tcol('black',0.7), cex=cex)
	points.err(rsp.bins$n, rsp.bins$mean, y.lower=rsp.bins$mean-rsp.bins$sd, y.upper=rsp.bins$mean+rsp.bins$sd, col='black', pch=pch.rsp, cex=cex.rsp)
	#points(L.noimp.ct$N, L.noimp.ct$estimate, pch=pch, col='darkblue', cex=big.cex)
	points(yres$n.genes, yres$model$psi['expr','abund'], pch=5, col='red', cex=big.cex)
	points.err(lg.rsp$N, lg.rsp$estimate, y.lower=lg.rsp$range[1], y.upper=lg.rsp$range[2], pch=pch.rsp, col='blue', cex=big.cex)
	if (file.out) dev.off()
	
}

if (fig.averaging.increases.corr) {

	x <- yres$raw
	sub.f <-c('prot.ghaem','prot.newman','prot.peng','prot.gygi')
	target.prot <- sub.f[[1]]
	prots <- sapply(sub.f, function(m){
		z <- match(na.omit(x[,c('orf',m)])$orf, x$orf)
		x[z,target.prot]
		})

	max.meas <- 8
	ys <- subset(yres$est, n.mrna>=max.meas & n.prot>=max.meas)
	x <- yres$Xc.ms.avg[zs <- match(ys$orf, yres$Xc.ms.avg$orf),]
	n <- 1:max.meas
	nreps <- 50
	corrs <- sapply(n, function(m) {
		sub.corrs <- sapply(1:nreps, function(r){
			mrna.inds <- sample(1:length(yres$ms.fields$mrna), m)
			prot.inds <- sample(1:length(yres$ms.fields$prot), m)
			if (m>1) {
				xm <- rowMeans(x[,yres$ms.fields$mrna[mrna.inds]],na.rm=T)
				xp <- rowMeans(x[,yres$ms.fields$prot[prot.inds]],na.rm=T)
			} else {
				xm <- x[,yres$ms.fields$mrna[mrna.inds]]
				xp <- x[,yres$ms.fields$prot[prot.inds]]
			}
			cortest(xm, xp, meth='p')$estimate
		})
		c(mean=mean(sub.corrs),sd=sd(sub.corrs))
		})

	slopes <- sapply(n, function(m) {
		sub.slopes <- sapply(1:nreps, function(r){
			mrna.inds <- sample(1:length(yres$ms.fields$mrna), m)
			prot.inds <- sample(1:length(yres$ms.fields$prot), m)
			if (m>1) {
				xm <- exp(rowMeans(x[,yres$ms.fields$mrna[mrna.inds]],na.rm=T))
				xp <- ys$prot
			} else {
				xm <- exp(x[,yres$ms.fields$mrna[mrna.inds]])
				xp <- ys$prot
			}
			# 4,3 = RMA, slope
			llmodel2(xp~xm)$regression.results[4,3]
		})
		c(mean=mean(sub.slopes),sd=sd(sub.slopes))
		})

	subset.prots <- list(yres$raw$prot.ghaem, yres$raw[zs <- match(ys$orf, yres$Xc.ms.avg$orf),]$prot.ghaem)

	mar = c(4,4,1,1)
	xlim <- c(10,1e7)

	if (file.out) dev.out("fig.averaging.increases.corr", width=9, height=3, output.type=output.type)
	split.screen(c(1,3))
	screen(1)
	par(mar=mar)
	cols <- myrainbow(length(sub.f))
	multidens(prots, log=T, xlim=xlim, xlab='Protein level (mol./cell, western blot)', fill=T, col=cols, rel='c', line.col='black')
	abline(v=sapply(prots, median, na.rm=T), lwd=2, col=cols)
	screen(2)
	par(mar=mar)
	multidens(subset.prots, log=T, xlim=xlim, rel='c', fill=T, col=c('gray30','gray80'), line.col='black', xlab='Protein level (mol./cell, western blot)')
	text(1e5, 0.5, label=paste('N =', nrow(ys)))
	screen(3)
	par(mar=mar)
	plot.err(n, corrs['mean',], y.lower=corrs['mean',]-corrs['sd',], y.upper=corrs['mean',]+corrs['sd',], pch=16, ylim=c(0,1), yaxs='i', xlab='Number of measurements averaged', ylab='mRNA--protein correlation', las=1)
	#plot.err(n, slopes['mean',], y.lower=slopes['mean',]-slopes['sd',], y.upper=slopes['mean',]+slopes['sd',], pch=16, ylim=c(0,1.5), yaxs='i', xlab='Number of measurements averaged', ylab='mRNA--protein correlation', las=1)
	axis(4,axTicks(4))
	close.screen(all=TRUE)
	if (file.out) dev.off()
}



if (fig.mrna.prot.scatter) {
	mrna.lim <- c(1e-2, 4e2)
	prot.lim <- c(2e-1, 4e6)
	mar=c(4,4,1,1)
	
	if (file.out) dev.out("fig.mrna.prot.scatter", width=8, height=8, output.type=output.type)
	split.screen(c(2,2))
	screen(1)
	par(mar=mar)
	# Reproducible subset
	x <- yres$est
	rep.sub <- subset(x, n.mrna>1 & n.prot>1)
	nrep.sub <- subset(x, n.mrna==1 | n.prot==1)
	det.sub <- subset(x, n.prot>0)
	noprot.sub <- subset(x, n.prot==0)
	cols <- tcol(c(myrainbow(3),'black','gray'),0.7)
	names(cols) <- c('noprot','nrep','det','rep','all')

	lplot(1, 1, type='n', xlim=mrna.lim, ylim=prot.lim, pch=16, cex=0.5, col='gray80', xaxs='i', yaxs='i', xlab='mRNA level (mol./cell)', ylab='Protein level (mol./cell)')
	points(noprot.sub$mrna, noprot.sub$prot, pch=16, col=cols[['noprot']], cex=0.5)
	points(nrep.sub$mrna, nrep.sub$prot, pch=16, col=cols[['nrep']], cex=0.5)
	points(det.sub$mrna, det.sub$prot, pch=16, col=cols[['det']], cex=0.7)
	points(rep.sub$mrna, rep.sub$prot, pch=16, col=cols[['rep']], cex=0.7)
	screen(2)
	par(mar=mar)
	plist <- list('all'=yres$est$mrna, 'det'=det.sub$mrna, 'noprot'=noprot.sub$mrna, 'rep'=rep.sub$mrna, 'nrep'=nrep.sub$mrna)
	multidens(plist, xlim=mrna.lim, log=T, col=unlist(cols[names(plist)]), fill=T, xaxs='i', rel='c')
	screen(3)
	par(mar=mar)
	plist <- list('all'=yres$est$prot, 'det'=det.sub$prot, 'noprot'=noprot.sub$prot, 'rep'=rep.sub$prot, 'nrep'=nrep.sub$prot)
	multidens(plist, xlim=prot.lim, log=T, col=unlist(cols[names(plist)]), fill=T, xaxs='i', rel='c')
	screen(4)
	par(mar=mar)
	x <- yres$est.sample
	rep.sub <- subset(x, n.mrna>1 & n.prot>1)
	nrep.sub <- subset(x, n.mrna==1 | n.prot==1)
	det.sub <- subset(x, n.prot>0)
	noprot.sub <- subset(x, n.prot==0)
	cols <- tcol(c(myrainbow(3),'black'),0.7)
	lplot(1, 1, type='n', xlim=mrna.lim, ylim=prot.lim, pch=16, cex=0.5, col='gray80', xaxs='i', yaxs='i', xlab='mRNA level (mol./cell)', ylab='Protein level (mol./cell)')
	points(noprot.sub$mrna, noprot.sub$prot, pch=16, col=cols[1], cex=0.5)
	points(nrep.sub$mrna, nrep.sub$prot, pch=16, col=cols[2], cex=0.5)
	points(det.sub$mrna, det.sub$prot, pch=16, col=cols[3], cex=0.7)
	points(rep.sub$mrna, rep.sub$prot, pch=16, col=cols[4], cex=0.7)

	close.screen(all=TRUE)
	if (file.out) dev.off()
}

if (fig.scm.abs.standard.comparison) {
	meth = 'p'
	pic <- read.table("../data/abs/picotti09-absolute-abundances.txt", header=T)
	z.pic <- match(pic$orf, yres$bg$orf)
	zk <- read.table("../data/abs/zenklusen08-fish-mrna.txt", header=T)
	z.zk <- match(zk$orf, yres$bg$orf)
	avg <- yres$est[z.pic,]
	cg <- coef(lm(log.nz(prot)~log.nz(mrna), data=yres$est))
	cg2 <- lmodel2(log.nz(prot)~log.nz(mrna), data=yres$est, range.x='interval', range.y='interval', nperm=1)$regression.results[4,]
	cat("# Slope of prot~mrna", as.numeric(cg2[3]), '\n')
	rib <- yres$est[yres$subset$ribosome,'prot']
	cat("# Median ribosomal protein level", med.rib <- median(rib,na.rm=T), "\n")
	
	# Accuracy assessment
	cat("# Average mRNA fold-change =",  format(exp(abs(mean(log(yres$est[z.zk,]$mrna/zk$mrna)))),2), "+/-", format(sd(yres$est[z.zk,]$mrna/zk$mrna),2), "\n")
	cat("# Average protein fold-change =",  format(exp(abs(mean(log(avg$prot/pic$srm)))),2), "+/-", format(sd(avg$prot/pic$srm),2), "\n")
	cat("# Average Ghaem fold-change =",  format(exp(abs(mean(log(yres$raw[z.pic,]$prot.ghaem/pic$srm)))),2), "+/-", format(sd(yres$raw[z.pic,]$prot.ghaem/pic$srm),2), "\n")
	prot.per.cell <- sum(yres$est$prot*yres$bg$mw/6.022e23, na.rm=T)
	cat("# Total protein mass per haploid cell =", format(prot.per.cell,2),'\n')
	cat("# Total protein molecules per haploid cell =", sum(yres$est$prot, na.rm=T),'\n')
	
	mar = c(4,4,1,1)
	
	if (file.out) dev.out('fig.scm.abs.standard.comparison', width=3.5, height=7, output.type=output.type)
	split.screen(c(2,1))
	screen(1)
	par(mar=mar)
	lim <- c(1,30)
	plot.err(yres$est[z.zk,]$mrna, zk$mrna, x.lower=yres$est[z.zk,]$mrna.lower.sd, x.upper=yres$est[z.zk,]$mrna.upper.sd, y.lower=zk$mrna-zk$mrna.se, y.upper=zk$mrna+zk$mrna.se, xlim=lim, ylim=lim, pch=16, las=1,
		xlab='Estimated absolute mRNA level (mol./cell)', ylab='Measured absolute mRNA level (mol./cell)')
	abline(0,1,lty='dotted')
	screen(2)
	par(mar=mar)
	prot.lim <- c(20,4e6)
	sd.const <- 1.96
	plot.err(avg$prot, pic$srm, x.lower=avg$prot.lower.sd, x.upper=avg$prot.upper.sd, y.upper=pic$srm+sd.const*pic$sd.srm, y.lower=pic$srm-sd.const*pic$sd.srm, log='xy', pch=16, xlim=prot.lim, ylim=prot.lim,
		xaxt='n', yaxt='n', ylab='Measured absolute protein level (mol./cell)', xlab='Estimated absolute protein level (mol./cell)')
	my.axis(1, prot.lim, log=T, las=1)
	my.axis(2, prot.lim, log=T, las=1)
	pcor(cortest(avg$prot, pic$srm, log=T, meth=meth))
	abline(0,1, lty='dotted')
	close.screen(all=TRUE)
	if (file.out) dev.off()
}


if (fig.undetected.prot) {
	# Isolate set of undetected proteins
	y <- subset(data.frame(yres$est, rdens=yres$te$data$rd.median), n.mrna>0)
	yesprot <- subset(y, n.prot>0)
	noprot <- subset(y, n.prot==0)

	z <- match(y$orf, yres$bg$orf)
	z.no <- match(noprot$orf, yres$bg$orf)
	cat("# Found", nrow(subset(y, n.mrna>0 & n.prot==0)), "genes with undetected protein but detected mRNA\n")
	cat("# Found", nrow(subset(noprot, !is.na(rdens) & rdens>0)), "genes with undetected protein but detected ribosome footprints\n")
	cat("# Found", nrow(subset(y, n.mrna>0 & n.prot==0 & is.na(rdens))), 
		"genes with undetected protein and no detected ribosome footprints.\n")
	pcor(ct <- cortest(yesprot$prot, yesprot$rdens, log=T, meth='p'))
	pcor(no.ct <- cortest(noprot$prot, noprot$rdens, log=T, meth='p'))
	
	cols <- tcol(c(myrainbow(3),'black','gray'),0.7)
	names(cols) <- c('noprot','nrep','det','rep','all')

	mar=c(4,4,1,1)
	if (file.out) dev.out("fig.undetected.prot", width=4, height=8, output.type=output.type)
	split.screen(c(2,1))
	screen(1)
	par(mar=mar)
	xlim <- c(1e-4, 1e3)
	lplot(noprot$rdens, noprot$prot, xlim=xlim, col=tcol('black',0.6), pch=16, xlab='Ribosome footprint density (rpkm)', ylab='Estimated protein level (mol./cell)')
	g <- llmodel2(prot~rdens, data=noprot)
	labline(g, x=xlim, method='RMA', lty='dotted')
	text(0.06, 5e2, label=paste("r =",round(no.ct$estimate, 2)), pos=2)
	screen(2)
	par(mar=mar)
	multidens(list(y$rdens, noprot$rdens, yesprot$rdens), col=unlist(cols[c('all','noprot','rep')]), lwd=2, log=T, fill=T, xlim=xlim, xlab='Ribosome footprint density (rpkm)', rel='c')
	close.screen(all=TRUE)
	if (file.out) dev.off()
}

if (fig.paxdb.vs.scm) {
	pax <- read.table("~/research/data/scerevisiae/paxdb-proteome/paxdb-scerevisiae-integrated.txt", header=T, sep='\t')
	z <- match(yres$bg$orf, pax$orf)
	pcor(cortest(pic$srm, pax[z,][z.pic,]$abundance, meth='p', log=T))
	pcor(cortest(pic$srm, yres$est[z.pic,]$prot, meth='p', log=T))
	pcor(cortest(pic$srm, pax[z,][z.pic,]$abundance, meth='p', log=F))
	pcor(cortest(pic$srm, yres$est[z.pic,]$prot, meth='p', log=F))
	pcor(cortest(yres$est$prot, pax[z,]$abundance, meth='s'))
}

if (fig.mrna.prot.overview.distribution) {
	if (file.out) dev.out("fig.mrna.prot.overview.distribution", width=8, height=4, output.type=output.type)
	multidens(yres$est[,c('mrna','prot')], fill=T, col=tcol(myrainbow(2),0.8), 
		line.col='black', log=T, xlab='Average cellular level (molecules/cell)')
	if (file.out) dev.off()
}

model2.results.fname <- '../data/model2-results.Rdata'

if (calc.exponent) {
	# Estimate the exponents by model-II regression
	m2f <- function(x, y, fn) {
		d <- na.omit(data.frame(x=x,y=y))
		g <- fn(y~x, data=d, nperm=1, range.x='interval', range.y='interval')
		c(g$regression.results[,3], nrow(d))
	}

	res <- NULL
	mflds <- yres$fields$mrna
	pflds <- yres$fields$prot
	d <- yres$raw
	fn <- llmodel2
	for (xm in 1:length(mflds)) {
		for (xp in 1:length(pflds)) {
			r <- m2f(d[,mflds[xm]], d[,pflds[xp]], fn)
			res <- rbind(res, c(xm, xp, r))
		}
	}
	m2res <- as.data.frame(res)

	colnames(m2res) <- c('xm','xp','OLS','MA','SMA','RMA','n')
	rownames(m2res) <- NULL
	save(m2res, file=model2.results.fname)
	cat("# Wrote Model II regression results to", model2.results.fname, '\n')
} else {
	# Load
	m2res <- get(load(model2.results.fname))
	cat("# Read Model II regression results from", model2.results.fname, '\n')
}


if (fig.exponent) {
	method = 'p'
	clog = TRUE

	d <- data.frame(yres$est[,c('mrna','prot')], samp=yres$est.sample[,c('mrna','prot')],
		yres$te$data[,c('ing1',yres$te$te.mrna.names, yres$te$te.rd.names,'mrna.median','rd.median')],prot1=yres$raw$prot.degodoy)
	db <- d

	#dcomp <- na.omit(log.nozero(data.frame(prot1=yres$raw$prot.degodoy, prot2=yres$raw$prot.ghaem, mrna1=exp(yres$raw.ms.avg$mrna.yassour), mrna2=exp(yres$raw.ms.avg$mrna.lipson), rdens1=yres$bg$rdens.ingolia, rdens2=yres$bg$rdens.gerashchenko)))
	meancorr <- mean(yres$corrs)
	confcorr95 <- 1.96*sd(yres$corrs)
	scm.est <- list(estimate=meancorr, conf.int=c(meancorr-confcorr95, meancorr+confcorr95))
	ct.list <- list(
		cortest(db$ing1, db$prot1, meth=method, log=clog),
		cortest(db$rd.ing1, db$prot1, meth=method, log=clog),
		cortest(db$rd.ger1, db$prot1, meth=method, log=clog),
		cortest(db$rd.ger2, db$prot1, meth=method, log=clog),
		cortest(db$rd.sub1, db$prot1, meth=method, log=clog),
		cortest(db$rd.mcm1, db$prot1, meth=method, log=clog),
		cortest(db$ing1, db$prot, meth=method, log=clog),
		cortest(db$rd.ing1, db$prot, meth=method, log=clog),
		cortest(db$rd.ger1, db$prot, meth=method, log=clog),
		cortest(db$rd.ger2, db$prot, meth=method, log=clog),
		cortest(db$rd.sub1, db$prot, meth=method, log=clog),
		cortest(db$rd.mcm1, db$prot, meth=method, log=clog),
		scm.est
		#cor.sp.matrix(yres$te$data[,yres$te$te.rd.names], yres$est$prot, meth=method, log=clog) #cortest(db$rd.median, db$prot, meth=method, log=clog)
		)

	# Estimate from the G factors
	prot.names <- c('futcher','gygi','newman','degodoy','lee','lu','nagaraj','peng','thakur','washburn','ghaem')
	mrna.names <- c('causton','holstege','lipson2','roth','miura','dudley','garcia','mackay','pelechano','ingolia','lipson','nagalakshmi','yassour','velc')
	G.prot <- yres$model$G[prot.names]
	G.mrna <- yres$model$G[mrna.names]
	medex <- median(G.prot)/median(G.mrna)

	# Subset of the data in which at least half the genes are detected
	l1.m2 <- subset(m2res, n>5854/2)

	# Barplot of slopes
	# Scatterplot comparison

	mar = c(4,4,1,1)
	if (file.out) dev.out("fig.exponent", width=9, height=6, output.type=output.type)
	split.screen(c(2,3))
	screen(1)
	par(mar=mar)
	dosq <- 1 # Square the results?
	#corr.means <- cortest()
	mids <- barplot.err(ct.est <- sapply(ct.list, function(x){x$estimate^dosq}), 
		sapply(ct.list, function(x){x$conf.int[1]^dosq}), 
		sapply(ct.list, function(x){x$conf.int[2]^dosq}), 
		las=1, ylim=c(0,1.0), space=c(0.2,0.2,0.2,0.2,0.2,0.2,1,0.2,0.2,0.2,0.2,0.2,1), 
		col=c(rep(c('orange','darkgreen','darkgray','darkgray','darkgray','darkgray'),2),'blue'), names.arg=NA, ylab='Correlation with protein levels')
	
	screen(2)
	par(mar=mar)
	reg.lim <- c(0,5)
	plot(m2res$OLS, m2res$RMA, cex=0.8, pch=16, col='gray', xlim=reg.lim, ylim=reg.lim, las=1, xlab='Ordinary least-squares regression slope', ylab='Reduced major-axis regression slope')
	points(l1.m2$OLS, l1.m2$RMA, cex=0.8, pch=16)
	abline(0,1, lty='dotted', lwd=2)
	abline(h=medex, lty='dashed', lwd=2)
	
	screen(3)
	par(mar=mar)
	slopes <- m2res[,c('OLS','RMA')]
	x <- apply(slopes, 2, ecdf)
	s <- seq(reg.lim[1],reg.lim[2],0.01)
	matplot(s, g <- sapply(x,function(m){m(s)}), type='l', lwd=2, lty='solid', col=myrainbow(4), las=1, ylab='Proportion', xlab='Regression slope')
	#x <- multi.ecdf(g <- m2res[,c('OLS','MA','SMA','RMA')], lwd=2, xlim=c(0,3), xlab='Slope (mRNA\u2013protein exponent)', plot=F)
	#s <- seq(0,3,0.1); matplot(s, sapply(x,function(m){m(s)}), type='l', lwd=2, col=myrainbow(4))
	abline(h=0.5, col='gray')
	abline(v=apply(slopes, 2, median), col=myrainbow(4), lwd=2, lty='dotted')
	abline(v=medex, lwd=2, lty='dashed')
	#plot(n.prot, rma.slopes, cex=0.8, pch=16, col='black', log='x', las=1, ylim=c(1,2), xlab='Number of proteins detected', ylab='Exponent')

	# Use subset of new ribosome-profiling data with measurements in all three studies
	#y <- subset(yres$te$data, rd.n>=3)
	y <- subset(yres$te$data, rd.n==5 & mrna.n==3)
	cat("RD vs mRNA:", nrow(y),"genes\n")

	screen(4)
	par(mar=mar)
	lplot(y$mrna, y$rd.median, col=tcol('black',0.4), xlab='mRNA level (mol./cell)', ylab='Relative ribosome density (AU)')
	labline(g <- llmodel2(rd.median~mrna, data=y), method='RMA', lwd=2, col='gray')
	labline(g, lwd=2, slope=1, lty='dotted',col='gray')

	cat("TE vs mRNA:", nrow(y),"genes\n")
	screen(5)
	par(mar=mar)
	lplot(y$mrna, y$te, col=tcol('black',0.4), xlab='mRNA level (mol./cell)', ylab='Relative translational efficiency (AU)')

	# Use subset of new mRNA data with measurements in all three studies
	ym <- subset(yres$te$data, rd.n>=3 & mrna.n==3)

	cat("RD and protein vs SCM and new mRNA:", nrow(ym),"genes\n")
	screen(6)
	par(mar=mar)
	lslxy <- function(x,y,...) {llmodel2(y~x,...)$regression.results[4,3]}
	lslyx <- function(y,x) {llmodel2(y~x)$regression.results[4,3]}
	lslxy.ci <- function(x,y) {as.matrix(unlist(llmodel2(y~x)$confidence.intervals[4,4:5]))[,drop=F]}
	lslyx.ci <- function(y,x) {unlist(llmodel2(y~x)$confidence.intervals[4,4:5])}
	m.slopes <- c(lslxy(ym$mrna, ym$rd.median), lslxy(ym$mrna, ym$prot))
	mm.slopes <- c(lslxy(ym$mrna.median, ym$rd.median), lslxy(ym$mrna.median, ym$prot))
	m.slopes.ci <- t(cbind(lslxy.ci(ym$mrna,ym$rd.median),lslxy.ci(ym$mrna,ym$prot)))
	mm.slopes.ci <- t(cbind(lslxy.ci(ym$mrna.median,ym$rd.median),lslxy.ci(ym$mrna.median,ym$prot)))
	barplot.err(cbind(m.slopes,mm.slopes), cbind(m.slopes.ci[,1], mm.slopes.ci[,1]), cbind(m.slopes.ci[,2], mm.slopes.ci[,2]), 
		beside=T, ylim=c(0,2), space=c(0,0.2), las=1)
	abline(h=1.69, lty='dashed')
	abline(h=1, lty='dotted')
	close.screen(all=TRUE)
	if (file.out) dev.off()
}

if (fig.relative.dynamic.range) {
	if (file.out) dev.out("fig_relative_dynamic_range", width=9, height=3, output.type=output.type)
	par(mar=mar)
	cols <- c('darkblue','darkblue','orange','purple')
	scaleto <- function(x, y) {
		med <- median(y, na.rm=T)
		xl <- log.nz(x)
		exp(xl-median(xl,na.rm=T)+log.nz(med))
	}
	# 85% of 200k ribosomes engaged in translation gives us how many per mRNA species?
	ribosomes <- yres$te$data$rd.median*(yres$bg$len-3)
	ribosomes <- ribosomes*200000*0.85/sum(ribosomes,na.rm=T)
	normexp <- data.frame(mrna=yres$est$mrna, 
		mrna.median=scaleto(yres$te$data$mrna.median, yres$est$mrna),
		rd.median=ribosomes, #scaleto(yres$te$data$rd.median, yres$est$mrna), 
		prot=yres$est$prot, rd.n=yres$te$data$rd.n)
	normexp.sub <- na.omit(subset(normexp, rd.n>4))
	multidens(normexp[,1:4], 
		xlab='Normalized level', log=T, lwd=2, col=cols, lty=c('solid','2222','1212','solid'))

	# The width of 80% of the data
	distr.width <- 0.8
	distr.diff <- (1-distr.width)/2
	rel.width <- function(w, x){
		# Get quantiles for fraction w of the data, symmetrically
		y <- apply(log.nz(x), 2, quantile, probs=c((1-w)/2,1-(1-w)/2),na.rm=T)
		# Divide the difference in quantiles by the first one
		# to get relative values
		as.numeric((y[2,]-y[1,])/(y[2,1]-y[1,1]))
	}
	rel.width2 <- function(w,x){
		# Get quantiles for fraction w of the data, symmetrically
		y <- apply(log.nz(x), 2, sd, na.rm=T)
		# Divide the difference in quantiles by the first one
		# to get relative values
		y/y[1]
	}
	sd.ratios <- rel.width(distr.width,na.omit(normexp[,1:4]))
	s <- seq(0.1,1,0.01)
	sdrlist <- as.data.frame(t(sapply(s, rel.width, x=normexp[,1:4])))
	colnames(sdrlist) <- colnames(normexp[,1:4])
	sdrlist$s <- s
	print(sd.ratios)
	if (file.out) dev.off()
}

if (fig.sma.ma.rma.slopes) {
	mar = c(4,4,1,1)
	if (file.out) dev.out("fig-exponent-supp", width=4, height=4, output.type=output.type)
	flds <- c('OLS','RMA','SMA','MA')
	slopes <- m2res[,flds]
	cols <- myrainbow(4)
	x <- apply(slopes, 2, ecdf)
	s <- seq(reg.lim[1],reg.lim[2],0.01)
	matplot(s, g <- sapply(x,function(m){m(s)}), type='l', lwd=2, lty='solid', col=myrainbow(4), las=1, ylab='Proportion', xlab='Regression slope')
	#x <- multi.ecdf(g <- m2res[,c('OLS','MA','SMA','RMA')], lwd=2, xlim=c(0,3), xlab='Slope (mRNA\u2013protein exponent)', plot=F)
	#s <- seq(0,3,0.1); matplot(s, sapply(x,function(m){m(s)}), type='l', lwd=2, col=myrainbow(4))
	abline(h=0.5, col='gray')
	abline(v=(meds <- apply(slopes, 2, median)), col=cols, lwd=2, lty='dotted')
	abline(v=medex, lwd=2, lty='dashed')
	text(meds, 0.5, label=flds, pos=3, col=cols)
	text(medex, 0.5, label='SCM', pos=3)
	if (file.out) dev.off()
}


if (fig.dynamic.range) {
#	d <- na.omit(data.frame(m=yres$raw$mrna.yassour.ypd0.1, p=yres$raw$prot.ghaem))#, sentinel=yres$raw$prot.lu))
	d <- na.omit(data.frame(m=yres$raw$mrna.holstege, p=yres$raw$prot.ghaem))#, sentinel=yres$raw$prot.lu))
	d <- data.frame(m=yres$est$mrna, p=yres$est$prot)#, sentinel=yres$raw$prot.lu))

	prop <- 0.95
	minp <- (1-prop)/2
	cat("Dynamic range for ", 100*prop, "% of the proteome\n", sep='')
	qm <- quantile(d$m, probs=c(minp, 1-minp), na.rm=T)
	qp <- quantile(d$p, probs=c(minp, 1-minp), na.rm=T)
	dyn.range.m <- qm[2]/qm[1]
	dyn.range.p <- qp[2]/qp[1]
	dyn.range <- log(dyn.range.p)/log(dyn.range.m)
	cat("mRNA =", qm[1], 'to', qm[2], 'mol./cell\n')
	cat("protein =", qp[1], 'to', qp[2], 'mol./cell\n')
	cat("N =", nrow(d), ": M =",format(dyn.range.m,dig=3),"/ P =",format(dyn.range.p,dig=3), ", total = ", format(dyn.range,digits=3), '\n')
	#g <- llmodel2(p~m, data=d)
	#print(g)

	#multidens(list(d$m, d$p), log=T, xlim=c(1e-2, 1e7), xlab='Cellular level (molecules/cell)', fill=T, col=tcol(myrainbow(2)), line.col='black')
}


if (fig.toy.model) {
	# Random number seed
	set.seed(115)

	# Number of genes
	n <- 5854
	# Exponent of empirical (evolved) relationship between steady-state
	# mRNA levels and translation rates
	gamma <- 0.56
	# Scaling factor, 1/time
	alpha <- 0.1
	# Degradation rate, 1/time
	delta <- 0.001
	# Standard deviation of mean-zero variation added to log mRNA levels to yield
	# unscaled log translation rates
	te.variation <- 1.1
	# Steady-state mRNA levels in molecules/cell (log-normal)
	# Mean and variance are equal to those of the SCM mean estimates
	log.m <- rnorm(n, mean=1.09, sd=1.25)
	m <- exp(log.m)
	# Translation rate -- add log-normal variation to, and scale, mRNA levels
	tau <- alpha*exp(log.m + rnorm(n,mean=0,sd=te.variation))^gamma
	# Steady-state protein levels in molecules/cell (log-normal)
	prot.variation <- 0.55
	p <- (tau/delta)*exp(log.m + rnorm(n, mean=0, sd=prot.variation))

	# Plot protein vs. mRNA
	plot(m, p, log='xy', pch=16, las=1, 
		xlab='mRNA level (mol./cell)', ylab='Protein level (mol./cell)')
	# Plot translation rate vs. mRNA
	plot(m, tau, log='xy', pch=16, las=1,
		xlab='mRNA level (mol./cell)', 
		ylab='Translation rate per mRNA (proteins/sec)')

	if (T) {

	# Use translational efficiency measurements without missing data
	y <- yres$te$data #subset(yres$te$data, mrna.n==3)
	trans.eff <- y$te

	# Plot the levels
	xlim <- c(1e-2,3e2)
	ylim <- c(1e-1,2e6)

	if (file.out) dev.out("fig.toy.model", width=8, height=8, output.type=output.type)
	split.screen(c(2,2))
	screen(1)
	par(mar=c(0,4,4,1))
	lplot(yres$est.sample$mrna, yres$est.sample$prot, col=tcol('black',0.3), xlim=xlim, ylim=ylim, xlab='mRNA level (mol./cell)', ylab='Protein level (mol./cell)')
	labline(g <- llmodel2(prot~mrna, data=yres$est.sample), x=xlim, lwd=2, lty='dotted', slope=1)
	screen(2)
	par(mar=c(0,4,4,1))
	lplot(m, p, xlim=xlim, ylim=ylim, col=tcol('black',0.3), xlab='mRNA level (mol./cell)', ylab='Protein level (mol./cell)')
	labline(g <- llmodel2(p~m), x=xlim, lwd=2, lty='dotted', slope=1)
	screen(3)
	par(mar=c(4,4,0,1))
	lplot(y$mrna, trans.eff, xlim=xlim, ylim=c(3e-3,3e1), col=tcol('black',0.3), 
		xlab='mRNA level (mol./cell)', 
		ylab='Normalized ribosome footprint counts per mRNA (AU)')
	screen(4)
	par(mar=c(4,4,0,1))
	lplot(m, tau, xlim=xlim, ylim=c(1e-3,1e1), col=tcol('black',0.3), 
		xlab='mRNA level (mol./cell)', 
		ylab='Translation rate per mRNA (proteins/mRNA/sec)')
	close.screen(all=TRUE)
	if (file.out) dev.off()
	
	cat("Toy protein vs. mRNA\n")
	print((g2t <- llmodel2(p~m))$regression.results[4,3])
	cat("SCM sample protein vs. mRNA\n")
	print((g2r <- llmodel2(prot~mrna, data=yres$est.sample))$regression.results[4,3])
	cat("Toy TE vs. mRNA\n")
	print((g2t.rd <- llmodel2(tau~m))$regression.results[4,3])
	cat("Real TE vs. SCM mRNA\n")
	print((g2r.rd <- llmodel2(trans.eff~yres$est$mrna))$regression.results[4,3])
	pcor(cortest(m, tau, meth='s', log=F))
	pcor(cortest(yres$est.sample$mrna, trans.eff,meth='s',log=F))
	pcor(cortest(m, p, meth='p', log=T))
	pcor(cortest(yres$est.sample$mrna, yres$est.sample$prot,meth='p',log=T))
	}
}

if (fig.new.vs.old.mrna) {
	d <- data.frame(yres$te$data[,c('mrna','mrna.median','mrna.n','prot')])
	ds <- subset(d, mrna.n>2)
	y <- d
	y$mrna.median <- y$mrna.median/median(y$mrna.median, na.rm=T)*median(y$mrna, na.rm=T)

	if (file.out) dev.out("fig_new_vs_old_mrna", width=6, height=3, output.type=output.type)
	split.screen(c(1,2))
	screen(1)
	par(mar=c(4,4,1,1))
	lplot(y$mrna, y$mrna.median, xlab='SCM mRNA', ylab='Recent RNA-seq mRNA', cex=0.6, col=tcol('black',0.3))
	abline(0,1,lwd=2, col='darkgray',lty='dotted')
	labline(g <- llmodel2(mrna.median~mrna, data=y), method='RMA', lwd=2, col='darkgray')
	screen(2)
	par(mar=c(4,4,1,1))
	multi.ecdf(list('SCM'=y$mrna, 'Recent'=y$mrna.median), xlab='mRNA level', log=T, lwd=2, legend.at=c(7e-3,1), legend.cex=0.6)
	close.screen(all=TRUE)
	if (file.out) dev.off()
	pstat(ks.test(y$mrna, y$mrna.median))
	pstat(wilcox.test(y$mrna, y$mrna.median))
}