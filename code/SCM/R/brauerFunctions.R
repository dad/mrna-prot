
get_calibration_genes <- function (this_table, thresh_slope=1e-5,thresh_r_sq=1e-5) {

	calibration_list <- as.character(this_table$"ORF"[(this_table$"Boot.Pval.Slope"<thresh_slope) & (this_table$"Boot.Pval.R.Sq."<thresh_r_sq)])
	return(calibration_list)
}

build_table <- function (obs_lms,boot_pvals,file_out_table=NULL,DEBUG=FALSE) {
	table_s1 <- NULL
	table_s1 <- data.frame( as.character(obs_lms$orfs), obs_lms$beta, obs_lms$alpha, obs_lms$r_sq, sqrt(obs_lms$gene_var), boot_pvals )
	colnames(table_s1) <- c("ORF","GR.Slope","GR.Int.","R.Sq.","Std.Dev.","Boot.Pval.Slope","Boot.Pval.Int.","Boot.Pval.R.Sq.","Boot.Pval.Std.Dev.")
	if (DEBUG) {
		OFFSET <- 5
		table_s1[1:10,]
### DEBUG
		par(mfrow=c(2,2))
### DEBUG
		for (i in (OFFSET+1:4)) { plot(sort(table_s1[,i]),type="l",ylab="P-value (bootstrap)",xlab="Gene ID") }
### DEBUG
		#for (i in (OFFSET+1:4)) { hist(table_s1[,i],breaks=100,prob=TRUE,col="grey") }
### DEBUG
	}
	if (!is.null(file_out_table)) {
		write.table(table_s1, file=file_out_table, sep=",", row.names=FALSE, col.names=TRUE)
	} else {
		return(table_s1)
	}
}

compute_boot_pvals <- function (obs_lms,boot_lms) {
	N <- obs_lms$num_genes
	B <- boot_lms$num_boot_samples
	boot_pvals <- matrix(NA,nrow=N,ncol=4)
	prev <- progressBar()
	for (g in 1:N) {
		boot_pvals[g,1] <- sum(abs(boot_lms$beta)>abs(obs_lms$beta[g]))/B
		boot_pvals[g,2] <- sum(abs(boot_lms$alpha)>abs(obs_lms$alpha[g]))/B
		boot_pvals[g,3] <- sum(boot_lms$r_sq>obs_lms$r_sq[g])/B
		boot_pvals[g,4] <- sum(boot_lms$gene_var>obs_lms$gene_var[g])/B
		prev <- progressBar(g/N,prev)
	}
	return(boot_pvals)
}

run_bootstrap_lms <- function (expression,dilution_rate,col_sets=NULL,B=100000,with_cov=NULL,cov=NULL) {
	if (is.null(with_cov)) {
		b_exp <- resample_matrix_by_col(expression,B)
	} else {
		b_exp <- resample_matrix_by_col_sets(expression,B,col_sets)
	}
	if (is.null(with_cov)) {
		b_beta <- numeric(B)
		b_alpha <- numeric(B)
		b_r_sq <- numeric(B)
		b_gene_var <- numeric(B)
		prev <- progressBar()
		for (g in 1:B) {
			this_fit <- lm(as.formula("b_exp[g,] ~ dilution_rate+1"))
			b_beta[g] <- summary(this_fit)$coef[2]
			b_alpha[g] <- summary(this_fit)$coef[1]
			b_r_sq[g] <- summary(this_fit)$r.sq
			b_gene_var[g] <- var(b_exp[g,])
			prev <- progressBar(g/B,prev)
		}
	} else if (with_cov=="slope") {
		# ...
	} else if (with_cov=="intercept") {
		b_beta <- numeric(B)
		k <- length(levels(cov))
		b_alpha <- matrix(NA,nrow=B,ncol=k)
		b_r_sq <- numeric(B)
		b_gene_var <- numeric(B)
		prev <- progressBar()
		for (g in 1:B) {
			this_fit <- lm(as.formula("b_exp[g,] ~ dilution_rate+1"))
			b_beta[g] <- summary(this_fit)$coef[2]
			b_alpha[g,] <- summary(this_fit)$coef[c(1,3:(k+1))] # baseline, other levels (fitted; contrast.treatment) ### DEBUG
			b_r_sq[g] <- summary(this_fit)$r.sq
			b_gene_var[g] <- var(b_exp[g,])
			prev <- progressBar(g/B,prev)
		}
	}
	return(list(num_boot_samples=B, alpha=b_alpha, beta=b_beta, r_sq=b_r_sq, gene_var=b_gene_var))
}

run_lms <- function (expression,dilution_rate,with_cov=NULL,cov=NULL) {
	N <- dim(expression)[1]
	if (is.null(with_cov)) {
		beta <- numeric(N)
		alpha <- numeric(N)
		r_sq <- numeric(N)
		gene_var <- numeric(N)
		prev <- progressBar()
		for (g in 1:N) {
			this_fit <- lm(as.formula("as.numeric(expression[g,]) ~ dilution_rate+1"))
			alpha[g] <- summary(this_fit)$coef[1]
			beta[g] <- summary(this_fit)$coef[2]
			r_sq[g] <- summary(this_fit)$r.sq
			gene_var[g] <- var(as.numeric(expression[g,]),na.rm=TRUE)
			prev <- progressBar(g/N,prev)
		}
	} else if (with_cov=="slope") {
		# ...
	} else if (with_cov=="intercept") {
		beta <- numeric(N)
		k <- length(levels(cov))
		alpha <- matrix(NA,nrow=N,ncol=k)
		r_sq <- numeric(N)
		gene_var <- numeric(N)
		prev <- progressBar()
		for (g in 1:N) {
			this_fit <- lm(as.formula("as.numeric(expression[g,]) ~ dilution_rate+cov+1"))
			alpha[g,] <- summary(this_fit)$coef[c(1,3:(k+1))] # baseline, other levels (fitted; contrast.treatment) ### DEBUG
			beta[g] <- summary(this_fit)$coef[2]
			r_sq[g] <- summary(this_fit)$r.sq
			gene_var[g] <- var(as.numeric(expression[g,]),na.rm=TRUE)
			prev <- progressBar(g/N,prev)
		}
	}
	return(list(orfs=rownames(expression),num_genes=N,alpha=alpha,beta=beta,r_sq=r_sq,gene_var=gene_var))
}

resample_matrix_by_col <- function (x,B) {
	x_boot <- matrix(NA,nrow=B,ncol=dim(x)[2])
	for (c in 1:dim(x)[2]) {
		x_boot[,c] <- sample(x[,c],B,replace=TRUE)
	}
	return(x_boot)
}

resample_matrix_by_col_sets <- function (x,B,col_sets) {
	x_boot <- matrix(NA,nrow=B,ncol=dim(x)[2])
	for (s in 1:length(col_sets)) {
		pool_s <- as.vector(x[,col_sets[[s]]])
		M <- length(col_sets[[s]])
		x_boot[,c(col_sets[[s]])] <- matrix(sample(pool_s,B*M,replace=TRUE),ncol=M)
	}
	return(x_boot)
}

est_gr <- function (file=NULL, microarray=NULL, condition=1, gold_std=gold_std_brauer,
	calibration_list=calibration_orf, baseline_x=NULL, baseline_gr=NULL, verbose=FALSE,
	known_rate=NULL, file_out_stats=NULL) {

	## 0. open an output file connection (append mode) if provided
	if(!is.null(file_out_stats)) { zz <- file(file_out_stats, "a") }

	## 1. load expression from file, unless otherwise specified
	try( {microarray <- read.table(file, header=FALSE); microarray <- as.matrix(microarray)}, silent=TRUE )
	microarray_list <- rownames(microarray)
	if (verbose) {
		cat("Number of genes & conditions in input microarray:\n")
		print(dim(microarray))
	}

	if(!is.null(file_out_stats)) {
		cat("Evaluating the following conditions:\n ", file=zz)
		cat(colnames(microarray),"\n", sep="\t", file=zz) }

	if(!is.null(baseline_x) && ((baseline_x < 1) || (baseline_x > dim(microarray)[2]))) {
		baseline_x <- NULL
		baseline_gr <- NULL }

	## 3. pull out relevant genes and their expression data
	relevant_list <- intersect(intersect(calibration_list,gold_std[,1]),microarray_list)

	if(!is.null(file_out_stats)) {
		cat("\nNote: this data contains ",as.character(length(relevant_list))," genes out of ",
			as.character(length(calibration_list))," in the calibration list:\n ", sep="", file=zz)
		cat(calibration_list,"\n", sep=" ", file=zz) }
	if( length( relevant_list ) < 1 ) {
		print( "No yeast calibration genes found in input file!" )
		return( NULL ) }

	relevant_list_ids_gstd <- match(relevant_list,gold_std[,1])
	beta_hat <- gold_std[relevant_list_ids_gstd,2]
	alpha_hat <- gold_std[relevant_list_ids_gstd,3]
	relevant_list_ids_mrry <- match(relevant_list,microarray_list)
	#try( x <- microarray[relevant_list_ids_mrry,], silent=TRUE ) ## for whole matrix
	x <- matrix(microarray[relevant_list_ids_mrry,condition],ncol=1)
	relevant_list_mrry <- rownames(microarray)[relevant_list_ids_mrry]
	rownames(x) <- relevant_list_mrry
	#print(cbind( as.character(relevant_list), as.character(rownames(microarray)[relevant_list_ids_mrry]), as.character(gold_std[relevant_list_ids_gstd,1]) )) ## DEBUG
	not_na_ids <- !is.na(x)
	beta_hat <- beta_hat[not_na_ids]
	alpha_hat <- alpha_hat[not_na_ids]
	x <- x[not_na_ids,1]

	dAbsoluteSkew <- 0
	if(!(is.null(baseline_x) || is.null(baseline_gr))) {
		if(!is.null(file_out_stats)) {
			cat( "\nGiven absolute condition and rate: ", as.character( baseline_x ), ", ",
				as.character( baseline_gr ), sep="", file=zz )
			cat( "\nCalculating absolute offset...\n", file=zz ) }

		mdTmp <- matrix( microarray[,baseline_x] )
		rownames( mdTmp ) <- rownames( microarray )
		lsRatesBaseline <- est_gr( known_rate = baseline_gr, microarray = mdTmp, gold_std = gold_std,
			calibration_list = calibration_list, verbose = verbose )
		if(is.list(lsRatesBaseline)) {
			dAbsoluteSkew <- as.numeric( lsRatesBaseline$g_rate ) } }

	x_prime <- x
	alpha_hat_prime <- alpha_hat
	beta_hat_prime <- beta_hat

	## 5. estimate the growth rate & compute residuals
	## 5.1. initialize
	rate <- compute_rate(x_prime, alpha_hat_prime,beta_hat_prime, 0.0, dAbsoluteSkew, known_rate)
	baseline_shift <- compute_baseline_shift(x_prime, alpha_hat_prime,beta_hat_prime, rate,
		dAbsoluteSkew, known_rate)
	epsilon_rate <- 1.0
	epsilon_baseline_shift <- 1.0
	tolerance <- 1e-7
	iter <- 0
	## 5.2. iterate till convergence
	while (epsilon_rate>tolerance && epsilon_baseline_shift>tolerance) {
		iter <- iter +1
		#print(iter); ## DEBUG
		#print(epsilon_rate); print(epsilon_baseline_shift); print(tolerance) ## DEBUG

		# it would be nice to do this every iteration, but that can prevent convergence
		if(iter == 1) {
			if(is.null(known_rate)) {
				vdPred <- baseline_shift + alpha_hat + beta_hat * (dAbsoluteSkew + rate) }
			else {
				vdPred <- baseline_shift + alpha_hat + beta_hat * (known_rate + rate) }
			lsFit <- lm( as.formula("vdPred ~ x" ))
			errors <- lsFit$residuals
			dLQ <- quantile( errors, probs = 0.25 )
			dUQ <- quantile( errors, probs = 0.75 )
			dIQR <- dUQ - dLQ
			dFence <- 1.5
			dLF <- dLQ - ( dFence * dIQR )
			dUF <- dUQ + ( dFence * dIQR )
			prime_ids <- which((errors > dLF) & (errors < dUF))
			x_prime <- matrix(x[prime_ids], ncol=1)
			names(x_prime) <- names(x[prime_ids])
			alpha_hat_prime <- alpha_hat[prime_ids]
			beta_hat_prime <- beta_hat[prime_ids] }

		## 5.2.1. update estimate of rate
		rate_old <- rate
		rate <- compute_rate(x_prime, alpha_hat_prime,beta_hat_prime, baseline_shift, dAbsoluteSkew,
			known_rate)
		## 5.2.2. update estimate of baseline_shift
		baseline_shift_old <- baseline_shift
		baseline_shift <- compute_baseline_shift(x_prime, alpha_hat_prime,beta_hat_prime, rate,
			dAbsoluteSkew, known_rate)
		## 5.2.3. check convergence
		epsilon_rate <- abs(rate - rate_old)
		epsilon_baseline_shift <- abs(baseline_shift - baseline_shift_old)
		if(is.na(epsilon_rate) || is.na(epsilon_baseline_shift)) {
			return( NULL ) }
	}
	## 5.3. compute residuals
	if(is.null(known_rate)) {
		errors <- x_prime - (baseline_shift + alpha_hat_prime + beta_hat_prime * (dAbsoluteSkew + rate)) }
	else {
		errors <- x_prime - (baseline_shift + alpha_hat_prime + beta_hat_prime * (known_rate + rate)) }

	if(!is.null(file_out_stats)) {
		cat( "\nTotal iterations: ", as.character( iter ), sep="", file=zz )
		cat( "\nPredicted rate and baseline: ", as.character( rate ), ", ",
			as.character( baseline_shift ), sep="", file=zz )
		cat( "\nRetained calibration genes:\n", file=zz )
		cat( names(x_prime), "\n", sep="\t", file=zz )
		cat( "\nResiduals:\n", file=zz )
		cat( errors, "\n", sep="\t", file=zz ) }

	if(!is.null(file_out_stats)) { close(zz) }

	## 6. return results
	return(list( g_rate=rate, b_shift=baseline_shift, orfs=names(x_prime), res=errors ))
}

compute_rate <- function (x, alpha_hat,beta_hat, baseline_shift, abs_skew = 0, known_rate = NULL) {
	num_genes <- length(alpha_hat)
	rate_i <- numeric(num_genes)
	for (i in 1:num_genes) {
		if(is.null(known_rate)) {
			rate_i[i] <- (x[i] - baseline_shift - alpha_hat[i]) / beta_hat[i] - abs_skew }
		else {
			rate_i[i] <- (x[i] - baseline_shift - alpha_hat[i]) / beta_hat[i] - known_rate }
	}
	#print(rate_i); print(median(rate_i)) ## DEBUG
	return( mean(rate_i) )
}

compute_baseline_shift <- function (x, alpha_hat,beta_hat, rate, abs_skew = 1, known_rate = NULL) {
	num_genes <- length(alpha_hat)
	baseline_shift_i <- numeric(num_genes)
	for (i in 1:num_genes) {
		if(is.null(known_rate)) {
			baseline_shift_i[i] <- x[i] - alpha_hat[i] - beta_hat[i] * (rate + abs_skew) }
		else {
			baseline_shift_i[i] <- x[i] - alpha_hat[i] - beta_hat[i] * (known_rate + rate) }
	}
	return( mean(baseline_shift_i) )
}

progressBar <- function (frac=0, prev=NULL) {
#	print( frac )
	return( prev )
}

cverror <- function (mdExpression, vdGrowthRates, frmeColSets, iN = 100) {

	lsResults = list( c(), c(), c(), c(), c(), c() )
	for( i in 1:iN ) {
		print( i )
		vdSubsample = runif( dim( mdExpression )[2] )
		viTrain = which( vdSubsample <= 0.66 )
		viTest = which( vdSubsample > 0.66 )
		mdTrain <- matrix( mdExpression[,viTrain], ncol = length( viTrain ) )
		rownames( mdTrain ) <- rownames( mdExpression )

		lsReal <- run_lms( mdTrain, vdGrowthRates[viTrain] )
		lsBootstrap <- run_bootstrap_lms( mdTrain, vdGrowthRates[viTrain], frmeColSets[viTrain] )
		mdBootstrapPValues <- compute_boot_pvals( lsReal, lsBootstrap )
		frmeGRParameters <- build_table( lsReal, mdBootstrapPValues )
		lsCalibration <- get_calibration_genes( frmeGRParameters )
		if( length( lsCalibration ) == 0 ) {
			next }
		for( j in viTest ) {
			mdTest <- matrix( mdExpression[,j], ncol = 1 )
			rownames( mdTest ) <- rownames( mdExpression )
			lsGR <- est_gr( microarray = mdTest, gold_std = frmeGRParameters,
				calibration_list = lsCalibration )
			if( is.na( lsGR ) || ( length( lsGR ) == 0 ) || !is.numeric( lsGR[[1]] ) ) {
				next }
			k <- ( ( j - 1 ) %% 6 ) + 1
			lsResults[[k]] <- append( lsResults[[k]], as.numeric( lsGR[[1]] ) ) } }

	return( lsResults )
}

calculateRates <- function( mdData, frmeGRParameters, lsCalibration, strFileLog = NULL,
	iAbsolute = FALSE, dAbsolute = FALSE ) {

	if( is.na( iAbsolute ) || ( iAbsolute == FALSE ) || ( iAbsolute < 1 ) ) {
		iAbsolute <- NULL
		dAbsolute <- NULL }

	iConditions <- dim( mdData )[2]
        vdRates <- numeric( iConditions )
	for( i in 1:iConditions ) {
		lsRate = est_gr( microarray = mdData, condition = i, gold_std = frmeGRParameters,
			calibration_list = lsCalibration, baseline_x = iAbsolute, baseline_gr = dAbsolute,
			file_out_stats = strFileLog )
		if( is.null( lsRate ) ) {
			print( "Rate inference failed" )
			return( NULL ) }
		vdRates[i] = lsRate$g_rate}

	return( as.numeric( vdRates ) )
}

calculateGrowthStats <- function( mdData, frmeGRParameters, lsCalibration, strFileLog = NULL,
	iAbsolute = FALSE, dAbsolute = FALSE ) {

	if( is.na( iAbsolute ) || ( iAbsolute == FALSE ) || ( iAbsolute < 1 ) ) {
		iAbsolute <- NULL
		dAbsolute <- NULL }

	iConditions <- dim( mdData )[2]
        growthStats <- list()
        growthStats$vdRates <- numeric(iConditions)
        growthStats$lenRes <- numeric(iConditions)
        growthStats$residuals <- list()
	for( i in 1:iConditions ) {
		lsRate = est_gr( microarray = mdData, condition = i, gold_std = frmeGRParameters,
			calibration_list = lsCalibration, baseline_x = iAbsolute, baseline_gr = dAbsolute,
			file_out_stats = strFileLog )
		if( is.null( lsRate ) ) {
			print( "Rate inference failed" )
			return( NULL ) }
                expName <- colnames(mdData)[i]
		growthStats$vdRates[i] = lsRate$g_rate
                resi <- as.numeric(lsRate$res)
                growthStats$lenRes[i] <- length(resi)
                names(resi) <- lsRate$orfs
                growthStats$residuals[[expName]] = resi
              }
        names(growthStats$vdRates) <- colnames(mdData)
        names(growthStats$lenRes) <- colnames(mdData)
	return( growthStats )
}

plotRates <- function( mdData, frmeGRParameters, lsCalibration, strFileRates = NULL,
	strFileLog = NULL, iAbsolute = FALSE, dAbsolute = FALSE ) {

	vdRates <- calculateRates( mdData, frmeGRParameters, lsCalibration, strFileLog, iAbsolute,
		dAbsolute )
	if( is.null( vdRates ) ) {
		return( FALSE ) }

	if( !is.na( iAbsolute ) && !is.na( dAbsolute ) && iAbsolute && dAbsolute ) {
		strType <- "Absolute" }
	else {
		strType <- "Relative" }
	strTitle <- paste( "Inferred", strType, "Growth Rates" )

	if( !is.null( strFileRates ) ) {
		mdOut <- matrix( c(colnames( mdData ), vdRates), ncol = 2 )
		colnames( mdOut ) <- c("Condition", strTitle)
		write.table( mdOut, file = strFileRates, quote = FALSE, sep = "\t", row.names = FALSE ) }

	plot( vdRates, type = "b", main = strTitle, ylab = "Growth Rate", xlab = "Conditions",
		xaxt = "n" )
	axis( 1, labels = colnames( mdData ), at = 1:length( vdRates ) )
	return( TRUE )
}

