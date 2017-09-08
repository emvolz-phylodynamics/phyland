library(phydynR)
library(bbmle)

.gen.island.model <- function( demes, xF = 1e4 )
{
	m <- length(demes)
	function(theta, x0, t0, t1, res = 10 ){
		times <- seq( t0, t1, l = res )
		
		# note Y = Ne * xF , births calibrated so co rate is 1/Ne
		
		b <- matrix( 0., nrow = m, ncol = m )
		rownames(b) = colnames(b) <- demes
		for (k in 1:m)
		{
			b[k,k] <- theta[paste0('Ne',k)] * xF * xF / 2.
		}
		
		f <- lapply( 1:res, function(i) b )
		
		y <- setNames( theta[ paste0('Ne', 1:m)  ]  , demes) * unname(xF)
		Y <- lapply( 1:res, function(i) y)
		
		g <- matrix( 0., nrow =m, ncol = m )
		rownames(g) = colnames(g) <- demes
		for (k in 1:m) for (l in 1:m){
			if (k!=l){
				g[k,l] <-  y[l] * theta[ paste0( 'mu', l, k ) ]
			}
		}
		G <- lapply( 1:res, function(i) g )
		
		list( rev(times), f, G, Y )
	}	
}


.tips2states <- function(tips, delimiter='_', index=NULL, regex=NULL)
{
	if (!is.null(delimiter) & !is.null(regex))
	 stop('At most one of the options *delimiter* or *regex* should be supplied for determining demes of sampled lineages.')
	
	if(!is.null(delimiter)){
	  if(is.null(index)){
		tipstates <- sapply( strsplit( tips, delimiter, fixed=TRUE) , function(x) tail(x,1) )
	  } else {
	    tipstates <- sapply( strsplit( tips, delimiter, fixed=TRUE) , function(x) x[index] )
	  }
	} else if( !is.null(regex)){
		tipstates <- regmatches( tips, regexpr( regex, tips ))
		if (length( tipstates)!= length(tips)) 
		  stop('Error with regex: number of matches does not equal number of tips.')
	}
	
	demes <- unique( tipstates )
	m <- length(demes)
	ssts <- matrix( 0., nrow = length(tips), ncol = m)
	colnames(ssts) <- demes
	rownames(ssts) <- tips
	for (i in 1:length(tips)){
		ssts[i, which( demes==tipstates[i] )] <- 1. 
	}
	list( ssts = ssts, demes = demes,  m = m, tipstates = tipstates )
}

#' Maximum likelihood inference of migration rates and effective population size between demes
#'
#' @param tree A time-scaled phylogeny in ape::phylo format
#' @param delimiter A character string which will be used to infer the deme of sampling from tip labels
#' @param index The integer position within tip labels of the sample deme information. If not provided will use the last position, e.g. <accession number>_locationOfSampling
#' @param regex As an alternative to delimiter/index, the location of sampling can be inferred by matching a regular expression string
#' @param design A matrix describing which migration rates to estimate and which migration rates are assumed to be equal. Row and column names should correspond to names of demes. Unique integers should specify rates to estimate. If not provided, all rates will be estimated
#' method Optimisation method to be used by optim
#' quiet If TRUE, will suppress likelihood printing to stdout
#' @param ... Additional parameters passed to optim
#' @return A fitted model
phylandml <- function( tree, delimiter= '_', index= NULL, regex = NULL, design=NULL
 , method = 'BFGS'
 , quiet = FALSE
 , ... )
{
	minLL = -Inf
	tts <- .tips2states ( tree$tip.label, delimiter, index, regex  )
	n <- length(tree$tip.label) 
	sts <- setNames(  node.depth.edgelength( tree )[1:n] , tree$tip.label )
	ssts <- tts$ssts 
	m <- tts$m
	demes <- tts$demes
	bdt <- DatedTree( tree, sts, ssts, tol = 1.)
	
	#log transform all vars
	dm <- .gen.island.model( demes )
	
	require(bbmle) 
	if (is.null( design)){
		design <- matrix( 1:(m*m), nrow =m ,ncol =m)
		rownames(design) = colnames(design) <- demes
	}
	diag(design) <- NA 
	design <- design[demes,demes]
	
	migrate_names <- paste0('migrate', na.omit(as.vector(unique(design) ) ) )
	muNames2mignames <- c() 
	muNames2niceNames <- c()
	for (k in 1:m) for (l in 1:m){
		if (k!=l){
			 migname <- paste0('migrate', design[k,l] )
			 muname <- paste0('mu', k, l )
			 #nicename <- paste(demes[k], demes[l], sep='->') 
			 nicename <- paste(demes[l], demes[k], sep='->') 
			 muNames2mignames <- c(muNames2mignames
			   , setNames( migname , muname )
			 )
			 muNames2niceNames <- c(muNames2niceNames
			   , setNames( nicename, muname )
			 )
		}
	}
	NeNames <- setNames( paste0( 'Ne', 1:m) , demes )
	pnames <- c( NeNames
	 , migrate_names
	)
	np <- length( pnames)
	
	Ne0 <- bdt$maxHeight 
	theta0 <- as.list( setNames( rep(1,np) , pnames) )
	for ( n in NeNames ) theta0[[n]] <- log( Ne0 )
	for ( n in migrate_names ) theta0[[n]] <- log( bdt$maxHeight / 20  )
	
	t1 <- bdt$maxSampleTime
	t0 <- t1 - bdt$maxHeight-1
	
	muNames <- names( muNames2mignames )
	theta0_dm <- setNames( c(rep(Ne0, m) ,  rep(bdt$maxHeight , length(muNames) ) ) 
	 , c( NeNames, muNames )
	)
	of0 <- function(...) {
		theta1 <- exp( unlist( as.list( match.call() )[pnames] ) )
		theta0 <- theta0_dm 
		theta0[ NeNames ] <- unname( theta1[NeNames] )
		theta0[ names(muNames2mignames) ] <- unname( theta1[ muNames2mignames ]  )
		if (any(theta0 < 0)) return(minLL)
		tfgy <- dm(theta0, x0=NA, t0, t1 )
		o <- colik.pik.fgy(bdt, tfgy, timeOfOriginBoundaryCondition=TRUE, maxHeight=Inf, forgiveAgtY=1, AgtY_penalty=0, returnTree=FALSE, step_size_res=10)
		if (!quiet) print( c( o, theta0) )
		if (is.na(o)) return(-minLL)
		-max(o , minLL)
	}
	
	formals( of0 ) <- theta0
	
	mlefit <- bbmle::mle2(of0 , theta0, method = method, optimizer='optim')
	theta1 <- exp( coef(mlefit) )
	
	#nice names: 
	estnames2niceNames <- setNames(c(demes, (muNames2niceNames) )
	  , c( NeNames, muNames2mignames ))
	dmnames2niceNames <- setNames(c(demes, (muNames2niceNames) )
	  , c( NeNames, names(muNames2niceNames) ))
	theta2 <- theta0_dm 
	theta2[ NeNames ] <- unname( theta1[NeNames] )
	theta2[ names(muNames2mignames) ] <- unname( theta1[ muNames2mignames ]  )
	names( theta2 ) <- dmnames2niceNames[ names(theta2) ]
	
	vcov_logcoef <- vcov(mlefit )
	rownames( vcov_logcoef ) = colnames( vcov_logcoef ) <- dmnames2niceNames[ rownames( vcov_logcoef ) ]
	
	o <- list( 
	  coef = theta2
	  , logcoef = coef(mlefit)
	  , fit = mlefit 
	  , vcov_logcoef = vcov_logcoef 
	  , muNames2mignames
	  , design=design
	  , demes = demes
	  , dmnames2niceNames = dmnames2niceNames
	  , estnames2niceNames = estnames2niceNames
	  , env = environment()
	  , of0 = of0
	)
	class(o) <- 'phylandml'
	o
}



#' Profile confidence intervals for fitted phylogeographic model 
#'
#' @param fit A fited model of class *phylandml*
#' @param whichparm Name of parameter to be profiled 
#' @param guess_se Initial guess of standard error of parameter
#' @param np Integer number of spline points to be used for interpolation. Increase this value to get more accurate profile
#' @param cntrl_bnd Tuning parameter for number of log likelihood units to search around the optimum when generating profile
#' @return Likelihood profile
confint.phylandml <- function(fit, whichparm, guess_se, np=6,  ... ) #cntrl_bnd = 20,
{
	nicenames2estnames <- with(environment(fit$of0),  setNames(names(estnames2niceNames), estnames2niceNames) )
	eparm <- nicenames2estnames[whichparm]
	
	.get.val <- function(y){
		environment(fit$of0)[['y']] <- y 
		environment(fit$of0)$whichparm <- whichparm
		environment(fit$of0)$eparm <- eparm
		environment(fit$of0)$x <- fit$logcoef
		
		with( environment(fit$of0), 
		{
			x[eparm] <- unname(y)
			do.call( of0, as.list(x))
		})
	}
	.of1 <- function(y) .get.val(y) - (-logLik(fit$fit)) - cntrl_bnd
	
	# search up 
	#~ 	fub <- uniroot(.of1, lower = fit$logcoef[eparm] , upper = 5 * guess_se + fit$logcoef[eparm] )
	#~ 	if (is.na(fub$estim.prec)) stop('Upper bound could not be found. Try increasing guess_se.')
	#~ 	ub <- fub$root 	
	ub <- fit$logcoef[eparm] + guess_se
	
	# search down
	#~ 	flb <- uniroot(.of1, upper = fit$logcoef[eparm] , lower = -5 * guess_se + fit$logcoef[eparm] )
	#~ 	if (is.na(flb$estim.prec)) stop('Lower bound could not be found. Try increasing guess_se.')
	#~ 	lb <- flb$root 
	lb <- fit$logcoef[eparm] - guess_se 
	
	# profile 
	grid <- seq( lb, ub, l = np )
print(grid)
print(exp(grid))
	.of2 <- function(y, pval){
		environment(fit$of0)$y <- y
		environment(fit$of0)$pval <- pval
		environment(fit$of0)$eparm <- eparm
		environment(fit$of0)$x <- fit$logcoef
		with( environment(fit$of0), 
		{
			x[names(y)] <- unname(y)
			x[eparm] <- pval
			do.call( of0, as.list(x))
		})
	}
	ll <- sapply( grid, function(pval){
		-optim(par = fit$logcoef[-which(names(fit$logcoef)==eparm) ]
		 , fn = .of2
		 , method = 'BFGS'
		 , pval = pval
		)$value
	})
	
	grid <- c( grid, fit$logcoef[eparm] )
	ll <- c( ll, logLik( fit$fit))
	ll <- ll[ order( grid ) ]
	grid <- sort(grid)
	prfun <- approxfun( grid, ll, rule = 2)
	
	# search up again
	.of3 <- function(y) prfun(y)  - (logLik( fit$fit ) - 1.96 )
	if( tail( ll,1) > (max( ll) - 1.96 ) )
	  stop('Upper bound could not be found. Try increasing cntrl_bnd')
	f_ci_ub <- uniroot( .of3, lower = fit$logcoef[eparm], upper = ub )
	if (is.na( f_ci_ub$estim.prec)) stop('Upper bound could not be found. Try increasing guess_se')
	ci_ub <- f_ci_ub$root 
	# search down again 
	if( ll[1] > (max( ll) - 1.96 ) )
	  stop('Lower bound could not be found. Try increasing cntrl_bnd')
	f_ci_lb <- uniroot( .of3, lower = lb, upper = fit$logcoef[eparm] )
	if (is.na( f_ci_lb$estim.prec)) stop('Lower bound could not be found. Try increasing guess_se')
	ci_lb <- f_ci_lb$root
	
	ci <- data.frame( exp(c(ci_lb, ci_ub) ) )
	colnames(ci) <- whichparm
	rownames(ci) <- c('2.5%', '97.5%' )
	rv <- list( 
	  ci = ci
	  , grid = grid
	  , ll = ll
	)
	class(rv) <- c('profile', 'profile.phylandml' )
	rv
}

print.profile.phylandml <- function(x, ...){
	stopifnot( 'profile.phylandml' %in% class(x))
	print( x$ci)
	invisible(x) 
}



summary.phylandml <- function(x, ... )
{
	stopifnot( 'phylandml' %in% class(x))
	cat('Summary of log transformed parameters:\n')
	X <- attr( bbmle::summary( x$fit  ), 'coef' )
	rownames( X ) <- x$estnames2niceNames[ rownames(X) ] 
	
	print( X  )
	cat('\n')
	cat( 'Design matrix\n')
	print( x$design )
	cat('\n') 
	cat( 'Estimated effective sizes and rates:\n')
	X <- as.data.frame( x$coef )
	colnames(X) <- ''
	print(X)
	cat('\n')
	invisible( x )
}

coef.phylandml <- function( x, ... )
{
	stopifnot( 'phylandml' %in% class(x))
	x$coef 
}

print.phylandml <- function(x, ... )
{
	summary(x)
}


