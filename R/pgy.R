library(phydynR)
library(bbmle)


.gen.island.model <- function( demes, xF = 1e4 )
{
	m <- length(demes)
	function(theta, x0, t0, t1, res = 10, integrationMethod=NA ){
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
		
		deaths <- lapply( 1:res, function(i) setNames(rep(0, m), demes) )
		
		list( rev(times), f, G, Y , NA, deaths=deaths)
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
#' @param method Optimisation method to be used by optim
#' @param quiet If TRUE, will suppress likelihood printing to stdout
#' @param start_migrate Optional starting value for optimisation of migration rates
#' @param start_Ne Optional starting value for optimisation of Ne 
#' @param Ne_logprior Optional log prior density for Ne
#' @param migrate_logprior Optional log prior density for migration rate
#' @param nstarts integer number of starting conditions to generate for optimisation
#' @param ncpu On multi-cpu computers, will conduct optimisation from different starting conditions in parallel on this many cpu's
#' @param minBL Any short branch lengths will be rounded up to this value
#' @param ... Additional parameters passed to optim
#' @return A fitted model with summary and coef methods.
phylandml <- function( tree, delimiter= '_', index= NULL, regex = NULL, design=NULL
 , method = 'BFGS'
 , quiet = FALSE
 , start_migrate = NA
 , start_Ne = NA
 , Ne_logprior = function(x) 0
 , migrate_logprior = function(x) 0
 , nstarts = 1
 , ncpu = 1
 , minBL = 0.
 , ... )
{
	require(phydynR)
	minLL = -Inf
	tts <- .tips2states ( tree$tip.label, delimiter, index, regex  )
	n <- length(tree$tip.label) 
	sts <- setNames(  node.depth.edgelength( tree )[1:n] , tree$tip.label )
	ssts <- tts$ssts 
	m <- tts$m
	demes <- tts$demes
	
	tree$edge.length <- pmax( tree$edge.length, minBL)
	if ( min ( tree$edge.length ) <= 0){
		stop('The provided tree has branch lengths < or == 0. Try setting the minBL option to a small value > 0.')
	}
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
	
	Ne0 <- ifelse( is.na(start_Ne), bdt$maxHeight / (2*length(demes) ), unname(start_Ne) ) 
	start_Ne <- Ne0 
	start_migrate <- ifelse( is.na( start_migrate), (bdt$maxHeight / 100) / length(demes)  , unname(start_migrate) )
	theta0 <- as.list( setNames( rep(1,np) , pnames) )
	for ( n in NeNames ) theta0[[n]] <- log( Ne0 )
	for ( n in migrate_names ) theta0[[n]] <- log( start_migrate )
	theta0s <- lapply( 1:nstarts, function(x) theta0)
	if (nstarts > 1 ){
		require(lhs)
		lmat <- improvedLHS( n = nstarts, k = np )
		Ne_lb <- start_Ne / 5
		Ne_ub <- start_Ne * 5
		mr_lb <- start_migrate / 5
		mr_ub <- start_migrate * 5
		theta_lb <- setNames( rep(1,np) , pnames)
		theta_lb[NeNames] <- Ne_lb
		theta_lb[migrate_names] <- mr_lb
		theta_ub <- setNames( rep(1,np) , pnames)
		theta_ub[NeNames] <- Ne_ub
		theta_ub[migrate_names] <- mr_ub
		theta0s <- lapply( 1:nstarts , function(i){
			x <- theta_lb + lmat[i,] * (theta_ub - theta_lb )
			as.list( log( x))
		})
	}
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
		lp <- sum( sapply( NeNames, function(x) Ne_logprior( theta1[x] )) ) + 
		  sum( sapply( muNames2mignames , function(x) migrate_logprior(theta1[x] ) ) )
		if (any(theta0 < 0)) return(minLL)
		tfgy <- dm(theta0, x0=NA, t0, t1 )
		suppressWarnings( {
			o <- colik.pik.fgy(bdt, tfgy, timeOfOriginBoundaryCondition=TRUE, maxHeight=Inf, forgiveAgtY=1, AgtY_penalty=0, returnTree=FALSE, step_size_res=10)
		})
		if (!quiet) print( c( o, lp, theta0) )
		if (is.na(o)) return(-minLL)
		-max(o + lp , minLL)
	}
	
	.trwithstates <- function(...)
	{
		theta1 <- exp( unlist( as.list( match.call() )[pnames] ) )
		theta0 <- theta0_dm 
		theta0[ NeNames ] <- unname( theta1[NeNames] )
		theta0[ names(muNames2mignames) ] <- unname( theta1[ muNames2mignames ]  )
		if (any(theta0 < 0)) return(minLL)
		#tfgy <- dm(theta0, x0=NA, t0, t1 )
		#colik.pik.fgy(bdt, tfgy, timeOfOriginBoundaryCondition=TRUE, maxHeight=Inf, forgiveAgtY=1, AgtY_penalty=0, returnTree=TRUE, step_size_res=10)$tree	
		phydynR::ace(bdt, theta0, demographic.process.model=dm, x0=NA, t0, res = 1e3
		  , integrationMethod='lsoda'
		  , timeOfOriginBoundaryCondition = TRUE
		  , maxHeight = Inf 
		  , forgiveAgtY = 1 #can be NA; if 0 returns -Inf if A > Y; if 1, allows A>Y everywhere
		  , AgtY_penalty = 0 # penalises likelihood if A > Y
		  , step_size_res = 10 # for adaptive ode solver, set to value < 1
		)

	}
	ace0 <- function(tr){
		#note indexed by node number, gives lineage ancestral to node
		ace.funcs <- lapply( 1:ncol( tr$acestates ), function(i){
			function(h) {
				ei <- which( tr$edge[,2] == i )
				if (length( ei )==0 ){
					# probably the root node 
					return( setNames( tr$acestates[,i] , colnames(tr$sampleStates )) )
				}
				u <- tr$parent[i ] #tr$edge[i, 1]
				v <- i
				h0 <- tr$heights[v]
				h1 <- tr$heights[u] 
				if ( h > h1 | h < h0 ) stop(paste( 'Lineage is not extant at given time before most recent sample. Lineage is extant: ', h0, ',', h1 , ' time units before most recent sample.' ) )
				sapply( 1:tr$m, function(k) {
					approx( 
					 c(h0, h1 )
					 , c(tr$acestates[k,i], tr$mstates[k,i]) # TODO should probably have a 'upper' acestates
					 , xout = h)$y
				}) -> x 
				setNames( x / sum(x) , colnames(tr$sampleStates) )
			}
		})
	}
	
	formals( of0 ) <- theta0
	formals( .trwithstates ) <- theta0
	
	if (nstarts > 1){
		# multiple fits, return best 
		if (ncpu > 1 ){
			require(parallel)
			mlefits <- mclapply( theta0s, function(x){
				tryCatch( bbmle::mle2(of0 , x, method = method, optimizer='optim', skip.hessian=TRUE, ...)
				 , error = function(e) NULL)
			}, mc.cores = ncpu)
		} else{
			mlefits <- lapply( theta0s, function(x){
				tryCatch( bbmle::mle2(of0 , x, method = method, optimizer='optim', skip.hessian=TRUE, ...)
				 , error = function(e) NULL )
			})
		}
		mlefits <- mlefits [ order(sapply( mlefits, function(f) ifelse(is.null(f), -Inf, bbmle::logLik(f)) ), decreasing=TRUE) ]
		mlefit <- mlefits[[1]] 
	} else{
		mlefit <- bbmle::mle2(of0 , theta0, method = method, optimizer='optim', skip.hessian=TRUE, ...)
		mlefits <- list( mlefit )
	}
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
	
	# state est 
	bdt <- do.call(.trwithstates, as.list( coef(mlefit )) )
	ace.funcs <- NULL
	acetab <- NULL
	{
		ace.funcs <- NA #ace0( bdt )
		acetab <- as.data.frame( t( bdt$acestates ) )
		colnames( acetab) <- colnames(bdt$sampleStates)
	}
	
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
	  , ace.funcs = ace.funcs 
	  , ace = acetab 
	  , bdt = bdt 
	  , mlefits = mlefits
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
#' @param ncpu Number of CPUs to use on multi-core machines
#' @return Likelihood profile
confint.phylandml <- function(fit, whichparm, guess_se, np=13, ncpu = 1,  ... ) #cntrl_bnd = 20,
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
	if (ncpu  > 1 ){
		require(parallel)
		ll <- unlist( mclapply( grid, function(pval){
			tryCatch( -optim(par = fit$logcoef[-which(names(fit$logcoef)==eparm) ]
			 , fn = .of2
			 , method = 'BFGS'
			 , pval = pval
			)$value , error = function(e) NA )
		} , mc.cores = round(ncpu) ) )
	} else{
		ll <- sapply( grid, function(pval){
			tryCatch( -optim(par = fit$logcoef[-which(names(fit$logcoef)==eparm) ]
			 , fn = .of2
			 , method = 'BFGS'
			 , pval = pval
			)$value, error = function(e) NA )
		})
	}
	
	grid <- c( grid, fit$logcoef[eparm] )
	ll <- c( ll, logLik( fit$fit))
	ll <- ll[ order( grid ) ]
	grid <- sort(grid)
	#prfun <- approxfun( grid, ll, rule = 2)
	
	sumdat <- data.frame( Parameter.value = exp(grid), Profile.log.likelihood = ll )
	print(sumdat)
	
	require(akima)
	adat <- suppressWarnings( aspline( grid, ll ) )
	prfun <- approxfun( adat$x, adat$y, rule=2)
	
	# search up again
	success <- TRUE
	.of3 <- function(y) prfun(y)  - (logLik( fit$fit ) - 1.96 )
	if( tail( ll,1) > (max( ll) - 1.96 ) ) {
	  cat("
*** Upper bound could not be found. Likelihood surface may be flat. Inspect $grid and $ll. Try increasing guess_se  ***
	  ")
	  success <- FALSE 
	}
	f_ci_ub <- tryCatch( uniroot( .of3, lower = fit$logcoef[eparm], upper = ub )
	 , error = function(e) list(root=NA, estim.prec = NA) )
	if (is.na( f_ci_ub$estim.prec)) {
		cat("
*** Upper bound could not be found. Likelihood surface may be flat. Inspect $grid and $ll. Try increasing guess_se  ***
	  ")
	  success <- FALSE 
	}
	ci_ub <- ifelse(success, f_ci_ub$root , NA )
	
	# search down again 
	if( ll[1] > (max( ll) - 1.96 ) ){
	  cat("
*** Lower bound could not be found. Rate may be indistinguishable from zero. Inspect $grid and $ll. Try increasing guess_se  ***
	  ")
	  success <- FALSE 
	}
	f_ci_lb <- tryCatch( uniroot( .of3, lower = lb, upper = fit$logcoef[eparm] )
	  , error = function(e) list(root=NA, estim.prec= NA))
	if (is.na( f_ci_lb$estim.prec)) {
		cat("
*** Lower bound could not be found. Rate may be indistinguishable from zero. Inspect $grid and $ll. Try increasing guess_se  ***
	   ")
	   success <- FALSE
	}
	ci_lb <- ifelse(success, f_ci_lb$root, NA )
	
	ci <- data.frame( exp(c(ci_lb, ci_ub) ) )
	colnames(ci) <- whichparm
	rownames(ci) <- c('2.5%', '97.5%' )
	rv <- list( 
	  ci = ci
	  , grid = grid
	  , ll = ll
	  , sumdat = sumdat 
	)
	class(rv) <- c('profile', 'profile.phylandml' )
	rv
}

plot.confint.phylandml <- function(x, ...){
	stopifnot( 'profile.phylandml' %in% class( x) )
	with( x$sumdat, plot( Parameter.value, Profile.log.likelihood, 'l', ... ))
	with( x$sumdat, points( Parameter.value, Profile.log.likelihood))
	invisible(x)
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


