#' best.match
#'
#' @name best.match
#' @description matches values to ones from a given nominal vector by their closest distance to it.
#' @param actual value or vector. the values to be matched
#' @param nominal value or vector. the values to be matched against
#'
#' @return vector. matched values. values are matched by their minimal absolut distance to its target.
#' @export
#'
#' @examples
#' plot( best.match( dnorm( .05 * -100:100 ), seq( 0, .5, by = .1 ) ) )
best.match <-
	function( actual, nominal ) {

		sapply(
			actual,
			function( a )
				nominal[ which.min( abs( a - nominal ) ) ] )
	}

#' which.best.match
#'
#' @name which.best.match
#' @description matches values to ones from a given nominal vector by their closest distance to it.
#' @param actual value or vector. the values to be matched
#' @param nominal value or vector. the values to be matched against
#'
#' @return vector of indices. indeices of matched values in nominal. values are matched by their minimal absolut distance to its target.
#' @export
#'
#' @examples
#' plot( which.best.match( dnorm( .05 * -100:100 ), seq( 0, 1, by = .1 ) ) )
which.best.match <-
	function( actual, nominal ) {

		sapply(
			actual,
			function( a )
				which.min( abs( a - nominal ) ) )
	}


#' new.Call
#' @name new.Call
#' @description  creates a new.Call object for using with do.Call
#' @param func function's name as string
#' @param ...  arguments of the function
#' @param canonical if TRUE a call with all neccessary arguments is called. Missing will be filled by args of the function's default call.
#'
#' @return a new.Call object
#' @export
#'
#' @examples
#' n.C <- new.Call( "print", "Hello World!" )
#' do.Call( n.C )
new.Call <-
	function( func = "Sys.time", ..., canonical = F, arg.len.cut.off = 25 ) {
		args.list <- list( ... )
		if( canonical ) {
			frmls           <- formals( func )
			frmls.names     <- names( frmls )
			args.list.names <- names( args.list )
			empty.names     <- c( )
			if( is.null( args.list.names ) && 0 < length( args.list ) )
				empty.names <- c( 1 : length( args.list ) )
			else
				empty.names <- which( args.list.names == "" )

			common.names                             <- setdiff( frmls.names, args.list.names )
			names( args.list )[ empty.names ]        <- common.names[ 1 : length( empty.names ) ]
			frmls[ names( args.list ) ]              <- args.list
			frmls.chars                              <- as.character( frmls )
			frmls.chars[ nchar( frmls.chars ) > arg.len.cut.off ] <- paste0( substr( frmls.chars[ nchar( frmls.chars ) > arg.len.cut.off ], 1, arg.len.cut.off - 5 ), " ··· " )
			return(
				list(
					FUNC = func,
					ARGS = frmls,
					CALL = paste0(
						func,
						"( ",
						paste( frmls.names, frmls.chars, sep = "=", collapse = ", " ),
						ifelse( length( frmls ) < 1, ")", " )" )
					)
				)
			)
		}
		list(
			FUNC = func,
			ARGS = args.list )
		#			ARGS = ifelse( is.null( args.list ), list( ), args.list ) )
	}

#' do.Call
#'
#' @name do.Call
#' @description  calls a function defined by new.Call
#' @param new.Call a new.Call object
#' @param verbose.doCall show whats going on beghind the scenes
#'
#' @return notting
#' @export
#'
#' @examples
#' n.C <- new.Call( "print", "Hello World!" )
#' do.Call( n.C )
do.Call <-
	function( new.Call = new.Call( canonical = T ), verbose.doCall = F ) {
		if( verbose.doCall == T  && ! is.null( new.Call$CALL ) )
			print( new.Call$CALL )
		if( is.null( new.Call$ARGS ) )
			do.call( new.Call$FUNC, list( ) )
		else
			do.call( new.Call$FUNC, new.Call$ARGS )
	}

#' dist.Call
#'
#' @name dist.Call
#' @description  a new.Call object for distribution's use
#' @param prefix a prefix like "", "d", "p", "q", "r"
#' @param family a gamlss family or "norm"
#' @param ...  arguments of the call
#' @param canonical if TRUE a call with all neccessary arguments is called. Missing will be filled by args of the function's default call.
#'
#' @return a new.Call object for calling with do.Call
#' @export
#'
#' @examples
#' d.C <- dist.Call( "r", "BCPE", 10, canonical = T )
#' do.Call( d.C, verbose.doCall = T )
dist.Call <-
	function( prefix = "p", family = "norm", ... ) {
		new.Call( paste0( prefix, family ), ... )
	}

#' gamlss.families
#'
#' @name gamlss.families
#' @description  lists all available gamlss families
#' @return a vector of gamlss family names
#' @export
#'
#' @examples
#' gamlss.families( )
gamlss.families <-
	function( ) {
		library( gamlss.dist )
		intersect(
			intersect(
				gsub( "^d", "", ls( "package:gamlss.dist" )[ grepl( "^d[A-Z]+", ls( "package:gamlss.dist" ) ) ] ),
				gsub( "^p", "", ls( "package:gamlss.dist" )[ grepl( "^p[A-Z]+", ls( "package:gamlss.dist" ) ) ] ) ),
			intersect(
				gsub( "^q", "", ls( "package:gamlss.dist" )[ grepl( "^q[A-Z]+", ls( "package:gamlss.dist" ) ) ] ),
				gsub( "^r", "", ls( "package:gamlss.dist" )[ grepl( "^r[A-Z]+", ls( "package:gamlss.dist" ) ) ] ) ) ) }

#' xdist
#'
#' @name xdist
#' @description  call a dist, ddist, pdist, qdist, rdist
#' @param x prefix "no", "d", "p", "q", "r"
#' @param family a gamlss family like BCTo
#' @param ... arguments for the constructed call
#' @param verbose tell what doing
#'
#' @return the expected data
#' @export
#'
#' @examples
#' xdist( "r", "NO", 10 )
xdist <-
	function ( x = "p", family = "NO", ..., verbose = F ) {
		if( verbose ) {
			print( paste0( "pdf of ", family, " with arguments:" ) )
			print( paste0( names( list( ... ) ),": ", list( ... ) ) ) }
		do.call( eval( parse( text =  paste0( "gamlss.dist::", x, family ) ) ), list( ... ) )
	}

#' ddist
#'
#' @name ddist
#' @description  call a density function of some gamlss family
#' @param family a gamlss family like BCTo
#' @param ... arguments for density function
#' @param verbose tell what doing
#'
#' @return the expected data
#' @export
#'
#' @examples
#' ddist( "NO", c( -1, 0, +1 ), mu = 0, sigma = 1 )
ddist <-
	function ( family = "NO", ..., verbose = F ) {
		if( verbose ) {
			print( paste0( "pdf of ", family, " with arguments:" ) )
			print( paste0( names( list( ... ) ),": ", list( ... ) ) ) }
		do.call( paste0( "d", family ), list( ... ) )
	}

#' pdist
#'
#' @name pdist
#' @description  call a pdf of some gamlss family
#' @param family a gamlss family like BCTo
#' @param ... arguments for the call
#' @param verbose tell what doing
#'
#' @return the expected data
#' @export
#'
#' @examples
#' pdist( "NO", c( -1, 0, +1 ), mu = 0, sigma = 1 )
pdist <-
	function ( family = "NO", ..., verbose = F ) {
		if( verbose ) {
			print( paste0( "cdf of ", family, " with arguments:" ) )
			print( paste0( names( list( ... ) ),": ", list( ... ) ) ) }
		do.call( paste0( "p", family ), list( ... ) )
	}

#' qdist
#'
#' @name qdist
#' @description  call a pdf of some gamlss family
#' @param family a gamlss family like BCTo
#' @param ... arguments for the call
#' @param verbose tell what doing
#'
#' @return the expected data
#' @export
#'
#' @examples
#' qdist( "NO", c( 0.1586553, 0.5000000, 0.8413447 ), mu = 0, sigma = 1 )
qdist <-
	function ( family = "NO", ..., verbose = F ) {
		if( verbose ) {
			print( paste0( "inv cdf of ", family, " with arguments:" ) )
			print( paste0( names( list( ... ) ),": ", list( ... ) ) ) }
		do.call( paste0( "q", family ), list( ... ) )
	}

#' rdist
#'
#' @name rdist
#' @description  create a random values of some gamlss family
#' @param family a gamlss family like BCTo
#' @param ... arguments for the call
#' @param verbose tell what doing
#'
#' @return the expected data
#' @export
#'
#' @examples
#' rdist( "NO", c( 0.1586553, 0.5, 0.8413447 ), mu = 0, sigma = 1 )
rdist <-
	function ( family = "NO", ..., verbose = F ) {
		if( verbose ) {
			print( paste0( "random values of ", family, " with arguments:" ) )
			print( paste0( names( list( ... ) ),": ", list( ... ) ) ) }
		do.call( paste0( "r", family ), list( ... ) )
	}


#' inverse.probability
#'
#' @name qDist
#' @description wrapper: is the quantile function for the respective family or prediction
#' @param cent value. centiles that should be calculated. 0.0 < cent & cent < 1.0
#' @param prediction a data.frame conraining distribution parameters (mu,sigma,nu,tau) for a certain gamlss model.
#' @param family the family of prediction, for the case it's not an attribute
#'
#' @return values relating to probabilities and prediction or family
#' @export
#'
#' @examples
#' plot(d<-data.frame(x=1:1000,y=rBCTo(1000,1000,100,10,1)+rnorm( 1000,0,.1*1:1000)+1:1000))
#' (m<-lms(y,x,data=d)) # somtimes BCPEo or BCTo
#' (f<-predictAll(m,data.frame(x=10*c(1:100))))
#' (p<-as.data.frame(sapply(cnt<-c(.025,.5,.975),qDist,f)))
qDist <-
    function( cent, prediction, family = NULL ) {

    	a <-
    		attr( prediction, "family" )[ 1 ]

    	fam <-
    		ifelse(
    			is.null( family ),
    			ifelse(
    				is.null( a ) || nchar( a ) < 1,
    				"NO",
    				a ),
    			family )

    	tryCatch(
    		switch(
            	as.character( fam ),
	            BCT   = { gamlss.dist::qBCT(   cent, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
	            BCTo  = { gamlss.dist::qBCTo(  cent, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
	            BCPE  = { gamlss.dist::qBCPE(  cent, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
	            BCPEo = { gamlss.dist::qBCPEo( cent, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
	            BCCG  = { gamlss.dist::qBCCG(  cent, prediction$mu, prediction$sigma, prediction$nu ) },
	            BCCGo = { gamlss.dist::qBCCGo( cent, prediction$mu, prediction$sigma, prediction$nu ) },
	            NO    = { gamlss.dist::qNO(    cent, prediction$mu, prediction$sigma ) },
	            PO    = { gamlss.dist::qPO(    cent, prediction$mu ) },
            	LNO   = { gamlss.dist::qLNO(   cent, prediction$mu, prediction$sigma, prediction$nu ) },
            	TF    = { gamlss.dist::qTF(    cent, prediction$mu, prediction$sigma, prediction$nu ) }
            ),
    		warning = function( msg ) {
    			message( paste0( "warning! something went wrong with ", fam ) )
    			message( msg )
    			return( NULL ) },
    		error = function( msg ) {
    			message( paste0( "error! something went wrong with ", fam ) )
    			message( msg )
    			return( NULL ) } )
    }

#' probability
#'
#' @name pDist
#' @description not for direct use
#' @param x independent values
#' @param prediction a gamlss prediction
#' @param family the family of prediction, for the case it's not an attribute
#'
#' @return values relating to probabilities and prediction or family
#' @export
#'
pDist <-
	function( x, prediction, family = NULL ) {

		a<-
			attr( prediction, "family" )[ 1 ]

		fam <-
			ifelse(
				is.null( family ),
				ifelse(
					is.null( a ),
					"NO",
					a ),
				family )

		switch(
			as.character( fam ),
			BCT   = { gamlss.dist::pBCT(   x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCTo  = { gamlss.dist::pBCTo(  x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCPE  = { gamlss.dist::pBCPE(  x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCPEo = { gamlss.dist::pBCPEo( x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCCG  = { gamlss.dist::pBCCG(  x, prediction$mu, prediction$sigma, prediction$nu ) },
			BCCGo = { gamlss.dist::pBCCGo( x, prediction$mu, prediction$sigma, prediction$nu ) },
			NO    = { gamlss.dist::pNO(    x, prediction$mu, prediction$sigma ) },
			PO    = { gamlss.dist::pPO(    x, prediction$mu ) },
			LNO   = { gamlss.dist::pLNO(   x, prediction$mu, prediction$sigma, prediction$nu ) },
			TF    = { gamlss.dist::pTF(    x, prediction$mu, prediction$sigma, prediction$nu ) }
		)
	}

#' density
#'
#' @name dDist
#' @description not for direct use
#' @param x independent values
#' @param prediction a gamlss prediction
#' @param family the family of prediction, for the case it's not an attribute
#'
#' @return values relating to probabilities and prediction or family
#' @export
#'
dDist <-
	function( x, prediction, family = NULL ) {

		a<-
			attr( prediction, "family" )[ 1 ]

		fam <-
			ifelse(
				is.null( family ),
				ifelse(
					is.null( a ),
					"NO",
					a ),
				family )

		switch(
			as.character( fam ),
			BCT   = { gamlss.dist::dBCT(   x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCTo  = { gamlss.dist::dBCTo(  x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCPE  = { gamlss.dist::dBCPE(  x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCPEo = { gamlss.dist::dBCPEo( x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCCG  = { gamlss.dist::dBCCG(  x, prediction$mu, prediction$sigma, prediction$nu ) },
			BCCGo = { gamlss.dist::dBCCGo( x, prediction$mu, prediction$sigma, prediction$nu ) },
			NO    = { gamlss.dist::dNO(    x, prediction$mu, prediction$sigma ) },
			PO    = { gamlss.dist::dPO(    x, prediction$mu ) },
			LNO   = { gamlss.dist::dLNO(   x, prediction$mu, prediction$sigma, prediction$nu ) },
			TF    = { gamlss.dist::dTF(    x, prediction$mu, prediction$sigma, prediction$nu ) }
		)
	}


#' random distribution
#'
#' @name rDist
#' @description not for direct use
#' @param x independent values
#' @param prediction a gamlss prediction
#' @param family the family of prediction, for the case it's not an attribute
#'
#' @return values relating to probabilities and prediction or family
#' @export
#'
rDist <-
	function( n = 1, mu = 1, sigma = 1, nu = 1, tau = 1, fam = "NO" ) {

		switch(
			as.character( fam ),
			BCT   = { gamlss.dist::rBCT(   n, mu, sigma, nu, tau ) },
			BCTo  = { gamlss.dist::rBCTo(  n, mu, sigma, nu, tau ) },
			BCPE  = { gamlss.dist::rBCPE(  n, mu, sigma, nu, tau ) },
			BCPEo = { gamlss.dist::rBCPEo( n, mu, sigma, nu, tau ) },
			BCCG  = { gamlss.dist::rBCCG(  n, mu, sigma, nu ) },
			BCCGo = { gamlss.dist::rBCCGo( n, mu, sigma, nu ) },
			NO    = { gamlss.dist::rNO(    n, mu, sigma ) },
			PO    = { gamlss.dist::rPO(    n, mu ) },
			LNO   = { gamlss.dist::rLNO(   n, mu, sigma, nu ) },
			TF    = { gamlss.dist::rTF(    n, mu, sigma, nu ) }
		)
	}

#' compute.model
#'
#' @name compute.model
#' @description computes the gamlss model
#' @param y y
#' @param x x
#' @param families families try to fit
#' @param n.cyc number of iterations
#' @param refit continue computation
#'
#' @return a gamlss model
#' @export
#'
compute.model <-
    function( y, x, families = c( "BCCG", "BCPE", "BCT", "BCCGo", "BCPEo", "BCTo", "LNO" ), n.cyc = 30, refit = F ) {

    	d <-
    		data.frame( x = x, y = y )

    	gamlss::lms(
    		y = y,
    		x = x,
    		data = d,
    		families = families,
    		n.cyc = n.cyc )
    }

#' compute.prediction
#'
#' @name compute.prediction
#' @description computes a prediction for a gamlss model and new x values
#' @param x the independent value
#' @param model a gamlss model
#'
#' @return a gamlss prediction of the model with new x values
#' @export
#'
compute.prediction <-
    function( x, model ) {

    	tryCatch(
    		gamlss::predictAll( model, data.frame( x = x ) ),
    		error = function( msg ) {
    			message( msg )
    			return( NA ) } )
    }

#' compute.percentiles
#'
#' @name compute.percentiles
#' @description computes percentiles for a gamlss prediction
#' @param cent the centiles
#' @param prediction a prediction
#' @param family the family of prediction, for the case it's not an attribute
#'
#' @return a data frame where every centile is a column along x
#' @export
#'
compute.percentiles <-
	function( cent = c( .025, .100, .500, .900, .975 ), prediction, family = NULL ) {

		l <-
			tryCatch(
				as.data.frame(
					lapply(
						cent,
						qDist,
						prediction,
						family ) ),
				error = function( msg ) {
					message( msg )
					return( NULL )
				} )

		if( ! is.null( l ) )

			names( l ) <-
				paste0( 100 * round( cent, 4 ), "%" )

		l
	}

#' compute.sds
#'
#' @name compute.sds
#' @description computes sds for a gamlss model
#' @param model a gamlss model
#' @param y the dependent variable
#' @param x the independent variable
#'
#' @return the residuals of the model
#' @export
#'
compute.sds <-
    function( model, y = model$y, x = model$xvar ) {

    	gamlss::z.scores( model, y, x )
    }

#' do.the.whole.thing
#'
#' @name do.the.whole.thing
#' @description does what it does
#' @param cent the centiles that should be calculated
#' @param y.col.name the dependent variable
#' @param x.col.name the independent variable
#' @param group.col.name group, the whole thing is done for every group seperately
#' @param data data frame with at least a dependent an indepenent and a group variable
#' @param x.pred x variable for prediction
#' @param n.cyc number of iterations for lms
#'
#' @return a list with everything
#' @export
#'
do.the.whole.thing <-
    function(
    	cent,
    	y.col.name = col.name,
    	x.col.name = "AGE",
    	group.col.name = "SEX",
    	data,
    	x.pred = NULL,
    	fam = c( "BCPE", "BCT", "BCCG", "BCPEo", "BCTo", "BCCGo" ),
    	n.cyc = 30, refit = F ) {

    	grps <-
            levels( data[ , group.col.name ] )

    	as.data.frame(
            sapply(
                grps,
                function( g ) {

                	d.g <-
                        data[ data[ , group.col.name ] == g, ]

                	d.g.mdl <-
                        compute.model( d.g[ , y.col.name ], d.g[ , x.col.name ], fam, n.cyc, refit = refit )

                	if( is.null( x.pred ) ) {

                		x.pred <-
                            d.g.mdl$xvar
                	}

                	d.g.prd <-
                        tryCatch(
                        	compute.prediction( x.pred, d.g.mdl ),
                        	error = function( msg ) {
                        		message( "Prediction could not be calcultated." )
                        		message( msg )
                        		return( NULL ) } )

                	if( ! is.null( d.g.prd ) ) {

                		d.g.prcntls <-
                			tryCatch(
                				compute.percentiles( cent, d.g.prd ),
                        		error = function( msg ) {
                        			message( "Percentiles could not be calcultated." )
                        			message( msg )
                        			return( NULL ) } ) }

                	l <-
                		list( model = NA, pred = NA, cent = NA, sds = NA )

                	if( ! is.null( d.g.prcntls ) ) {

                		d.g.prcntls[ , x.col.name ] <-
                			x.pred

                		d.g.sds <-
                			compute.sds( d.g.mdl )

                		l <-
                			list(
                        		model = d.g.mdl,
                        		pred  = d.g.prd,
                        		cent  = d.g.prcntls,
                        		sds   = d.g.sds ) }
                	l
                }
            )
        )
    }
