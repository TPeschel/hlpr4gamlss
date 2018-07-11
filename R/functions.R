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


#' inverse.probability
#'
#' @name inverse.probability
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
#' (p<-as.data.frame(sapply(cnt<-c(.025,.5,.975),inverse.probability,f)))
inverse.probability <-
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
	            BCT = { gamlss.dist::qBCT( cent, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
	            BCTo = { gamlss.dist::qBCTo( cent, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
	            BCPE = { gamlss.dist::qBCPE( cent, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
	            BCPEo = { gamlss.dist::qBCPEo( cent, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
	            BCCG = { gamlss.dist::qBCCG( cent, prediction$mu, prediction$sigma, prediction$nu ) },
	            BCCGo = { gamlss.dist::qBCCGo( cent, prediction$mu, prediction$sigma, prediction$nu ) },
	            NO = { gamlss.dist::qNO( cent, prediction$mu, prediction$sigma ) },
	            PO = { gamlss.dist::qPO( cent, prediction$mu ) } ),
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
#' @name probability
#' @description not for direct use
#' @param x independent values
#' @param prediction a gamlss prediction
#' @param family the family of prediction, for the case it's not an attribute
#'
#' @return z-scores for
#' @export
#'
probability <-
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
			BCT = { gamlss.dist::pBCT( x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCTo = { gamlss.dist::pBCTo( x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCPE = { gamlss.dist::pBCPE( x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCPEo = { gamlss.dist::pBCPEo( x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCCG = { gamlss.dist::pBCCG( x, prediction$mu, prediction$sigma, prediction$nu ) },
			BCCGo = { gamlss.dist::pBCCGo( x, prediction$mu, prediction$sigma, prediction$nu ) },
			NO = { gamlss.dist::pNO( x, prediction$mu, prediction$sigma ) },
			PO = { gamlss.dist::pPO( x, prediction$mu ) } )
	}

#' density
#'
#' @name density
#' @description not for direct use
#' @param x independent values
#' @param prediction a gamlss prediction
#' @param family the family of prediction, for the case it's not an attribute
#'
#' @return z-scores for
#' @export
#'
density <-
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
			BCT = { gamlss.dist::dBCT( x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCTo = { gamlss.dist::dBCTo( x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCPE = { gamlss.dist::dBCPE( x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCPEo = { gamlss.dist::dBCPEo( x, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
			BCCG = { gamlss.dist::dBCCG( x, prediction$mu, prediction$sigma, prediction$nu ) },
			BCCGo = { gamlss.dist::dBCCGo( x, prediction$mu, prediction$sigma, prediction$nu ) },
			NO = { gamlss.dist::dNO( x, prediction$mu, prediction$sigma ) },
			PO = { gamlss.dist::dPO( x, prediction$mu ) } )
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
    function( y, x, families = c( "BCCG", "BCPE", "BCT", "BCCGo", "BCPEo", "BCTo" ), n.cyc = 30, refit = F ) {

    	d <-
    		data.frame( x = x, y = y )

    	d.lms <-
    		gamlss::lms(
    			y = y,
    			x = x,
    			data = d,
    			families = families,
    			n.cyc = n.cyc )

    	if( ! d.lms$converged && refit ) {

    		d.lms. <-
    			refit( d.lms )

    		a <-
    			setdiff( names( d.lms ), names( d.lms. ) )
    		
    		for( i in a ) {
    			d.lms.[[ a[ i ] ]] <-
    				d.lms[[ a[ i ] ]] }
    		
    		d.lms <-
    			d.lms.
    	}

    	d.lms
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
					inverse.probability,
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
