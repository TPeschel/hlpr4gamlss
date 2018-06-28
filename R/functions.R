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


#' calc.quant
#'
#' @name calc.quant
#' @description not for direct use
#' @param cent value. centiles that should be calculated. 0.0 < cent & cent < 1.0
#' @param prediction a data.frame conraining distribution parameters (mu,sigma,nu,tau) for a certain gamlss model.
#' @param family the family of prediction, for the case it's not an attribute
#'
#' @return gives the values of a gamlls model prediction for a certain centile
#' @export
#'
#' @examples
#' plot(d<-data.frame(x=1:1000,y=rBCTo(1000,1000,100,10,1)+rnorm( 1000,0,.1*1:1000)+1:1000))
#' (m<-lms(y,x,data=d)) # somtimes BCPEo or BCTo
#' (f<-predictAll(m,data.frame(x=10*c(1:100))))
#' (p<-as.data.frame(sapply(cnt<-c(.025,.5,.975),calc.quant,f)))
calc.quant <-
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

    	print( fam )

    	switch(
            as.character( fam ),
            BCT = { gamlss.dist::qBCT( cent, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
            BCTo = { gamlss.dist::qBCTo( cent, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
            BCPE = { gamlss.dist::qBCPE( cent, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
            BCPEo = { gamlss.dist::qBCPEo( cent, prediction$mu, prediction$sigma, prediction$nu, prediction$tau ) },
            BCCG = { gamlss.dist::qBCCG( cent, prediction$mu, prediction$sigma, prediction$nu ) },
            BCCGo = { gamlss.dist::qBCCGo( cent, prediction$mu, prediction$sigma, prediction$nu ) },
            NO = { gamlss.dist::qNO( cent, prediction$mu, prediction$sigma ) },
            PO = { gamlss.dist::qPO( cent, prediction$mu ) }
        )
    }

#' calc.prob
#'
#' @name calc.prob
#' @description not for direct use
#' @param x
#' @param prediction
#' @param family the family of prediction, for the case it's not an attribute
#'
#' @return
#' @export
#'
#' @examples
calc.prob <-
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
            PO = { gamlss.dist::pPO( x, prediction$mu ) } ) }

#' compute.model
#'
#' @name compute.model
#' @description computes the gamlss model
#' @param y y
#' @param x x
#' @param families families try to fit
#'
#' @return a gamlss model
#' @export
#'
#' @examples
compute.model <-
    function( y, x, families = c( "BCCG", "BCPE", "BCT", "BCCGo", "BCPEo", "BCTo" ) ) {
        gamlss::lms( y = y, x = x, data = data.frame( x, y ), families = families )
    }

#' compute.prediction
#'
#' @name compute.prediction
#' @description computes a prediction for a gamlss model and new x values
#' @param x
#' @param model
#'
#' @return
#' @export
#'
#' @examples
compute.prediction <-
    function( x, model ) {
        gamlss::predictAll( model, data.frame( x = x ) )
    }

#' compute.percentiles
#'
#' @name compute.percentiles
#' @description computes percentiles for a gamlss prediction
#' @param cent the centiles
#' @param prediction a prediction
#' @param family the family of prediction, for the case it's not an attribute
#'
#' @return
#' @export
#'
#' @examples
compute.percentiles <-
    function( cent = c( .025, .100, .500, .900, .975 ), prediction, family = NULL ) {
        l <-
            as.data.frame(
                lapply(
                    cent,
                    calc.quant,
                    prediction,
                    family ) )
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
#' @return
#' @export
#'
#' @examples
compute.sds <-
    function( model, y = model$y, x = model$xvar ) {
        gamlss::z.scores( model, y, x )
    }

#' compute.references
#'
#' @name compute.references
#' @description computes reference values for a gamlss prediction
#' @param cent
#' @param prediction a gamlss prediction
#' @param family the family of prediction, for the case it's not an attribute
#'
#' @return
#' @export
#'
#' @examples
compute.references <-
    function( cent, prediction, family = NULL ) {
        s <-
            as.data.frame(
                sapply(
                    cent,
                    calc.quant,
                    prediction,
                    family ) )
        names( s ) <-
            paste0( round( 100 * cent, 2 ), "%" )
        s
    }

#' do.the.whole.thing
#'
#' @name do.the.whole.thing
#' @description does what it does
#' @param cent
#' @param y.col.name
#' @param x.col.name
#' @param group.col.name
#' @param data
#' @param x.pred
#'
#' @return
#' @export
#'
#' @examples
do.the.whole.thing <-
    function( cent, y.col.name = col.name, x.col.name = "AGE", group.col.name = "SEX", data, x.pred = NULL, fam = c( "BCPE", "BCT", "BCCG", "BCPEo", "BCTo", "BCCGo" ) ) {
        grps <-
            levels( data[ , group.col.name ] )
        as.data.frame(
            sapply(
                grps,
                function( g ) {
                    d.g <-
                        data[ data[ , group.col.name ] == g, ]
                    d.g.mdl <-
                        compute.model( d.g[ , y.col.name ], d.g[ , x.col.name ], fam )
                    if( is.null( x.pred ) ) {
                        x.pred <-
                            d.g.mdl$xvar }
                    d.g.prd <-
                        compute.prediction( x.pred, d.g.mdl )
                    d.g.prcntls <-
                        compute.percentiles( cent = cent, d.g.prd )
                    d.g.prcntls[ , x.col.name ] <-
                        x.pred
                    d.g.sds <-
                        compute.sds( d.g.mdl )
                    d.g.ref <-
                        compute.references( cent, d.g.prd )
                    list(
                        model = d.g.mdl,
                        pred  = d.g.prd,
                        cent  = d.g.prcntls,
                        sds   = d.g.sds,
                        ref   = d.g.ref
                    )
                }
            )
        )
    }
