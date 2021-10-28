
##
## Utility functions for weak null simulations
##
## Mainly this generates datasets with different sorts of potential outcome
## structures.
##

##### Utility function ######


scat = function(str, ...) {
    cat(sprintf(str, ...))
}


##### Theoretical power calculation functions #####

# Take a dataset and calculate asymptotic power from the neyman approximation.
analyze.power = function( dat, N = nrow(dat), p.tx=0.5 ) {
    # Calc summary statistics
    dd = data.frame( mu0 =  mean( dat$Y.0 ),
                     mu1 = mean( dat$Y.1 ),
                     sd0 = sd( dat$Y.0 ),
                     sd1 = sd( dat$Y.1 ) )

    # Get power
    dd = mutate( dd,
                 tau.true = mu1 - mu0,
                 SE = sqrt( sd0^2 / (N*(1-p.tx)) + sd1^2/(N*p.tx) ),
                 pow = 1 - pnorm( qnorm(1-alpha), mean = tau.true / SE, sd=1 ) )
    ##dd$pow
    dd
}



# Calculate power of t-test by getting charactaristics of the DGP (by generating
# huge sample) and calculating some summary statistics.
calc.power = function( ..., DGP, p.tx, alpha = 0.10 ) {

    # make a fake dataset to get moments from
    dots = list( ... )
    N = dots$N
    dots$N = NULL
    dat = do.call( make.obs.data, c( N = 10000, dots, DGP=DGP, p=0.5 ) )

    analyze.power( dat, N=N, p.tx=p.tx )
}



# Calculate the theoretical power of a t-test for a given DGP by generating a
# huge sample, pulling moments off it and then doing some math.
#
# This one is taylored for the "extreme" DGP (make.data.mix)
#calc.power.extreme = function( tau, extreme.p, N, p.tx, alpha = 0.10 ) {
#    dat = make.data.mix( 100000, tau=tau, p.extreme=extreme.p )
#    analyze.power( dat, N=N, p.tx=p.tx )
#}



##### Data generating functions #####


control.side.distribution = function( N, distribution = "normal" ) {
    if ( distribution == "logistic" ) {
        rlogis( N, scale = sqrt(3) / pi )
    } else {
        rnorm( N )
    }
}

control.side.cdf = function( y, distribution = "normal" ) {
    if ( distribution == "logistic" ) {
        plogis( y )
    } else {
        pnorm( y )
    }
}

control.side.cdf.inv = function( p, distribution = "normal" ) {
    if ( distribution == "logistic" ) {
        qlogis( p )
    } else {
        qnorm( p )
    }
}



if ( FALSE ) {
    A = data.frame( type = "logistic", Y = rlogis( 10000 ) * sqrt(3) / pi, stringsAsFactors = FALSE )
    B = data.frame( type = "normal", Y = rnorm( 10000 ), stringsAsFactors = FALSE )
    C = bind_rows( A, B )
    table(C$type )
    ggplot( C, aes( x = Y, col=type ) ) +
        geom_density()
    plot( ecdf( A$Y ) )
    plot( ecdf( B$Y ), add=TRUE, col="green" )
}


#' Generate data where Y0 is N( 0, 1 ) and Y1 has a treatment impact of an
#' exponential added on to it with parameter tau.lambda
#' @return Dataframe with N rows and Y.0, Y.1, tau for each row.
make.data.exp = function( N, tau.lambda = 1, distribution = "normal" ) {
    Y0 = control.side.distribution( N, distribution )
    Y1 = Y0 + rexp( N, tau.lambda )
    data.frame( Y.0 = Y0, Y.1 = Y1, tau = Y1 - Y0 )
}



#' Generate data where Y0 is N( 0, 1 ) and we add an exponential distributed
#' treatment impact to a small proportion of the units (the rest have 0 impact).
#' @param tau The average tx impact of the units that have a treatment effect.
#' @return Dataframe with N rows and Y.0, Y.1, tau for each row.
make.data.mix = function( N, tau = 0.2, p.extreme = 0.20, distribution = "normal"  ) {
    Y0 = control.side.distribution( N, distribution )
    t.extr = tau # / p.extreme
    Y1 = Y0 + ifelse( runif( N ) <= p.extreme, rexp( N, 1 / t.extr ), 0 )
    data.frame( Y.0 = Y0, Y.1 = Y1, tau = Y1 - Y0 )
}


#' Generate data where Y0 is N( 0, 1 ) and we add 0 or a constant positive impact with
#' probability (1-p.extreme) or p.extreme
#' @param tau The average tx impact of the units that have a treatment effect.
#' @return Dataframe with N rows and Y.0, Y.1, tau for each row.
make.data.const.rare = function( N, tau = 0.2, p.extreme = 0.20, distribution = "normal"  ) {
    Y0 = control.side.distribution( N, distribution )
    t.extr = tau# / p.extreme
    Y1 = Y0 + ifelse( runif( N ) <= p.extreme, t.extr, 0 )
    data.frame( Y.0 = Y0, Y.1 = Y1, tau = Y1 - Y0 )
}



#' @param tau The target average treatment impact
#' @param p.extreme The proportion of units with impacts.
#'
make.data.uni.uni = function( N, tau = 0.2, p.extreme = 0.20 ) {
    Y0 = runif( N )
    mu0 = 0.5
    mu1 = mu0 + tau
    #mx = (2/p.extreme) *( mu1 - (1-p.extreme)^2/2) + 1 - p.extreme
    mx = (2 * mu1 - 1 + p) / p
    #scat( "Mx = %.3f\n", mx )
    tau = ifelse( Y0 <= 1 - p.extreme, 0, (Y0 - (1-p.extreme)) * (mx - 1)/p.extreme )
    data.frame( Y.0 = Y0, Y.1 = Y0 + tau, tau = tau )
}


# Take the high control values and add treatment impacts to them.  (So an impact
# on the high end.)
make.data.rare.tail = function( N, scale = 2, p.extreme = 0.20, distribution = "normal" ) {
    Y0 = control.side.distribution( N, distribution )
    cut = control.side.cdf.inv( 1 - p.extreme )
    Y1 = Y0 + ifelse( Y0 <= cut, 0, (Y0 - cut) * scale )
    data.frame( Y.0 = Y0, Y.1 = Y1, tau = Y1 - Y0 )

}

#' Generate data where Y0 is N( 0, 1 ) and tau is tau + t() * scale / (df/df-2)
#'
#' @return Dataframe with N rows and Y.0, Y.1, tau for each row.
make.data.t = function( N, tau = 0.2, scale = 0.2, df = 3, distribution = "normal"  ) {
    Y0 = control.side.distribution( N, distribution )
    taus = tau + rt( N, df=df ) * scale / sqrt( df/(df-2))
    Y1 = Y0 + taus
    data.frame( Y.0 = Y0, Y.1 = Y1, tau = taus )
}

# This implements, I think, the DGP analyzed by Rosenbaum's work.
# The core idea is this: generate the quantiles of Y0 and Y1 with a uniform.
# To get Y0, simply invert the CDF to get the corresponding y.
# For Y1, we flip a coin, and then if we have an effect we calculate
# y = F^{-1}( p^(1/m) )
# The CDF of this will be the inverse of y, i.e., F(y)^m
make.data.exp.cdf = function( N, p.extreme = 0.2, m = 2, distribution = "normal" ) {
    p0 = runif( N )
    p = ifelse( runif( N ) <= p.extreme, p0^(1/m), p0 )
    Y0 = control.side.cdf.inv( p0, distribution )
    Y1 = control.side.cdf.inv( p, distribution )

    data.frame( Y.0 = Y0, Y.1 = Y1, tau = Y1 - Y0 )
}


#' Generate data where Y0 is N( 0, 1 ) and Y1 is a constant shift
#'
#' @return Dataframe with N rows and Y.0, Y.1, tau for each row.
make.data.const = function( N, tau = 0.2 ) {
    Y0 = rnorm( N )
    Y1 = Y0 + tau
    data.frame( Y.0 = Y0, Y.1 = Y1, tau = Y1 - Y0 )
}





#' Given potential outcomes schedule, randomize within block and generate
#' observed potential outcomes and add them to the schedule.
#'
#' @param dat Dataframe with a named Y0, Y1, and block column
#' @param p Proportion of units treated.  Can be vector and will treat that
#'   proportion in each block, rounded to nearest and with at least 1 treatment
#'   and control unit.
#' @param Y0 name of Y0 column
#' @param Y1 name of Y1 column
#' @param blockvar name of blocking column.  This column will be converted to a
#'   factor, if it is not already, and the order of p corresponds to the levels
#'   of this factor.
#'
#' @return augmented `dat` with Z and Yobs columns
#' @export
add.obs.data = function(dat,
                        p = 0.5,
                        Y0 = "Y.0",
                        Y1 = "Y.1" ) {
    N = nrow(dat)

    Z = sample( nrow( dat ) ) <= p * nrow(dat)

    dat$Z = as.numeric( Z )

    dat$Yobs = ifelse(Z, dat[[Y1]], dat[[Y0]])

    dat
}




#' @title Simulated data generator
#'
#' @description Make data with randomized treatment
#' @param N Number of observations
#' @param p proportion units treated
#' @param DGP function to make the response schedule
#' @param ... parameters to pass to the DGP.
#' @return data.frame
#' @rdname make.obs.data
#' @export
make.obs.data = function( N, p, DGP = make.data.mix, ... ) {
    dat = DGP(N = N, ... )
    add.obs.data( dat, p = p )
}



if ( FALSE ) {
    library( tidyverse )

    debugonce(make.data.rare.tail )
    dat = make.obs.data( N=10, p=0.5, DGP = make.data.rare.tail )
    dat



    check.dat = function( p ) {

        dat = make.obs.data( 10000, p=0.5, DGP = make.data.rare.tail, p.extreme=p )

        ggplot( dat, aes( Y.0, Y.1 ) ) + geom_point() +
            geom_abline( intercept=0, slope=1, col="red")

        data.frame( p.extreme = p,
                    mu0 =  mean( dat$Y.0 ),
                    sd0 = sd( dat$Y.0 ),
                    sd1 = sd( dat$Y.1 ),
                    ATE = mean( dat$tau ) )
    }

    rs = map_df( c( 0.5, 0.2, 0.05 ), check.dat )
    rs
    rs = mutate( rs, SE = sqrt( sd0^2 / 25 + sd1^2/25 ),
                 pow = 1 - pnorm( 2, mean= 0.2 / SE, sd=1 ) )
    rs

}


