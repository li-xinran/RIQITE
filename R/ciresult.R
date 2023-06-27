



#' @title Result object for results of ci_quantile call
#'
#' @name ciresult
#'
#' @description
#' The ciresult object is an S3 class that holds the results from `ci()`.
#'
#' It has several methods that pull different information from this object, and
#' some printing methods for getting nicely formatted results.
#'
#'
#' @param x a ciresult object
#' (except for is.ciresult, where it is a generic object to check).
#'
#' @rdname ciresult
NULL


#' @return The parameters of the call
#'
#' @rdname ciresult
#' @export
arguments <- function( x ) {
    stopifnot( is.ciresult(x) )

    attr( x, "args" )
}



#' @return is.ciresult: TRUE if object is a ciresult object.
#'
#' @rdname ciresult
#'
#' @export
is.ciresult <- function(x) {
    inherits(x, "ciresult")
}



#' @title Get top few rows of simulation
#'
#' @param n Number of rows to take
#' @param ... extra options passed
#' @rdname ciresult
#'
#' @return Data.frame of top rows
#'
#' @export
head.ciresult <- function( x, n = 6L, ... ) {
    head( as.data.frame( x, n = n, ... ) )
}






scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}



#' @title Pretty print ci result
#'
#' @param ... extra options passed.
#' @rdname ciresult
#'
#' @return print: No return value; prints results.
#'
#' @export
print.ciresult <- function(x, ...) {

    ss = arguments(x)
    if ( is.null(ss$k.vec) ) {
        ss$k.vec = 1:ss$n
    }
    ss$min_k = min(ss$k.vec)
    ss$max_k = max(ss$k.vec)
    ss$method_name = paste0( ss$method.list, collapse = "-" )
    strg = glue::glue( "CI quantile on {ss$n} units with {ss$m} treated:
                        Method: {ss$method_name}
                        Alternative: {ss$alternative} with alpha={ss$alpha}
                        Quantiles covered from {ss$min_k} to {ss$max_k} ({length(ss$k.vec)} quantiles)
                        {ss$nperm} permutations with tolerance of {ss$tol}")
    scat( "%s\n", strg )

    invisible( ss )
}




sum_snip <- function( object, ss, c ) {

    n = 0
    np = 0
    nstr = "unknown"
    if ( ss$alternative == "greater" ) {
        n = sum( object$lower > c )
        np = round( 100 * n / ss$n )
        nstr = glue::glue( "{n} or more units ({np}%) with effects larger than {c}" )
        if ( n > 0 ) {
            nstr = glue::glue( "{nstr}
                           ..quantile {1+ss$n-n} and above have lower CI bounds > {c}" )
        }

    } else if (ss$alternative == "less" ) {
        n = sum( object$upper < c )
        np = round( 100 * n / ss$n )
        nstr = glue::glue( "{n} or more units ({np}%}) with effects smaller than {c}" )
        if ( n > 0 ) {
            nstr = glue::glue( "{nstr}\n..quantile {n-1} and below have higher CI bounds < {c}" )
        }
    }

    cat( nstr )
    cat( "\n" )

    tibble( c = c, n = n, np = np )
}



#' @title Pretty print ci result
#'
#' @description
#' Calculate some summary statistics on the result object.
#'
#' @param object object to summarize.
#' @param c Threshold for testing.
#' @param k Quantiles to report on.
#' @param ... extra options passed to print.ciresult
#' @rdname ciresult
#'
#' @return summary: No return value; prints results.
#'
#' @export
summary.ciresult <- function(object, c = 0, k = NULL, ...) {

    ss = arguments(object)
    if ( is.null( ss$k.vec ) ) {
        ss$k.vec = sort( unique( object$k ) )
    }

    n = 0
    nstr = "unknown"
    if ( ss$alternative == "greater" ) {
        fin = min( which( is.finite(object$lower ) ) )
        fin_per = round( 100 * fin / ss$n )
        fin_n = 1 + ss$n - fin
        if ( max(ss$k.vec) < ss$n ) {
            warning( "Summary with missing quantiles for maximal effects may be suspect.  Fix k.vec" )
        }
        nstr = glue::glue( "Quantile {fin} and above ({fin_per}%, or {fin_n} units) have finite lower CIs" )
    } else if (ss$alternative == "less" ) {
        if ( min(ss$k.vec) > 1 ) {
            warning( "Summary with missing quantiles for minimal effects may be suspect.  Fix k.vec." )
        }
        fin = max( which( is.finite(object$upper ) ) )
        fin_per = round( 100 * fin / ss$n )
        fin_n = fin
        nstr = glue::glue( "Quantile {fin} and below ({fin_per}%, or {fin_n} units) have finite lower CIs" )
    } else {
        warning( glue::glue( "Alternative '{alternative}' not fully implemented for summary" ) )
    }

    print.ciresult(object)
    cat( "\n" )
    cat( nstr )
    cat( "\n" )

    if ( ss$alternative == "greater" ) {
        fbnd = object %>%
            dplyr::filter( is.finite( lower ) )
        ss$bnd_low = min( fbnd$lower )
        ss$bnd_high = max( fbnd$lower )
        if ( nrow(fbnd) > 0 ) {
            scat( "Finite bounds range from %.2f to %.2f\n", min( fbnd$lower ), max( fbnd$lower ) )
            scat( "CI for maximal effect is %.2f and above.\n", max( fbnd$lower ) )
        }
    } else {
        fbnd = object %>%
            dplyr::filter( is.finite( upper ) )
        ss$bnd_low = min( fbnd$upper )
        ss$bnd_high = max( fbnd$upper )
        if ( nrow(fbnd) > 0 ) {
            scat( "Finite bounds range from %.2f to %.2f\n", min( fbnd$upper ), max( fbnd$upper ) )
            scat( "CI for minimal effect is %.2f and below.\n", min( fbnd$upper ) )
        }
    }

    cat( "\n" )

    mp = purrr::map_df(c, sum_snip, object = object, ss = ss )

    ss$c = mp

    if ( !is.null( k ) ) {
        kvec = k
        if ( length(k) == 1 && k == "distinct" ) {
            #pst = paste0( round( object$lower, digits=k_step ), "-", round( object$upper, digits=k_step ) )
            pst = paste0( object$lower, "-", object$upper )
            kvec = object$k[ !duplicated( pst ) ]
        }
        flt = object %>%
            dplyr::filter( k %in% kvec ) %>%
            knitr::kable( format = "simple",
                          caption = "Selected quantile confidence intervals") %>%
            print()

    }


    ss$k_cut = fin

    invisible( ss )
}


