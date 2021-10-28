
# Utility to help printing out nicely formatted stuff.
scat = function( str, ... ) {
    cat( sprintf( str, ... ) )
}



#' Run the stephenson test.
#'
#' This method runs coin's independence test with the stephenson rank.  It
#' basically is a shorthand call to coin.
#'
#' @param data The dataframe with the variables defined in 'formula'
#' @param formula Defines the testing relationship.  See coin package's
#'   independence_test
#' @param distribution Passed to independence_test.  approximate(K) will do K
#'   permutations.  "asymptotic" will calculate an asymptotic approximation.
#'   "exact" will do exact.
#' @param alternative Direction of test.  E.g., "greater" will test whether the
#'   first factor in the list of factors of treatment is higher than the other.
#'   So the first factor level is the assumed treatment.
#' @param tail Do we focus on the right tail or left tail of the outcomes (i.e.,
#'   do we flip the sign of the outcomes).
#' @param subset_size Look at all subsets of this size and score 1 if treatment
#'   group has largest value.
#' @param ... Extra arguments to pass to the coin's `independence_test()` method.
#'
#' @return The result of coin's `independence_test()`
#' @export
stephenson_test <- function (data, formula, distribution = approximate(1000),
                             alternative = "less", tail = c( "right", "left" ),
                             subset_size = 5, ...) {
    tail = match.arg( tail )

    if (tail == "right") {
        sign <- 1
        alt <- alternative
    } else {
        sign <- -1
        alt <- ifelse(alternative == "less", "greater", "less")
    }
    res = coin::independence_test(
        formula = formula, data = data, distribution = distribution,
        alternative = alt, ...,
        ytrafo = function (data) {
            stephenson_trafo(data = data, sign = sign,
                             subset_size = subset_size)
        }
    )
    #if ( tail == "left" ) {
    #    warning( "Need to flip signs on outcomes for results" )
    #}
    res
}


#' Bundle the stephenson rank function into a `transformation` that allows coin package to use it in its independence tests.
#' @param data Dataframe
#' @param subset_size Size of subsets for seeing which observation is largest
#' @param sign Do we flip sign of all y to look for extreme low values rather than extreme high?
#' @export
stephenson_trafo  <- function (data, subset_size = 6, sign = 1) {
    coin::trafo(data = data, numeric_trafo = function (y)
        stephenson_rank(y = sign * y, subset_size = subset_size))
}



#' Calculate stephenson ranks
#'
#' For each observation, calculate number of subsets of size `subset_size` where
#' this observation is the largest in the set.
#'
#' @param y Vector of numerical values.
#' @param subset_size Look at all subsets of this size
#' @return Vector of stephenson ranks for the given input vector of values.
#' @export
stephenson_rank <- function (y, subset_size = 6) {
    r <- rank(y, na.last = "keep")
    q_tilde <- ifelse(r < subset_size, 0, choose(r - 1, subset_size - 1))
    return(q_tilde)
}


#' Subtract off constant treatment effect from data
#'
#' @return adjusted dataframe with new "Yadj" column.
adjust_outcomes = function( data, Y, Z, tau0 ) {
    data$Yadj = data[[Y]] - ifelse( as.numeric( data[[Z]] ) == 1, tau0, 0 )
    data
}

#' Find confidence interval threshold by looking for p-values crossing alpha.
#'
#' Utility function for actually doing the line search for the CI methods.
#' @param step Initial step size for searching for threshold value.
find_ci_internal <- function ( Y, Z, S, coin_test, ci_level = .9,
                               precision = 0.001,
                               step = precision * 2^5,
                               alternative,
                               verbosity = 1,
                               round_digits,
                               distribution = "asymptotic", init_tau0 = NULL, tau_grid = NULL,
                               ...) {
    if (verbosity >= 2) {
        scat( "\nBeginning binary search at init_tau0 = %.4f, initial step size of %.4f, and target precision of %.6f\n",
              init_tau0, step, precision)
    }

    test.tau = function( tau0 ) {
        if (verbosity >= 3) {
            scat( "Testing tau0=%f\n", tau0, "\n" )
        }
        Y_adj <- ifelse(Z == "treated", Y - round(tau0, round_digits), Y)
        test_data <- data.frame(Z, S, Y_adj)
        test_args <- list(
            formula = Y_adj ~ Z | S,
            data = test_data,
            alternative = alternative,
            distribution = distribution,
            ...
        )
        test_res <- do.call(coin_test, args = test_args)
        if (verbosity >= 3) print(test_res)
        p <- coin::pvalue(test_res)
        p
    }


    #browser()
    tau.list = rep( NA, 50 )
    pvalue.list = rep( NA, 50 )
    max_tx = diff(range(Y))
    alpha = 1-ci_level
    if ( verbosity >= 2 ) {
        scat( "Testing with alpha=%f\n\tMax possible tx = %f\n", alpha, max_tx )
    }

    # step 1: bracket our tau0
    p.orig = test.tau( init_tau0 )
    tau.list[1] = init_tau0
    pvalue.list[1] = p.orig
    cntr = 1
    if ( verbosity >= 2 ) {
        scat( "Pvalue = %f for initial tau=%f\n", p.orig, init_tau0 )
    }
    seek_low_pvalue = !( p.orig <= alpha )
    done = FALSE
    failed = FALSE
    double = 1

    # The direction of a step (go up or go down) for the search to bracket the
    # threshold value
    alt_sign <- if (identical(alternative, "greater")) -1 else 1
    if ( seek_low_pvalue ) {
        alt_sign = alt_sign * -1
    }
    p.low = NA
    p.high = NA

    tau0 = init_tau0
    while( !done ) {
        double = double*2
        init_tau0 = tau0
        tau0 <- tau0 - alt_sign * step * double
        p = test.tau( tau0 )
        tau.list[cntr] = tau0
        pvalue.list[cntr] = p
        cntr = cntr+1
        if ( verbosity >= 2 ) {
            scat( "Pvalue = %f for tau=%f\n", p, tau0 )
        }
        if ( ( seek_low_pvalue && (p <= alpha)) || (!seek_low_pvalue && (p > alpha) ) ) {
            done = TRUE
        } else if ( abs(tau0) > max_tx ) {
            failed = TRUE
            done = TRUE
        }
    }

    if ( failed ) {
        stop( "Failed to bracket p-values" )
    }

    # the above search gives us fenceposts.   Do a binary search within these
    # posts to find our threshold tau.
    if ( verbosity >= 1 ) {
        scat( "Search range of %f (pv=%f) - %f (pv=%f)\n", init_tau0, p.orig, tau0, p )
        scat( "\tWidth of range = %f\n", abs( init_tau0 - tau0 ) )
    }
    if ( seek_low_pvalue ) {
        tau0.low = tau0
        p.low = p
        p.high = p.orig
        tau0.high = init_tau0
    } else {
        tau0.low = init_tau0
        p.low = p.orig
        p.high = p
        tau0.high = tau0
    }

    # Zero in on threshold
    done = FALSE
    while ( abs( tau0.low - tau0.high ) > precision ) {
        tau0 = (tau0.low + tau0.high) / 2

        p = test.tau( tau0 )
        tau.list[cntr] = tau0
        pvalue.list[cntr] = p
        cntr = cntr + 1

        if ( p <= alpha ) {
            tau0.low = tau0
            p.low = p
        } else {
            tau0.high = tau0
            p.high = p
        }
        if ( verbosity >= 2 ) {
            scat( "Fenceposts are now %f - %f (delta=%f)\n", tau0.low, tau0.high, tau0.high - tau0.low )
        }
    }


    pvalues = data.frame( tau0 = tau.list, pvalue = pvalue.list )
    pvalues = pvalues[1:(cntr-1),]

    # add in extra taus, if any
    if ( !is.null( tau_grid ) ) {
        pvs = sapply( tau_grid, test.tau )
        pvalues = rbind( pvalues,
                         data.frame( tau0 = tau_grid, pvalue = pvs ) )
    }

    pvalues = pvalues[ order( pvalues$tau0 ), ]

    ci_bound <- (tau0.low+tau0.high)/2

    out <- list(ci_bound = ci_bound, pvalues = pvalues)

    #out$is_approximate = grepl("Approx", class(test_res@distribution))
    #out$is_approximate = TRUE #distribution != "asymptotic"
    out$tau0 = (tau0.low+tau0.high)/2
    out$tau0.low = tau0.low
    out$tau0.high = tau0.high
    out$pvalue.low = p.low
    out$pvalue.high = p.high
    out$range = abs( tau0.low - tau0.high )
    out$step = step
    out$precision = precision

    if ( verbosity >= 2 ) {
        scat( "Total number of grid search steps = %d\n", cntr - 1 )
    }
    out
} # find_ci_internal




#' @title Confidence intervals for maximal (minimal) effects
#'
#' @description Find confidence interval with test inversion
#'
#'   Search a range of values and return a confidence interval for maximal
#'   (minimal) treatment impact.
#'
#'   If there is a stratification variable, permute within strata/blocks defined
#'   by this categorical variable in data, Default: NULL
#'
#' @param data Dataframe of data. This MUST contain the variables refered to in formula.
#' @param formula Formula of outcome on treatment potentiall blocked by some
#'   stratifying factor, e.g., Y ~ Z | D.  Treatment and strat var must be
#'   factors. It is assumed these variables are all inside the passed dataframe.
#' @param coin_test Function to calculate test statistic.
#' @param ci_level Level of confidence, Default: 0.9
#' @param step_fraction Resolution of grid search, Default: 0.001 (i.e., each
#'   step is this portion of range of Y).
#' @param alternative What direction (max or minimal effects), Default: 'less'
#' @param tr_level Name of tx level to consider "treatment"  Other levels are
#'   control., Default: NULL (and will take first of factor).
#' @param verbosity Print out stuff as we go, Default: 1
#' @param round_digits Round grid search to this resolution, Default: 3
#' @param precision How fine should the final grid be for finding the threshold
#'   tau0 for the CI?
#' @param distribution What kind of ref distribution to calculate (asymp or
#'   approximate), Default: 'asymptotic'
#' @param init_tau0 Starting point for search, Default: NULL
#' @param tau_grid List of tau values to test in addition to the search.
#' @param pvalue_ci Return a confidence interval for the calculated pvalue of
#'   the crossing tau0 (if using approximate MC simulation).
#' @param ... Passed to the coin_test function
#' @return Confidence interval (one-sided)
#' @seealso \code{\link[coin]{pvalue}}
#' @rdname find_ci
#' @export
#' @importFrom coin pvalue
find_ci <- function (data,
                     coin_test,
                     formula,
                     tr_level = NULL,
                     alternative = c( "greater", "less" ),
                     distribution = "asymptotic",
                     ci_level = .9,
                     step_fraction = .001,
                     precision = 0.001,
                     verbosity = 0,
                     round_digits = 3,
                     init_tau0 = NULL,
                     tau_grid = NULL,
                     pvalue_ci = TRUE,
                     #                     guess_and_refine = is.null( init_tau0 ),
                     ...) {
    alternative = match.arg(alternative)
    stopifnot(alternative %in% c("greater", "less"))
    stopifnot(verbosity %in% 0L:3L)
    if (verbosity >= 1) {
        if (identical(alternative, "less")) {
            cat("\n\nTesting Null Hypothesis:",
                "Y_i(1) - Y_i(0) >= tau0 for all i ....\n")
        }
        if (identical(alternative, "greater")) {
            cat("\n\nTesting Null Hypothesis:",
                "Y_i(1) - Y_i(0) <= tau0 for all i ....\n")
        }
    }
    terms = all.vars( formula )
    stopifnot( length( terms ) == 2 || length( terms ) == 3 )
    tr_var = terms[[2]]
    stopifnot( !is.null( data[[tr_var]] ) )

    y_var = terms[[1]]
    stopifnot( !is.null( data[[y_var]] ) )
    if ( length( terms ) == 3 ) {
        strat_var = terms[[3]]
        stopifnot( !is.null( data[[strat_var]] ) )

    } else {
        strat_var = NULL
    }

    if (is.null(tr_level)) {
        tr_level <- levels(factor(data[[tr_var]]))[1]
    } else {
        stopifnot(tr_level %in% levels(factor(data[[tr_var]])))
    }
    Z <- ifelse(data[[tr_var]] == tr_level, "treated", "control")
    Z = factor( Z, levels = c("treated","control" ) )
    Y <- data[[y_var]]
    S <- if (is.null(strat_var)) gl(1, nrow(data)) else data[[strat_var]]

    step <-  diff(range(Y)) * step_fraction

    # If we want a two-step process where we do a rough detection and then the
    # final, actual permutation detection, we can put the "guess and refine"
    # code here (if we put it back in package) See below for removed code.
    if (is.null(init_tau0)) {
        init_tau0 <- 0
    }

    out = find_ci_internal( Y, Z, S,
                            coin_test=coin_test, ci_level=ci_level,
                            precision=precision, step=step, alternative=alternative,
                            verbosity=verbosity,
                            round_digits=round_digits, distribution=distribution,
                            init_tau0 = init_tau0, tau_grid = tau_grid, ... )

    # Clean up results and send back.
    if ( alternative == "less" ) {
        out$ci = c( -Inf, out$ci_bound )
    } else {
        out$ci = c( out$ci_bound, Inf )
    }
    names( out$ci ) = c( "low", "high" )
    out$alternative = alternative
    out$ci_level = ci_level
    out$round_digits = round_digits
    out$pvalue_ci = pvalue_ci

    class( out ) = "permconfint"
    invisible(out)

}


# Doc for this flag (now removed)
# #' @param guess_and_refine TRUE if method should conduct a two-step process,
# #'   first using approximate and rough grid search to get limits, and then to
# #'   find more exact tau0.


# if ( guess_and_refine ) {
#     if (is.null(init_tau0)) {
#         init_tau0 <- diff(range(Y)) * alt_sign
#     }
#
#     stopifnot( step_fraction < 0.01 )
#     rstep = step * 100
#     rdistribution = "asymptotic"
#     if ( verbosity > 1 ) {
#         scat( "Making initial asymptotic pass with step %.3f\n", rstep )
#     }
#
#     out = find_ci_internal( Y, Z, S,
#                             coin_test=coin_test, ci_level=ci_level,
#                             step=rstep, alt_sign=alt_sign,
#                             alternative=alternative, verbosity=verbosity,
#                             round_digits=round_digits, distribution=rdistribution,
#                             init_tau0 = init_tau0, ... )
#
#     init_tau0 = out$ci_bound + 5 * rstep * alt_sign
#     if ( verbosity > 1 ) {
#         scat( "Found rough initial tau0 of %.3f\n", init_tau0 )
#     }
# } else {
#     if (is.null(init_tau0)) {
#         init_tau0 <- diff(range(Y)) * alt_sign
#     }
# }



#' @title Print method for results from confidence interval call.
#'
#' @description This is a utility function.
#'
#' @param x Main argument.
#' @param ... Additional arguments.
#'
#' @export
print.permconfint = function( x, ... ) {
    # [1] "ci_bound"     "pvalues"      "tau0"         "tau0.low"     "tau0.high"    "pvalue.low"   "pvalue.high"
    #    [8] "range"        "step"         "precision"    "ci"           "alternative"  "ci_level"     "round_digits"
    #    [15] "pvalue_ci"
    with( x, {

        cat( "Permutation confidence interval for extreme effects\n" )

        if (alternative == "greater") {
            cat("\n", ci_level * 100,
                "% confidence interval for maximum treatment effect: [",
                round(ci_bound, round_digits), ", Inf)\n\n", sep = "")
        } else if (alternative == "less") {
            cat("\n", ci_level * 100,
                "% confidence interval for minimum treatment effect: (-Inf, ",
                round(ci_bound, round_digits), "]\n\n", sep = "")
        }

        #        scat( "Resolution of grid: tau0 = %f - %f (delta = %f)\n",
        #              round(tau0.low, round_digits),
        #              round(tau0.high, round_digits),
        #              round(range, round_digits))
        #       if ( pvalue_ci && !is.null( attr(p, "conf.int") ) ) {
        #           cat("\n99% confidence interval for the p-value of test of final tau0 = ",
        #               round(tau0, round_digits), ": [", round(attr(p, "conf.int")[1], 4),
        #                ", ", round(attr(p, "conf.int")[2], 4), "]", sep = "")
        #       }

    } )

}


#' Wrapper for find_ci to use the Stephenson Rank statistic.
#'
#' @rdname find_ci
#' @export
find_ci_stephenson <- function (data, formula, tr_level = NULL,
                                distribution = "asymptotic",
                                alternative = "less", tail = c( "right", "left" ),
                                subset_size = 5,
                                ci_level = .9,
                                precision = 0.001,
                                step_fraction = .001,
                                verbosity = 0, round_digits = 3,
                                init_tau0 = NULL,
                                tau_grid = NULL,
                                pvalue_ci = TRUE ) {

    tail = match.arg( tail )

    if (tail == "right") {
        alt <- alternative
    } else {
        terms = all.vars( formula )
        stopifnot( length( terms ) == 2 || length( terms ) == 3 )
        y_var = terms[[1]]

        data[ y_var ] = data[ y_var ] * -1
        alt <- ifelse(alternative == "less", "greater", "less")
    }

    res = find_ci( data, formula, tr_level = tr_level,
                   coin_test = "stephenson_test",
                   distribution=distribution,
                   alternative=alt,
                   ci_level=ci_level, step_fraction=step_fraction, precision = precision,
                   verbosity=verbosity, round_digits=round_digits,
                   init_tau0=init_tau0,
                   tau_grid = tau_grid,
                   pvalue_ci = pvalue_ci,
                   #             guess_and_refine = guess_and_refine,
                   subset_size=subset_size )

    if ( tail == "left" ) {
        res$ci = -1 * rev( res$ci )
        res$tau0 = -1 * res$tau0
        res$pvalues$tau0 = -1 * res$pvalues$tau0
        res$ci_bound = -1 * res$ci_bound
    }
    res
}



#' Generate the four possible stephenson rank based confidence intervals.
#'
#' @param data The data
#' @param formula The formula
#' @param ... Parameters passed to `stephenson_test` and `find_ci_stephenson`
#' @return Table of the four confidence intervals determined by which tail and
#'   which alternative.
#' @export
generate_four_stephenson_cis = function( data, formula, tr_level = NULL, ... ) {

    CIs.raw = expand.grid( alternative=c("less","greater"),
                           tail = c("right","left" ),
                           stringsAsFactors = FALSE)

    process.one = function( alternative, tail ) {
        cat( "Processing alt=", alternative, "  -  tail=", tail, "\n", sep="" )

        st1 <- stephenson_test(data = data,
                               formula = formula,
                               alternative = alternative, tail = tail,
                               ...)


        ci1 <- find_ci_stephenson( data, formula,
                                   tr_level = tr_level,
                                   alternative = alternative,
                                   tail = tail,
                                   ... )

        tibble( pvalue = as.numeric( pvalue( st1 ) ),
                CI.low = ci1$ci[[1]],
                CI.high = ci1$ci[[2]],
                pvaluelist = list( ci1$pvalues ) )
    }

    CIs.raw = as_tibble( CIs.raw )
    CIs.raw$res = map2( CIs.raw$alternative, CIs.raw$tail, process.one )

    # TODO The following line unnest() makes a warning.  Not sure why, or what
    # to do about it.
    CIs = CIs.raw %>%
        unnest( cols=res )

    CIs
}



