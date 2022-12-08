
##
## Testing the code
##

library( testthat )
library( RIQITE )
library( tidyverse )

context("Checking power calculators")

test_that("Tx function generator works", {

    Y0s = 1:10
    ff = tx_function_factory( "constant", ATE = 1 )
    expect_true( all( ff( Y0s ) == 1 ) )

    ff = tx_function_factory("scale", rho=1.5 )
    imps <- ff( Y0s )
    expect_true( all( imps == ( (1:10) * 1.5 ) ) )

    ff = tx_function_factory("rexp", ATE = 0 )
    expect_true( length( ff( Y0s ) ) == 10 )

    ff = tx_function_factory("rexp", ATE = 0.5, tx_scale = 0 )
    expect_true( all.equal( ff( Y0s ),  rep( 0.5, 10 ) ) )
} )



test_that( "data generation works", {
    dat = generate_finite_data( 350,
                                tx_function = tx_function_factory( "constant", ATE  = 0.2 ) )
    expect_true( nrow( dat ) == 350 )
    expect_true( ncol( dat ) == 2 )

    expect_true( all.equal( dat$tau, rep( 0.2, 350 ) ) )
})




test_that( "CI methods work", {
    dat = generate_finite_data( 150,
                                tx_function = tx_function_factory( "rnorm", ATE  = 0.2 ) )

    expect_true( sd( dat$tau) > 0.90 )

    res <- calc_power_finite( dat$Y0, dat$tau, 0.5, R=5, summarise = FALSE, nperm = 200, targeted_power = FALSE )
    expect_true( is.matrix(res) )
    expect_equal( dim( res ), c( 5, 150 ) )
    expect_true( !is.null( colnames( res ) ) )

    rs <- RIQITE:::summary_power_results( res )
    expect_true( nrow( rs ) == 150 )
    expect_true( !is.null( rs$k ) )
    expect_true( all( rs$power >= 0 ) )
    expect_true( is.data.frame(rs) )

    res3 <- calc_power_finite( dat$Y0, dat$tau, 0.5, R=5, summarise = TRUE, nperm = 200, targeted_power = FALSE,
                                      c = 0, quantile_n_CI = 0.10,
                                      k.vec = 140:150 )
    expect_true( is.data.frame(res3) )
    expect_true( is.numeric( res3$k ) )
    expect_equal( colnames(res3), colnames(rs) )
})


test_that( "explore_stephenson_s_finite CI works", {

    dat = generate_finite_data( 150,
                                tx_function = tx_function_factory( "rnorm", ATE  = 0.2 ) )

    res <- explore_stephenson_s_finite( s = c(2, 5), Y0 = dat$Y0, tau = dat$tau, R = 10, nperm = 100,
                                        k.vec = 140:150,
                                        p_tx = 0.5, c = -3, targeted_power = FALSE )
    res
    table( res$s )
    expect_true( nrow(res) == 22 )

    res2 <- explore_stephenson_s_finite( s = c(2, 5), Y0 = dat$Y0, tau = dat$tau, R = 10,
                                         p_tx = 0.5, c = -3, targeted_power = TRUE, percentile = 0.95 )

    res2
    expect_equal( res2$k,  c( 142, 142 ) )
})



test_that( "finite power calc works", {
    set.seed( 40404 )
    dat = generate_finite_data( 500,
                                tx_function = tx_function_factory( "rexp" ) )
    pow = calc_power_finite( dat$Y0, dat$tau, p_tx = 0.5, targeted_power = TRUE )
    pow
    pow2 = calc_power_finite( dat$Y0, dat$tau, p_tx = 0.5, percentile = 0.99, targeted_power = TRUE )
    pow2
    expect_true( pow2$k == 495 )
    expect_true( pow2$power < pow$power )


    # Get raw CI results as a matrix
    nrow(dat)
    dat$tau = 4 #max( dat$Y0 ) - min( dat$Y0 )
    pow_raw = calc_power_finite( dat$Y0, dat$tau, p_tx = 0.5, targeted_power = FALSE, summarise = FALSE, R = 10, nperm = 100 )
    expect_equal( dim( pow_raw ), c( 10, 500 ) )

})




test_that( "raw results for superpop power calc works", {
    set.seed( 404044 )

    # Get raw CI results as a matrix
    pow_raw = calc_power( n = 15,
                          p_tx = 0.5, targeted_power = FALSE, summarise = FALSE, R = 10, iter_per_set = 5, nperm = 100,
                          k.vec = 10:15 )
    pow_raw
    expect_equal( dim( pow_raw ), c( 10, 6 ) )
    expect_true( is.matrix( pow_raw ) )

    sps <- summary_power_results( pow_raw, quantile_n_CI = NA, c = -2 )
    sps
    expect_true( is.data.frame(sps) )
    expect_true( nrow(sps) == 6 )
    expect_true( all( sps$n > 0 ) )
    expect_true( sd( sps$n ) == 0 )
})





test_that( "Explore s works", {
    esses = c( 2, 4, 5, 7, 10, 30 )
    ess = explore_stephenson_s( s = esses,
                                n = 500,
                                tx_function = "rexp",
                                p_tx = 0.5,
                                R = 30 )
    ess
    expect_true( nrow( ess ) == 6 )
    expect_true( all( ess$s == esses ) )
})




test_that( "Explore s with CI approach works", {
    esses = c( 2, 5, 10 )
    ess = explore_stephenson_s( s = esses,
                                n = 100,
                                tx_function = "rexp",
                                p_tx = 0.5, c = -1,
                                nperm = 200,
                                R = 100, k.vec = 95:100, verbose = FALSE,
                                targeted_power = FALSE, calc_ICC = TRUE )

    ess
    expect_true( nrow( ess ) == 6 * 3 )
    expect_true( all( ess$s == rep( esses, each=6 ) ) )
})





test_that( "Explore s for emperical data works", {
    esses = c( 2, 4, 30 )
    Y0 = c( rep( 0, 19 ), 2 )
    Y0
    sd( Y0 )

    smpdat = generate_finite_data(100, Y0,
                                  tx_function_factory("rexp",
                                                      ATE = 0.50) )
    smpdat
    table( smpdat$Y0 )

    calc_power_finite( smpdat$Y0, smpdat$tau, p_tx = 0.5, targeted_power = TRUE )

    calc_power( n = 100, Y0,
                tx_function_factory("rexp", ATE = 0.25),
                p_tx = 0.5, calc_ICC = TRUE )

    ess = explore_stephenson_s( s = esses,
                                n = 50,
                                Y0_distribution = Y0,
                                tx_function = "rexp",
                                ATE = 0.50,
                                p_tx = 0.5,
                                R = 1000, iter_per_set = 50,
                                calc_ICC = TRUE )
    ess
    expect_true( nrow( ess ) == 3 )
    expect_true( ncol( ess ) == 10 )

})






test_that( "Fully superpop simulation works (no reusing data)", {
    esses = c( 3, 10 )

    raw_res <- calc_power( n = 90, Y0_distribution = rnorm, R = 7, nperm = 100, iter_per_set = 1,
                tx_function_factory("rexp", ATE = 0.25),
                p_tx = 0.5, calc_ICC = TRUE,
                summarise = FALSE, targeted = FALSE )
    dim( raw_res )
    expect_equal( dim( raw_res ), c( 7, 90 ) )

    ss = summary_power_results(raw_res)
    expect_true( nrow( ss ) == 90 )

})



test_that( "Fully superpop stephenson-s works (no reusing data)", {

    sres <- explore_stephenson_s( s = c(4, 30),
                                     n = 90, Y0_distribution = rexp, R = 7, nperm = 100, iter_per_set = 1,
                           tx_function_factory("rexp", ATE = 0.25),
                           p_tx = 0.5, targeted = FALSE )
    expect_true( nrow(sres) == 90 * 2 )
    expect_true( length( unique( sres$n ) ) == 2 )
})





