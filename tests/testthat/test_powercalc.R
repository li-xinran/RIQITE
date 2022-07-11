
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

  ff = tx_function_factory("scale", scale=1.5 )
  imps <- ff( Y0s )
  expect_true( all( imps == ( (1:10) * 0.5 ) ) )

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


test_that( "finite power calc works", {
  set.seed( 40404 )
  dat = generate_finite_data( 500,
                              tx_function = tx_function_factory( "rexp" ) )
  pow = calc_power_finite( dat$Y0, dat$tau, p_tx = 0.5 )
  pow
  pow2 = calc_power_finite( dat$Y0, dat$tau, p_tx = 0.5, percentile = 0.99 )
  pow2
  expect_true( pow2$k == 495 )
  expect_true( pow2$power < pow$power )
})



test_that( "Explore s works", {
  esses = c( 2, 4, 5, 7, 10, 30 )
  ess = explore_stephenson_s( s = esses,
                              n = 500,
                              tx_function = tx_function_factory("rexp"),
                              p_tx = 0.5,
                              R = 30 )
  expect_true( nrow( ess ) == 6 )
  expect_true( all( ess$s == esses ) )
})




test_that( "Explore s for emperical works", {
  esses = c( 2, 4, 30 )
  Y0 = c( rep( 0, 19 ), 2 )
  Y0
  sd( Y0 )

  smpdat = generate_finite_data(100, Y0,
                                 tx_function_factory("rexp",
                                                            ATE = 0.50) )
  smpdat
  table( smpdat$Y0 )

  calc_power_finite( smpdat$Y0, smpdat$tau, p_tx = 0.5 )

  calc_power( n = 100, Y0,
              tx_function_factory("rexp", ATE = 0.25),
              p_tx = 0.5, calc_ICC = TRUE )

  ess = explore_stephenson_s( s = esses,
                              n = 50,
                              Y0_distribution = Y0,
                              tx_function = tx_function_factory("rexp",
                                                                ATE = 0.50),
                              p_tx = 0.5,
                              R = 1000, iter_per_set = 50,
                              calc_ICC = TRUE )
  ess
  expect_true( nrow( ess ) == 3 )
  expect_true( ncol( ess ) == 5 )

})


