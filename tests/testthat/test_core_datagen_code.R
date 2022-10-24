
##
## Testing the data generation code for the power simulations
##

library( testthat )
library( RIQITE )
library( tidyverse )

test_that("correlation of tx impact and Y0 setting works", {

    Y = rnorm( 100 )
    tau = rnorm( 100 )
    tau2 = RIQITE:::correlate_Y0_and_tx( Y, tau, rho = -1 )

    expect_true( cor( rank(Y), rank(tau2) ) == -1 )

    tau3 = RIQITE:::correlate_Y0_and_tx( Y, tau, rho = 0 )
    expect_true( all( tau == tau3 ) )


    tfunc = tx_function_factory( "rexp", rho = 1 )
    tt = tfunc( Y )
    expect_true( cor( rank(Y), rank(tt) ) == 1 )
    expect_true( cor( Y, tt ) > 0 )

    tfunc = tx_function_factory( "rexp", rho = 0.2 )
    tt = tfunc( Y )
    cor( tt, Y )
    expect_true( cor( rank(Y), rank(tt) ) > 0 )
    expect_true( cor( Y, tt ) > 0 )

})








