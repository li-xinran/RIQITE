


# Comparing different linear models with same marginal variance


library( tidyverse )
library( RIQITE )

# rho = 1, omega = 1
# Y(1) variation is (1+1)^2 + 1 = 5

N = 150

# single giant dataset
dat = generate_finite_data( 20000,
                            tx_function = tx_function_factory("linear", ATE = 0.0, rho = 1, tx_scale = 1) )
dat
skimr::skim( dat )
var( dat$tau )
var( dat$Y0 + dat$tau )
cor( dat$Y0, dat$tau )

explore_stephenson_s_finite(  Y0 = dat$Y0[1:N], tau = dat$tau[1:N],
                              s = c( 2, 5, 20 ),
                              p_tx = 1/3 )

calc_power_finite( Y0 = dat$Y0, tau = dat$tau, p_tx = 1/3 )

calc_power( n = N,
            tx_function = tx_function_factory("linear", ATE = 0.0, rho = 1, tx_scale = 1),
            p_tx = 1/3 )


# Repeatedly do this across collection of datasets
cat( "Starting first exploration\n" )

essA = explore_stephenson_s( s = c( 2, 5, 20 ),
                            n = N,
                            tx_function = tx_function_factory("linear", ATE = 0.0, rho = 1, tx_scale = 1),
                            p_tx = 1/3,
                            calc_ICC = TRUE,
                            parallel = TRUE, R = 2000, nperm=10000, verbose=TRUE )
essA



## Example 2: Negative rho

# rho = -1, omega = sqrt(5)
# Y(1) variation is (1+-1)^2 + 5 = 5

dat2 = generate_finite_data( 20000,
                            tx_function = tx_function_factory("linear", ATE = 0.0, rho = -1, tx_scale = sqrt(5)) )
dat2
var( dat2$tau )
var( dat2$Y0 + dat2$tau )
cor( dat2$Y0, dat2$tau )


# Look at power across scenarios for this alternative DGP
essB = explore_stephenson_s( s = c( 2, 5, 20 ),
                             n = N,
                             tx_function = tx_function_factory("linear", ATE = 0.0, rho = -1, tx_scale = sqrt(5)),
                             p_tx = 1/3,
                             calc_ICC = TRUE,
                             R = 2000, nperm = 10000,
                             parallel = TRUE, verbose = TRUE )
essB


## Example 3: 0 rho

# rho = 0, omega = sqrt(4) = 2
dat3 = generate_finite_data( 20000,
                             tx_function = tx_function_factory("linear", ATE = 0.0, rho = 0, tx_scale = sqrt(4)) )
dat3
var( dat3$tau )
var( dat3$Y0 + dat3$tau )
cor( dat3$Y0, dat3$tau )


# Look at power across scenarios for this alternative DGP
essC = explore_stephenson_s( s = c( 2, 5, 20 ),
                             n = N,
                             tx_function = tx_function_factory("linear", ATE = 0.0, rho = 0, tx_scale = sqrt(4)),
                             p_tx = 1/3,
                             calc_ICC = TRUE,
                             R = 2000, nperm = 10000,
                             parallel = TRUE, verbose = TRUE )
essC


## Table of results

bind_rows( A = essA, B = essB, C = essC, .id = "scenario" ) %>%
    filter( s == 5 ) %>%
    dplyr::select( -s ) %>%
    knitr::kable( digits = 2)
