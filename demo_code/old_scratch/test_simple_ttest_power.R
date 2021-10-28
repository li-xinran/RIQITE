
# Simulation checking the DGP code and the power of the simple t-test as
# compared to the asymptotic approximation of that power.
#
# Finding:
# The asymp. approx can be fairly poor with outliers and moderate sample sizes.


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library( tidyverse )
source( "data_simulators.R" )





# Run a simulated experiment and following t-test 1000 times and calculate
# rejection, etc., of the simulated experiments.
run_single = function( N, tau, p.extreme, p.tx = 0.5 ) {
    rps = plyr::rdply( 1000, {
        dat = make.obs.data( N=N, p=p.tx, tau=tau, p.extreme=p.extreme )
        #dat = make.obs.data( N=N, p=0.5, DGP = make.data.norm, tau=tau )

        dd = dat %>% group_by( Z ) %>%
            summarize( sd = sd( Yobs ),
                       mean = mean( Yobs ),
                       n = n())

        tau.hat = diff( dd$mean )
        SE = with( dd, sqrt(  sum( sd^2 / n ) ) )

        tt = t.test(Yobs ~ Z, data = dat, alternative = "less")

        data.frame( tau.hat = tau.hat, z = tau.hat / SE, pv.ttest = tt$p.value )
    })
    rps
}


##### The Simulation #####
if ( FALSE ) {
    # Do simple check of power vs. theoretical power
    N = 100
    p.extreme = 0.2
    tau = 0.6
    p.tx = 0.3

    # Run the experiment
    rps = run_single( N=N, tau=tau, p.extreme=p.extreme, p.tx = p.tx )

    mean( rps$z >= qnorm( 0.95 ) )

    mean( rps$pv.ttest <= 0.05 )

    calc.power.extreme( tau, p.extreme, N, p.tx )

    calc.power( DGP = make.data.mix, tau=tau, p.extreme=p.extreme, N=N, p.tx=p.tx )

}