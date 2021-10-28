##
## Simulation to examine the coverage and performance of confidence intervals
## using different test statistics.
##
##
## This file looks at the coverage of Confidence Intervals for the both the ATE
## and the maximal effect.  We should have invalid converage for the ATE, and
## serious overcoverage for the maximal effect.
##
##
## Results:
## So far, coverage for the ATE is good for everything, if not a bit
## conservative.  Sad.
##

library(tidyverse)
library(coin)
library(perminf)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("data_simulators.R")


# Type I error rate
alpha = 0.10

RUN_SIMULATION_CI_MIX = TRUE

R = 100
#R = 1000
p = 0.5



# Make some data, analyze it and get some pvalues from the various testing
# approaches.
#
# @return pvalues and some test statistics as a data.frame (With 1 row).
one.run = function( N, ..., DGP, how.run =  approximate(10000)) {
    dat = make.obs.data( N, p=0.5, ..., DGP=DGP )

    dat$Z = factor(dat$Z, levels=c(1,0), labels=c("T","C") )

    steph <- find_ci_stephenson(
        data = dat,
        tr_var = "Z", y_var = "Yobs",
        alternative = "greater",
        subset_size = 6,
        distribution = how.run,
        precision = 0.01
    )

    wil = find_ci(
        data = dat,
        tr_var = "Z", y_var = "Yobs",
        coin_test = wilcox_test,
        alternative = "greater",
        distribution = how.run,
        precision = 0.01
    )

    dmn = find_ci(
        data = dat,
        tr_var = "Z", y_var = "Yobs",
        coin_test = independence_test,
        alternative = "greater",
        distribution = "asymptotic"
    )

    tt = t.test(Yobs ~ Z, data = dat, alternative = "greater")

    data.frame(
        max.tau = max( dat$tau ),
        ATE = mean( dat$tau ),
        Ybar.0 = tt$estimate[[1]],
        Ybar.1 = tt$estimate[[2]],
        CI.steph = steph$ci[[1]],
        CI.wilcox = wil$ci[[1]],
        CI.diff = dmn$ci[[1]],
        CI.t = tt$conf.int[[1]],
        prop.above.steph = mean( dat$tau > steph$ci[[1]] ),
        prop.above.wilcox = mean( dat$tau > wil$ci[[1]] ),
        prop.above.diff = mean( dat$tau > dmn$ci[[1]] ),
        prop.above.t = mean( dat$tau > tt$conf.int[[1]] )
    )
}

if ( FALSE ) {
    one.run( N = 50, DGP = make.data.mix)
}


scat = function(str, ...) {
    cat(sprintf(str, ...))
}


one.simulation = function(N, ..., DGP, R = 100, how.run = approximate(10000) )  {
    scat("Sim: N=%d, R=%d\n", N, R)
    scat( "\tParam: %s\n", paste( names(list(...)), list(...), sep="=", collapse="; " ) )

    rps = plyr::rdply(R, one.run(N = N, ..., DGP=DGP, how.run=how.run ))

    rps
}


# Run simulation for each of the multifactor combinations listed in levels (so
# it will run nrow(levels) separate experiments).
#
# Results will be saved in file 'filename' with a time stamp.
run_a_simulation = function( levels, DGP, filename, R=10, how.run = approximate(10000) ) {

    theo = mutate( levels,
                   power = pmap( levels, calc.power, DGP=DGP, p.tx=0.5 ) ) %>%
        unnest()
    theo$method = "theoretical"
    theo = rename( theo, power = pow )
    theo

    # Run the simulation
    res = mutate( levels,
                  results = pmap(levels, one.simulation, DGP = DGP, R = R, how.run = how.run ) )

    res
    r2 = unnest(res)
    r2

    cat( "Saving simulation results\n" )
    FILE_NAME = paste( "simulation_study_results_", filename, "_", format(Sys.time(), "%Y_%m_%d_%H_%M"), ".RData", sep="" )
    save.image( file=FILE_NAME )

    # name of the factors of the experiment
    group = names( levels )

    # Aggregate results
    pvs = grep( "prop.above.", names( r2 ) )
    taus = grep( "CI.", names( r2 ) )

    names = gsub( "CI.", "", names(r2)[taus] )
    r2b = reshape( as.data.frame( r2 ),  idvar=c(".n", group),
                    varying=list( taus, pvs ),
                    timevar = "method",
                    times = names,
                    v.names=c("CI","prop.above"), direction="long")
    head( r2b )

    p3<- r2b %>% purrrlyr::slice_rows(group)
    r2.sum = p3 %>% group_by( method, add=TRUE ) %>%
        summarise(power = mean( 0 < CI ),
                  coverage.ATE = mean( CI <= ATE ),
                  coverage.max = mean( CI <= max.tau ),
                  mean.max.tau = mean( max.tau ),
                  mean.CI = mean( CI ),
                  mean.prop.above = mean( prop.above ) )
   # head( r2.sum )
  #  head( theo )
    #r2.sum = merge( r2.sum, select( theo, -mu0, -mu1, -sd0, -sd1, -tau1, -SE, -test),
    #                by = group,
    #                all.x = TRUE )
    r2.sum = bind_rows( r2.sum, select( theo, -mu0, -mu1, -sd0, -sd1, -tau1, -SE ) )
   # head( r2.sum )
   # tail( r2.sum )
    invisible( list( results = res, summary=r2.sum ) )
}


# Testing one.simulation() and Looking at DGP
if (FALSE) {
    #dat = make.obs.data( 1000, p=0.5, p.extreme=0.05, tau=0.4 )
    dat = make.obs.data( 1000, p=0.5, p.extreme=0.05, tau=0.4, DGP=make.data.mix )
    dat = make.obs.data( 1000, p=0.5, tau.lambda=1, DGP=make.data.exp )
    dat = make.obs.data( 1000, p=0.5, tau = 0.1, df=3, scale=1, DGP=make.data.t )

    head( dat )
    ggplot(dat, aes(x = Yobs, group = Z, col=as.factor(Z))) + geom_density()
    ggplot(dat, aes(x = Y.0, y=Y.1) ) + geom_point()

    analyze.power( dat, N=100, p.tx=0.5 )

    # Proportion of obs where tx P.O. is larger than _all_ Y(0)
    mean( dat$Y.1 > max( dat$Y.0 ) )
    sum( dat$Y.1 > max( dat$Y.0 ) )
    nrow(dat)

    # Hack check
    dat = make.obs.data( 200, p=0.5, p.extreme=0.05, tau=0.4, DGP=make.data.mix )
    steph <- find_ci_stephenson(
        data = dat,
        tr_var = "Z", y_var = "Yobs",
        alternative = "greater",
        subset_size = 6,
        distribution =approximate(10000)
    )
    steph$ci
    steph

    # Check CI generation
    steph <- find_ci_stephenson(
        data = dat,
        tr_var = "Z", y_var = "Yobs",
        alternative = "greater",
        subset_size = 6,
        distribution =approximate(10000)
    )

    wil = find_ci(
        data = dat,
        tr_var = "Z", y_var = "Yobs",
        coin_test = wilcox_test,
        alternative = "greater",
        distribution = approximate(10000)
    )


    # Check simulation code
    one.run(N = 100, DGP=make.data.exp, tau.lambda = 1)

    one.run(N = 100, DGP=make.data.mix, p.extreme=0.05)

    one.simulation(N = 10, DGP=make.data.mix, p.extreme=0.05, R=10)

    # How many extreme ones do we get and what do they look like?
    dat = make.obs.data( 50, p=0.5, p.extreme=0.10, tau=0.4 )
    boxplot( Yobs ~ Z, data=dat )
    ggplot( dat, aes( x=Yobs, group=Z ) ) +
        facet_wrap( ~ Z, ncol=1 ) +geom_dotplot()
    dat$tau
    qplot( dat$tau )
    dat$tau[ dat$tau != 0 ]
    table( dat$tau > 0, dat$Z )


    # check run_simulation
    one.simulation( N=20, tau=1.0, p.extreme=0.10, DGP=make.data.mix, R= 10, how.run="asymptotic" )
}


#### Mix (positive rare effect) simulation ####

if (RUN_SIMULATION_CI_MIX) {
    levels = expand.grid(N = c(80, 160, 320),
                         tau = c(0.5, 1.0),
                         p.extreme = c(0.10, 0.40) )
    DGP = make.data.const.rare

    #DGP = make.data.mix


   # levels = expand.grid(N = c(20, 100),
    #                     tau = c( 1.0),
   #                      p.extreme = c(0.10, 0.20))
    levels
    levels = as.tibble(levels)

    scat( "Running %d experiments with %d reps each\n", nrow( levels ), R )

    # Run the simulation
    res = run_a_simulation( levels, DGP=DGP, filename = "CI_mix", R=R )

    res$summary

    # Plot power
    ggplot(res$summary, aes(N, power, group = method, col = method)) +
        facet_grid( tau ~ p.extreme, labeller=label_both) +
        geom_line() + geom_point() +
        geom_hline( yintercept = 0.90, lty=2, alpha=0.7 )

    # Plot coverage of ATE
    ggplot(res$summary, aes(N, coverage.ATE, group = method, col = method)) +
        facet_grid( tau ~ p.extreme, labeller=label_both) +
        geom_line() + geom_point() +
        geom_hline( yintercept = 0.90, lty=2, alpha=0.7 )

    # how tight the CI
    ggplot(res$summary, aes(N, mean.CI, group = method, col = method)) +
        facet_grid( tau ~ p.extreme, labeller=label_both) +
        geom_line() + geom_point()


}

