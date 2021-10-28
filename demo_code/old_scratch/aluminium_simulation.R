##
## Conduct a simulation using the aluminium data
##
## As the other simulations, this will compare different tx-co testing
## approaches for detecting impacts when there are rare but large impacts.
##

library(tidyverse)
library(coin)
library(perminf)
library(directlabels)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

source("data_simulators.R")
source( "aluminium_example.R" )

RUN_SIMULATION = TRUE

R = 1000

p = 0.1

                                        # Type I error rate
alpha = 0.10


make.data.aluminium = function( N, shift = 0 ) {
    smp = sample_n( simdat, N, replace=TRUE )
    tweaks = rnorm( N, mean=0, sd=0.02 )
    smp = transmute( smp,
                    Y.0 = Y0 + tweaks,
                    Y.1 = Y1 + tweaks - shift,
                    tau = Y.1 - Y.0 )
    smp
}

if ( FALSE ) {
    make.data.aluminium( 5 )
    make.data.exp( 5 )
}



                                        # Make some data, analyze it and get some pvalues from the various testing
                                        # approaches.
                                        #
                                        # @return pvalues and some test statistics as a data.frame (With 1 row).
one.run = function( N, ..., DGP) {
    dat = make.obs.data( N, p=0.5, ..., DGP=DGP )

    dat$Z = factor(dat$Z, levels=c(1,0), labels=c("T","C") )

    steph <- stephenson_test(
        formula = Yobs ~ Z,
        data = dat,
        alternative = "greater",
        subset_size = 6,
        ## distribution = "asymptotic"
        approximate(B = 1000)
    )

                                        # steph.left <- stephenson_test(
                                        #     formula = Yobs ~ Z,
                                        #     data = dat,
                                        #     alternative = "greater",
                                        #     tail = "left",
                                        #     subset_size = 6,
                                        #     distribution = "asymptotic"
                                        #     ## approximate(B = 1000)
                                        # )


    wil = wilcox_test(
        formula = Yobs ~ Z,
        data = dat,
        alternative = "greater",
        ## distribution = "asymptotic"
        ## distribution = "exact"
        approximate(B = 1000)
    )

    dmn = independence_test(
        formula = Yobs ~ Z,
        data = dat,
        alternative = "greater",
        ## distribution = "asymptotic"
        ## distribution = "exact"
        approximate(B = 1000)
    )

    tt = t.test(Yobs ~ Z, data = dat, alternative = "greater")

    data.frame(
        big.Y1 = sum( dat$Y.1 > max( dat$Y.0 ) ),
        Ybar.0 = tt$estimate[[1]],
        Ybar.1 = tt$estimate[[2]],
        tau.hat = as.numeric(diff(tt$estimate)),
        t.test = tt$p.value,
        stephenson = as.numeric(pvalue(steph)),
                                        #     stephenson.left = as.numeric(pvalue(steph.left)),
        perm.mean = as.numeric(pvalue(dmn)),
        wilcox = as.numeric(pvalue(wil))
    )
}

                                        # Call one.run() multiple times and tidy up the results
one.simulation = function(N, ..., DGP, R = 100)  {
    scat("Sim: N=%d, R=%d\n", N, R)
    scat( "\tParam: %s\n", paste( names(list(...)), list(...), sep="=", collapse="; " ) )

    rps = plyr::rdply(R, one.run(N = N, ..., DGP=DGP))
    head(rps)

    rs = gather(rps,
                stephenson,
                                        #stephenson.left,
                wilcox,
                perm.mean,
                t.test,
                key = "test",
                value = "p.value")
    rs
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
}


                                        # Run simulation for each of the multifactor combinations listed in levels (so
                                        # it will run nrow(levels) separate experiments).
                                        #
                                        # Results will be saved in file 'filename' with a time stamp.
run_a_simulation = function( levels, DGP, filename, R=10 ) {

    scat( "Running simulation with %d levels\n", nrow( levels ) )

    theo = mutate( levels,
                  power = pmap( levels, calc.power, DGP=DGP, p.tx=0.5 ) ) %>%
        unnest()
    theo$test = "theoretical"
    theo = rename( theo, power = pow )
    theo

                                        # Run the simulation
    res = mutate( levels,
                 results = pmap(levels, one.simulation, DGP = DGP, R = R) )

    res
    r2 = unnest(res)
    r2


                                        # Aggregate results
    group = names( levels )
    p3<- r2 %>% purrrlyr::slice_rows(group)
    r2.sum = p3 %>% group_by( test, add=TRUE ) %>%
        summarise(power = mean(p.value <= alpha),
                  SE.pow = sqrt( power * (1-power) / n() ),
                  mean.big.Y1 = mean( big.Y1) )

    head( r2.sum )


    r2.sum = bind_rows( r2.sum, select( theo, -mu0, -mu1, -sd0, -sd1, -tau.true, -SE ) )

    FILE_NAME = paste( "simulation_study_results_", filename, "_", format(Sys.time(), "%Y_%m_%d_%H_%M"), ".RData", sep="" )
    scat( "Saving simulation results to '%s'\n", FILE_NAME )
    save.image( file=FILE_NAME )

    invisible( list( results = res, summary=r2.sum ) )
}



if ( RUN_SIMULATION ) {
    scat( "Running the Aluminium simulation\n" )

                                        #levels = expand.grid( N = c( 10, 20, 30, 40, 50 ), tau.lambda = c( 0.5, 1, 2 ) )
    levels = expand.grid( N = c(25, 35, 50, 75, 100, 150 ), shift=c(0, 0.3, 0.6) )
                                        #levels = expand.grid( N = c(25, 50, 100 ), shift=c(0, -0.3, -0.6) )
    levels
    levels = as.tibble(levels)
    scat( "Aluminium Experiment: Running %d experiments with %d reps each\n",
         nrow( levels ), R  )

                                        #theo = mutate( levels, power = pmap_dbl( list( tau, p.extreme, N ), calc.power ) )
    sim.res = run_a_simulation( levels, DGP=make.data.aluminium,
                               filename="aluminium", R=R )
    sim.res

    r2.sum = sim.res$summary
    r2.sum

    r2.sum = filter( r2.sum, test != "theoretical" )

    library( ggthemes )


    my_t = theme_tufte() + theme( legend.position="bottom",
                                 legend.direction="horizontal", legend.key.width=unit(1,"cm"),
                                 panel.border = element_blank() )
    theme_set( my_t )
    r2.sum$shift = as.factor( r2.sum$shift )

                                        # Plot power curves
    head( r2.sum )

    revalues <- c(t.test = "Neyman t",
                  perm.mean = "permutation mean diff.",
                  wilcox = "Wilcoxon",
                  stephenson = "Stephenson")

    plt <- r2.sum %>%
        mutate(Test = plyr::revalue(test, revalues),
               Test = factor(Test, revalues),
               ATE = paste("ATE =", .61 - as.numeric(as.character(shift)))) %>%
        ggplot() +
        aes(x = N, y = power, group = Test, col = Test, linetype = Test,
            shape = Test) +
        facet_wrap( ~ ATE) +
        scale_y_continuous(breaks = seq(0, 1, .1)) +
    geom_line() +
        geom_point() +
        ## geom_errorbar( aes( ymin = power - 2*SE.pow, ymax = power + 2*SE.pow ),
        ##               width=0.1) +
        labs( x = "Sample Size", y = "Power", color = NULL, linetype = NULL,
             shape = NULL) +
        geom_hline( yintercept=c(0 ), lty=c(1) )
    plt


    p2 <- plt +
        geom_dl(aes(label = test),
                method = list(dl.combine("first.points", "last.points"),
                              cex = 0.8)) +
        coord_cartesian( xlim=c(1, 250 ) ) +
                                        # scale_x_log10()   +
        labs( title="Power as a function of sample size" )

    print( p2 )

    ggsave( "aluminium_simulation.pdf", plt, width = 6, height = 4 )

}
