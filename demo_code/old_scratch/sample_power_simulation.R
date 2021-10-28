##
## Conduct a simulation to compare different tx-co testing approaches for
## detecting impacts when there are rare but large impacts.
##

library(tidyverse)
library(coin)
library(perminf)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

source("data_simulators.R")

RUN_SIMULATION_T = FALSE
RUN_SIMULATION_MIX = FALSE
RUN_SIMULATION_EXP_CDF = FALSE
RUN_SIMULATION_BENIN = TRUE

SAVE_PLOT = TRUE

CONTROL_DISTRIBUTION = "normal"  # or logistic

R = 1000
p = 0.1

# Type I error rate
alpha = 0.10


# Make some data, analyze it and get some pvalues from the various testing
# approaches.
#
# @return pvalues and some test statistics as a data.frame (With 1 row).
one.run = function( N, ..., distribution, DGP) {
    dat = make.obs.data( N, p=0.5, ..., distribution = distribution, DGP=DGP )

    dat$Z = factor(dat$Z, levels=c(1,0), labels=c("T","C") )

    steph <- stephenson_test(
        formula = Yobs ~ Z,
        data = dat,
        alternative = "greater",
        subset_size = 6,
        ## distribution = "asymptotic"
        ## distribution = "exact"
        distribution = approximate(B = 1000)
    )

    steph.left <- stephenson_test(
        formula = Yobs ~ Z,
        data = dat,
        alternative = "less",
        tail = "left",
        subset_size = 6,
        ## distribution = "asymptotic"
        ## distribution = "exact"
        distribution = approximate(B = 1000)
    )


    wil = wilcox_test(
        formula = Yobs ~ Z,
        data = dat,
        alternative = "greater",
        ## distribution = "asymptotic"
        ## distribution = "exact"
        distribution = approximate(B = 1000)
    )

    dmn = independence_test(
        formula = Yobs ~ Z,
        data = dat,
        alternative = "greater",
        ## distribution = "asymptotic"
        ## distribution = "exact"
        distribution = approximate(B = 1000)
    )

    tt = t.test(Yobs ~ Z, data = dat, alternative = "greater")

    data.frame(
        big.Y1 = sum( dat$Y.1 > max( dat$Y.0 ) ),
        Ybar.0 = tt$estimate[[1]],
        Ybar.1 = tt$estimate[[2]],
        tau.hat = as.numeric(diff(tt$estimate)),
        t.test = tt$p.value,
        stephenson = as.numeric(pvalue(steph)),
        stephenson.left = as.numeric(pvalue(steph.left)),
        perm.mean = as.numeric(pvalue(dmn)),
        wilcox = as.numeric(pvalue(wil))
    )
}

# Call one.run() multiple times and tidy up the results
one.simulation = function(N, ..., distribution, DGP, R = 100)  {
    scat("Sim: N=%d, R=%d\n", N, R)
    scat( "\tParam: %s\n", paste( names(list(...)), list(...), sep="=", collapse="; " ) )

    rps = plyr::rdply(R, one.run(N = N, ..., distribution=distribution, DGP=DGP))
    head(rps)

    rs = gather(rps,
                stephenson,
                stephenson.left,
                wilcox,
                perm.mean,
                t.test,
                key = "test",
                value = "p.value")
    return(rs)
    gc()
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
                   power = pmap( levels, calc.power, DGP=DGP,
                   distribution=CONTROL_DISTRIBUTION, p.tx=0.5 ) ) %>%
        unnest()
    theo$test = "theoretical"
    theo = rename( theo, power = pow )
    theo

    # Run the simulation
    res = mutate( levels,
                 results = pmap(levels, one.simulation, DGP = DGP,
                                distribution=CONTROL_DISTRIBUTION,
                                R = R) )

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

    FILE_NAME = paste( "simulation_study_results_", filename, "_", CONTROL_DISTRIBUTION, "_", format(Sys.time(), "%Y_%m_%d_%H_%M"), ".RData", sep="" )
    scat( "Saving simulation results to '%s'\n", FILE_NAME )
    save.image( file=FILE_NAME )

    invisible( list( results = res, summary=r2.sum ) )
}

####Power sim with the t-distribution on the treatment impacts. ####

if (RUN_SIMULATION_T) {

    #levels = expand.grid( N = c( 10, 20, 30, 40, 50 ), tau.lambda = c( 0.5, 1, 2 ) )
    levels = expand.grid(N = c(20, 40, 80, 160, 320),
                         tau = c( -0.10, 0.50 ),
                         scale = c( 0, 0.5, 2.0 ),
                         df = c( 3 ) )
    levels
    levels = as.tibble(levels)
    scat( "Running %d experiments with %d reps each\n", nrow( levels ), R  )


    #theo = mutate( levels, power = pmap_dbl( list( tau, p.extreme, N ), calc.power ) )
    sim.res = run_a_simulation( levels, DGP=make.data.t, filename="t", R=R )
    r2.sum = sim.res$summary

    head( r2.sum )

    r2.sum = mutate( r2.sum, sd.ratio = round( sqrt( 1 + scale ), digits=1 ) )

    # Plot power curves
    plt <- ggplot(r2.sum, aes(N, power, group = test, col = test)) +
        facet_grid( sd.ratio ~ tau + df, labeller=label_both) +
        geom_line() + geom_point() +
        geom_errorbar( aes( ymin = power - 2*SE.pow, ymax = power + 2*SE.pow ), width=0.1, alpha=0.25) +
        geom_hline( yintercept=c(0 ) ) +
        geom_hline( yintercept=alpha, lty=2, alpha=0.7 ) +
        scale_x_log10()
    print( plt )

    if (SAVE_PLOT) ggsave( "power_simulation_t.pdf", plt, width = 6, height = 6)
}



#### Mix (positive rare effect) simulation ####

# Looking at Mix DGP
if (FALSE) {
    p = 0.5
    p.extreme = 0.20
    tau = 1.5
    #DGP = make.data.mix
    DGP = make.data.const.rare
    #DGP = make.data.uni.uni

    #dat = make.obs.data( 1000, p=0.5, p.extreme=0.05, tau=0.4 )
    dat = make.obs.data( 10000, p=p, p.extreme=p.extreme, tau=tau, DGP=DGP )

    mean( dat$Y.1 - dat$Y.0 )

    head( dat )
    ggplot(dat, aes(x = Yobs, group = Z, col=as.factor(Z))) + geom_density()
    ggplot(dat, aes(x = Y.0, y=Y.1) ) + geom_point()

    analyze.power( dat, N=100, p.tx=0.5 )

    # Proportion of obs where tx P.O. is larger than _all_ Y(0)
    mean( dat$Y.1 > max( dat$Y.0 ) )
    sum( dat$Y.1 > max( dat$Y.0 ) )
    nrow(dat)


    # How many extreme ones do we get and what do they look like?
    dat = make.obs.data( 50, p=p, p.extreme=p.extreme, tau=tau, DGP = DGP )
    boxplot( Yobs ~ Z, data=dat )
    ggplot( dat, aes( x=Yobs, group=Z ) ) +
        facet_wrap( ~ Z, ncol=1 ) +geom_dotplot( size=0.25)
    dat$tau
    qplot( dat$tau )
    sort( dat$tau[ dat$tau != 0 ] )
    table( dat$tau > 0, dat$Z )
}

if (RUN_SIMULATION_MIX) {

    levels = expand.grid(N = c( 20, 80, 320 ), #N = c(20, 40, 80, 160, 320),
                         tau = c(1.0, 1.5, 3, 6),
                         p.extreme = c(0.05, 0.10, 0.20, 0.40) )
    levels
    levels = as.tibble(levels)
    DGP = make.data.const.rare

    #levels$tau = levels$tau * levels$p.extreme

    # Run the simulation
    res = run_a_simulation( levels, DGP=DGP, filename = "mix", R=R )
    res

    head( res$summary )

    sum = res$summary %>% ungroup()

    # undo tx rescaling (no longer needed)
    # sum <- sum %>% mutate( tau = tau / p.extreme )

    # Plot power curves
    plt <- ggplot( sum, aes(N, power, group = test, col = test)) +
        facet_grid( tau ~ p.extreme, labeller=label_both) +
        geom_line() + geom_point() +
        geom_errorbar( aes( ymin = power - 2*SE.pow, ymax = power + 2*SE.pow ), width=0.1, alpha=0.25) +
        geom_hline( yintercept=c(0,1), alpha=0.7, size=0.5 ) +
        geom_hline( yintercept=c(alpha, 0.8), lty=2, alpha=0.7, size=0.5 ) +
        scale_x_log10( breaks=unique( res$summary$N ) )

    library( ggthemes )

    my_t = theme_tufte() + theme( legend.position="bottom",
                                  legend.direction="horizontal", legend.key.width=unit(1,"cm"),
                                  panel.border = element_blank() )
    theme_set( my_t )
    print( plt )
    if (SAVE_PLOT) {
        ggsave( "power_simulation_mix_normal.pdf", plt, width = 6, height = 4 )
    }
}

#### Exponentiated CDF a la Rosenbaum ####


if (FALSE) {
    p = 0.5
    p.extreme = 0.20
    power = 6 #DGP = make.data.mix
    DGP = make.data.exp.cdf

    #dat = make.obs.data( 1000, p=0.5, p.extreme=0.05, tau=0.4 )
    dat = make.obs.data( 10000, p=p, p.extreme=p.extreme, m = power, DGP=DGP )

    mean( dat$Y.1 - dat$Y.0 )

    head( dat )
    ggplot(dat, aes(x = Yobs, group = Z, col=as.factor(Z))) + geom_density()
    ggplot(dat, aes(x = Y.0, y=Y.1) ) + geom_point()
    ggplot(dat, aes(x = Y.0, y=tau) ) + geom_point()

    analyze.power( dat, N=100, p.tx=0.5 )

    # Proportion of obs where tx P.O. is larger than _all_ Y(0)
    mean( dat$Y.1 > max( dat$Y.0 ) )
    sum( dat$Y.1 > max( dat$Y.0 ) )
    nrow(dat)


    # How many extreme ones do we get and what do they look like?
    dat = make.obs.data( 50, p=p, p.extreme=p.extreme, tau=tau, DGP = DGP )
    boxplot( Yobs ~ Z, data=dat )
    ggplot( dat, aes( x=Yobs, group=Z ) ) +
        facet_wrap( ~ Z, ncol=1 ) +geom_dotplot( size=0.25)
    dat$tau
    qplot( dat$tau )
    sort( dat$tau[ dat$tau != 0 ] )
    table( dat$tau > 0, dat$Z )
}

if (RUN_SIMULATION_EXP_CDF) {

    levels = expand.grid(N = c( 80, 320 ), #N = c(20, 40, 80, 160, 320),
                         p.extreme = c(0.10, 0.40),
                         m = c( 2, 4, 6) )
    levels
    levels = as.tibble(levels)
    DGP = make.data.exp.cdf
    nrow( levels )
    #levels$tau = levels$tau * levels$p.extreme

    # Run the simulation
    res = run_a_simulation( levels, DGP=DGP, filename = "expCDF", R=R )
    res

    head( res$summary )

    sum = res$summary %>% ungroup()

    # undo tx rescaling (no longer needed)
    # sum <- sum %>% mutate( tau = tau / p.extreme )

    # Plot power curves
    plt <- ggplot( sum, aes(N, power, group = test, col = test)) +
        facet_grid( m ~ p.extreme, labeller=label_both) +
        geom_line() + geom_point() +
        geom_errorbar( aes( ymin = power - 2*SE.pow, ymax = power + 2*SE.pow ), width=0.1, alpha=0.25) +
        geom_hline( yintercept=c(0,1), alpha=0.7, size=0.5 ) +
        geom_hline( yintercept=c(alpha, 0.8), lty=2, alpha=0.7, size=0.5 ) +
        scale_x_log10( breaks=unique( res$summary$N ) )

    library( ggthemes )

    my_t = theme_tufte() + theme(legend.position="bottom",
                                 legend.direction="horizontal",
                                 legend.key.width=unit(1,"cm"),
                                 panel.border = element_blank() )
    theme_set( my_t )
    print( plt )
    if (SAVE_PLOT) {
        ggsave( "exp_cdf_simulation_mix.pdf", plt, width = 6, height = 4 )
    }
}




#### Benin Simulation #####

if ( RUN_SIMULATION_BENIN ) {
    scat( "Running the Benin simulation\n" )

    # Get some parameters from benin data
    benin <- perminf::benin %>%
        select(Village, District, Candidate, Treatment, VoteShare)
    benin$District = factor( benin$District )
    benin_sub = benin %>% filter(Treatment != "Clientelist")
    benin_sub

    rs = benin_sub %>% group_by( Treatment ) %>%
        summarize( mean.vote = mean( VoteShare ),
                   sd.vote = sd( VoteShare ) )

    rs
    diff( rs$mean.vote )
    v = rs$sd.vote^2
    benin.var.ratio = v[[2]] / v[[1]]
    benin.var.ratio


    #levels = expand.grid( N = c( 10, 20, 30, 40, 50 ), tau.lambda = c( 0.5, 1, 2 ) )
    levels = expand.grid(N = c(16, 160),
                         tau = c(0), ## c( -0.06 )
                         scale = sqrt( seq( 0, 9, by=0.5 ) ),
                         df =  c(.Machine$double.xmax
                                        #3
                                 )
                         )
    levels
    levels = as.tibble(levels)
    scat( "Benin Experiment: Running %d experiments with %d reps each\n",
         nrow( levels ), R  )

    #theo = mutate( levels, power = pmap_dbl( list( tau, p.extreme, N ), calc.power ) )
    sim.res = run_a_simulation( levels, DGP=make.data.t, filename="benin", R=R )
    r2.sum = sim.res$summary

    # calculate ratio of treatment outcomes and control outcomes (given
    # population model)
    # This is a measure of the heterogeniety.
    r2.sum$sd.ratio = sqrt( 1 + r2.sum$scale^2 )


    names( sim.res )
    res = sim.res$results
    head( res )
    head( res$results[[1]] )

    r2 = unnest(res)
    head( r2 )

    r2 = select( r2, N, tau, scale, .n, test, p.value )
    r2 = spread( r2, test, p.value )
    head( r2 )

    double.pow = r2 %>%
        group_by( N, tau, scale ) %>%
        summarise(stephenson_right = mean( stephenson <= alpha ),
                  stephenson_left = mean( stephenson.left <= alpha ),
                  ## I assume that alpha is divided by two to control the FWER
                  stephenson_both =
                      mean( stephenson <= alpha & stephenson.left <= alpha ),
                  stephenson_either =  mean( stephenson <= alpha/2 |
                                  stephenson.left <= alpha/2 ),
                  t = mean(t.test <= alpha),
                  wilcoxon = mean(wilcox <= alpha)
                  ) %>%
        mutate(sd_ratio = sqrt( 1 + scale^2 ))
    head( double.pow )
    unique( double.pow$N )
    unique( double.pow$tau )
    max( double.pow$stephenson_both )

    double.pow = double.pow %>%
        gather(stephenson_right:wilcoxon, key="measure", value="power" )
    head( double.pow )

    revalues <- c(t = "Neyman t",
                  wilcoxon = "Wilcoxon",
                  stephenson_right = "Stephenson (positive)",
                  stephenson_left = "Stephenson (negative)",
                  stephenson_both = "Stephenson (both)")

    library( ggthemes )

    my_t = theme_tufte() +
        theme( legend.position="bottom",
              legend.direction="horizontal",
              legend.key.width=unit(1,"cm"),
              panel.border = element_blank() )
    theme_set( my_t )

    double.pow %>%
        filter(measure %in% names(revalues)) %>%
        mutate(sample_size = paste("N =", N),
               Test = plyr::revalue(measure, revalues),
               Test = factor(Test, revalues)) %>%
        ggplot() +
        aes( x = sd_ratio, y = power, col = Test, linetype = Test) +
        geom_line() +
        geom_hline( yintercept=c(0), lty=c(1) ) +
        geom_vline( xintercept = sqrt( benin.var.ratio ), col="grey",
                   linetype = "dashed") +
        annotate("text", x = sqrt( benin.var.ratio ) + .01, y = .6,
                 label = "observed SD ratio", hjust = 0, size = rel(2)) +
        facet_wrap(~sample_size) +
        labs( x = "Ratio of Standard Deviations", y = "Power",
             color = NULL, linetype = NULL) +
        scale_y_continuous(breaks = seq(0, 1, .1))
    if (SAVE_PLOT) {
        ggsave( "benin_simulation2.pdf", width = 7.5, height = 3.5)
    }


    double.pow %>%
        filter(sd_ratio == 1)

    double.pow %>%
        filter(sd_ratio >= sqrt( benin.var.ratio ) - .1 &
               sd_ratio <= sqrt( benin.var.ratio ) + .1) %>%
        group_by(N, measure) %>%
        summarise(power = mean(power))

    r2.sum = filter( r2.sum, test != "theoretical", test != "perm.mean" )

    # Plot power curves
    plt <- r2.sum %>%
        ungroup() %>%
        mutate(Test = factor(test, labels = c("Stephenson (right tail)",
                                              "Stephenson (left tail)",
                                              "Neyman t", "Wilcoxon")),
               sample_size = paste("N =", N)) %>%
        ggplot() +
        aes(x=sd.ratio, y=power, group = Test, col = Test, linetype = Test,
            shape = Test) +
        facet_wrap(~sample_size) +
        geom_line() +
        geom_point() +
        ## geom_errorbar( aes( ymin = power - 2*SE.pow, ymax = power + 2*SE.pow ),
        ##                width=0.025) +
        geom_hline( yintercept=c(0,0.10 ), lty=c(1,2,1,2) ) +
        geom_vline( xintercept = sqrt( benin.var.ratio ), col="red" ) +
        annotate("text", x = sqrt( benin.var.ratio ) + .01, y = .3,
                 label = "observed SD ratio", hjust = 0) +
        labs( x = "Ratio of Standard Deviations (Tr/Co)", y = "Power") +
        scale_y_continuous(breaks = seq(0, 1, .1))
    print( plt )

    if (SAVE_PLOT) {
        ggsave( "benin_simulation.pdf", plt, width = 6, height = 4 )
    }
}
