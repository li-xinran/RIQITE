#
# Code that implements teacher professional development analysis in
# the main paper.
#
# This code is SOMEWHAT OUT OF DATE.  See the coe in the RIllustration
# file for the final analysis presented in the actual paper. (This
# file does have some misc. extensions and detours that may be of
# interest, however.)
#


## LIBRARIES
library(tidyverse)
library( RIQITE )

## DATA
data( "electric_teachers" )
dat = electric_teachers

set.seed( 1010101 )

# scramble order of data
dat = sample_n( dat, nrow(dat) )


ggplot( dat, aes( gain ) ) +
    facet_wrap( ~ Tx ) +
    geom_histogram() +
    geom_vline( xintercept = 0, col="red" )


# Classic OLS analysis with fixed effects for site: What can we say
# about the average effect?
mod = lm( gain ~ 0 + TxAny + as.factor(Site), data=dat )
summary( mod )
confint( mod )
confint( mod, level = 0.90 )


#### Looking at variation of scores, and standardizing ####

nrow(dat)

# control side data
cdat = filter( dat, TxAny == 0 )
sd( dat$gain )
dat$std_gain = dat$gain / sd( cdat$gain )
sd( dat$std_gain )

dat %>% group_by( TxAny ) %>%
    summarise( sdGstd = sd( std_gain ),
               sdG = sd( gain ) )
summary( lm( std_gain ~ TxAny, data=dat ) )



#### Check power and performance for specific quantile ####

rho = 0
R = 100
cdat = filter( dat, TxAny == 0 )
EstTx = 0.65
TxVar = 0.5

rs <- explore_stephenson_s( s = c(2, 5, 10, 20),
                      n = nrow( dat ),
                      Y0_distribution = cdat$std_gain,
                      tx_function = "rexp",
                      scale_tx = 0.5, ATE = EstTx, rho = -0.45,
                      R = R, calc_ICC = TRUE,
                      targeted_power = TRUE, k.vec = 200 )

rs


#### Select best test s-value for test statistic  ####

cdat = filter( dat, TxAny == 0 ) # Get standardized outcome
EstTx = 0.65
TxVar = 0.5
R = 100

s_list = c( 2, 3, 4, 5, 6, 8, 10, 15, 20 )

checks = expand_grid( TxVar = c( 0.5, 1.0 ),
                      tx_dist = c( "constant", "rexp", "rnorm" ),
                      rho = c( -0.5, 0, 0.5 ) )
checks = filter( checks, (TxVar == 0.5 & rho == 0) | tx_dist != "constant" )
checks$TxVar[ checks$tx_dist == "constant" ] = 0
cat( "Number of scenarios:", nrow(checks), "\n" )

# Calculate power across all listed scenarios. Each check will use
# parallel, so no need to call parallel for this outer loop.
checks$data = pmap( checks, function( TxVar, tx_dist, rho ) {
    cat( glue::glue("Running {TxVar} {tx_dist} rho={rho}\n" ) )
    explore_stephenson_s( s = s_list,
                          n = nrow( dat ),
                          Y0_distribution = cdat$std_gain,
                          tx_function = tx_function_factory(tx_dist,
                                                            ATE = EstTx,
                                                            tx_scale=TxVar,
                                                            rho = rho ),
                          R = R, calc_ICC = TRUE, parallel = TRUE,
                          targeted_power = FALSE, k.vec = (233-100+1):233 )
} )
checks
s_selector = unnest( checks, cols = "data" )

s_selector
saveRDS( s_selector, here::here( "demo_code/heller_s_check_results.rds" ) )


#### Visualizing the results of selecting s ####

# Make plot of power curves under all the simulated scenarios to see
# what s is best.

s_selector = readRDS( here::here( "demo_code/heller_s_check_results.rds" ) )

s_selector <- s_selector %>%
    mutate( range = round( pow_h - pow_l, digits = 2 ),
            txd = as.numeric(as.factor(tx_dist)) - 2,
            s_adj = s * 1.05^txd ) %>%
    mutate( tx_dist = fct_recode(tx_dist,
                                 exponential = "rexp",
                                 normal = "rnorm" ) )

s_selector

# Peek at structure of single unit
s_selector %>%
    filter( s == 2, k == 220 )

n_units = nrow(dat)
n_units

s_selector = s_selector %>%
    filter( k > 233 - 35 ) %>%
    mutate( group = case_when( k > n_units - 5 ~ "high",
                               k > n_units - 15 ~ "med",
                               TRUE ~ "low" ) ) %>%
    mutate( group = factor( group, levels = c( "high", "med", "low" ) ) )

# The number of units in each of the high, med, and low groups
s_selector %>%
    filter( tx_dist == "constant", s==2, TxVar == 1 ) %>%
    pull( "group" ) %>%
    table()

avg = s_selector %>%
    filter( tx_dist != "constant" ) %>%
    group_by( s, group, rho ) %>%
    summarise( power = mean( power ),
               n = n() ) %>%
    rename( s_adj = s )
avg




# make plot shown in current work
n_units = nrow(electric_teachers)
s_selector = mutate( s_selector,
                     group = cut( k,
                                  breaks = c(0, 170, 200, 220, 233),
                                  labels = c( "v low", "low", "med", "high" ) ) )
summary( s_selector$k )
table( s_selector$group )

s_selector$rho.f = factor( s_selector$rho,
                           levels = c( -0.5, 0, 0.5 ),
                           labels = c( "rho = -0.5", "rho = 0", "rho = 0.5" ) )


s_selector %>% filter( rho == 0, tx_dist =="constant", TxVar == 0.5, k == 200 )

s_selector %>% filter( rho == 0, tx_dist =="constant", TxVar == 0.5, s==2 ) %>%
    pull( group ) %>% table()

avg = s_selector %>% group_by( s, group, rho.f ) %>%
    summarise( q_ci = median( q_ci ), .groups="drop" ) %>%
    filter( group != "v low", is.finite( q_ci ) )

s_selector %>%
    filter( group != "v low" ) %>%
    ggplot( aes( s, q_ci ) ) +
    facet_grid( rho.f ~ group ) +
    geom_hline( yintercept = 0 ) +
    geom_smooth( aes( col=tx_dist ), method = "loess", se = FALSE,
                 span = 1, lwd= 0.5 ) +
    labs( x = "s", y = "Median lower CI bound", col = "tx distribution:" ) +
    scale_x_log10( breaks = unique( s_selector$s ) ) +
    theme_minimal() +
    theme( panel.grid.minor = element_blank() ) +
    geom_point( data = avg, col="black", size=2 ) +
    coord_cartesian( ylim = c( -1, 2 ) ) +
    theme( #legend.position="bottom",
           #legend.direction="horizontal",
           legend.key.width=unit(1,"cm"),
           panel.border = element_rect(colour = "grey", fill=NA, linewidth=1) )
ggsave( filename = "demo_code/s_selector.pdf", width = 8, height = 3.5 )



# Make a separate plot not currently shown in the published work
library( ggthemes )
my_t = theme_tufte() +
    theme( legend.position="bottom",
           legend.direction="horizontal",
           legend.key.width=unit(1,"cm"),
           panel.border = element_blank() )
theme_set( my_t )

ggplot( s_selector, aes( s_adj, power, col=tx_dist ) ) +
    #facet_grid(  group ~ rho, labeller = label_both ) +
    facet_grid( rho ~ group ) +
    geom_hline(yintercept = 0.80, lty = 2 ) +
    # geom_point() +
    # geom_errorbar( aes( ymax = power + 2*SEpower,
    #                     ymin = power - 2*SEpower ),
    #                width = 0 ) +
    geom_smooth( se = FALSE, method = "loess", span = 1 ) +
    #geom_line( ) +
    labs( x = "s", y = "Power" ) +
    scale_x_log10( breaks = unique( s_selector$s ) ) +
    geom_line( data = avg, col="black" ) +
    coord_cartesian( ylim=c(0,1) ) +
    theme_minimal() +
    labs( col = "distribution:" ) +
    scale_y_continuous( breaks = c(0,0.2,0.4,0.6,0.8,1) ) +
    theme_tufte()
ggsave( filename = "demo_code/s_selector.pdf", width = 8, height = 3.5 )


# Find s with max average power for each group
avg %>% ungroup() %>%
    group_by( group ) %>%
    filter( power == max(power) )


# Looking at number of units with significant effects.
s_num_n <- s_selector %>%
    group_by( s, TxVar, tx_dist ) %>%
    summarise( n = mean( n ) )

ggplot( s_num_n, aes( s, n, col=as.factor(TxVar) ) ) +
    geom_smooth( se=FALSE )

# For number of significant units...
s_num_n %>%
    group_by( TxVar,tx_dist ) %>%
    filter( n == max(n) )



#
# We note for many quantiles the median lower bound to the CI is $-\infty$ across many if not all scenarios.
# To allow for aggregation, we impute a lower bound for these quantiles of -3.79, based on our calculation above:
#     ```{r}
# s_selector$q_ci[ !is.finite( s_selector$q_ci ) ] = -3.79
# ```
# Alternatively, we can just focus on quantiles where we always found non-infinite bounds.
#
#

# Check number of units in each group
# s_selector %>% dplyr::filter( s==2, tx_dist=="constant", TxVar == 0 ) %>%
#     pull( group ) %>% table()

#
# We can also identify which $s$ generally gave the most informative confidence intervals across the different groups as follows:
#     ```{r, warning=FALSE}
# avg %>% ungroup() %>%
#     group_by( rho.f, group ) %>%
#     filter( q_ci == max(q_ci) ) %>%
#     dplyr::select(-q_ci, -n) %>%
#     pivot_wider( names_from="rho.f", values_from="s" ) %>%
#     knitr::kable()
# ```



#### Analysis with the chosen test statistic of s=5  ####

method.list = list( name = "Stephenson", s = 5 )

n = nrow( dat )
n
n * 0.95

# test the null hypothesis that the 70% quantile of individual effect is
# less than or equal to 0
#
# Note: If k = n, then you are testing the max effect (sharp null).
pval = pval_quantile( Z = dat$TxAny, Y = dat$gain,
                      k = ceiling(n*0.70),
                      c = 0,
                      alternative = "greater",
                      method.list = method.list, nperm = 10^5 )
pval


##### construct simultaneous confidence intervals for all quantiles of individual effects #####

ci = ci_quantile( Z = dat$TxAny, Y = dat$gain,
                  alternative = "greater",
                  method.list = method.list,
                  nperm = 10^5,
                  alpha = 0.10 )

ci = mutate( ci, per = (k+0.5) / (nrow(ci)+1) )
filter( ci, sign(lower) > 0 )

plot( ci$lower )

# Quantile effects of different levels of effect.  E.g., the 65% is bounded
# below by 0.  The 84% is an impact of 10 or above.
ci2 = ci %>%
    filter( lower > 0 ) %>%
    mutate( lower = round( lower, digits=1 ) ) %>%
    group_by( lower ) %>%
    arrange( n ) %>%
    slice_head( n=1 ) %>%
    relocate( n, per, lower )
ci2

# NOTE: We are not following the randomization within site in the above.



# Try smaller set size.
method.list2 = list( name = "Stephenson", s = 5 )
ci.s5 = ci_quantile( Z = dat$TxAny, Y = dat$gain,
                     alternative = "greater", method.list = method.list2,
                     nperm = 10^5,
                     alpha = 0.10 )
ci.s5$n = 1:nrow(ci.s5)
ci.s5 = mutate( ci.s5, per = (n+0.5) / (nrow(ci)+1) )
ci2.s5 = ci.s5 %>%
    filter( lower > 0 ) %>%
    mutate( lower = round( lower, digits=1 ) ) %>%
    group_by( lower ) %>%
    arrange( n ) %>%
    slice_head( n=1 ) %>%
    relocate( n, per, lower )
ci2.s5
ci2





#### Analysis with blocking structure ####

# Here we randomize within site using a different implementation of
# the test for maximal effects that takes into account the blocking
# structure of the RCT.

dat = mutate( dat,
              TxAny.f = factor( TxAny, levels=c(1,0), labels=c("Tx","Co") ),
              Site = as.factor( Site ) )
table( dat$TxAny.f )
table( dat$TxAny.f, dat$Site )
levels( dat$TxAny.f )


# Testing null
steph_nonsup_test <-
    stephenson_test( dat, formula = gain ~ TxAny.f | Site,
                     subset_size = 6,
                     alternative = "greater",
                     tail = "right",
                     distribution = coin::approximate(1e6))
steph_nonsup_test

# NOTE: The factors of TxAny.f have to have Tx as the first factor.  If we
# reverse, we then need to flip the tail to get our low p-value.

# CI for extreme effects
steph_pos_ci <-
    find_ci_stephenson(dat, formula = gain ~ TxAny.f | Site,
                       alternative = "greater",
                       tail = "right",
                       tr_level = "Tx",
                       verbosity = 0,
                       distribution = coin::approximate(1e4),
                       subset_size = 6,
                       ci_level = 0.9)
steph_pos_ci
steph_pos_ci$ci



# Four versions of the CI depending on tail and sign flip on data.
CIs = generate_four_stephenson_cis( dat,
                                    gain ~ TxAny.f | Site,
                                    subset_size = 6,
                                    tr_level = "Tx",
                                    distribution=coin::approximate(100000) )
CIs
