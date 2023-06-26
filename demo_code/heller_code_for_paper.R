## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache.lazy = FALSE)
library(dplyr)
library(ggplot2)
library( RIQITE )
library( tidyverse )
set.seed(1)
nperm = 1000



## ----get the professional development data, message=FALSE----------------------------------------------------------------------------
data( electric_teachers )

mod <- lm(gain~TxAny, data = electric_teachers)
est_ATE <- round(mod$coefficients[2], digits = 2)
est_ATE


## ----electric_teacher_histograms, fig.height = 4, fig.width = 3.2, fig.align = "center"----------------------------------------------
hist(electric_teachers$gain[electric_teachers$TxAny==1],
     xlim = range(electric_teachers$gain), freq = FALSE, ylim = c(0, 0.065),
     col = "grey", border=FALSE, breaks = 20,
     main = NULL, xlab = "gain score")
hist(electric_teachers$gain[electric_teachers$TxAny==0],
     freq = FALSE, add = TRUE, breaks = 20, col = NULL)


## ----standardize_outcomes------------------------------------------------------------------------------------------------------------
sd_0 = sd( electric_teachers$Per.1[ electric_teachers$TxAny == 0 ] )
electric_teachers$std_gain = electric_teachers$gain / sd_0
electric_teachers %>% group_by( TxAny ) %>%
    summarise( raw_mean = mean( gain ),
               mean = mean( std_gain ),
               var = var( std_gain ),
               max = max( std_gain ),
               min = min( std_gain ),
               raw_sd = sd( gain ),
               sd = sd( std_gain ) )


## ------------------------------------------------------------------------------------------------------------------------------------
tx_het <- sqrt( 0.852^2 - 0.806^2 )
tx_het


## ----run_single_s_explore, cache=TRUE------------------------------------------------------------------------------------------------
set.seed(303020)
cdat = filter( electric_teachers, TxAny == 0 )
ate = est_ATE / sd_0
R = 1000
nperm = 10^3
s_list = c( 2, 3, 4, 5, 6, 8, 10, 15, 20 )
res <- explore_stephenson_s( s = s_list,
                      n = nrow( electric_teachers ),
                      Y0_distribution = cdat$std_gain,
                      tx_function = "rnorm", ATE = ate,
                      tx_scale = tx_het, rho = 0,
                      nperm = nperm,
                      parallel = TRUE, n_workers = 6,
                      R = R, calc_ICC = FALSE,
                      targeted_power = FALSE, k.vec = (233-99):233 )


## ------------------------------------------------------------------------------------------------------------------------------------
filter( res, k == 170 )


## ------------------------------------------------------------------------------------------------------------------------------------
checks = expand_grid( TxVar = c( 0.25, 0.75 ),
                      tx_dist = c( "rexp", "rnorm" ),
                      rho = c( -0.5, 0, 0.5 ) )

# There is no treatment variation, nor correlation with Y0, when
# the tx impact is a constant.
checks = bind_rows( checks,
                    tibble( TxVar = 0, tx_dist = "constant", rho = 0 ) )


## ----secret_simulation_block, include=FALSE------------------------------------------------------------------------------------------
if ( !file.exists( here::here( "demo_code/heller_s_check_results.rds" ) ) ) {
    checks$data = pmap( checks, function( TxVar, tx_dist, rho ) {
        cat( glue::glue("Running {TxVar} {tx_dist} rho={rho}\n" ) )
        explore_stephenson_s( s = s_list,
                              n = nrow( electric_teachers ),
                              Y0_distribution = cdat$std_gain,
                              tx_function = tx_dist,
                              nperm = nperm,
                              ATE = ate, tx_scale = TxVar, rho = rho,
                              R = R, calc_ICC = FALSE,
                              parallel = TRUE, n_workers = 6,
                              verbose = TRUE,
                              targeted_power = FALSE, k.vec = (233-99):233 )
    } )
    s_selector = unnest( checks, cols = "data" )

    saveRDS( s_selector, here::here( "demo_code/heller_s_check_results.rds" ) )
} else {
    s_selector = readRDS( here::here( "demo_code/heller_s_check_results.rds" ) )
}


## ----demo_simulation_code_block, eval=FALSE------------------------------------------------------------------------------------------
## checks$data = pmap( checks, function( TxVar, tx_dist, rho ) {
##     cat( glue::glue("Running {TxVar} {tx_dist} rho={rho}\n" ) )
##     explore_stephenson_s( s = s_list,
##                           n = nrow( dat ),
##                           Y0_distribution = cdat$std_gain,
##                           tx_function = tx_dist,
##                           ATE = ate, tx_scale = TxVar, rho = rho,
##                           nperm = nperm, R = R,
##                           calc_ICC = FALSE,
##                           parallel = TRUE, n_workers = 5,
##                           targeted_power = FALSE, k.vec = (233-99):233 )
## } )
## checks
## s_selector = unnest( checks, cols = "data" )


## ----plot_s_selector, warning=FALSE--------------------------------------------------------------------------------------------------
n_units = nrow(electric_teachers)
s_selector = mutate( s_selector,
                     group = cut( k,
                                  breaks = c(0, 233-60, 233-30, 233-10, 233),
                                  labels = c( "v low", "low", "med", "high" ) ) )

avg = s_selector %>% group_by( s, group, rho ) %>%
    summarise( q_ci = median( q_ci ),
               n = n(), .groups="drop" ) %>%
    filter( group != "v low" )
avg

s_selector$rho.f = factor( s_selector$rho,
                          levels = c( -0.5, 0, 0.5 ),
                          labels = c( "rho = -0.5", "rho = 0", "rho = 0.5" ) )

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
    coord_cartesian( ylim = c( -1.5, 2 ) ) +
    theme( legend.position="bottom",
                                legend.direction="horizontal", legend.key.width=unit(1,"cm"),
           panel.border = element_rect(colour = "grey", fill=NA, size=1) )



## ----some_checking_code, include=FALSE, eval=FALSE-----------------------------------------------------------------------------------
## if ( FALSE ) {
##     s_selector %>% dplyr::filter( s==2, tx_dist=="constant", TxVar == 0.5 ) %>%
##         pull( group ) %>% table()
## }
##
## # Save plot for main paper.
## ggsave( filename = "demo_code/s_selector.pdf", width = 8, height = 4 )


## ------------------------------------------------------------------------------------------------------------------------------------
qis <- s_selector %>%
    group_by( TxVar, tx_dist, rho, group, k ) %>%
    filter( q_ci == max(q_ci) )
table( group = qis$group, s = qis$s )


## ------------------------------------------------------------------------------------------------------------------------------------
n_score <- s_selector %>%
    group_by( rho, tx_dist, s ) %>%
    summarize( n = mean( n ), .groups = "drop")
ggplot( n_score, aes( s, n, col=tx_dist, pch=as.factor(rho) ) ) +
    facet_wrap( ~ rho, nrow=1 ) +
    geom_point( size=3) +
    labs( pch="rho" ) +
    theme_minimal()


## ---- cache=TRUE---------------------------------------------------------------------------------------------------------------------
nperm = 10^6
pval.steph8 <- pval_quantile(Z=electric_teachers$TxAny, Y=electric_teachers$gain,
                              k=n_units, c=0, alternative = "greater",
                              method.list=list(name="Stephenson", s=8), nperm=nperm )
pval.steph8 # final p-value


## ---- fig.height = 4, fig.width = 6.4, fig.align = "center"--------------------------------------------------------------------------
ci6 = ci_quantile( Z=electric_teachers$TxAny, Y=electric_teachers$gain, alternative="greater",
                  method.list=list( name="Stephenson", s=8 ), nperm=nperm, alpha=0.10 )
ci2 = ci_quantile( Z=electric_teachers$TxAny, Y=electric_teachers$gain, alternative="greater",
                  method.list=list( name="Stephenson", s=2 ), nperm=nperm, alpha=0.10 )
par(mfrow = c(1, 2))
plot_quantile_CIs(ci6, k_start = 102, main = "s=8" )
plot_quantile_CIs(ci2, k_start = 102, main = "s=2" )


## ---- cache=TRUE---------------------------------------------------------------------------------------------------------------------
tau.max.lower = ci_quantile( Z=electric_teachers$TxAny, Y=electric_teachers$gain, k.vec=n_units,
                             alternative="greater", method.list=list( name="Stephenson", s=8 ),
                             nperm=nperm, alpha=0.05 )$lower
tau.min.upper = ci_quantile( Z=electric_teachers$TxAny, Y=electric_teachers$gain, k.vec=1,
                             alternative="less", method.list=list( name="Stephenson", s=8 ),
                             nperm=nperm, alpha=0.05 )$upper
max(0, tau.max.lower - tau.min.upper)




