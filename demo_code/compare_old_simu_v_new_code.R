


#
# Comparing xinran's code to luke's code to make it line up
#
# Modified version of Xinran's simulations script used on the cluster.
#


jobid = 96

library( tidyverse )
library(RIQITE)
library(utils)

iter.max = 500

nperm = 10^4
sign.level = 0.1


tau0.grid = c(-1, 0, 1)
rho.grid = seq(from = -1, to = 1, by = 0.1)
omega.grid = c(0.5, 1)
para.table = expand.grid(tau0.grid, omega.grid, rho.grid)
colnames(para.table) = c("tau0", "omega", "rho")

nrow( para.table )

# Get the parameters for this simulation run
tau0 = para.table$tau0[jobid]
omega = para.table$omega[jobid]
rho = para.table$rho[jobid]

para.table[jobid, ]



#### Step 1: Getting mean CIs on number of significant units ####


## Old version:

# 1 + omega^2 + 2 * rho * omega
normal_gen_PO <- function(n){
    # Y0tau = rmvnorm( n, mean = c(0, tau0),
    #         sigma = matrix(c(1, rho*omega, rho*omega, omega^2), nrow = 2) )
    # Y0 = Y0tau[, 1]
    # Y1 = Y0tau[, 1] + Y0tau[, 2]

    Y0 = rnorm(n)
    epsilon = rnorm(n)
    tau = tau0 + omega * (rho * Y0 +  sqrt(1-rho^2) * epsilon)
    Y1 = Y0 + tau

    # tau = tau0 + omega * rho Y0 + omega * sqrt(1-rho^2) epsilon

    return(data.frame(Y1=Y1, Y0=Y0))
}


result = simu_power( gen = normal_gen_PO, n = 120, treat.prop = 0.5, iter.max = iter.max,
                     Steph.s.vec = c(2, 4, 6, 10, 30, 60),
                     alternative = "greater", alpha = sign.level,
                     nperm = nperm, tol = 10^(-2), switch = TRUE,
                     keep_data = TRUE )

perform = summary_power(result, measure = "CI_n(c)", alternative = "greater", cutoff = 0)
perform

## New version:
alt <- explore_stephenson_s( s = c(2, 4, 6, 10, 30, 60),
                             n = 120,
                             Y0_distribution = rnorm,
                             tx_function = "linear",
                             ATE = tau0,
                             tx_scale = omega * sqrt(1 - rho^2),
                             rho = omega * rho,
                             p_tx = 0.5,
                             R = iter.max, iter_per_set = 10,
                             alpha = sign.level,
                             c = 0,
                             quantile_n_CI = NA,
                             targeted_power = FALSE,
                             alternative = "greater",
                             nperm = 1000,
                             k.vec = NULL,
                             parallel = TRUE,
                             n_workers = 5,
                             calc_ICC = FALSE,
                             verbose = TRUE )


# Stats on each quantile.  The stats on n is the last column, same
# across all quantiles for a given s.

# Subset to the unique CIs on n (rather than the individual stats on individual units)
aa <- alt %>% dplyr::select( s, n ) %>%
    unique()
aa

dd = tibble( s = perform$s,
             old = perform$ci.mean ) %>%
    left_join( aa ) %>%
    mutate( per = n / old )


# These are in the ballpark of each other, I think
dd


dd %>% pivot_longer( cols = c( "old", "n" ) )  %>%
    ggplot( aes( s, value, col=name ) ) +
    geom_point() +
    geom_line()




# Step 2: Looking at DiM for baseline, I think?


# Old version:

## ALERT!!!!! I changed c = 0.5 to have non-unity power.  Remove this line for the real simulation!!!!!
c = 0.5
warning( "c changed to 0.5; update script to change it back" )

pval.dim = rep(NA, iter.max)

for(iter in 1:iter.max){
    dat = result$dat.all[[iter]]
    Z = result$Z.CRE[, iter]
    Y = Z * dat$Y1 + (1-Z) * dat$Y0

    pval.dim[iter] = pval_bound( Z = Z, Y = Y, c = c,
                                 method.list = list(name = "DIM"),
                                 alternative = "greater",
                                 Z.perm = result$Z.perm,
                                 impute = "control" )

    if ( iter %% 10 == 0 ) {
        cat( glue::glue( "{iter}/{iter.max}" ), "\n" )
    }
}
power.dim = mean( pval.dim <= sign.level )
power.dim

# New version:
alt_DIM <- calc_power( n = 120,
                       Y0_distribution = rnorm,
                       tx_function = "linear",
                       ATE = tau0,
                       tx_scale = omega * sqrt(1 - rho^2),
                       rho = rho * omega,
                       p_tx = 0.5,
                       R = iter.max,
                       percentile = 1,
                       alpha = sign.level,
                       targeted_power = TRUE,
                       use_pval_bound = TRUE,
                       c = c,
                       alternative = "greater",
                       method.list = list(name = "DIM") )


# These are actually fairly close:
alt_DIM
power.dim





#### Bundle and save simulation results ####




results <- list(tau0 = tau0,
                omega = omega,
                rho = rho,
                result = result,
                perform = perform,
                pval.dim = pval.dim,
                power.dim = power.dim )
results

saveRDS( results, file = paste0("comp_n120_s_n0_", jobid, ".rds") )




