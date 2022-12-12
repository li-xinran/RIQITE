
#
# Run simulation in the paper using the RIQITE package.
#

# This tells us what simulation scenario to run (indexed by jobid)
jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
setwd( "~/scratch/Simulation_Quantile_Normal_random_20221017" )

set.seed(jobid * 13)

# jobid = 28
# setwd("~/Documents/Quantile_Simulation_2022Sep/Simulation_Quantile_Normal_random_20221017")

library(tidyverse)
library(RIQITE)
library(utils)
# library(mvtnorm)

iter.max = 500
nperm = 10^5
sign.level = 0.1


tau0.grid = c(-1, 0, 1)
rho.grid = seq(from = -1, to = 1, by = 0.1)
omega.grid = c(0.5, 1)
para.table = expand.grid(tau0.grid, omega.grid, rho.grid)
colnames(para.table) = c("tau0", "omega", "rho")

tau0 = para.table$tau0[jobid]
omega = para.table$omega[jobid]
rho = para.table$rho[jobid]




#### Step 1: Getting mean CIs on number of significant units ####

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
                             parallel = FALSE,
                             calc_ICC = FALSE,
                             verbose = TRUE )


# The above gives stats for each quantile.  The stats on n is the last
# column, the same will be listed across all quantiles for a given s
# as the CI for n is for the entire analysis (all the quantiles
# aggregated).


# Get the CI on n as a function of s by subsetting to the unique CIs
# on n (rather than the individual stats on individual units)
n_CIs <- alt %>% dplyr::select( s, n ) %>%
    unique()

n_CIs



# Step 2: Looking at DiM for baseline, I think?

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
                       c = 0,
                       alternative = "greater",
                       method.list = list(name = "DIM") )


alt_DIM


# Step 3: Save the results


saveRDS( list(tau0 = tau0,
              omega = omega,
              rho = rho,
              result = alt,
              perform = n_CIs,
              pval.dim = alt_DIM,
              power.dim = alt_DIM$power
              ),
         file = paste0("comp_n120_s_n0_", jobid, ".rds")
)




