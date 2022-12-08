
#
# Run simulation in the paper using the RIQITE package.
#

jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(jobid)
setwd( "~/scratch/Simulation_Quantile_Normal_random_20221017" )

# jobid = 28
# setwd("~/Documents/Quantile_Simulation_2022Sep/Simulation_Quantile_Normal_random_20221017")

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


# 1 + omega^2 + 2 * rho * omega

normal_gen_PO <- function(n){
  # Y0tau = rmvnorm( n, mean = c(0, tau0),
  #         sigma = matrix(c(1, rho*omega, rho*omega, omega^2), nrow = 2) )
  # Y0 = Y0tau[, 1]
  # Y1 = Y0tau[, 1] + Y0tau[, 2]

  Y0 = rnorm(n)
  epsilon = rnorm(n)
  tau = omega * ( rho * Y0 + sqrt(1-rho^2) * epsilon ) + tau0
  Y1 = Y0 + tau

  return(data.frame(Y1=Y1, Y0=Y0))
}


result = simu_power( gen = normal_gen_PO, n = 120, treat.prop = 0.5, iter.max = iter.max,
                     Steph.s.vec = c(2, 4, 6, 10, 30, 60),
                     alternative = "greater", alpha = sign.level,
                     nperm = nperm, tol = 10^(-3), switch = TRUE )

perform = summary_power(result, measure = "CI_n(c)", alternative = "greater", cutoff = 0)

pval.dim = rep(NA, iter.max)

for(iter in 1:iter.max){
  dat = result$dat.all[[iter]]
  Z = result$Z.CRE[, iter]
  Y = Z * dat$Y1 + (1-Z) * dat$Y0

  pval.dim[iter] = pval_bound(Z = Z, Y = Y, c = 0, method.list = list(name = "DIM"),
                              alternative = "greater", Z.perm = result$Z.perm, impute = "control")

  print(iter)
}

power.dim = mean( pval.dim <= sign.level )

saveRDS( list(tau0 = tau0,
              omega = omega,
              rho = rho,
              result = result,
              perform = perform,
              pval.dim = pval.dim,
              power.dim = power.dim
              ),
         file = paste0("comp_n120_s_n0_", jobid, ".rds")
)




