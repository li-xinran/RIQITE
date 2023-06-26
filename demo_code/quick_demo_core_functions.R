

dat
p_tx = 0.5
n = 120
head(dat)
dat$tau = dat$Y1 - dat$Y0
dat <- mutate( dat,
               Z = as.numeric( sample(n) <= p_tx * n ),
               Yobs = Y0 + ifelse( Z, tau, 0 ) )

a = ci_quantile( dat$Z, dat$Yobs, nperm = 1000 )
class(a)

pval_bound( dat$Z, dat$Yobs, c = -0.5 )

pval_quantile( dat$Z, dat$Yobs, c = -0.5, k = 120 )
