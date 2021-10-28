

## LIBRARIES
library(tidyverse)
library( RIQITE )

## DATA
dat = make_fake_data("A")
dat = bind_rows( dat, dat )


#### choose the test statistic
method.list = list( name = "Stephenson", s = 10 )

n = nrow( dat )
n

#### test the null hypothesis that the 95\% quantile of individual effect is less than or equal to 0
pval = pval_quantile( Z = dat$Z, Y = dat$Yobs,
                      k = ceiling(n*0.95),
                      c = 0, alternative = "greater",
                      method.list = method.list, nperm = 10^4 )
pval


#### construct simultaneous confidence intervals for all quantiles of individual effects
ci = ci_quantile( Z = dat$Z, Y = dat$Yobs,
                  alternative = "two.sided", method.list = method.list, nperm = 10^5,
                  alpha = 0.10 )
ci
ci$tau = sort( dat$tau )
ci <- relocate( ci, lower, tau, upper )
ci


## Non-superiority null (all Policy <= Control)
steph_nonsup_test <-
  stephenson_test( dat, formula = Yobs ~ Z.f,
                  subset_size = 6,
                  alternative = "greater",
                  tail = "right",
                  distribution = coin::approximate(1e6))
steph_nonsup_test



### 90% ONE-SIDED CONFIDENCE INTERVALS
steph_pos_ci <-
  find_ci_stephenson(dat,
                     formula = Yobs ~ Z,
                     alternative = "greater",
                     tail = "right",
                     verbosity = 0,
                     distribution = coin::approximate(1e6),
                     subset_size = 6,
                     ci_level = 0.9)
steph_pos_ci
steph_pos_ci$ci


# Four versions
CIs = generate_four_stephenson_cis( dat, Yobs ~ Z, subset_size = 6,
                                    distribution=coin::approximate(100000) )
CIs
