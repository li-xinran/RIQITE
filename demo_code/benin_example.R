

## LIBRARIES
library(tidyverse)
library( RIQITE )

## DATA
benin <- RIQITE::benin %>%
  select(Village, District, Candidate, Treatment, VoteShare) %>%
  mutate(District = factor(District),
         Treatment = factor(Treatment, c("Policy", "Control"))) %>%
  filter(!is.na(Treatment)) %>%
  mutate( Z = as.numeric( Treatment == "Policy" ) )
benin




#### choose the test statistic
method.list = list( name = "Stephenson", s = 10 )

n = nrow( benin )

#### test the null hypothesis that the 95\% quantile of individual effect is less than or equal to 0
pval = pval_quantile( Z = benin$Z, Y = benin$VoteShare,
                      k = ceiling(n*0.95),
                      c = 0, alternative = "greater",
                      method.list = method.list, nperm = 10^4 )
pval


#### construct simultaneous confidence intervals for all quantiles of individual effects
ci = ci_quantile( Z = benin$Z, Y = benin$VoteShare,
                  alternative = "two.sided", method.list = method.list, nperm = 10^5,
                  alpha = 0.10 )
ci




##### Old weak null version #####

## Non-superiority null (all Policy <= Control)
steph_nonsup_test <- benin %>%
  stephenson_test(formula = VoteShare ~ Treatment | District,
                           subset_size = 6,
                           alternative = "greater",
                           tail = "right",
                           distribution = coin::approximate(1e6))
coin::pvalue(steph_nonsup_test)         # p = 0.086


### 90% ONE-SIDED CONFIDENCE INTERVALS
steph_pos_ci <- benin %>%
  find_ci_stephenson(formula = VoteShare ~ Treatment | District,
                              alternative = "greater",
                              tail = "right",
                              verbosity = 0,
                              distribution = coin::approximate(1e6),
                              subset_size = 6,
                              ci_level = 0.9)
steph_pos_ci$ci                         # at least one effect > 0.009

