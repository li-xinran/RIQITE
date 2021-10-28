set.seed(1)

## LIBRARIES
library(coin)
library(tidyverse)
library(perminf)

## DATA
benin <- perminf::benin %>%
    select(Village, District, Candidate, Treatment, VoteShare) %>%
    mutate(District = factor(District),
           Treatment = factor(Treatment, c("Policy", "Control"))) %>%
    filter(!is.na(Treatment))

### TESTS
## Paired t
t.test(VoteShare ~ Treatment, paired = TRUE, data = benin) # ATE est. = -0.0575, p = 0.45 (two-sided)

## Non-superiority null (all Policy </= Control)
steph_nonsup_test <- benin %>%
    perminf::stephenson_test(formula = VoteShare ~ Treatment | District,
                             subset_size = 6,
                             alternative = "greater",
                             tail = "right",
                             distribution = approximate(1e6))
coin::pvalue(steph_nonsup_test)         # p = 0.086

## Non-inferiority null (all Policy >/= Control)
steph_noninf_test <- benin %>%
    perminf::stephenson_test(formula = VoteShare ~ Treatment | District,
                             subset_size = 6,
                             alternative = "less",
                             tail = "left",
                             distribution = approximate(1e6))
coin::pvalue(steph_noninf_test)         # p = 0.125

### 90% ONE-SIDED CONFIDENCE INTERVALS
## Most-positive effect
steph_pos_ci <- benin %>%
    perminf::find_ci_stephenson(formula = VoteShare ~ Treatment | District,
                                alternative = "greater",
                                tail = "right",
                                verbosity = 0,
                                distribution = approximate(1e6),
                                subset_size = 6,
                                ci_level = 0.9)
steph_pos_ci$ci                         # at least one effect > 0.009

## Most-negative effect
steph_neg_ci <- benin %>%
    perminf::find_ci_stephenson(formula = VoteShare ~ Treatment | District,
                                alternative = "less",
                                tail = "left",
                                verbosity = 0,
                                distribution = approximate(1e6),
                                subset_size = 6,
                                ci_level = 0.9)
steph_neg_ci$ci                         # at least one effect < 0.051
