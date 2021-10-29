

## LIBRARIES
library(tidyverse)
library( RIQITE )

## DATA
dat = read_csv( "demo_code/cleanTeacher.csv",na = "." )
head( dat )
dat = dplyr::select( dat, T.ID, Site, Tx, Per.1, Per.2, gain )
head( dat )
table( dat$Tx, dat$Site, is.na( dat$gain ) )
filter( dat, is.na( gain ), Site =="S1" )
dat = mutate( dat,
              TxAny = ifelse( Tx == "D", 0, 1 ) )

# Site 1 apparently has no control scores; not sure why.
dat = filter( dat, !is.na( gain ) & Site != "S1" )
nrow( dat )
dat$gain

set.seed( 1010101 )
dat = sample_n( dat, nrow(dat) )

# Three testing occasions: Baseline (before Tx), Post treatment, and then long
# term follow up.  (1, 2, and 3 in the data).

# The gain (post score - pre score) for all teachers in percentage point correct
# on the test of electric circuit knowledge (we call this "content knowledge").

# For more on the data see
#
# Heller, J. I., Shinohara, M., Miratrix, L., Hesketh, S. R., & Daehler, K. R.
# (2010). Learning Science for Teaching: Effects of Professional Development on
# Elementary Teachers, Classrooms, and Students. Society for Research on
# Educational Effectiveness.
#
# https://files.eric.ed.gov/fulltext/ED514193.pdf

ggplot( dat, aes( gain ) ) +
  facet_wrap( ~ Tx ) +
  geom_histogram() +
  geom_vline( xintercept = 0, col="red" )


# Classic OLS: What can we say about the average effect?
mod = lm( gain ~ 0 + TxAny + as.factor(Site), data=dat )
summary( mod )
confint( mod )
confint( mod, level = 0.90 )



#### choose the test statistic
method.list = list( name = "Stephenson", s = 10 )

n = nrow( dat )
n
n * 0.95

# ## test the null hypothesis that the 70\% quantile of individual effect is
# less than or equal to 0
#
# Note: If k = n, then you are testing the max effect (sharp null) --- like the
# original paper of Alan, Devin and Luke.
pval = pval_quantile( Z = dat$TxAny, Y = dat$gain,
                      k = ceiling(n*0.70),
                      c = 0, alternative = "greater",
                      method.list = method.list, nperm = 10^5 )
pval


#### construct simultaneous confidence intervals for all quantiles of individual effects
ci = ci_quantile( Z = dat$TxAny, Y = dat$gain,
                  alternative = "greater", method.list = method.list, nperm = 10^5,
                  alpha = 0.10 )
ci$n = 1:nrow(ci)
ci = mutate( ci, per = (n+0.5) / (nrow(ci)+1) )
filter( ci, sign(lower) == sign(upper) )


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
                  alternative = "greater", method.list = method.list2, nperm = 10^5,
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





##### Old code analysis ######

# Here we randomize within site using the old code.
#
# This code takes into account the blocking structure of the RCT.

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
