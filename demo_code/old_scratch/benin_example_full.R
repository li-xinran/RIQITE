##
## This file illustrates the perminf package using the Benin data.
##
## TO DO: It should be converted to an Rmd file and then a vignette to be put in
## the package.
##

## Libraries
library(coin)
library(tidyverse)
library(perminf)

#data( benin, package = "perminf" )

## Data
benin <- perminf::benin %>%
    select(Village, District, Candidate, Treatment, VoteShare)


head( benin )

# Make district a factor for the permutation inference
benin$District = factor( benin$District )

benin_sub = benin %>% filter(Treatment != "Clientelist")
benin_sub


## Illustrate Stephenson ranks
data.frame( y = benin_sub$VoteShare,
            rank = rank(benin_sub$VoteShare),
            tx = benin_sub$Treatment,
            s_rank = stephenson_rank(benin_sub$VoteShare, subset_size = 6)) %>%
  arrange(y) %>%
  print(digits = 2)


## Use the coin package's independence_test with Stephenson ranks as scores for
## Y (with default subset size of 6) (Interpretation: Is there evidence of any
## positive treatment effects?)
independence_test(formula = VoteShare ~ Treatment | District,
                  data = benin_sub,
                  alternative = "less", # control "less" than treatment (Policy)
                  distribution = "asymptotic",
                  ytrafo = stephenson_trafo )


## Same as above, but using our standalone function (specifying subset size = 6).
stephenson_test(formula = VoteShare ~ Treatment | District,
                data = benin_sub,
                alternative = "less", subset_size = 6,
                tail = "right", distribution = "asymptotic")



## Now specifying subset size = 2
stephenson_test(formula = VoteShare ~ Treatment | District,
                data = benin_sub,
                alternative = "less", subset_size = 2,
                tail = "right", distribution = "asymptotic")




## Note similarity of subset of size 2 to Wilcoxon rank sum
wilcox_test(formula = VoteShare ~ Treatment | District,
            data = benin_sub,
            alternative = "less",
                  distribution = "asymptotic")


## Contrast with a permutation test using the difference of means as a test
## statistic (no score transformation)
independence_test(formula = VoteShare ~ Treatment | District,
                  data = benin_sub,
                  alternative = "less",
                  distribution = "asymptotic")


## Now we test the left tail
## (Interpretation: Is there evidence of any negative treatment effects?)
stephenson_test(formula = VoteShare ~ Treatment | District,
                data = benin_sub,
                subset_size = 6, distribution = "asymptotic",
                alternative = "greater", tail = "left")

## Same as this:
stephenson_test(formula = I(-VoteShare) ~ Treatment | District,
                data = benin_sub,
                subset_size = 6,
                distribution = approximate(100000),
                alternative = "less", tail = "right")



#####  Confidence Intervals ######

# Get a confidence interval for maximal effects using Stephenson Rank
benin_ci_precise.steph <- benin_sub %>%
    find_ci_stephenson(## formula = VoteShare ~ Treatment |District,
                       tr_var = "Treatment", 
                       tr_level = "Policy",
                       y_var = "VoteShare",
                       strat_var = "District",
                       alternative = "less",
                       verbosity = 0,
                       distribution = approximate(100000),
                       subset_size = 6 )
benin_ci_precise.steph$ci




# The following more general method first gets a starting range of values by
# looking for low and high p-values, and then zeros in on the threshold p-value
# for the CI with a binary search.
benin_ci_precise.onestep <- find_ci( benin_sub,
                             coin_test = "stephenson_test",
                             tr_var = "Treatment",
                             tr_level = "Policy",
                             y_var = "VoteShare",
                             strat_var = "District",
                             alternative = "less",
                             verbosity = 2,
                             ci_level = .9,
                             step_fraction = .0001,
                             distribution = approximate(100000),
                             round_digits = 5,
                             subset_size = 6
)

benin_ci_precise.onestep$ci





##
## Looking at impact of blocking by district.
##

## Our "significant" test:
stephenson_test(formula = VoteShare ~ Treatment | District,
                data = benin_sub,
                distribution = "asymptotic", subset_size = 6,
                alternative = "greater", tail = "left")
# It suggests that there were positively impacted units that were at the low end
# of the outcome range.

## Original experiment was blocked, but if we ignore the blocking we lose power:
stephenson_test(formula = VoteShare ~ Treatment,
                data = benin_sub,
                distribution = "asymptotic", subset_size = 6,
                alternative = "greater", tail = "left")


# Looking at pairwise differences of our blocks
b_w = benin_sub %>% select( -Village, -Candidate ) %>%
    spread( Treatment, VoteShare ) %>%
    mutate( delta.per = Policy - Control ) %>%
    arrange( delta.per )
b_w


b_w2 = benin_sub %>% mutate( steph_rank = stephenson_rank( -1 * VoteShare, subset_size=6) ) %>%
    select( -VoteShare, -Candidate, -Village ) %>%
    spread( Treatment, steph_rank ) %>%
    mutate( delta = Control - Policy )
b_w2

merge( b_w, b_w2, by="District", all=TRUE ) %>%
    arrange( delta )



b2 <- benin_sub %>% group_by( District ) %>%
    mutate( Vcent = VoteShare - mean(VoteShare),
            Vavg = mean(VoteShare) ) %>%
    ungroup()
b2
head(b2)
ggplot( b2, aes( x=Treatment, y=Vcent ) ) +
    geom_point()



###### Looking at the precise vs. approximate calculation of p-values


### Grid search to find confidence bound
## This gives a quick asymptotic approximation
benin_ci_init <- find_ci( benin_sub,
                          coin_test = "stephenson_test",
                          tr_var = "Treatment",
                          tr_level = "Policy",
                          y_var = "VoteShare",
                          strat_var = "District",
                          alternative = "less",
                          verbosity = 1,
                          ci_level = .9,
                          step_fraction = .01,
                          distribution = "asymptotic",
                          subset_size = 6 )
benin_ci_init$ci


## Now give starting value based on above, reduce step, and use resampling distribution for
## more precise answers

# the value
benin_ci_init$ci_bound + .03

benin_ci_precise <- find_ci( benin_sub,
                             coin_test = "stephenson_test",
                             tr_var = "Treatment",
                             tr_level = "Policy",
                             y_var = "VoteShare",
                             strat_var = "District",
                             alternative = "less",
                             verbosity = 0,
                             ci_level = .9,
                             step_fraction = .0001,
                             distribution = approximate(100000),
                             #init_tau0 = benin_ci_init$ci_bound + .03,
                             round_digits = 5,
                             subset_size = 6
)
benin_ci_precise

# The sequence of examined tests
benin_ci_precise$pvalues

# The final bound
benin_ci_precise$ci_bound

# Looking at the search pattern of what tau0 are tested with each step
pvs = benin_ci_precise$pvalues
pvs = pvs[ order( as.numeric( rownames( pvs ) ) ), ]
plot( pvs$tau0, type="b" )
abline( h = benin_ci_precise$ci_bound, lty=2 )


# Compare to approximate
pv = bind_rows( approximate = benin_ci_init$pvalues,
                precise = benin_ci_precise$pvalues,
                .id = "test")
head( pv )

ggplot( pv, aes( x=tau0, y=pvalue, col=test ) ) +
    geom_line() + geom_point() +
    #    coord_cartesian( xlim=range( -0.1, 0.1 ) ) +
    #benin_ci_precise$pvalues$tau0 - 0.05,  benin_ci_precise$pvalues$tau0 + 0.10  ) ) +
    geom_hline(yintercept=0.1) +
    geom_vline(xintercept=benin_ci_precise$ci_bound, col="grey")


# Compare to original method call above (they should basically match)
benin_ci_precise.onestep$ci_bound - benin_ci_precise$ci_bound




