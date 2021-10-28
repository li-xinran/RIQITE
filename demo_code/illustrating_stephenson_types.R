##
## This script is to help understand the four possible stephenson rank tests
## (based on the four directions) we might use.
##
## The given scenarios give illustrations of the separated CI case where we have
## evidence that tx both hurts and helps some people.
##

library( tidyverse )
library( RIQITE )
library( coin )
library( xtable )



#####  Make the Data #####

dat = make_fake_data("D")

# dat$Yobs = dat$Yobs * -1



#### Examine the scenario and associated ranks #####

# Our scenario
library( ggthemes )

my_t = theme_tufte() + theme( legend.position="bottom",
                              legend.direction="horizontal", legend.key.width=unit(1,"cm"),
                              panel.border = element_blank(),
                              axis.ticks.y = element_blank(),
                              axis.text.y = element_blank())
theme_set( my_t )

plt <- ggplot( dat, aes( x=Yobs ) ) +
    facet_wrap( ~ Z, ncol=1 ) +
    geom_dotplot( binwidth = 1) +
    labs( y="" )

print( plt )
ggsave( paste0( "scenario", SCENARIO, ".pdf" ), plt, height = 2.5, width = 6 )

dat %>% group_by( Z ) %>% summarise( Ybar = mean( Yobs ),
                                     sd = sd( Yobs ) )


levels( dat$Z )

# Look at the rankings depending on ordering of Y
dat = mutate( dat,
             rank = rank(Yobs),
             rnk = stephenson_rank( dat$Yobs ),
             ng_rnk = stephenson_rank( -dat$Yobs ) )

dat %>%
    select(Tr = Z,
           Y = Yobs,
           `ordinary rank` = rank,
           `Steph rank (pos)` = rnk,
           `Steph rank (neg)` = ng_rnk) %>%
    xtable::xtable(digits = 1) %>%
    print(floating = FALSE, include.rownames = FALSE)


dat

#### Baseline, classical tests of the null hypothesis ####

# Can we detect impact with t test?

t.test( Yobs ~ Z, data=dat, alternative="greater" )
t.test( Yobs ~ Z, data=dat, alternative="less" )


# Or wilcox?
wilcox_test(formula = Yobs ~ Z,
            data = dat,
            alternative = "greater",
            distribution = "exact")

wilcox_test(formula = Yobs ~ Z,
            data = dat,
            alternative = "less",
            distribution = "exact")



#### Generate the four CIs and four tests using Stephenson Rank ####

# The four different confidence intervals and tests against 0 treatment (strong null)
CIs = generate_four_stephenson_cis( dat, Yobs ~ Z, subset_size = 6, distribution=approximate(100000) )
CIs

cat( "Four tests give 4 CIs.\nRoughly speaking, we believe there are units that have treatment impacts in the 4 CIs\n")
print( as.data.frame( CIs[-6] ) )

print( xtable(  as.data.frame( CIs[-6] ), digits=c(0,0,2,2,2,2) ), include.rownames = FALSE )




if ( FALSE ) {
    ## Look at p-value curves
    head( CIs )
    CIl = CIs %>%
        rename( pval = pvalue ) %>%
                unnest( cols=pvaluelist,
                          names_repair="unique" )
    head( CIl )
    ggplot( CIl, aes( tau0, pvalue ) ) +
        facet_grid( alternative ~ tail ) +
        geom_line() +
        geom_hline( yintercept = 0.10, col="red", lty=2) +
        geom_hline( yintercept = 0.90, col="red", lty=2 )


    ggplot( filter( CIl, alternative=="greater", tail=="right"), aes( tau0, pvalue ) ) +
        geom_line() + geom_point() +
        geom_hline( yintercept = 0.10, col="red", lty=2) +
        geom_hline( yintercept = 0.90, col="red", lty=2 )

    CIs
    CIs$pvaluelist[[1]]$tau0
    1 + min( abs( diff( CIs$pvaluelist[[3]]$tau0 ) ) )
}


if ( FALSE ) {
    # look at shifted distributions
    dat = adjust_outcomes(dat, "Yobs", "Z", 9 )
    dat
    arrange( dat, Yadj )
    plt <- ggplot( dat, aes( x=Yadj ) ) +
        facet_wrap( ~ Z, ncol=1 ) +
        geom_dotplot( binwidth = 1)
    print( plt )

}

if  (FALSE ) {

    # Compare the pvalue curves to each other (for diagnostic purposes).
    CIl$mpvalue = 1 - CIl$pvalue1
    CIl$pp = ifelse( CIl$alternative == "less", CIl$mpvalue, CIl$pvalue1 )
    ggplot( CIl, aes( tau0, pvalue1, col=alternative, lty=tail ) ) +
        geom_line() + geom_point() +
        geom_hline( yintercept = 0.10, col="red" )
    ggplot( CIl, aes( tau0, pp, col=alternative, lty=tail ) ) +
        facet_wrap( ~ tail ) +
        geom_line() + geom_point() +
        geom_hline( yintercept = 0.10, col="red" ) +
        geom_hline( yintercept = 0.90, col="red" )
}
