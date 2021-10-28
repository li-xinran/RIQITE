#
# Example and data taken from Rosenbaum 2007 confidence interval paper
#

library(coin)
library(tidyverse)
library(perminf)

#### Make the data ####

dat = data.frame(
    W.ID = c(8, 31, 5, 40, 10, 18, 24, 23, 13, 4, 45, 20, 15, 11, 7, 33, 2,
             16, 43, 26, 38, 12, 27, 29, 39),
    C.ID = c(104, 101, 110, 103, 102, 106, 108, 105, 111, 112, 113, 125,
             126, 123, 118, 119, 120, 124, 121, 122, 116, 114, 115, 109,
             117),
    W.smoke = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 10, 10, 15, 20, 20, 20, 20,
                20, 20, 20, 20, 30, 30, 30),
    C.smoke = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 10, 10, 15, 20, 20, 20, 20,
                20, 20, 20, 20, 30, 30, 30),
    W.age = c(26, 33, 34, 43, 45, 45, 46, 49, 51, 54, 28, 31, 49, 36, 32,
              33, 36, 37, 39, 43, 49, 50, 37, 38, 43),
    C.age = c(27, 29, 34, 44, 45, 46, 47, 50, 53, 54, 24, 23, 34, 38, 32,
              36, 36, 37, 40, 43, 51, 58, 42, 46, 50),
    W.adducts = c(1.83, 1.06, 0.82, 1.05, 0.87, 0.99, 1.15, 0.31, 1.91, 1.05,
                  2.05, 1.8, 0.74, 3.02, 1.75, 6.37, 0.8, 4.12, 2.08, 2.95,
                  2.32, 1.11, 3.9, 2.94, 0.36),
    C.adducts = c(0.65, 1.38, 1.48, 1.39, 1.26, 0.23, 0.4, 0.71, 1.32, 1.24,
                  1.12, 2.22, 1.31, 1.59, 0.79, 1.2, 0.62, 2.03, 1.3, 1.3, 1.55,
                  2.42, 1.42, 2.26, 0.93)
)



head( dat )
dat$pairID = 1:nrow(dat)


tx = dat[c(1,3,5,7,9)]
co = dat[c(2,4,6,8,9)]
names(co) = names(tx) = c("ID","smoke","age","adducts", "pairID")

datl = bind_rows( Tx = tx, Co=co, .id="Z" )
head( datl )

datl$pairID = as.factor( datl$pairID )
datl$Z = factor( datl$Z, levels= c("Tx","Co" ) )


ggplot( datl, aes( adducts ) ) +
    facet_wrap( ~ Z, ncol=1) +
    geom_dotplot( binwidth=0.1 )

ggplot( datl, aes( adducts ) ) +
    facet_wrap( ~ Z, ncol=1) +
    geom_histogram( aes( y = ..density.. ), binwidth = 0.25 )


datl %>% group_by( Z ) %>% summarise( mean = mean( adducts ),
                                      sd = sd( adducts ) )
nrow( datl )
head( dat )

rm( tx, co )



#### Make a dataset with all potential outcomes for simulation study ####
simdat = dat
simdat = transform( simdat, ID = paste( W.ID, C.ID, sep="-" ),
                    smoke = (W.smoke + C.smoke)/2,
                    age = (W.age + C.age)/2,
                    Y0 = C.adducts,
                    Y1 = W.adducts )

qplot( Y0, Y1, data=simdat ) +
    geom_abline( intercept = 0, slope=1 ) +
    geom_smooth( method="lm" )


# tweak to make the impacts monotonic (but preserving the margin)
# simdat$Y0 = sort(simdat$Y0)[ order( simdat$Y1 ) ]
# (Doesn't quite work as intended)
qplot( Y0, Y1, data=simdat ) +
    geom_abline( intercept = 0, slope=1 ) +
    geom_smooth( method="lm" )

simdat = mutate( simdat, tau = Y1 - Y0 )

qplot( simdat$tau )

# How far apart are the units when sorted
median( diff( sort( simdat$Y0 ) ) )
mean( diff( sort( simdat$Y1 ) ) )

qqnorm( simdat$Y0 )





#### Conduct tests on oringal data ####
if  (FALSE ) {

    head( datl )

    stephenson_test(formula = adducts ~ Z | pairID,
                    data = datl,
                    alternative = "less", subset_size = 6,
                    tail = "right", distribution = "asymptotic")

    stephenson_test(formula = adducts ~ Z,
                    data = datl,
                    alternative = "less", subset_size = 6,
                    tail = "right", distribution = "asymptotic")


    # Look at ranks and how they load by tx and co
    dd = data.frame( y = datl$adducts,
                     rank = rank(datl$adducts),
                     tx = datl$Z,
                     s_rank = stephenson_rank(datl$adducts, subset_size = 6)) %>%
        arrange(-y) %>%
        print(digits = 2)
    head( dd, 20 )

    qplot( tx, s_rank, data=dd, geom="boxplot" )
    qplot( tx, rank, data=dd, geom="boxplot" )

    levels( datl$Z )

    CIs = generate_four_stephenson_cis( datl, adducts ~ Z | pairID, tr_level = "Tx", subset_size = 6 )

    print( CIs )


    alum.lr <- find_ci_stephenson( datl,
                                   formula = adducts ~ Z | pairID,
                                   tr_level = "Tx",
                                   alternative = "less",
                                   tail = "right",
                                   verbosity = 0,
                                   subset_size = 6 )
    alum.lr
    alum.lr$ci



    # Ignoring blocking
    CIs.ig = generate_four_stephenson_cis( datl, adducts ~ Z, tr_level = "Tx", subset_size = 6 )
    CIs.ig

}



