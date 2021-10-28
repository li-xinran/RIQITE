
##
## Testing the code
##

library( testthat )
library( perminf )
library( tidyverse )

context("Do things run")

test_that("we can call functions without crashing", {

                                        # data( benin )
    benin <- perminf::benin
    ## Make district a factor for the permutation inference
    benin$District = factor( benin$District )

    benin_sub = benin %>% filter(Treatment != "Clientelist")

    ci <- find_ci(data = benin_sub,
                  coin_test = "stephenson_test",
                  formula = VoteShare ~ Treatment | District,
                  tr_level = "Policy",
                  alternative = "less",
                  step_fraction = .0001,
                  subset_size = 6
                  )
    names( ci )

    expect_equal( ci$ci[[1]], -Inf )

    ci <- find_ci_stephenson( benin_sub,
                             formula = VoteShare ~ Treatment | District,
                             alternative = "greater",
                             subset_size = 6
                             )
    ci$ci
    expect_equal( ci$ci[[2]], Inf )

    ci <- find_ci_stephenson( benin_sub,
                             formula = VoteShare ~ Treatment | District,
                             tr_level = "Policy",
                             tail="left",
                             alternative = "greater",
                             subset_size = 6
                             )
    ci$ci
    expect_equal( ci$ci[[2]], Inf )

    ci <- find_ci_stephenson( benin_sub,
                             formula = VoteShare ~ Treatment | District,
                             tr_level = "Policy",
                             tail="left",
                             alternative = "less",
                             subset_size = 6
                             )
    ci$ci
    expect_equal( ci$ci[[1]], -Inf )
})








