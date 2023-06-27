test_that("CI_bound works", {


    # Possible negative constant treatment effect situation
    dd = make_fake_data(scenario = "B")


    pp = pval_bound( dd$Z, dd$Yobs, nperm=1000, alternative = "less" )
    expect_true( pp < 0.1 )

    pp = ci_bound( dd$Z, dd$Yobs, nperm=1000, alternative = "less" )
    expect_true( pp < 0 )

    # We have evidence that tx effects can't be smaller than -25!
    pp = ci_bound( dd$Z, dd$Yobs, nperm=1000, alternative = "greater" )
    expect_true( pp < 0 )


    dd = make_fake_data(scenario = "A")
    pp1 = ci_bound( dd$Z, dd$Yobs, nperm=1000, alternative = "greater",
                    method.list = list(name = "Stephenson", s = 10) )
    pp1
    pp2 = ci_bound( dd$Z, dd$Yobs, nperm=1000, alternative = "less",
                    method.list = list(name = "Stephenson", s = 10) )
    pp2

    # We have tx heterogeniety
    expect_true( pp2 < pp1 )


})
