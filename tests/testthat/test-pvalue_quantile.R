


test_that("pvalue_quantile works", {

    set.seed(10103)
    dd = make_fake_data(scenario = "C")
    nrow(dd)
    pv = pval_quantile( Z = dd$Z, Y = dd$Yobs, c = 3, nperm = 100, k = 20 )
    pv2 = pval_quantile( Z = dd$Z, Y = dd$Yobs, c = 3, nperm = 100, k = 24 )
    expect_true( pv2 <= pv )


    pv2 = pval_quantile( Z = dd$Z, Y = dd$Yobs, c = 0, nperm = 100, k = 1, alternative = "less" )
    expect_true( pv2 > 0.1 )



})



test_that("automatic label switching of ci_quantile works", {

    dd = make_fake_data(scenario = "C")
    dd = dd %>% filter( Z == 1 ) %>%
        bind_rows( dd )
    dd = dd %>% filter( Z == 1 ) %>%
        bind_rows( dd )
    table(dd$Z )


    set.seed(10444)
    n = nrow(dd)
    pv = ci_quantile( Z = dd$Z, Y = dd$Yobs, nperm = 200, k.vec = 1:n, switch = TRUE )
    pv2 = ci_quantile( Z = dd$Z, Y = dd$Yobs, nperm = 200, k.vec = 1:n, switch = FALSE )
    tt <- left_join( pv, pv2, by="k" ) %>%
        dplyr::select( -upper.x, -upper.y ) %>%
        mutate( tau = sort( dd$Y1 - dd$Y0 ) ) %>%
        dplyr::filter( k > 20 )
    tt

    expect_true( all( tt$lower.x <= tt$tau ) )
    expect_true( all( tt$lower.y <= tt$tau ) )

    expect_true( sum( is.finite(tt$lower.x) ) < sum( is.finite(tt$lower.y) ) )



})
