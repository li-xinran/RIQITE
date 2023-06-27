

#' Summarize result of a ci_quantile call
#'
#' This will print out the confidence interval on number of significant units and so forth.
#'
#' @param ci Result object from the ci_quantile call.
#' @param c Threshold to test confidence intervals against
#'
#' @return Invisible object of some summary stats, which are also printed out.
#' @export
summarize_quantile_CIs <- function( ci, c ) {

    npos = sum( ci6$lower >= c )


}



test_that("plot and summary doesn't crash", {

    nperm = 100
    data( "electric_teachers")
    n_units = nrow(electric_teachers)
    pval.steph8 <- pval_quantile(Z=electric_teachers$TxAny, Y=electric_teachers$gain,
                                 k=n_units, c=0, alternative = "greater",
                                 method.list=list(name="Stephenson", s=8), nperm=nperm )
    expect_true( pval.steph8 < 0.01 )


    ci6 = ci_quantile( Z=electric_teachers$TxAny, Y=electric_teachers$gain, alternative="greater",
                       method.list=list( name="Stephenson", s=8 ), nperm=nperm, alpha=0.10 )
    ci2 = ci_quantile( Z=electric_teachers$TxAny, Y=electric_teachers$gain, alternative="greater",
                       method.list=list( name="Wilcoxon" ),
                       k.vec = 200:220,
                       nperm=nperm, alpha=0.10 )

    pr <- capture_output( print( ci6 ) )
    expect_true( stringr::str_detect(pr, "CI quantile on 233 units with 164 treated" ) )

    op <- capture_output( summary( ci6 ) )
    expect_true( stringr::str_detect( op, "greater" ) )

    cc <- capture_output( expect_warning( summary( ci2 ) ) )

    plot_quantile_CIs(ci6, k_start = 102, main = "s=8" )


    ci6 = ci_quantile( Z=electric_teachers$TxAny, Y=electric_teachers$gain, alternative="less",
                       method.list=list( name="Stephenson", s=8 ), nperm=nperm, alpha=0.10 )

    aa = capture_output( ss <- summary( ci6 ) )
    expect_true( ss$c$np == 0 )

})



