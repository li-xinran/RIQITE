

#' Plot quantile CIs
#'
#' Given the results of a quantile CI call, plot the results.
#'
#' @param result Result of a
#' @param k_start Quantile to start with (if NULL, plot all informative CIs).
#' @param main Title of plot
#'
#' @return The result object, again (invisibly).
#'
#' @export
plot_quantile_CIs <- function( result, k_start = NULL, main = NULL ) {
  ci.limit = result$lower

  n = length( ci.limit )

  if ( is.null( k_start ) ) {
    k_start = n - sum( !is.nan( ci.limit ) & !is.infinite(ci.limit) )
  }

  ylim = c( k_start, n + 1)
  xlim = range(ci.limit[ci.limit>-Inf]) * 1.1
  plot(NA, ylab="k", xlab=expression( "lower"~"confidence"~"limit"~"for" ~tau[(k)] ),
       ylim=ylim,
       xlim=xlim,
       main = main )

  for(k in 1:length(ci.limit)){
    lines( c( max( ci.limit[k], min(ci.limit[ci.limit>-Inf]) - 100 ), max(ci.limit)+10),
           rep(k,2), col = "grey" )
  }
  points( ci.limit, c(1:length(ci.limit)), pch = 20 )
  abline(v = 0, lty = 2)

  invisible( result )
}

