#' An example function for generating the potential outcomes
#'
#' A example function for generating the potential outcomes, which can be used for in the function simu_power for power comparison.
#'
#' @param n A numerical object specifying the sample size.
#'
#' @export
normal_gen_PO <- function(n){
  mu0 = 0; sd0 = 1
  mu1 = 1; sd1 = 1

  Y1 = rnorm(n, mean = mu1, sd = sd1)
  Y0 = rnorm(n, mean = mu0, sd = sd0)

  return(data.frame(Y1=Y1, Y0=Y0))
}


#' Simulation for various Stephenson rank statistics
#'
#' Conduct simulation for inferring quantiles of individual treatment effects using various
#' Stephenson rank statistics.
#'
#' @param gen A function that generates the potential outcomes.
#' @param n A numerical object specifying the sample size.
#' @param treat.prop A numerical object specifying the proportion of treated units.
#' @param iter.max A numerical object specifying the number of simulated data sets.
#' @param Steph.s.vec A vector specifying the parameters for a sequence of Stephenson rank statistics under investigation.
#' @param alternative A character takes value "greater" and "less", indicating the
#' direction of alternative hypotheses.
#' @param alpha A numerical object specifying the significance level.
#' @param nperm A positive integer representing the number of permutations for
#'   approximating the randomization distribution of the rank sum statistic.
#' @param tol A numerical object specifying the precision of the obtained
#'   confidence intervals. For example, if tol = 10^(-3), then the confidence
#'   limits are precise up to 3 digits.
#' @param switch A logical object indicating whether performing treatment label
#'   switching to make the treated group have larger size.
#' @param print.iter An integer determining how frequently we report the progress of the simulation.
#' @param keep_data TRUE means keep the datasets themselves.  FALSE means do not keep.
#'
#' @return A list with the elements
#'   \item{Steph.s.vec}{An vector specifying the parameters for a sequence of Stephenson rank
#'   statistics under investigation.}
#'   \item{dat.all}{A list containing all simulated data sets.}
#'   \item{Z.CRE}{A \eqn{n} by iter.max matrix containing all simulated assignments.}
#'   \item{method.list.all}{A list describing the Stephenson rank statistics under investigation.}
#'   \item{ci.upper.all}{A list containing the inferred upper confidence limits for all quantiles of
#'   individual effects for all simulated data sets using the specified Stephenson rank statistics.}
#'   \item{ci.lower.all}{A list containing the inferred lower confidence limits for all quantiles of
#'   individual effects for all simulated data sets using the specified Stephenson rank statistics.}
#'
#' @export
simu_power <- function(gen = normal_gen_PO, n = 100, treat.prop = 0.5, iter.max = 100,
                        Steph.s.vec = c(2, 5, 10),
                        alternative = "greater", alpha = 0.05,
                        nperm = 10^3, tol = 10^(-3), switch = TRUE, print.iter = 10,
                       keep_data = FALSE
                        ){

  m = floor( n * treat.prop )
  Z.perm = assign_CRE(n, m, nperm)
  stat.null.all = list()
  method.list.all = list()
  for( s.ind in 1:length(Steph.s.vec) ){
    method.list.all[[s.ind]] = list( name = "Stephenson", s = Steph.s.vec[s.ind] )
    stat.null.all[[s.ind]] = null_dist(n, m, method.list = method.list.all[[s.ind]], Z.perm = Z.perm )
  }
  cat("Monte Carlo approximation for null distributions completed...\n")

  ci.upper.all = list()
  ci.lower.all = list()
  for( s.ind in 1:length(Steph.s.vec) ){
    ci.upper.all[[s.ind]] = matrix(NA, nrow = n, ncol = iter.max)
    ci.lower.all[[s.ind]] = matrix(NA, nrow = n, ncol = iter.max)
  }
  Z.CRE = assign_CRE(n, m, iter.max)

  dat.all = list()
  for(iter in 1:iter.max){
    dat = gen(n)
    dat.all[[iter]] = dat

    Z = Z.CRE[, iter]
    Y = Z * dat$Y1 + (1-Z) * dat$Y0

    # pval_max(Z = Z, Y = Y, c = 0, method.list = method.list.all[[s.ind]], Z.perm = Z.perm, impute = "control")

    for( s.ind in 1:length(Steph.s.vec) ){
      ci.iter.s = ci_quantile(Z = Z, Y = Y,
                            alternative = alternative, method.list = method.list.all[[s.ind]],
                            stat.null = stat.null.all[[s.ind]], alpha = alpha,
                            tol = tol, switch = switch)
      ci.upper.all[[s.ind]][, iter] = ci.iter.s$upper
      ci.lower.all[[s.ind]][, iter] = ci.iter.s$lower
    }

    if( iter %% print.iter == 0 ){
      cat( paste0( "Simulation for ", iter, "th iteration (of ", iter.max, ") completed...\n" ) )
    }
  }

  # for( s.ind in 1:length(Steph.s.vec) ){
  #   print( mean( colSums(ci.lower.all[[s.ind]] > 0) ) )
  # }

  res <- list(
    Steph.s.vec = Steph.s.vec,
    Z.CRE = Z.CRE,
    Z.perm = Z.perm,
    method.list.all = method.list.all,
    ci.upper.all = ci.upper.all,
    ci.lower.all = ci.lower.all
    )
  if ( keep_data ) {
    res$dat.all = dat.all
  }
  return( res )
}




#' Comparing the performance of various Stephenson rank statistics
#'
#' Summarize the simulation to compare the performance of various Stephenson rank statistics.
#'
#' @param result A list containing the result from the simu_power function.
#' @param measure A character specifying how we compare different Stephenson rank statistics.
#' If it equals "CI_n(c)", then we compare the average confidence limits for the number of units with effects passing the threshold cutoff.
#' If it equals "CI_tau(k)", then we compare the average confidence limits for the individual treatment effect at rank k.
#' If it equals "test_tau(k)", then we compare the power for testing whether the individual treatment effect at rank k passes the threshold cutoff.
#' @param alternative A character takes value "greater" and "less", indicating the
#' direction of alternative hypotheses.
#' @param cutoff A numerical object specifying the cutoff.
#' @param k An integer that specifies the quantile of individual effects under investigation.
#'
#' @export
summary_power <- function(result, measure = c( "CI_n(c)", "test_tau(k)", "CI_tau(k)" ),
                          alternative = "greater", cutoff = 0, k = NULL){

  measure = match.arg(measure)

  if(measure == "CI_n(c)" & alternative == "greater"){
    tab = data.frame(s = result$Steph.s.vec, ci.mean = rep(NA, length(result$Steph.s.vec)))
    for(i in 1:nrow(tab)){
      tab$ci.mean[i] = mean( colSums( result$ci.lower.all[[i]] > cutoff ) )
    }
    return(tab)
  }

  if(measure == "CI_n(c)" & alternative == "less"){
    tab = data.frame(s = result$Steph.s.vec, ci.mean = rep(NA, length(result$Steph.s.vec)))
    for(i in 1:nrow(tab)){
      tab$ci.mean[i] = mean( colSums( result$ci.upper.all[[i]] < cutoff ) )
    }
    return(tab)
  }

  if ( is.null( k ) || !is.numeric( k ) ) {
    stop( sprintf( "Need to supply k (from 1 to %d) for quantile to target for summary calculation",
                   nrow( result$Z.CRE ) ) )
  }

  if(measure == "test_tau(k)" & alternative == "greater"){
    tab = data.frame(s = result$Steph.s.vec, power.mean = rep(NA, length(result$Steph.s.vec)))
    for(i in 1:nrow(tab)){
      tab$power.mean[i] = mean( result$ci.lower.all[[i]][k, ] > cutoff )
    }
    return(tab)
  }

  if(measure == "test_tau(k)" & alternative == "less"){
    tab = data.frame(s = result$Steph.s.vec, power.mean = rep(NA, length(result$Steph.s.vec)))
    for(i in 1:nrow(tab)){
      tab$power.mean[i] = mean( result$ci.upper.all[[i]][k, ] < cutoff )
    }
    return(tab)
  }

  if(measure == "CI_tau(k)" & alternative == "greater"){
    tab = data.frame(s = result$Steph.s.vec, ci.mean = rep(NA, length(result$Steph.s.vec)))
    for(i in 1:nrow(tab)){
      tab$ci.mean[i] = median( result$ci.lower.all[[i]][k, ] )
    }
    return(tab)
  }

  if(measure == "CI_tau(k)" & alternative == "less"){
    tab = data.frame(s = result$Steph.s.vec, ci.mean = rep(NA, length(result$Steph.s.vec)))
    for(i in 1:nrow(tab)){
      tab$ci.mean[i] = median( result$ci.upper.all[[i]][k, ] )
    }
    return(tab)
  }
}

#### Testing Code ####


if(FALSE){

  library( RIQITE )

  normal_gen_PO <- function(n){
    mu0 = 0; sd0 = 1
    mu1 = 1; sd1 = 0.5

    Y1 = rnorm(n, mean = mu1, sd = sd1)
    Y0 = rnorm(n, mean = mu0, sd = sd0)

    return(data.frame(Y1=Y1, Y0=Y0))
  }

  # Small number of iters for easier exploration
  result = simu_power( gen = normal_gen_PO, n = 120, treat.prop = 0.5, iter.max = 100,
                       Steph.s.vec = c(2, 10, 60),
                       alternative = "greater", alpha = 0.05,
                       nperm = 10^3, tol = 10^(-3), switch = TRUE )

  summary_power(result, measure = "CI_n(c)", alternative = "greater", cutoff = 0)
  debugonce( summary_power )
  summary_power(result, measure = "test_tau(k)", k = 120 )
  summary_power(result, measure = "CI_tau(k)", k = 120 )

  names( result )
  dim( result$Z.CRE )
  dim( result$Z.perm )  # null distribution
  result$method.list.all

  cil <- result$ci.lower.all
  length( cil )
  dim( cil[[1]] )
  cil[[1]][,1]

  # one-sided, so all upper is infinity.
  result$ci.upper.all[[1]]



  ### NOT SURE WHAT FOLLOWING CODE DOES ###
  library(utils)

  # para.table = rbind( expand.grid( c(-1, 0, 1, 2), 1, c(0:9)/2 ),
  #        expand.grid( c(-1, 0, 1, 2), c(0:9)/2, 1 ) )

  sigma1.grid = c( c(1:4)/5, 1, c(3:9)/2 )
  para.table = expand.grid( c(-1, 0, 1, 2), 1, sigma1.grid )
  colnames(para.table) = c("tau0", "sigma0", "sigma1")
  para.table

  jobid = 1

  normal_gen_PO <- function(n){
    mu0 = 0; sd0 = para.table$sigma0[jobid]
    mu1 = para.table$tau0[jobid]; sd1 = para.table$sigma1[jobid]

    Y1 = rnorm(n, mean = mu1, sd = sd1)
    Y0 = rnorm(n, mean = mu0, sd = sd0)

    return(data.frame(Y1=Y1, Y0=Y0))
  }

}





