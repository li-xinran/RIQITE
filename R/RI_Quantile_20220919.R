###################################################################
####### basics for testing quantiles of indivial effects ##########
###################################################################

#' score function of the ranks
#'
#' Calculate the scores of n units (ordered from low to high score) given a
#' specified test statistic, not allowing for ties.
rank_score <- function( n, method.list = list(name = "Wilcoxon") ){

  if(method.list$name == "Wilcoxon"){
    score = c(1:n)
    score = score/max(score)
    return(score)
  }

  if(method.list$name == "Stephenson"){
    score = choose( c(1:n) - 1, method.list$s - 1 )
    score = score/max(score)
    return(score)
  }

}


#' Generate matrix of complete randomization assignments
#'
#' @param n Sample size
#' @param m Number units treated
#' @param nperm number of permutations.  Inf will generate all possible
#'   permutations.
assign_CRE <- function(n, m, nperm){
  if(is.finite(nperm)){
    Z.perm = matrix(0, nrow = n, ncol = nperm)
    for(iter in 1:nperm){
      Z.perm[sample( c(1:n), m, replace = FALSE ), iter] = 1
    }
  }

  if(is.infinite(nperm)){
    comb.all = combn(n, m)
    nperm = ncol(comb.all)
    Z.perm = matrix(0, nrow = n, ncol = nperm)
    for(iter in 1:nperm){
      Z.perm[comb.all[, iter], iter] = 1
    }
  }

  return(Z.perm)
}


#' rand dist of the rank score stat
#'
#'   Generate the null distribution of the given test statistic for an
#'   experiment with m treated out of n units.
#'
#' @param Z.perm Matrix of possible treatment assignments.  NULL will generate
#'   complete rand with nperm combinations.
#'
#' @param score The list of possible scores of the units (under a null,
#'   invariant under treatment assignment).  If null will generate from the
#'   passed method.list and given n, m.
#'
#'
null_dist <- function(n, m, method.list = NULL, score = NULL, nperm = 10^5, Z.perm = NULL){
  if(is.null(score)){
    score = rank_score( n, method.list )
  }

  if(is.null(Z.perm)){
    Z.perm = assign_CRE(n, m, nperm)
  } else {
    stopifnot( is.matrix(Z.perm) )
  }

  stat_null = rep(NA, ncol(Z.perm))
  for(iter in 1:ncol(Z.perm)){
    stat_null[iter] = sum( score[ Z.perm[, iter] == 1 ] )
  }

  return(stat_null)
}




#' @title sort the treated units using the "first" method
sort_treat <- function(Y, Z){
  r = rank(Y, ties.method = "first")
  ind.sort = sort.int(r, index.return = TRUE)$ix
  ind.sort.treat = ind.sort[Z[ind.sort] == 1]
  return(ind.sort.treat)
}



#' @title minimum stat value in H_{k,c}
min_stat <- function(Z, Y, k, c, method.list = NULL, score = NULL, ind.sort.treat = NULL){
  n = length(Y)
  m = sum(Z)
  if(is.null(score)){
    score = rank_score( n, method.list )
  }

  # sort the treated units
  if(is.null(ind.sort.treat)){
    ind.sort.treat = sort_treat(Y, Z)
  }

  # get xi vector
  xi = rep(c, n)
  if(k < n){
    xi[ ind.sort.treat[ ( m + 1 - min(m, n-k) ):m ] ] = Inf
  }

  stat.min = sum( score[ rank( Y - Z * xi, ties.method = "first" )[Z==1] ] )
  stat.min

  return(stat.min)
}



### p-val for testing tau_{(k)} <= c ###
pval_H_k_c <- function(Z, Y, k, c, method.list = NULL, score = NULL, stat.null = NULL, nperm = 10^5, Z.perm = NULL, ind.sort.treat = NULL){
  n = length(Z)
  m = sum(Z)

  # get score if score is null #
  if(is.null(score)){
    score = rank_score( n, method.list )
  }

  # emp null dist #
  if(is.null(stat.null)){
    stat.null = null_dist(n, m, score = score, nperm = nperm, Z.perm = Z.perm)
  }

  # min stat value under H_{k,c} #
  stat.min = min_stat(Z, Y, k, c, score = score, ind.sort.treat = ind.sort.treat)

  # p-value #
  pval = mean( stat.null >= stat.min )

  return(pval)
}


#' @title one-sided confidence interval (c, Inf) for all quantiles
#'
#' @description generate one-sided confidence intervals
#'
#' @param Z treatment assignment (vector)
#' @param Y outcome (vector)
conf_quant_larger <- function( Z, Y, k.vec = NULL, method.list = NULL, score = NULL, stat.null = NULL, nperm = 10^5,
                               Z.perm = NULL, alpha = 0.05, tol = 10^(-3), ind.sort.treat = NULL ){
  n = length(Z)
  m = sum(Z)

  # get score if score is null #
  if(is.null(score)){
    score = rank_score( n, method.list )
  }

  # emp null dist #
  if(is.null(stat.null)){
    stat.null = null_dist(n, m, score = score, nperm = nperm, Z.perm = Z.perm)
    nperm = length(stat.null)
  }else{
    nperm = length(stat.null)
  }

  # find threshold such that p-value <= alpha <===> test stat > threshold #
  thres = sort(stat.null, decreasing = TRUE)[ floor(nperm * alpha) + 1 ]

  # range of c #
  Y1.max = max(Y[Z==1])
  Y1.min = min(Y[Z==1])
  Y0.max = max(Y[Z==0])
  Y0.min = min(Y[Z==0])

  c.max = Y1.max - Y0.min + tol
  c.min = Y1.min - Y0.max - tol

  

  # sort the treated units
  if(is.null(ind.sort.treat)){
    ind.sort.treat = sort_treat(Y, Z)
  }
  
  if( is.null(k.vec) ){
    
    # conf interval for all quantiles tau_{(k)} #
    c.limit = rep(NA, n)
    
    for(k in n:(n-m)){
      # define the target fun #
      # f > 0 <==> p-value <= alpha
      # f decreases in c, p value increases in c
      if(k < n){
        c.max = c.limit[k+1]
      }
      f <- function(c){
        stat.min = min_stat(Z, Y, k, c, score = score, ind.sort.treat)
        return(stat.min - thres)
      }
      # check whether f(-Inf) = f(c.min) > 0 #
      if( f(c.min) <= 0 ){
        c.sol = -Inf
      }
      if( f(c.max) > 0){
        c.sol = c.max
      }
      if( f(c.min) > 0 & f(c.max) <= 0 ){
        c.sol = uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        # find the min c st p-value > alpha <==> f <= 0 #
        c.sol = round(c.sol, digits = -log10(tol))
        if( f(c.sol) <= 0 ){
          while( f(c.sol) <= 0 ){
            c.sol = c.sol - tol
          }
          c.sol = c.sol + tol
        }else{
          while(f(c.sol) > 0 ){
            c.sol = c.sol + tol
          }
        }
      }
      c.limit[k] = c.sol
    }
    
    if( n-m > 1){
      c.limit[1:(n-m-1)] = c.limit[n-m]
    }
    
    c.limit[c.limit > (Y1.max - Y0.min) + tol/2] = Inf
  }else{
    
    k.vec.sort = sort(k.vec, decreasing = FALSE)
    j.max = length(k.vec.sort)
    j.min = max( sum(k.vec <= (n-m)), 1) 
    
    # conf interval for tau_{(k)} with k in k.vec#
    c.limit = rep(NA, j.max)
    
    for(j in j.max:j.min){
      # define the target fun #
      # f > 0 <==> p-value <= alpha
      # f decreases in c, p value increases in c
      k = k.vec.sort[j]
      if(j < j.max){
        c.max = c.limit[j+1]
      }
      f <- function(c){
        stat.min = min_stat(Z, Y, k, c, score = score, ind.sort.treat)
        return(stat.min - thres)
      }
      # check whether f(-Inf) = f(c.min) > 0 #
      if( f(c.min) <= 0 ){
        c.sol = -Inf
      }
      if( f(c.max) > 0){
        c.sol = c.max
      }
      if( f(c.min) > 0 & f(c.max) <= 0 ){
        c.sol = uniroot(f, interval = c(c.min, c.max), extendInt = "downX", tol = tol)$root
        # find the min c st p-value > alpha <==> f <= 0 #
        c.sol = round(c.sol, digits = -log10(tol))
        if( f(c.sol) <= 0 ){
          while( f(c.sol) <= 0 ){
            c.sol = c.sol - tol
          }
          c.sol = c.sol + tol
        }else{
          while(f(c.sol) > 0 ){
            c.sol = c.sol + tol
          }
        }
      }
      c.limit[j] = c.sol
    }
    
    
    if(j.min > 1){
      c.limit[1:(j.min-1)] = c.limit[j.min]
    }
    
    c.limit[c.limit > (Y1.max - Y0.min) + tol/2] = Inf
  }
  
  return( c.limit )
}

###################################################################
#####  functions for testing quantiles of individual effects ######
###################################################################
### p value for testing tau_{(k)} <= or >= c ###


#' Randomization test for quantiles of individual treatment effects
#'
#' Obtain the p-value for testing the null hypothesis H0: tau_{(k)} <= c, or H0:
#' tau_{(k)} >= c, or H0: tau_{(k)} = c, where tau_{(k)} denotes individual
#' treatment effect at rank k.
#'
#' @param Z An \eqn{n} dimensional treatment assignment vector.
#' @param Y An \eqn{n} dimensional observed outcome vector.
#' @param k An integer between 1 and n specifying which quantile of individual
#'   effect is of interest.
#' @param c A numerical object specifying the threshold for the null hypothesis.
#' @param alternative A character takes value "greater", "less" and "two.sided",
#'   indicating the alternative hypothesis.
#' @param method.list A list specifies the choice of the rank sum test
#'   statistic. For example, list(name="Wilcoxon") means the Wilcoxon rank sum
#'   statistic, and list(name = "Stephenson", s = 10) means the Stephenson rank
#'   sum statistic with parameter s=10.
#' @param score An \eqn{n} dimensional transformed ranks, i.e., (phi(1), phi(2),
#'   ..., phi(n)), where phi() denotes the rank transformation function.
#' @param stat.null An vector whose empirical distribution approximates the
#'   randomization distribution of the rank sum statistic.
#' @param nperm A positive integer representing the number of permutations for
#'   approximating the randomization distribution of the rank sum statistic.
#' @param Z.perm A \eqn{n \times nperm} matrix that specifies the permutated assignments for approximating the null distribution of the test statistic.
#' @param switch A logical object indicating whether performing treatment label
#'   switching to make the treated group have larger size.
#'
#' @return The p-value for testing the specified null hypothesis of interest.
#' @export
pval_quantile <- function(Z, Y, k, c, alternative = "greater",
                          method.list = list( name = "Stephenson", s = 10 ),
                          score = NULL, stat.null = NULL, nperm = 10^6, Z.perm = NULL, 
                          switch = TRUE ){

  # Check inputs
  if ( is.null( Y ) || !is.numeric( Y ) ) {
    stop( "Need to pass Y as numeric value for outcome" )
  }
  if ( !(length(unique(Z)) == 2 &&
               all( sort( unique( Z ) ) == c(0,1) ) ) ) {
    stop( "Need to pass Z as vector of 0/1 values for control/treatment assignment" )
  }
  stopifnot( length(Y) == length(Z) )

  n = length(Z)

  # switching for larger treated group #
  if( sum(Z==1) < n/2 & switch ){
    Z = 1-Z
    Y = -1 * Y
    if(!is.null(stat.null)){
      if(is.null(score)){
        score = rank_score( n, method.list )
      }
      stat.null = sum(score) - stat.null
    }
  }

  if(alternative == "greater"){
    pval = pval_H_k_c(Z = Z, Y = Y, k = k, c = c, method.list = method.list, score = score,
                      stat.null = stat.null, nperm = nperm, Z.perm = Z.perm)
    return(pval)
  }

  if(alternative == "less"){
    pval = pval_H_k_c(Z = Z, Y = -1 * Y, k = n+1-k, c = -1 * c, method.list = method.list, score = score,
                      stat.null = stat.null, nperm = nperm, Z.perm = Z.perm)
    return(pval)
  }

  if(alternative == "two.sided"){
    pval.greater = pval_H_k_c(Z = Z, Y = Y, k = k, c = c, method.list = method.list, score = score,
                              stat.null = stat.null, nperm = nperm, Z.perm = Z.perm)
    pval.less = pval_H_k_c(Z = Z, Y = -1 * Y, k = n+1-k, c = -1 * c, method.list = method.list, score = score,
                           stat.null = stat.null, nperm = nperm, Z.perm = Z.perm)
    pval = 2 * min(pval.greater, pval.less)
    return(pval)
  }
}



#' Randomization inference for quantiles of individual treatment effects
#'
#' Obtain one-sided or two-sided confidence intervals for all quantiles of
#' individual treatment effects
#'
#' @param Z An \eqn{n} dimensional treatment assignment vector.
#' @param Y An \eqn{n} dimensional observed outcome vector.
#' @param alternative A character takes value "greater", "less" and "two.sided",
#'   indicating whether the confidence intervals are one-sided with lower or
#'   upper confidence limits or two-sided.
#' @param method.list A list specifies the choice of the rank sum test
#'   statistic. For example, list(name="Wilcoxon") means the Wilcoxon rank sum
#'   statistic, and list(name = "Stephenson", s = 10) means the Stephenson rank
#'   sum statistic with parameter s=10.
#' @param score An \eqn{n} dimensional transformed ranks, i.e., (phi(1), phi(2),
#'   ..., phi(n)), where phi() denotes the rank transformation function.
#' @param stat.null An vector whose empirical distribution approximates the
#'   randomization distribution of the rank sum statistic.
#' @param nperm A positive integer representing the number of permutations for
#'   approximating the randomization distribution of the rank sum statistic.
#' @param Z.perm A \eqn{n \times nperm} matrix that specifies the permutated assignments for approximating the null distribution of the test statistic.
#' @param alpha A numerical object, where 1-alpha indicates the confidence
#'   level.
#' @param tol A numerical object specifying the precision of the obtained
#'   confidence intervals. For example, if tol = 10^(-3), then the confidence
#'   limits are precise up to 3 digits.
#' @param switch A logical object indicating whether performing treatment label
#'   switching to make the treated group have larger size.
#'
#' @return A list with the elements \item{lower}{An \eqn{n} dimensional vector
#'   consisting of lower confidence limits for individual effects sorted
#'   increasingly.} \item{upper}{An \eqn{n} dimensional vector consisting of
#'   upper confidence limits for individual effects sorted increasingly.}
#'
#' @export
ci_quantile <- function(Z, Y, k.vec = NULL, alternative = "greater", method.list = list( name = "Stephenson", s = 10 ),
                        score = NULL, stat.null = NULL, nperm = 10^6, Z.perm = NULL, alpha = 0.05,
                        tol = 10^(-3), switch = TRUE){
  # Check inputs
  if ( is.null( Y ) || !is.numeric( Y ) ) {
    stop( "Need to pass Y as numeric value for outcome" )
  }
  if ( !(length(unique(Z)) == 2 &&
         all( sort( unique( Z ) ) == c(0,1) ) ) ) {
    stop( "Need to pass Z as vector of 0/1 values for control/treatment assignment" )
  }
  stopifnot( length(Y) == length(Z) )

  n = length(Z)

  # switching for larger treated group #
  if( sum(Z==1) < n/2 & switch ){
    Z = 1-Z
    Y = -1 * Y
    if(!is.null(stat.null)){
      if(is.null(score)){
        score = rank_score( n, method.list )
      }
      stat.null = sum(score) - stat.null
    }
  }
  
  
  
  if(is.null(k.vec)){
    conf.int = data.frame(k = c(1:n), lower = rep(NA, n), upper = rep(NA, n))
  }else{
    k.vec = sort(k.vec, decreasing = FALSE)
    conf.int = data.frame(k = k.vec, lower = rep(NA, length(k.vec)), upper = rep(NA, length(k.vec)))
  }
  

  if( alternative == "greater" ){
    c.lower = conf_quant_larger(Z = Z, Y = Y, k.vec = k.vec, method.list = method.list, score = score,
                                stat.null = stat.null, nperm = nperm, Z.perm = Z.perm, alpha = alpha, tol = tol )
    conf.int$upper = Inf
    conf.int$lower = c.lower
    return(conf.int)
  }

  if( alternative == "less" ){
    if(is.null(k.vec)){
      k.vec.rev = NULL
    }else{
      k.vec.rev = n+1-k.vec
    }
    c.upper = -1 * rev( conf_quant_larger(Z = Z, Y = -1 * Y, k.vec = k.vec.rev, method.list = method.list, score = score,
                                          stat.null = stat.null, nperm = nperm, Z.perm = Z.perm,  alpha = alpha, tol = tol ) )
    conf.int$upper = c.upper
    conf.int$lower = -Inf
    return(conf.int)
  }

  if(alternative == "two.sided"){
    c.lower = conf_quant_larger(Z = Z, Y = Y, k.vec = k.vec, method.list = method.list, score = score,
                                stat.null = stat.null, nperm = nperm, Z.perm = Z.perm, alpha = alpha/2, tol = tol )
    
    if(is.null(k.vec)){
      k.vec.rev = NULL
    }else{
      k.vec.rev = n+1-k.vec
    }
    c.upper = -1 * rev( conf_quant_larger(Z = Z, Y = -1 * Y, k.vec = k.vec.rev, method.list = method.list, score = score,
                                          stat.null = stat.null, nperm = nperm, Z.perm = Z.perm,  alpha = alpha/2, tol = tol ) )
    conf.int$upper = c.upper
    conf.int$lower = c.lower
    return(conf.int)
  }
}


if( FALSE ){
  n = 100
  m = 50
  Z = sample( c(rep(1, m), rep(0, n-m)) )
  Y = rnorm(n) + 2 * Z
  method.list = list( name = "Stephenson", s = 6 )
  
  k.vec = c(1:n)
  alternative = "two.sided"
  
  stat.null = null_dist( n, m, method.list = method.list, nperm = 10^3 )
  
  a = ci_quantile(Z, Y, alternative = alternative, method.list = method.list,
                         stat.null = stat.null, alpha = 0.05,
                          tol = 10^(-3), switch = FALSE)
  
  b = ci_quantile(Z, Y, k.vec = k.vec, alternative = alternative, method.list = method.list,
                  stat.null = stat.null, alpha = 0.05,
                  tol = 10^(-3), switch = FALSE)
  
  cbind(a[k.vec, ], b)
}
