
test_stat <- function(z, y, method.list = list(name = "Wilcoxon")){
  if(method.list$name == "DIM"){
    stat = mean(y[z==1]) - mean(y[z==0])
    return(stat)
  }

  if(method.list$name == "Wilcoxon"){
    score = rank_score( n, method.list = method.list )
    stat = sum( score[rank(y, ties.method = "first")[z==1]] )
    return(stat)
  }

  if(method.list$name == "Stephenson"){
    n = length(y)
    score = rank_score( n, method.list = method.list )
    stat = sum( score[rank(y, ties.method = "first")[z==1]] )
    return(stat)
  }
}

pval_imp_both <- function(Z, Y, c, method.list, Z.perm = NULL, iter.max = 10^5){
  Y1.imp = Y + (1-Z) * c
  Y0.imp = Y - Z * c
  m =  sum(Z)
  n = length( Z )
  stopifnot( length( Y ) == n )

  if( is.null(Z.perm) ){
    Z.perm = matrix(0, nrow = n, ncol = iter.max)
    for(iter in 1:iter.max){
      Z.perm[sample(n, m), iter] = 1
    }
  }

  iter.max = ncol(Z.perm)
  stat.perm = rep(NA, iter.max)
  for(iter in 1:iter.max){
    Z.iter = Z.perm[, iter]
    Y.imp = Z.iter * Y1.imp + (1-Z.iter) * Y0.imp
    stat.perm[iter] = test_stat(z = Z.iter, y = Y.imp, method.list)
  }

  stat.obs = test_stat(z = Z, y = Y, method.list)

  pval <- mean( stat.perm >= stat.obs )
  return(pval)
}


pval_imp_control <- function(Z, Y, c, method.list, Z.perm = NULL, iter.max = 10^5){
  Y1.imp = Y + (1-Z) * c
  Y0.imp = Y - Z * c
  m =  sum(Z)
  n = length( Z )
  stopifnot( length( Y ) == n )

  if( is.null(Z.perm) ){
    Z.perm = matrix(0, nrow = n, ncol = iter.max)
    for(iter in 1:iter.max){
      Z.perm[sample(n, m), iter] = 1
    }
  }

  iter.max = ncol(Z.perm)
  stat.perm = rep(NA, iter.max)
  for(iter in 1:iter.max){
    Z.iter = Z.perm[, iter]
    stat.perm[iter] = test_stat(z = Z.iter, y = Y0.imp, method.list)
  }

  stat.obs = test_stat(z = Z, y = Y0.imp, method.list)

  pval <- mean( stat.perm >= stat.obs )
  return(pval)
}


pval_max <- function(Z, Y, c, method.list, Z.perm = NULL, iter.max = 10^5, impute = "both"){

  if(impute == "both"){
    pval = pval_imp_both(Z = Z, Y = Y, c = c, method.list = method.list, Z.perm = Z.perm, iter.max = iter.max)
    return(pval)
  } else if(impute == "control"){
    pval = pval_imp_control(Z = Z, Y = Y, c = c, method.list = method.list, Z.perm = Z.perm, iter.max = iter.max)
    return(pval)
  } else {
      stop( glue::glue( "Unrecognized imputation procedure of '{impute}'" ) )
  }
}






ci_max_effect <- function(Z, Y, method.list, Z.perm = NULL, iter.max = 10^5, impute = "both", alpha = 0.05, tol = 10^(-3) ){

  if( is.null(Z.perm) ){
    Z.perm = matrix(0, nrow = n, ncol = iter.max)
    for(iter in 1:iter.max){
      Z.perm[sample(n, m), iter] = 1
    }
  }

  f <- function(c){
    pval = pval_max(Z = Z, Y = Y, c = c, method.list = method.list, Z.perm = Z.perm, impute = impute)
    return(pval - alpha)
  }

  c_sol <- uniroot(f, interval = c(-max(data$Y) + min(data$Y), max(data$Y) - min(data$Y)), extendInt = "upX", tol = tol)$root
  c_sol <- round(c_sol, digits = -log10(tol))
  if(f(c_sol) > 0){
    while( f(c_sol) > 0 ){
      c_sol = c_sol - tol
    }
    c_sol = c_sol + tol
  }else{
    while( f(c_sol) <= 0 ){
      c_sol = c_sol + tol
    }
  }
  return(c_sol)
}



#' Randomization test for bounded null hypotheses
#'
#' Obtain the p-value for testing the bounded null hypothesis H0: tau
#' <= c or H0: tau >= c
#'
#' @param Z An \eqn{n} dimensional treatment assignment vector.
#' @param Y An \eqn{n} dimensional observed outcome vector.
#' @param c A scalar or vector specifying the bounded null hypothesis.
#' @param alternative A character takes value "greater" and "less"
#'   indicating the direction of the alternative hypothesis.
#' @param method.list A list specifies the choice of the rank sum test
#'   statistic. For example, list(name = "DIM") means the
#'   difference-in-means statistic, list(name="Wilcoxon") means the
#'   Wilcoxon rank sum statistic, and list(name = "Stephenson", s =
#'   10) means the Stephenson rank sum statistic with parameter s=10.
#' @param Z.perm A \eqn{n \times nperm} matrix that specifies the
#'   permutated assignments for approximating the null distribution of
#'   the test statistic.
#' @param nperm A positive integer representing the number of
#'   permutations for approximating the randomization distribution of
#'   the test statistic
#' @param impute A character specifies the form of the test statistic.
#'   "control" means we use the imputed control potential outcomes in
#'   the test statistic, and "both" means we use the imputed observed
#'   outcome in the test statistic.
#' @return The p-value for testing the specified null hypothesis of
#'   interest.
#' @export
pval_bound <- function(Z, Y, c=0,
                       method.list = list(name = "DIM"),
                       alternative = "greater", Z.perm = NULL, nperm = 10^5, impute = "control"){

  if(alternative == "greater"){
    pval = pval_max(Z=Z, Y=Y, c=c, method.list=method.list, Z.perm = Z.perm, iter.max = nperm, impute = impute)
    return(pval)
  } else if(alternative == "less"){
    pval_max(Z=Z, Y=-1*Y, c=-1*c, method.list=method.list, Z.perm = Z.perm, iter.max = nperm, impute = impute)
    return(pval)
  } else {
      stop( glue::glue( "Unrecognized alternative of {alternative}" ) )
  }
}



#' Confidence interval for the maximum or minimum individual effect
#'
#' Obtain confidence interval for the maximum or minimum individual
#' effect
#'
#' @inheritParams pval_bound
#'
#' @param alpha A numeric object specifies the level of the confidence
#'   interval. In particular, the confidence level is \eqn{1-\alpha}.
#' @param tol A numerical object specifying the precision of the
#'   obtained confidence intervals. For example, if tol = 10^(-3),
#'   then the confidence limits are precise up to 3 digits.
#'
#' @return The p-value for testing the specified null hypothesis of
#'   interest.
#' @export
ci_bound <- function(Z, Y, alternative = "greater", method.list = list(name = "DIM"),
                     Z.perm = NULL, nperm = 10^3, impute = "control", alpha = 0.05, tol = 10^(-3) ){
  if(alternative == "greater"){
    ci.bound = ci_max_effect(Z=Z, Y=Y, method.list=method.list, Z.perm = Z.perm, iter.max = nperm, impute = impute, alpha = alpha, tol = tol )
    return(ci.bound)
  }
  if(alternative == "less"){
    ci.bound = -1 * ci_max_effect(Z=Z, Y=-1*Y, method.list=method.list, Z.perm = Z.perm, iter.max = nperm, impute = impute, alpha = alpha, tol = tol )
    return(ci.bound)
  }
}




