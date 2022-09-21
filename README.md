## This is an R package for conducting randomization inference for quantiles of individual treatment effects


### Installation

```
devtools::install_github("li-xinran/RIQITE")
```

### load the package

```
library(RIQITE)
```

### explanation for the main function

```
?pval_bound # get the p value for testing bounded null

?ci_bound # get the confidence interval for the maximum (or minimum) individual effect

?pval_quantile # get the p value for testing null hypotheses on quantiles of individual effects

?ci_quantile # get the confidence interval for quantiles of individual effects

?simu_power # conduct simulation to compare tests using different Stephenson rank sum statistics 

?summary_power # summarize the simulation to compare performance of different Stephenson rank sum statistics 
```


### A simple example 
#### generate the observed data

```
n = 200

m = n * 0.5

Z = sample( c( rep(1, m), rep(0, n-m) ) )

Y = rnorm(n) + Z * rnorm(n, mean = 0, sd = 5)

#### choose the test statistic
method.list = list( name = "Stephenson", s = 10 )

#### test the null hypothesis that the 95\% quantile of individual effect is less than or equal to 0
pval = pval_quantile(Z, Y, k = ceiling(n*0.95), c = 0, alternative = "greater", method.list = method.list, nperm = 10^5 )
pval

#### construct simultaneous confidence intervals for all quantiles of individual effects
ci = ci_quantile( Z = Z, Y = Y, alternative = "two.sided", method.list = method.list, nperm = 10^6,  alpha = 0.05 )
ci
```


