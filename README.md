# The *bbreg* package

### CRAN badges

[![](https://www.r-pkg.org/badges/version/bbreg?color=green)](https://cran.r-project.org/package=bbreg)
[![](http://cranlogs.r-pkg.org/badges/grand-total/bbreg?color=black)](https://cran.r-project.org/package=bbreg)
[![](http://cranlogs.r-pkg.org/badges/last-month/bbreg?color=orange)](https://cran.r-project.org/package=bbreg)
---
### Dev badges:
[![](https://img.shields.io/badge/devel%20version-2.0.2-blue.svg)](https://github.com/vpnsctl/bbreg/main)
[![R build status](https://github.com/vpnsctl/bbreg/workflows/R-CMD-check/badge.svg)](https://github.com/vpnsctl/bbreg/actions)
[![Build Status](https://travis-ci.com/vpnsctl/bbreg.svg?branch=main)](https://travis-ci.com/vpnsctl/bbreg)
---

The *bbreg* package deals with regression models with response variables being continuous
and bounded. It currently provides implementation of two regression models:  bessel regression <https://arxiv.org/abs/2003.05157> and beta
regression <https://doi.org/10.1080/0266476042000214501>. For both of these models, the estimation is 
done with the EM-algorithm. The EM-algorithm approach for beta regression
was developed in <https://doi.org/10.1080/00949655.2017.1350679>.

## Installing the *bbreg* package from CRAN:

```{r}
install.packages("bbreg")
```

## Installing the *bbreg* package from this repository:

To install the *bbreg* package from this repository, just run the command:

```{r}
#install.packages("devtools")
devtools::install_github("vpnsctl/bbreg")
```

To install the *bbreg* package from this repository **with vignette**, run the command:
```{r}
#install.packages("devtools")
devtools::install_github("vpnsctl/bbreg", build_vignettes = TRUE)
```

This repository will always contain the most recent version of the *bbreg* package.

## Vignette

The *bbreg* package has a vignette. Check the most recent version of the vignette at <https://rpubs.com/alexandrebsimas/intro-bbreg>

## Basic usage

The usage of the *bbreg* package is analogous to the usage of standard regression functions and packages in *R*:

```{r}
library(bbreg)
fit <- bbreg(agreement ~ priming + eliciting, data = WT)
fit
#> 
#> Bessel regression via EM - Model selected via Discrimination test (DBB)
#> 
#> Call: 
#> bbreg(agreement ~ priming + eliciting | 1)
#> 
#> Coefficients modeling the mean (with logit link):
#> (intercept)     priming   eliciting 
#>  -1.1537851  -0.2548633   0.3392742 
#> Coefficients modeling the precision (with identity link):
#> (intercept) 
#>    4.924825

summary(fit)
#> 
#> Bessel regression via EM - Model selected via Discrimination test (DBB):
#> Call:
#> bbreg(agreement ~ priming + eliciting | 1)
#> Number of iterations of the EM algorithm = 297
#> 
#>  Results of the discrimination test DBB:
#>     sum(z2/n) sum(quasi_mu)    |D_bessel|      |D_beta| 
#>        0.0853       52.6201        0.0004        0.0030 
#> 
#>  Pearson residuals:
#>      RSS      Min       1Q   Median       3Q      Max 
#> 331.9432  -1.7387  -0.6606  -0.3847   0.5562   4.5786 
#> 
#>  Coefficients modeling the mean (with logit link):
#>             Estimate Std.error z-value Pr(>|z|)    
#> (intercept) -1.15379   0.05251 -21.971  < 2e-16 ***
#> priming     -0.25486   0.05733  -4.446 8.77e-06 ***
#> eliciting    0.33927   0.05846   5.804 6.49e-09 ***
#> 
#>  Coefficients modeling the precision (with identity link):
#>             Estimate Std.error z-value Pr(>|z|)    
#> (intercept)   4.9248    0.4493   10.96   <2e-16 ***
#> g(phi) = 0.1316
#> ---
#>  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

For further details we refer the reader to the **vignette** whose link can be found above.
