
<!-- README.md is generated from README.Rmd. Please edit that file -->

The `MultiHorizonSPA` package allows R users to run the Multi Horizon
Superior Predictive Ability (SPA) test proposed by Quaedvlieg (2021):
compare the predictive performance of two distinct models when jointly
considering all horizons of a forecast path.

## Installation

You can install `MultiHorizonSPA` from CRAN as follows:

``` r
install.packages("MultiHorizonSPA")
```

or from GitHub:

``` r
install.packages("devtools")
devtools::install_github("lucabarbaglia/MultiHorizonSPA")
```

## A start-up example

Test for uniform SPA (uSPA).

``` r
library(MultiHorizonSPA)
data(LossDiff_uSPA)
Test_uSPA(LossDiff=LossDiff_uSPA, L=3)
#> $p_value
#> [1] 0.003003003
#> 
#> $t_uSPA
#> [1] 0.0804403
```

The output of the `Test_uSPA` function is a list containing two objects:

  - p-value: the p-value for uSPA;

  - t\_uSPA: the statistics for uSPA;

Now test for average SPA (aSPA).

``` r
library(MultiHorizonSPA)
data(LossDiff_aSPA)
weights <- t(as.matrix(rep(1, ncol(LossDiff_aSPA))/ncol(LossDiff_aSPA)))
Test_aSPA(LossDiff=LossDiff_aSPA, weights=weights, L=3)
#> $p_value
#> [1] 0.002002002
#> 
#> $t_aSPA
#> [1] 2.439664
```

The output of the `Test_aSPA` function is a list containing two objects:

  - p-value: the p-value for aSPA;

  - t\_aSPA: the statistics for aSPA.

## References:

  - Quaedvlieg, Rogier. “Multi-horizon forecast comparison.” Journal of
    Business & Economic Statistics 39.1 (2021): 40-53.
