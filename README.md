
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
Trow <- 200 
H <- 12
Mmethods <- 5
Losses <- matrix(rnorm(Trow*H, mean = 0), nrow = Trow, ncol = H)

Test_uSPA(LossDiff=Losses, L=3, B=5)
```

The output of the `Test_uSPA` function is a list containing two objects:

  - p-value: the p-value for uSPA;

  - t\_uSPA: the statistics for uSPA;

Now test for average SPA (aSPA).

``` r
library(MultiHorizonSPA)
Trow <- 200 
H <- 12
Mmethods <- 5
weights <- rep(1/H,H)
Losses <- matrix(rnorm(Trow*H, mean = 0), nrow = Trow, ncol = H)

Test_aSPA(LossDiff=Losses, weights=weights, L=3, B=5)
```

The output of the `Test_aSPA` function is a list containing two objects:

  - p-value: the p-value for aSPA;

  - t\_aSPA: the statistics for aSPA.


## Fast SPA tests


``` r
library(MultiHorizonSPA)
Trow <- 200 
H <- 12
Mmethods <- 5
Losses <- matrix(rnorm(Trow*H, mean = 0), nrow = Trow, ncol = H)

Fast_Test_uSPA(LossDiff=Losses, L=3, B=10)
```


``` r
library(MultiHorizonSPA)
Trow <- 200 
H <- 12
Mmethods <- 5
weights <- rep(1/H,H)
Losses <- matrix(rnorm(Trow*H, mean = 0), nrow = Trow, ncol = H)

Fast_Test_aSPA(LossDiff=Losses, weights=weights, L=3, B=10)
```


## Multiple Horizon Model Confidence Set

Note: This test can be computationally expensive. The Fast MCS test is reccommended.

``` r
library(MultiHorizonSPA)
Trow <- 20 
H <- 12
Mmethods <- 9
weights <- rep(1/H,H)

loss_list <- vector(mode = "list", length = Mmethods)

loss_list[[1]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
loss_list[[2]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
loss_list[[3]] <- matrix(rnorm(Trow*H, mean = 3), nrow = Trow, ncol = H)
loss_list[[4]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
loss_list[[5]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
loss_list[[6]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
loss_list[[7]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
loss_list[[8]] <- matrix(rnorm(Trow*H, mean = 3), nrow = Trow, ncol = H)
loss_list[[9]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
loss_list[[10]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)


MultiHorizonMCS(loss_list, L=3,B=5,unif_or_average = 'u')
#'
```


## Fast Multiple Horizon Model Confidence Set


``` r
library(MultiHorizonSPA)
Trow <- 20 
H <- 12
Mmethods <- 9
weights <- rep(1/H,H)

loss_list <- vector(mode = "list", length = Mmethods)

loss_list[[1]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
loss_list[[2]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
loss_list[[3]] <- matrix(rnorm(Trow*H, mean = 3), nrow = Trow, ncol = H)
loss_list[[4]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
loss_list[[5]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
loss_list[[6]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)
loss_list[[7]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
loss_list[[8]] <- matrix(rnorm(Trow*H, mean = 3), nrow = Trow, ncol = H)
loss_list[[9]] <- matrix(rnorm(Trow*H, mean = 2), nrow = Trow, ncol = H)
loss_list[[10]] <- matrix(rnorm(Trow*H, mean = 1), nrow = Trow, ncol = H)


num_cores <- 1


MultiHorizonSPA:::MultiHorizonMCS_cpp(loss_list, #
                                      0.05, # alpha_t
                                      0.05, # alpha_mcs
                                      weights, #
                                      3,#l
                                      5,#b
                                      1, # 1 means uniform
                                      num_cores,
                                      seed
)
```






## References:

  - Quaedvlieg, Rogier. “Multi-horizon forecast comparison.” Journal of
    Business & Economic Statistics 39.1 (2021): 40-53.
