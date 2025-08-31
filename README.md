
<!-- README.md is generated from README.Rmd. Please edit that file -->
simtte
======

<!-- badges: start -->
<!-- badges: end -->
The goal of simtte is to fast simulate simple and complex time-to-event models also known as survival analysis. It also allows user-defined time-to-event models via the awesome package [**mrgsolve**](https://github.com/metrumresearchgroup/mrgsolve).

Installation
------------

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("csetraynor/simtte")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(simtte)
# Load M-splines data example
data("ms_data")
mu <- ms_data$mu
basis <- ms_data$basis
coefs <- ms_data$coefs
time <- ms_data$time
# Simulate prognostic index, linear predictor
lp <- matrix(runif(nrow(basis)),  nrow = nrow(basis))

# simulate M-splines model
sim_tte(pi = lp, mu = mu, basis = basis, coefs = coefs, time = time, type = "ms", end_time = 100)
#> Compiling ms ... done.
#> # A tibble: 299 x 4
#>    sim_time sim_status    ID    lp
#>       <dbl>      <dbl> <int> <dbl>
#>  1    1.03           1     1 0.406
#>  2    2.59           1     2 0.696
#>  3    6.73           0     3 0.840
#>  4    3.50           1     4 0.973
#>  5    6.73           0     5 0.935
#>  6    0.745          1     6 0.290
#>  7    0.562          1     7 0.594
#>  8    2.24           1     8 0.299
#>  9    6.73           0     9 0.706
#> 10    6.73           0    10 0.735
#> # b