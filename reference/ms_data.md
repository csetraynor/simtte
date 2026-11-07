# Simulated Data for M-Splines Model Demonstration

A list containing parameters for an M-spline survival model, suitable
for use with
[`sim_tte`](https://csetraynor.github.io/simtte/reference/sim_tte.md).

## Usage

``` r
ms_data
```

## Format

A list with 4 elements:

- mu:

  Numeric scalar. Intercept parameter.

- basis:

  Numeric matrix. M-spline basis matrix with rows corresponding to time
  points and columns to basis functions.

- coefs:

  Numeric vector. Coefficients for each basis function.

- time:

  Numeric vector. Time points corresponding to the rows of the basis
  matrix.

## Source

Simulated example data.
