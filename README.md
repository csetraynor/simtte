# simtte

<!-- badges: start -->
<!-- badges: end -->

**simtte** simulates time-to-event (survival) datasets for clinical trial
design and analysis. It supports Weibull and flexible M-spline baseline
hazard models, using [mrgsolve](https://github.com/metrumresearchgroup/mrgsolve)
as the ODE solver backend and inverse transform sampling to generate event
times.

## Installation

### Prerequisites

**simtte** depends on **mrgsolve**, which requires a working C++ compiler:

- **Windows**: Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- **macOS**: Install Xcode command line tools (`xcode-select --install`)
- **Linux**: Install `g++` via your package manager

### From GitHub

```r
# install.packages("devtools")
devtools::install_github("csetraynor/simtte")
```

## Quick Start

```r
library(simtte)

# Weibull simulation
set.seed(42)
lp <- matrix(rnorm(20, 0, 0.5), nrow = 20)
result <- sim_tte(
  pi = lp, mu = -1, coefs = 1.1,
  time = seq(0.1, 50, by = 0.1),
  type = "weibull", end_time = 50
)
head(result)

# M-splines simulation
data("ms_data")
lp <- matrix(runif(nrow(ms_data$basis)), nrow = nrow(ms_data$basis))
result <- sim_tte(
  pi = lp, mu = ms_data$mu, basis = ms_data$basis,
  coefs = ms_data$coefs, time = ms_data$time, type = "ms"
)
head(result)
```

## License

GPL (>= 2)
