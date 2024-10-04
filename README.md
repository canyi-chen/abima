
# abmediation

<!-- badges: start -->
<!-- badges: end -->

The goal of abmediation is to ...

## Installation

You can install the development version of abmediation like so:

``` r
# Install abmediation from GitHub:
# install.packages("pak")
devtools::install_github("canyi-chen/abmediation")
```

## Example

This is a basic example which shows you how to apply the adaptive bootstrap for testing no mediation effect under the classical linear structral equation model:

``` r
library(abmediation)

## Set up parameters
M.family <- gaussian()
Y.family <- gaussian()
alpha_S <- beta_M <- 1/8

set.seed(2)
data <- generate_all_data(
  n = 500,
  alpha_S = alpha_S,
  beta_M = beta_M,
  M.family = M.family,
  Y.family = Y.family
)
S <- data$S
M <- data$M
Y <- data$Y
X <- data$X

ab.mediation(
  S,
  M,
  Y,
  X,
  M.family = M.family,
  Y.family = Y.family,
  lambda = 2,
  B = 199
)
```

