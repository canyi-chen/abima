\donttest{
  ## Load libraries
  library(abima)

  set.seed(2)

  ## Set up parameters
  M.family <- poisson()

  simulation <- function(alpha_S = 0, beta_M = 0) {
    data <- generate_all_data(
      n = 500,
      alpha_S = alpha_S,
      beta_M = beta_M,
      M.family = M.family
    )
    S <- data$S
    M <- data$M
    Y <- data$Y
    X <- data$X

    out <- abYlm.Mglm(S,
                      M,
                      Y,
                      X,
                      M.family = M.family,
                      lambda = 2,
                      B = 199)
    out
  }


  simulation(1 / 8, 1 / 8)

  simulation(0, 0)

}
