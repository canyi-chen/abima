test_that("multiplication works", {
  ### Test 1
  ## Set up parameters
  M.family <- poisson()
  alpha_S <- beta_M <- 1/8

  set.seed(2)
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

  out <- abYlm.Mglm(
    S,
    M,
    Y,
    X,
    M.family = M.family,
    B = 199
  )
  expect_equal(out$p_value, 0, tolerance = 1e-2)
  expect_equal(out$mediation_effect, 0.05440099, tolerance = 1e-2)

})
