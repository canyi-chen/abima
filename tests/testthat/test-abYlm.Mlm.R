test_that("multiplication works", {
  ### Test 1
  ## Set up parameters
  alpha_S <- beta_M <- 1/8

  set.seed(2)
  data <- generate_all_data(
    n = 200,
    alpha_S = alpha_S,
    beta_M = beta_M
  )
  S <- data$S
  M <- data$M
  Y <- data$Y
  X <- data$X

  out <- abYlm.Mlm(
    S,
    M,
    Y,
    X,
    B = 199
  )
  expect_equal(out$p_value_NIE, 0.0201005, tolerance = 1e-2)
  expect_equal(out$NIE, 0.01697224, tolerance = 1e-2)

  ### Test 2
  ## Set up parameters
  alpha_S <- beta_M <- rep(1/8, 2)

  set.seed(2)
  data <- generate_all_data(
    n = 200,
    alpha_S = alpha_S,
    beta_M = beta_M
  )
  S <- data$S
  M <- data$M
  Y <- data$Y
  X <- data$X

  out <- abYlm.Mlm(
    S,
    M,
    Y,
    X,
    B = 199
  )
  expect_equal(out$p_value_NIE, 0.05527638, tolerance = 1e-2)
  expect_equal(out$NIE, 0.01673136, tolerance = 1e-2)
})
