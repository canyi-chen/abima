# Function to compute coefficient matrix Q for multiple Y using the least squares approach
Q_function_multi <- function(X, Y) {
  stopifnot(nrow(X) == nrow(Y))
  return(solve(t(X) %*% X, t(X) %*% Y))
}

# Function to compute the residual (perpendicular) component of Y on X for multiple mediators
perp_function_multi <- function(X, Y) {
  return(Y - X %*% Q_function_multi(X, Y))
}

# Main PoC adaptive bootstrap function for multiple mediators
#' @importFrom stats coef resid as.formula lm
PoC_AB_multi <- function(S, M, Y, X, B = 500, lambda = 2) {
  data <- data.frame(S = S, M = M, Y = Y, X = X)
  J <- ncol(M)  # Number of mediators
  M_variables <- paste0("M", 1:J)

  p <- ncol(X)
  colnames(data) <- c("S", M_variables, "Y", paste0("X", 0:(p-1)))

  # Define the regression formulas
  med.fit.formula <- as.formula(paste0("cbind(", paste(M_variables, collapse = ","), ") ~ S + ", paste0("X", 1:(p-1), collapse = "+")))
  out.fit.formula <- as.formula("Y ~ . - 1")

  n <- nrow(data)
  lambda_n <- lambda * sqrt(n) / log(n)  # Adaptive lambda parameter

  # Helper function for bootstrap resampling
  PoC_AB_help_function <- function(data, ind) {
    # Original data
    S <- data.matrix(data[, "S"])
    M <- data.matrix(data[, M_variables])
    Y <- data.matrix(data[, "Y"])
    X <- data.matrix(data[, (J + 3):ncol(data)])

    # Fit models on the original data
    med.fit <- lm(med.fit.formula, data = data)
    out.fit <- lm(out.fit.formula, data = data)

    eps_M_hat <- resid(med.fit)
    eps_Y_hat <- resid(out.fit)

    S_perp <- perp_function_multi(X, S)
    M_perp <- perp_function_multi(cbind(X, S), M)

    # Estimate standard deviations
    sigma_alpha_S_hat <- sqrt(colMeans(eps_M_hat^2) / mean(S_perp^2))
    sigma_beta_M_hat <- sqrt(mean(eps_Y_hat^2) / colMeans(M_perp^2))

    # Bootstrap resampling
    data_ast <- data[ind,]
    S_ast <- data.matrix(data_ast[, "S"])
    M_ast <- data.matrix(data_ast[, M_variables])
    X_ast <- data.matrix(data_ast[, (J + 3):ncol(data_ast)])

    # Refit models on the bootstrap sample
    med.fit_ast <- lm(med.fit.formula, data = data_ast)
    out.fit_ast <- lm(out.fit.formula, data = data_ast)

    eps_M_hat_ast <- resid(med.fit_ast)
    eps_Y_hat_ast <- resid(out.fit_ast)

    S_perp_ast <- perp_function_multi(X_ast, S_ast)
    M_perp_ast <- perp_function_multi(cbind(X_ast, S_ast), M_ast)

    sigma_alpha_S_hat_ast <- sqrt(colMeans(eps_M_hat_ast^2) / mean(S_perp_ast^2))
    sigma_beta_M_hat_ast <- sqrt(mean(eps_Y_hat_ast^2) / colMeans(M_perp_ast^2))

    Q_S_ast <- Q_function_multi(X_ast, S_ast)
    Q_M_ast <- Q_function_multi(cbind(X_ast, S_ast), M_ast)

    V_S_ast <- mean(S_perp_ast^2)
    V_M_ast <- colMeans(M_perp_ast^2)

    # Calculate residual-based bootstrap statistics
    Z_S_ast <- sqrt(n) * ((t(eps_M_hat[ind,]) %*% (S_ast - X_ast %*% Q_S_ast) / n -
                             (t(eps_M_hat) %*% (S - X %*% Q_S_ast)) / n)) / V_S_ast
    Z_M_ast <- sqrt(n) * t((
      t(eps_Y_hat[ind]) %*% (M_ast - cbind(X_ast, S_ast) %*% Q_M_ast) / n -
        (t(eps_Y_hat) %*% (M - cbind(X, S) %*% Q_M_ast)) / n
    )) / V_M_ast

    # Residual-based interaction
    R_ast <- crossprod(Z_S_ast, Z_M_ast)

    # Compute T statistics
    T_alpha_hat <- sqrt(n) * coef(med.fit)["S", ] / sigma_alpha_S_hat
    T_beta_hat <- sqrt(n) * coef(out.fit)[M_variables] / sigma_beta_M_hat
    T_alpha_hat_ast <- sqrt(n) * coef(med.fit_ast)["S", ] / sigma_alpha_S_hat_ast
    T_beta_hat_ast <- sqrt(n) * coef(out.fit_ast)[M_variables] / sigma_beta_M_hat_ast

    # Indicator functions for adaptive adjustments
    I_alpha_ast <- (max(abs(T_alpha_hat_ast)) <= lambda_n) * (max(abs(T_alpha_hat)) <= lambda_n)
    I_beta_ast <- (max(abs(T_beta_hat_ast)) <= lambda_n) * (max(abs(T_beta_hat)) <= lambda_n)

    # Adjusted bootstrap statistic
    U_ast <- (crossprod(coef(med.fit_ast)["S", ], coef(out.fit_ast)[M_variables]) -
                crossprod(coef(med.fit)["S", ], coef(out.fit)[M_variables])) * (1 - I_alpha_ast * I_beta_ast) +
      n^(-1) * R_ast * I_alpha_ast * I_beta_ast

    return(U_ast)
  }

  # Perform the bootstrap
  b <- boot::boot(data, PoC_AB_help_function, R = B)

  # Estimate the mediation effect from the original data
  med.fit <- lm(med.fit.formula, data = data)
  out.fit <- lm(out.fit.formula, data = data)
  alpha_S_hat <- coef(med.fit)["S", ]
  beta_M_hat <- coef(out.fit)[M_variables]
  NIE_hat <- crossprod(alpha_S_hat, beta_M_hat)

  NDE_hat <- stats::coef(out.fit)['S']
  p_value_NDE <- summary(out.fit)$coefficients['S', 4]

  NTE.fit.formula <- stats::as.formula(paste0("Y~S+", paste(paste("X", 1:(
    p - 1
  ), sep = ""), collapse = "+")))
  NTE.fit <- stats::lm(NTE.fit.formula, data = data)
  NTE_hat <- stats::coef(NTE.fit)['S']
  p_value_NTE <- summary(NTE.fit)$coefficients['S', 4]

  # Calculate the p-value
  p_value <- colMeans(sweep(abs(b$t), 2, abs(NIE_hat), ">"))

  return(list(NIE = NIE_hat, p_value_NIE = p_value,
              NDE = NDE_hat, p_value_NDE = p_value_NDE,
              NTE = NTE_hat, p_value_NTE = p_value_NTE))
}
