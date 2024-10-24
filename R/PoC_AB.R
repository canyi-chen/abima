# Function to compute coefficient matrix Q using the least squares approach
Q_function <- function(X, Y) {
  stopifnot(nrow(X) == nrow(Y))
  return(solve(t(X) %*% X, t(X) %*% Y))
}

# Function to compute the residual (perpendicular) component of Y on X
perp_function <- function(X, Y) {
  return(Y - X %*% Q_function(X, Y))
}

# Main PoC adaptive bootstrap function for mediation analysis
PoC_AB <- function(S, M, Y, X, B = 500, lambda = 2) {
  data <- data.frame(S = S,
                     M = M,
                     Y = Y,
                     X = X)
  n <- nrow(data)
  p <- ncol(X)

  colnames(data) <- c("S", "M", "Y", paste(paste("X", 0:(p - 1), sep = "")))

  # Define the regression formulas
  med.fit.formula <- stats::as.formula(paste0("M~S+", paste(paste("X", 1:(
    p - 1
  ), sep = ""), collapse = "+")))
  out.fit.formula <- stats::as.formula("Y~.-1")


  lambda_n <- lambda * sqrt(n) / log(n)  # Adaptive lambda parameter

  # Helper function for bootstrap resampling
  PoC_AB_help_function <- function(data, ind) {
    n <- nrow(data)
    # Original data
    S <- data.matrix(data[, 1])
    M <- data.matrix(data[, 2])
    X <- data.matrix(data[, 4:ncol(data)])
    X_tilde <- cbind(X, S)

    # Fit models on the original data
    med.fit <- stats::lm(med.fit.formula, data = data)
    out.fit <- stats::lm(out.fit.formula, data = data)

    eps_M_hat <- matrix(stats::resid(med.fit))
    eps_Y_hat <- matrix(stats::resid(out.fit))

    S_perp <- perp_function(X, S)
    M_perp <- perp_function(cbind(X, S), M)

    # Estimate standard deviations
    sigma_alpha_S_hat <- sqrt(mean(eps_M_hat ^ 2) / mean(S_perp ^ 2))
    sigma_beta_M_hat <- sqrt(mean(eps_Y_hat ^ 2) / mean(M_perp ^ 2))


    # Bootstrap resampling
    data_ast <- data[ind, ]
    S_ast <- data.matrix(data_ast[, 1])
    M_ast <- data.matrix(data_ast[, 2])
    X_ast <- data.matrix(data_ast[, 4:ncol(data_ast)])
    X_tilde_ast <- cbind(X_ast, S_ast)

    # Refit models on the bootstrap sample
    med.fit_ast <- stats::lm(med.fit.formula, data = data_ast)
    out.fit_ast <- stats::lm(out.fit.formula, data = data_ast)

    eps_M_hat_ast <- matrix(stats::resid(med.fit_ast))
    eps_Y_hat_ast <- matrix(stats::resid(out.fit_ast))

    S_perp_ast <- perp_function(X_ast, S_ast)
    M_perp_ast <- perp_function(cbind(X_ast, S_ast), M_ast)

    sigma_alpha_S_hat_ast <- sqrt(mean(eps_M_hat_ast ^ 2) / mean(S_perp_ast ^
                                                                   2))
    sigma_beta_M_hat_ast <- sqrt(mean(eps_Y_hat_ast ^ 2) / mean(M_perp_ast ^
                                                                  2))

    Q_S_ast <- Q_function(X_ast, S_ast)
    Q_M_ast <- Q_function(cbind(X_ast, S_ast), M_ast)

    V_S_ast <- mean(S_perp_ast ^ 2)
    V_M_ast <- mean(M_perp_ast ^ 2)

    # Calculate residual-based bootstrap statistics
    Z_S_ast <- sqrt(n) * (mean(eps_M_hat[ind] * (S_ast - X_ast %*% Q_S_ast) -
                                 mean(eps_M_hat * (S - X %*% Q_S_ast)))) / V_S_ast
    Z_M_ast <- sqrt(n) * (mean(
      eps_Y_hat[ind] * (M_ast - X_tilde_ast %*% Q_M_ast) -
        mean(eps_Y_hat * (M - X_tilde %*% Q_M_ast))
    )) / V_M_ast

    # Residual-based interaction
    R_ast <- Z_S_ast * Z_M_ast

    # Compute T statistics
    T_alpha_hat <- sqrt(n) * stats::coef(med.fit)['S'] / sigma_alpha_S_hat
    T_beta_hat <- sqrt(n) * stats::coef(out.fit)['M'] / sigma_beta_M_hat
    T_alpha_hat_ast <- sqrt(n) * stats::coef(med.fit_ast)['S'] / sigma_alpha_S_hat_ast
    T_beta_hat_ast <- sqrt(n) * stats::coef(out.fit_ast)['M'] / sigma_beta_M_hat_ast


    # Indicator functions for adaptive adjustments
    I_alpha_ast <- (abs(T_alpha_hat_ast) <= lambda_n) * (abs(T_alpha_hat) <= lambda_n)
    I_beta_ast <- (abs(T_beta_hat_ast) <= lambda_n) * (abs(T_beta_hat) <= lambda_n)

    # Adjusted bootstrap statistic
    U_ast <- (
      stats::coef(med.fit_ast)['S'] * stats::coef(out.fit_ast)['M'] -
        stats::coef(med.fit)['S'] * stats::coef(out.fit)['M']
    ) * (1 - I_alpha_ast * I_beta_ast) +
      n ^ (-1) * R_ast * I_alpha_ast * I_beta_ast
    return(U_ast)
  }

  # Perform the bootstrap
  b <- boot::boot(data, PoC_AB_help_function, R = B)

  # Estimate the mediation effect from the original data
  med.fit <- stats::lm(med.fit.formula, data = data)
  out.fit <- stats::lm(out.fit.formula, data = data)
  alpha_S_hat <- stats::coef(med.fit)['S']
  beta_M_hat <- stats::coef(out.fit)['M']
  NIE_hat <- alpha_S_hat * beta_M_hat

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
