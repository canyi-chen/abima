# Fit GLM on the mediator and LM on the outcome, then compute mediation estimate
#
#' @importFrom stats glm
fit_glm_lm_mediation <- function(boot_data, M.family, s = 1, s_star = 0, covariates, covariates_new) {
  # Fit GLM for the mediator
  mediator_reg_formula <- as.formula(paste0("mediator ~ exposure + ", paste(covariates, collapse = "+")))
  glm_alpha_fit <- glm(mediator_reg_formula, family = M.family, data = boot_data)

  # Extract coefficients from GLM fit
  alpha_hat <- coef(glm_alpha_fit)["exposure"]
  alpha_int_hat <- coef(glm_alpha_fit)["(Intercept)"]
  alpha_vec_hat <- coef(glm_alpha_fit)[-c(1, 2)] # Coefficients for covariates

  linkinv_func <- M.family$linkinv
  g_alpha_hat <- linkinv_func(alpha_hat * s + alpha_int_hat + crossprod(covariates_new, alpha_vec_hat))
  g_alpha_hat_0 <- linkinv_func(alpha_hat * s_star + alpha_int_hat + crossprod(covariates_new, alpha_vec_hat))

  d_g_alpha_hat <- g_alpha_hat - g_alpha_hat_0
  alpha_residual <- boot_data$mediator - glm_alpha_fit$fitted.values

  # Fit LM for the outcome
  outcome_reg_formula <- as.formula(paste0("outcome ~ exposure + mediator + ", paste(covariates, collapse = "+")))
  lm_beta_fit <- lm(outcome_reg_formula, data = boot_data)
  beta_hat <- coef(lm_beta_fit)["mediator"]

  beta_residual <- lm_beta_fit$residuals

  # Compute mediation effect estimate
  mediation_estimate <- sqrt(nrow(boot_data)) * beta_hat * d_g_alpha_hat

  return(list(
    mediation_estimate = mediation_estimate,
    alpha_residual = alpha_residual,
    beta_residual = beta_residual,
    glm_alpha_fit = glm_alpha_fit,
    lm_beta_fit = lm_beta_fit
  ))
}

# Compute the local expansion bootstrap statistic
compute_local_stat_GLM_LM <- function(boot_data, M.family, alpha_residual, beta_residual,
                                      glm_alpha_fit, s = 1, s_star = 0, covariates, covariates_new) {

  alpha_hat <- coef(glm_alpha_fit)["exposure"]
  alpha_int_hat <- coef(glm_alpha_fit)["(Intercept)"]
  alpha_vec_hat <- coef(glm_alpha_fit)[-c(1, 2)] # Coefficients for covariates
  g_alpha_fit_values <- glm_alpha_fit$fitted.values
  dispersion_hat <- summary(glm_alpha_fit)$dispersion

  mu.eta_func <- M.family$mu.eta

  W_alpha_t <- mu.eta_func(alpha_hat * s + alpha_int_hat + crossprod(covariates_new, alpha_vec_hat)) * c(s, 1, covariates_new) -
    mu.eta_func(alpha_hat * s_star + alpha_int_hat + crossprod(covariates_new, alpha_vec_hat)) * c(s_star, 1, covariates_new)

  D_boot_mat <- cbind(boot_data$exposure, 1, boot_data[, covariates])
  V_alpha <- t(D_boot_mat) %*% (D_boot_mat / (dispersion_hat * M.family$variance(g_alpha_fit_values) * mu.eta_func(g_alpha_fit_values)^2))

  alpha_stat_res <- W_alpha_t %*% solve(V_alpha) %*% (t(D_boot_mat) %*% diag(1 / (dispersion_hat * M.family$variance(g_alpha_fit_values))) %*% alpha_residual)
  beta_stat_res <- sum(beta_residual) / sum(alpha_residual^2)

  n <- nrow(boot_data)
  residual_stat <- sqrt(n) * alpha_stat_res * beta_stat_res

  return(residual_stat)
}

# Perform adaptive bootstrap on the dataset
adaptive_bootstrap_GLM_LM <- function(my_data, M.family, B_num, test_stat, alpha_residual, beta_residual,
                                      glm_alpha_fit, lm_beta_fit, lambda, s = 1, s_star = 0, covariates, covariates_new) {

  bootstrap_stat <- function(boot_data, boot_index) {
    # Bootstrap data
    boot_sample <- boot_data[boot_index,]

    # Refit models on bootstrap sample
    fit_res <- fit_glm_lm_mediation(boot_sample, M.family, s, s_star, covariates, covariates_new)
    mediation_est_boot <- fit_res$mediation_estimate

    z_alpha <- summary(fit_res$glm_alpha_fit)$coefficients["exposure", "z value"]
    z_beta <- summary(fit_res$lm_beta_fit)$coefficients["mediator", "t value"]

    # Check thresholds
    if (abs(z_alpha) <= lambda && abs(z_beta) <= lambda) {
      ab_stat <- compute_local_stat_GLM_LM(boot_sample, M.family, fit_res$alpha_residual, fit_res$beta_residual,
                                           fit_res$glm_alpha_fit, s, s_star, covariates, covariates_new) - test_stat
    } else {
      ab_stat <- mediation_est_boot
    }

    return(c(mediation_est_boot, ab_stat))
  }

  # Perform bootstrapping using boot package
  boot_results <- boot::boot(my_data, bootstrap_stat, R = B_num)
  return(boot_results)
}

# Main function to calculate mediation effect and p-value using adaptive bootstrap
PoC_AB_GLM_LM <- function(S, M, Y, X, M.family = stats::binomial(link = 'logit'),
                          s = 1, s_star = 0, covariates_new = rep(0, ncol(X) - 1), B = 500, lambda = 2) {

  # Prepare data frame
  data <- data.frame(S, M, Y, X)
  colnames(data) <- c("exposure", "mediator", "outcome", paste0("X", 0:(ncol(X) - 1)))

  # Set covariate names
  covariates <- paste0("X", 1:(ncol(X) - 1))

  # Fit initial models and compute test statistic
  fit_res <- fit_glm_lm_mediation(data, M.family, s, s_star, covariates, covariates_new)
  test_stat <- fit_res$mediation_estimate

  # Perform adaptive bootstrap
  lambda_alpha <- lambda * sqrt(nrow(data)) / log(nrow(data))
  lambda_beta <- lambda * sqrt(nrow(data)) / log(nrow(data))
  bootstrap_results <- adaptive_bootstrap_GLM_LM(data, M.family, B, test_stat, fit_res$alpha_residual,
                                                 fit_res$beta_residual, fit_res$glm_alpha_fit, fit_res$lm_beta_fit,
                                                 lambda_alpha, s, s_star, covariates, covariates_new)

  # Calculate p-value
  p_value <- 2 * min(mean(bootstrap_results$t[, 2] > 0), 1 - mean(bootstrap_results$t[, 2] > 0))

  # Return results as S3 object
  result <- list(mediation_effect = test_stat / sqrt(nrow(data)), p_value = p_value)
  class(result) <- "mediationTestResult"

  return(result)
}

# Print method for the mediationTestResult object
print.mediationTestResult <- function(x) {
  cat("Mediation Effect Estimate: ", x$mediation_effect, "\n")
  cat("P-Value: ", x$p_value, "\n")
}
