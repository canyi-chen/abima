################################################################
#fit alpha logistic regression, beta linear regression, and ME
################################################################
compute_stat_GLM_LM <- function(boot_data, M.family, s = 1, s_star = 0, covariates, covariates_new) {
  ########################################
  #1. fit glm: alpha coefficient
  ########################################
  ########################################
  linkinv_func <- M.family$linkinv

  outcome_reg_formula <- stats::as.formula(paste0("outcome", "~", paste(c("exposure", "mediator", covariates),
                                                                        collapse = "+"), "-1"))
  mediator_reg_formula <- stats::as.formula(paste0("mediator", "~", paste(c("exposure", covariates), collapse = "+"),
                                                   "-1"))

  glm_alpha_fit <-
    stats::glm(mediator_reg_formula,
               family = M.family,
               data = boot_data)
  glm_alpha_summary <- summary(glm_alpha_fit)
  result_alpha <- glm_alpha_summary$coefficients
  alpha_hat <- result_alpha["exposure", "Estimate"]
  # alpha_int_hat <- result_alpha["(Intercept)", "Estimate"]
  z_alpha <- result_alpha["exposure", 3]
  alpha_vec_hat <- result_alpha[-which(colnames(boot_data) == "exposure"), "Estimate"]
  dispersion_hat <- summary(glm_alpha_fit)$dispersion

  g_alpha_hat <- linkinv_func(alpha_hat * s +
                                # alpha_int_hat +
                                crossprod(covariates_new, alpha_vec_hat))
  g_alpha_hat_0 <- linkinv_func(alpha_hat * s_star +
                                  # alpha_int_hat +
                                  crossprod(covariates_new, alpha_vec_hat))
  d_g_alpha_hat <- g_alpha_hat - g_alpha_hat_0
  g_alpha_fit_values <- glm_alpha_fit$fitted.values
  alpha_residual <- boot_data$mediator - glm_alpha_fit$fitted.values

  ########################################
  #2. fit lm: beta coefficient
  ########################################
  lm_beta_fit <- stats::lm(outcome_reg_formula, data = boot_data)
  lm_beta_summary <-
    summary(lm_beta_fit)
  result_beta <- lm_beta_summary$coefficients
  beta_hat <- result_beta["mediator", "Estimate"]
  z_beta <- result_beta["mediator", "t value"]
  beta_residual <- lm_beta_summary$residuals

  # get projected M after bootstrap
  M_proj_boot <-
    stats::lm(mediator_reg_formula, data = boot_data)$residuals

  ########################################
  #3. compute test statistic
  ########################################
  n <- nrow(boot_data)
  mediation_estimate <-
    sqrt(n) * as.numeric(beta_hat * d_g_alpha_hat)

  NDE_hat <- stats::coef(lm_beta_fit)['exposure']*(s - s_star)
  p_value_NDE <- lm_beta_summary$coefficients['exposure', 4]


  return(
    list(
      me = mediation_estimate,
      NDE_hat = NDE_hat,
      p_value_NDE = p_value_NDE,
      NTE_hat = NDE_hat + mediation_estimate/sqrt(n),
      alpha_residual = alpha_residual,
      beta_residual = beta_residual,
      z_alpha = z_alpha,
      z_beta = z_beta,
      M_proj_boot = M_proj_boot,
      g_alpha_fit_values = g_alpha_fit_values,
      alpha_hat = alpha_hat,
      # alpha_int_hat = alpha_int_hat,
      alpha_vec_hat = alpha_vec_hat,
      dispersion_hat = dispersion_hat
    )
  )
}


################################################################
#compute the local bootstrap at (0,0)
################################################################
compute_local_stat_GLM_LM <-
  function(boot_data,
           M.family,
           M_proj_boot,
           alpha_residual,
           beta_residual,
           alpha_hat,
           # alpha_int_hat,
           g_alpha_fit_values,
           dispersion_hat,
           alpha_vec_hat,
           s = 1, s_star = 0, covariates, covariates_new) {

    mu.eta_func <- function(x) M.family$mu.eta(x)

    W_alpha_t <- c(mu.eta_func(alpha_hat * s +
                                 # alpha_int_hat +
                                 crossprod(covariates_new, alpha_vec_hat))) * c(s,
                                                                                # 1,
                                                                                covariates_new) -
      c(mu.eta_func(alpha_hat * s_star +
                      # alpha_int_hat +
                      crossprod(covariates_new, alpha_vec_hat))) * c(s_star,
                                                                     # 1,
                                                                     covariates_new)

    D_boot_mat <-
      as.matrix(cbind(boot_data[, "exposure"], # 1,
                      boot_data[, covariates]))

    V_alpha <-
      t(D_boot_mat) %*% sweep(D_boot_mat,
                              MARGIN = 1,
                              1 / (dispersion_hat * M.family$variance(g_alpha_fit_values)) * M.family$mu.eta(M.family$linkfun(g_alpha_fit_values))^2,
                              '*')

    alpha_stat_res <-
      W_alpha_t %*%
      solve(V_alpha,
            t(D_boot_mat) %*%
              diag(as.numeric(1 / (dispersion_hat * M.family$variance(g_alpha_fit_values)) * M.family$mu.eta(M.family$linkfun(g_alpha_fit_values)))) %*%
              alpha_residual)

    beta_stat_res <-
      sum(M_proj_boot * beta_residual) / sum(M_proj_boot^2)

    n <- length(alpha_residual)
    residual_stat <- sqrt(n) * alpha_stat_res * beta_stat_res
    return(residual_stat)
  }

################################################################
#compute adaptive bootstrap
################################################################
one_ab_bootstrap_GLM_LM <- function(my_data,
                                    M.family,
                                    B_num,
                                    test_stat,
                                    alpha_residual,
                                    beta_residual,
                                    z_alpha,
                                    z_beta,
                                    lambda_alpha,
                                    lambda_beta,
                                    g_alpha_fit_values,
                                    dispersion_hat,
                                    s = 1, s_star = 0, covariates, covariates_new) {

  ab_single <- function(my_data, boot_index) {

    #1. bootstrap indexes and then data
    boot_data <- my_data[boot_index,]
    boot_res_alpha <- alpha_residual[boot_index]
    boot_res_beta <- beta_residual[boot_index]
    boot_g_alpha_fit_values <- g_alpha_fit_values[boot_index]

    #2. compute classical bootstrap statistic
    tmp_boot_res <- compute_stat_GLM_LM(boot_data, M.family = M.family, s = s, s_star = s_star,
                                        covariates = covariates,
                                        covariates_new = covariates_new)
    mediation_est_boot <-
      tmp_boot_res$me #classical bootstrap estimate
    z_alpha_boot <- tmp_boot_res$z_alpha
    z_beta_boot <- tmp_boot_res$z_beta
    M_proj_boot <- tmp_boot_res$M_proj_boot
    # g_alpha_fit_values <- tmp_boot_res$g_alpha_fit_values
    alpha_hat_boot <- tmp_boot_res$alpha_hat
    # alpha_int_hat_boot <- tmp_boot_res$alpha_int_hat
    alpha_vec_hat_boot <- tmp_boot_res$alpha_vec_hat
    dispersion_hat_boot <- tmp_boot_res$dispersion_hat

    #3. depends on passing threshold or not
    t_alpha <-
      (abs(z_alpha) <= lambda_alpha) &
      (abs(z_alpha_boot) <= lambda_alpha)
    t_beta <-
      (abs(z_beta) <= lambda_beta) & (abs(z_beta_boot) <= lambda_beta)

    #4. compute local expansion bootstrap
    if (t_alpha & t_beta) {
      ab_stat <-
        compute_local_stat_GLM_LM(
          boot_data,
          M.family = M.family,
          M_proj_boot,
          boot_res_alpha,
          boot_res_beta,
          alpha_hat_boot,
          # alpha_int_hat_boot,
          boot_g_alpha_fit_values,
          # dispersion_hat = dispersion_hat,
          dispersion_hat = dispersion_hat_boot,
          alpha_vec_hat = alpha_vec_hat_boot, s = s, s_star = s_star, covariates = covariates,
          covariates_new =
            covariates_new
        ) - test_stat
    } else {
      ab_stat <- mediation_est_boot
    }

    return(c(mediation_est_boot, ab_stat, tmp_boot_res$NTE_hat))
  }

  one_boot_res <- boot::boot(my_data, ab_single,
                             R = B_num)

  return(one_boot_res)
}


### Main function
PoC_AB_GLM_LM <- function(S, M, Y, X, M.family = stats::binomial(link = 'logit'),
                          s = 1, s_star = 0, covariates_new = rep(0, ncol(X) - 1), B = 500, lambda = 2) {

  data <- data.frame(
    S = S,
    M = M,
    Y = Y,
    X = X
  )
  p <- ncol(X)
  n <- nrow(X)

  colnames(data) <- c("exposure", "mediator", "outcome", paste(paste("X", 0:(p - 1), sep = "")))
  covariates <- paste("X", 0:(p - 1), sep = "")


  ## AB TEST
  # get indirect effects
  tmp_res <- compute_stat_GLM_LM(data,
                                 M.family = M.family,
                                 s = s,
                                 s_star = s_star,
                                 covariates = covariates,
                                 covariates_new = covariates_new)

  test_stat <- tmp_res$me
  mediation_estimate <- test_stat / sqrt(n)


  alpha_residual <- tmp_res$alpha_residual
  beta_residual <- tmp_res$beta_residual
  z_alpha <- tmp_res$z_alpha
  z_beta <- tmp_res$z_beta
  g_alpha_fit_values <- tmp_res$g_alpha_fit_values
  dispersion_hat <- tmp_res$dispersion_hat

  lambda_alpha <-
    lambda * sqrt(n) / log(n)
  lambda_beta <- lambda * sqrt(n) / log(n)
  one_boot_res <- one_ab_bootstrap_GLM_LM(
    data,
    M.family = M.family,
    B,
    test_stat,
    alpha_residual,
    beta_residual,
    z_alpha,
    z_beta,
    lambda_alpha,
    lambda_beta,
    g_alpha_fit_values,
    dispersion_hat,
    s = s,
    s_star = s_star,
    covariates = covariates,
    covariates_new = covariates_new
  )

  tmp_p_ab <- mean(one_boot_res$t[, 2] > 0)
  p_value <- 2 * min(tmp_p_ab, 1 - tmp_p_ab)

  # ref: https://www.statmodel.com/download/FAQ-Bootstrap%20-%20Pvalue.pdf
  tmp_p_b_NTE <- mean(one_boot_res$t[, 3] > 0)
  p_value_NTE <- 2 * min(tmp_p_b_NTE, 1 - tmp_p_b_NTE)
  return(list(
    NIE = mediation_estimate,
    p_value_NIE = p_value,
    NDE = tmp_res$NDE,
    p_value_NDE = tmp_res$p_value_NDE,
    NTE = tmp_res$NTE,
    p_value_NTE = p_value_NTE
  ))
}
