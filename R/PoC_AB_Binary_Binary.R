################################################################
#fit alpha logistic regression, beta linear regression, and ME
################################################################
compute_stat_Binary_Binary <- function(boot_data, M.family, s = 1, s_star = 0, covariates_new = rep(0, length
(covariates))) {
  logit_transform_fun <- M.family$linkinv
  # logit_deriv <- function(x) deriv_inverse.link(M.family$linkinv(x),M.family$link)
  logit_deriv <- function(x) M.family$mu.eta(x)

  ########################################
  #1. fit logistic: alpha coefficient
  ########################################
  ########################################
  #equivalent closed-form
  # S_boot <- boot_data$exposure
  # M_boot <- boot_data$mediator
  outcome_reg_formula <- as.formula(paste0("outcome", "~",
                                           paste(c("exposure", "mediator", covariates), collapse = "+")))
  mediator_reg_formula <- as.formula(paste0("mediator", "~",
                                            paste(c("exposure", covariates), collapse = "+")))

  lm_alpha_tmp <- glm(mediator_reg_formula, family = M.family, data = boot_data) #glm(mediator ~ exposure + covariate1, family=binomial(link='logit'),data=boot_data)
  lm_alpha <- summary(lm_alpha_tmp)
  result_alpha <- lm_alpha$coefficients
  alpha_hat <- result_alpha["exposure", "Estimate"]
  alpha_int_hat <- result_alpha["(Intercept)", "Estimate"]
  z_alpha <- result_alpha["exposure", "z value"]
  alpha_vec_hat <- result_alpha[-c(1, 2), "Estimate"] # added by ccy

  # fit_cov_1 <- as.matrix(my_data[, c("covariate0","exposure","covariate1")])
  # g_alpha_hat <- logit_transform_fun( alpha_hat + alpha_int_hat)
  g_alpha_hat <- logit_transform_fun(alpha_hat * s +
                                       alpha_int_hat +
                                       crossprod(covariates_new, alpha_vec_hat))
  # g_alpha_hat_0 <- logit_transform_fun( alpha_int_hat)
  g_alpha_hat_0 <- logit_transform_fun(alpha_hat * s_star +
                                         alpha_int_hat +
                                         crossprod(covariates_new, alpha_vec_hat))
  d_g_alpha_hat <- g_alpha_hat - g_alpha_hat_0
  g_alpha_fit_values <- lm_alpha_tmp$fitted.values

  alpha_residual <- boot_data$mediator - g_alpha_fit_values #lm_alpha_tmp$fitted.values

  ########################################
  #2. fit linear: beta coefficient
  ########################################
  lm_beta_tmp <- glm(outcome_reg_formula, family = M.family, data = boot_data) #glm(outcome ~ mediator + exposure + covariate1, family=binomial(link='logit'),data=boot_data)
  lm_beta <- summary(lm_beta_tmp)
  result_beta <- lm_beta$coefficients
  beta_hat <- result_beta["mediator", "Estimate"]
  beta_int_hat <- result_beta["(Intercept)", "Estimate"]
  tau_hat <- result_beta["exposure", "Estimate"]
  z_beta <- result_beta["mediator", "z value"]
  beta_vec_hat <- result_beta[-c(1, 2, 3), "Estimate"] # added by ccy
  beta_residual <- lm_beta$residuals

  # g_beta_hat <- logit_transform_fun(beta_hat+beta_int_hat)
  g_beta_hat <- logit_transform_fun(beta_hat +
                                      beta_int_hat +
                                      crossprod(covariates_new, beta_vec_hat) +
                                      tau_hat * s)
  # g_beta_hat_0 <- logit_transform_fun(beta_int_hat)
  g_beta_hat_0 <- logit_transform_fun(beta_int_hat +
                                        crossprod(covariates_new, beta_vec_hat) +
                                        tau_hat * s)
  d_g_beta_hat <- g_beta_hat - g_beta_hat_0
  g_beta_fit_values <- lm_beta_tmp$fitted.values
  beta_residual <- boot_data$outcome - g_beta_fit_values

  ########################################
  #3. compute test statistic
  ########################################
  P0_hat <- d_g_beta_hat * g_alpha_hat_0 + g_beta_hat_0
  P1_hat <- d_g_beta_hat * g_alpha_hat + g_beta_hat_0

  n <- nrow(boot_data)

  logit_inverse_fun <- M.family$linkfun
  mediation_estimate <- sqrt(n) * (logit_inverse_fun(P0_hat) - logit_inverse_fun(P1_hat))
  # as.numeric(beta_hat * d_g_alpha_hat)

  return(list(me = mediation_estimate, alpha_residual = alpha_residual,
              beta_residual = beta_residual, z_alpha = z_alpha, z_beta = z_beta,
              g_alpha_fit_values = g_alpha_fit_values,
              g_beta_fit_values = g_beta_fit_values,
              alpha_hat = alpha_hat, alpha_int_hat = alpha_int_hat,
              alpha_vec_hat = alpha_vec_hat,
              beta_hat = beta_hat, beta_int_hat = beta_int_hat,
              beta_vec_hat = beta_vec_hat,
              tau_hat = tau_hat,
              P0_hat = P0_hat))
}


################################################################
#compute the local bootstrap at (0,0)
################################################################
compute_local_stat_Binary_Binary <- function(boot_data, M.family, alpha_residual, beta_residual,
                                             alpha_hat, alpha_int_hat, g_alpha_fit_values,
                                             beta_hat, beta_int_hat, g_beta_fit_values,
                                             P0_hat_boot,
                                             alpha_vec_hat,
                                             beta_vec_hat, tau_hat,
                                             s = 1, s_star = 0, covariates_new = rep(0, length(covariates))) {
  logit_transform_fun <- M.family$linkinv
  # logit_deriv <- function(x) deriv_inverse.link(M.family$linkinv(x),M.family$link)
  logit_deriv <- function(x) M.family$mu.eta(x)

  # W_alpha_t <- c( logit_deriv(alpha_hat + alpha_int_hat),  logit_deriv(alpha_hat + alpha_int_hat)-logit_deriv(alpha_int_hat), 0 )
  # W_alpha_t <- c( logit_deriv(alpha_int_hat),  0, 0 )

  W_alpha_t <- c(logit_deriv(alpha_hat * s +
                               alpha_int_hat +
                               crossprod(covariates_new, alpha_vec_hat))) * c(s,
                                                                              1,
                                                                              covariates_new) -
    c(logit_deriv(alpha_hat * s_star +
                    alpha_int_hat +
                    crossprod(covariates_new, alpha_vec_hat))) * c(s_star,
                                                                   1,
                                                                   covariates_new) # ccy
  # W_alpha_t <- c( logit_deriv(alpha_hat + alpha_int_hat),  0, 0 )

  # D_boot_mat <- as.matrix(boot_data[, c("exposure", "covariate0", "covariate1")])
  D_boot_mat <-
    as.matrix(cbind(boot_data[, "exposure"], 1, boot_data[, covariates]))
  V_alpha <- t(D_boot_mat) %*% sweep(D_boot_mat, MARGIN = 1, g_alpha_fit_values * (1 - g_alpha_fit_values), '*')
  # same as #V_alpha <- t(D_boot_mat) %*% diag(g_alpha_fit_values * (1-g_alpha_fit_values) ) %*% D_boot_mat
  alpha_stat_res <- W_alpha_t %*%
    solve(V_alpha) %*%
    t(D_boot_mat) %*%
    alpha_residual

  # W_beta_t <- c(logit_deriv(beta_int_hat), 0, 0, 0)
  W_beta_t <- c(logit_deriv(beta_hat +
                              beta_int_hat +
                              crossprod(covariates_new, beta_vec_hat) +
                              tau_hat * s)) * c(1,
                                                1,
                                                covariates_new, s) -
    c(logit_deriv(beta_int_hat +
                    crossprod(covariates_new, beta_vec_hat) +
                    tau_hat * s)) * c(0,
                                      1,
                                      covariates_new, s) # ccy

    #   W_beta_t <- c(logit_deriv(beta_hat +
    #                           beta_int_hat +
    #                           crossprod(covariates_new, beta_vec_hat) +
    #                           tau_hat * s)) * c(1,
    #                                             1,
    #                                             covariates_new, 0) -
    # c(logit_deriv(beta_int_hat +
    #                 crossprod(covariates_new, beta_vec_hat) +
    #                 tau_hat * s)) * c(0,
    #                                   1,
    #                                   covariates_new, 0) # ccy

  D_boot_mat_beta <- as.matrix(cbind(boot_data[, "mediator"], 1, boot_data[, covariates], boot_data[, "exposure"])) #as.matrix( boot_data[, c("mediator","covariate0", "covariate1", "exposure")] )
  V_beta <- t(D_boot_mat_beta) %*% sweep(D_boot_mat_beta, MARGIN = 1, g_beta_fit_values * (1 - g_beta_fit_values), '*')
  beta_stat_res <- W_beta_t %*%
    solve(V_beta) %*%
    t(D_boot_mat_beta) %*%
    beta_residual

  l_prime_P_0 <- 1 / (P0_hat_boot * (1 - P0_hat_boot))
  n <- length(alpha_residual)
  residual_stat <- sqrt(n) *
    alpha_stat_res *
    beta_stat_res *
    l_prime_P_0
  return(residual_stat)
}

################################################################
#compute adaptive bootstrap
################################################################
one_ab_bootstrap_Binary_Binary <- function(my_data, M.family,
                                           B_num, test_stat,
                                           alpha_residual, beta_residual,
                                           z_alpha, z_beta,
                                           lambda_alpha, lambda_beta,
                                           g_alpha_fit_values, g_beta_fit_values,
                                           s = 1, s_star = 0, covariates_new = rep(0, length(covariates))) {
  n <- nrow(my_data)
  one_boot_res <- sapply(1:B_num, function(b, my_data, n) {
    #1. bootstrap indexes and then data
    boot_index <- sample(n, replace = T)
    boot_data <- my_data[boot_index,]
    boot_res_alpha <- alpha_residual[boot_index]
    boot_res_beta <- beta_residual[boot_index]
    boot_g_alpha_fit_values <- g_alpha_fit_values[boot_index]
    boot_g_beta_fit_values <- g_beta_fit_values[boot_index]

    #2. compute classical bootstrap statistic
    tmp_boot_res <- compute_stat_Binary_Binary(boot_data, M.family = M.family,
                                               s = s, s_star = s_star, covariates_new =
                                                 covariates_new)
    mediation_est_boot <- tmp_boot_res$me #classical bootstrap estimate
    z_alpha_boot <- tmp_boot_res$z_alpha
    z_beta_boot <- tmp_boot_res$z_beta
    # M_proj_boot <- tmp_boot_res$M_proj_boot
    # g_alpha_fit_values <- tmp_boot_res$g_alpha_fit_values
    alpha_hat_boot <- tmp_boot_res$alpha_hat
    alpha_int_hat_boot <- tmp_boot_res$alpha_int_hat
    alpha_vec_hat_boot <- tmp_boot_res$alpha_vec_hat

    beta_hat_boot <- tmp_boot_res$beta_hat
    beta_int_hat_boot <- tmp_boot_res$beta_int_hat
    beta_vec_hat_boot <- tmp_boot_res$beta_vec_hat
    tau_hat_boot <- tmp_boot_res$tau_hat

    P0_hat_boot <- tmp_boot_res$P0_hat

    #3. depends on passing threshold or not
    t_alpha <- (abs(z_alpha) <= lambda_alpha) & (abs(z_alpha_boot) <= lambda_alpha)
    t_beta <- (abs(z_beta) <= lambda_beta) & (abs(z_beta_boot) <= lambda_beta)

    #4. compute local expansion bootstrap
    if (t_alpha & t_beta) {
      # S_boot <- boot_data$exposure; M_boot <- boot_data$mediator
      # ab_stat <- compute_local_stat_Binary_Binary(boot_data, M_proj_boot, alpha_residual, beta_residual,
      #                               alpha_hat_boot, alpha_int_hat_boot,
      #                               boot_g_alpha_fit_values) - test_stat
      ab_stat <- compute_local_stat_Binary_Binary(boot_data, M.family = M.family,
                                                  boot_res_alpha, boot_res_beta,
                                                  alpha_hat_boot, alpha_int_hat_boot,
                                                  boot_g_alpha_fit_values,
                                                  beta_hat_boot, beta_int_hat_boot,
                                                  boot_g_beta_fit_values,
                                                  P0_hat_boot,
                                                  alpha_vec_hat_boot,
                                                  beta_vec_hat_boot, tau_hat_boot,
                                                  s = s, s_star = s_star, covariates_new = covariates_new) - test_stat
    }else {
      # ab_stat <- mediation_est_boot - test_stat
      ab_stat <- mediation_est_boot
    }

    return(c(mediation_est_boot, ab_stat))
  }, my_data = my_data, n = n)

  return(one_boot_res)
}

#' @export
PoC_AB_Binary_Binary <- function(S, M, Y, X, M.family = binomial(), s = 1, s_star = 0, covariates_new = rep(0, ncol
(X) - 1), B = 500, lambda = 2) {
  my_data <- data.frame(
    S = S,
    M = M,
    Y = Y,
    X = X
  )
  p <- ncol(X)
  n <- nrow(X)

  lambda_alpha <- lambda_beta <- lambda * sqrt(n) / log(n)

  colnames(my_data) <- c("exposure", "mediator", "outcome", paste(paste("X", 0:(p - 1), sep = ""))) # <<<-
  covariates <<- paste("X", 1:(p - 1), sep = "") # make it global


  ## AB TEST
  tmp_res <- compute_stat_Binary_Binary(my_data, M.family = M.family,
                                        s = s,
                                        s_star = s_star,
                                        covariates_new = covariates_new)

  test_stat <- tmp_res$me
  alpha_residual <- tmp_res$alpha_residual
  beta_residual <- tmp_res$beta_residual
  z_alpha <- tmp_res$z_alpha
  z_beta <- tmp_res$z_beta
  g_alpha_fit_values <- tmp_res$g_alpha_fit_values
  g_beta_fit_values <- tmp_res$g_beta_fit_values


  one_boot_res <- one_ab_bootstrap_Binary_Binary(my_data, M.family = M.family,
                                                 B, test_stat,
                                                 alpha_residual, beta_residual,
                                                 z_alpha, z_beta,
                                                 lambda_alpha, lambda_beta,
                                                 g_alpha_fit_values, g_beta_fit_values,
                                                 s = s, s_star = s_star, covariates_new = covariates_new)
  #Result 1: classical bootstrap
  tmp_p_class <- mean(one_boot_res[1,] > 0)
  p_val_classical_boot <- 2 * min(tmp_p_class, 1 - tmp_p_class)

  #Result 2: adaptive bootstrap
  # tmp_p_ab <- mean( one_boot_res[2, ] > test_stat)
  tmp_p_ab <- mean(one_boot_res[2,] > 0)
  p_val_ab_boot <- 2 * min(tmp_p_ab, 1 - tmp_p_ab)
  return(list(
    mediation_effect = test_stat / sqrt(n),
    p_value = p_val_ab_boot
  ))
}

################################################################
#generate data
################################################################
data_generate_model_binary_binary <- function(n, alpha_coef, beta_coef) {

  v0 <- 1
  alpha_int <- -1; beta_int <- -1
  alpha_cov <- v0; beta_cov <- 1
  tauS <- 1

  #1. generate exposure
  S <- rbinom(n, 1, 0.5)
  # S <- rnorm(n)
  X <- rbinom(n, 1, 0.5)

  #2. generate mediator
  mu_vec <- S * alpha_coef + alpha_int + alpha_cov * X
  logit_transform_fun <- binomial()$linkinv
  mu_logit <- logit_transform_fun(mu_vec)
  M <- rbinom(n, size = 1, prob = mu_logit)

  #3. generate outcome
  mean_y <- M * beta_coef +
    beta_int +
    beta_cov * X +
    tauS * S
  mean_y_logit <- logit_transform_fun(mean_y)
  Y <- rbinom(n, size = 1, prob = mean_y_logit)

  my_data = data.frame(exposure = S, mediator = M, outcome = Y, covariate1 = X, covariate0 = rep(1, n))

  return(my_data)
}


