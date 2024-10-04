################################################################
#fit alpha logistic regression, beta linear regression, and ME
################################################################
compute_stat_Binary_LM <- function(boot_data, s = 1, s_star = 0, covariates_new = rep(0, length(covariates))) {
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

  lm_alpha_tmp <-
    glm(mediator_reg_formula,
        family = binomial(link = 'logit'),
        data = boot_data)
  lm_alpha <- summary(lm_alpha_tmp)
  result_alpha <- lm_alpha$coefficients
  alpha_hat <- result_alpha["exposure", "Estimate"]
  alpha_int_hat <- result_alpha["(Intercept)", "Estimate"]
  z_alpha <- result_alpha["exposure", "z value"]
  alpha_vec_hat <- result_alpha[-c(1, 2), "Estimate"]

  g_alpha_hat <- logit_transform_fun(alpha_hat*s + alpha_int_hat + crossprod(covariates_new, alpha_vec_hat))
  g_alpha_hat_0 <- logit_transform_fun(alpha_hat*s_star + alpha_int_hat + crossprod(covariates_new, alpha_vec_hat))
  d_g_alpha_hat <- g_alpha_hat - g_alpha_hat_0
  g_alpha_fit_values <- lm_alpha_tmp$fitted.values

  alpha_residual <- boot_data$mediator - lm_alpha_tmp$fitted.values

  ########################################
  #2. fit linear: beta coefficient
  ########################################
  lm_beta <-
    summary(lm(outcome_reg_formula, data = boot_data))
  result_beta <- lm_beta$coefficients
  beta_hat <- result_beta["mediator", "Estimate"]
  z_beta <- result_beta["mediator", "t value"]
  beta_residual <- lm_beta$residuals

  #get projected M after bootstrap
  M_proj_boot <-
    lm(mediator_reg_formula, data = boot_data)$residuals

  ########################################
  #3. compute test statistic
  ########################################
  n <- nrow(boot_data)
  mediation_estimate <-
    sqrt(n) * as.numeric(beta_hat * d_g_alpha_hat)

  return(
    list(
      me = mediation_estimate,
      alpha_residual = alpha_residual,
      beta_residual = beta_residual,
      z_alpha = z_alpha,
      z_beta = z_beta,
      M_proj_boot = M_proj_boot,
      g_alpha_fit_values = g_alpha_fit_values,
      alpha_hat = alpha_hat,
      alpha_int_hat = alpha_int_hat,
      alpha_vec_hat = alpha_vec_hat
    )
  )
}

################################################################
#compute classical bootstrap
################################################################
# one_bootstrap_Binary_LM <- function(my_data, B_num) {
#   n <- nrow(my_data)
#   one_boot_res <- sapply(1:B_num, function(b, my_data, n) {
#     boot_index <- sample(n, replace = T)
#     boot_data <- my_data[boot_index,]
#     tmp_estiamte <- compute_stat_Binary_LM(boot_data)
#     mediation_estiamte <- tmp_estiamte$me
#     return(mediation_estiamte)
#   }, my_data = my_data, n = n)
#
#   return(one_boot_res)
# }

################################################################
#compute the local bootstrap at (0,0)
################################################################
compute_local_stat_Binary_LM <-
  function(boot_data,
           M_proj_boot,
           alpha_residual,
           beta_residual,
           alpha_hat,
           alpha_int_hat,
           g_alpha_fit_values,
           alpha_vec_hat,
           s = 1, s_star = 0, covariates_new = rep(0, length(covariates))) {

        logit_transform_fun <- M.family$linkinv
    # logit_deriv <- function(x) deriv_inverse.link(M.family$linkinv(x),M.family$link)
    logit_deriv <- function(x) M.family$mu.eta(x)

    W_alpha_t <- c(logit_deriv(alpha_hat*s + alpha_int_hat + crossprod(covariates_new, alpha_vec_hat))) * c(s,
                                                                                                            1,
                                                                                                            covariates_new) -
      c(logit_deriv(alpha_hat*s_star + alpha_int_hat + crossprod(covariates_new, alpha_vec_hat))) * c(s_star,
                                                                                                      1,
                                                                                                      covariates_new)

    D_boot_mat <-
      as.matrix(cbind(boot_data[,"exposure"], 1, boot_data[,covariates]))
    V_alpha <-
      t(D_boot_mat) %*% sweep(D_boot_mat,
                              MARGIN = 1,
                              g_alpha_fit_values * (1 - g_alpha_fit_values),
                              '*')

    alpha_stat_res <-
      W_alpha_t %*%
      solve(V_alpha) %*%
      t(D_boot_mat) %*%
      alpha_residual
    beta_stat_res <-
      sum(M_proj_boot * beta_residual) / sum(M_proj_boot^2)

    n <- length(alpha_residual)
    residual_stat <- sqrt(n) * alpha_stat_res * beta_stat_res
    return(residual_stat)
  }

################################################################
#compute adaptive bootstrap
################################################################
one_ab_bootstrap_Binary_LM <- function(my_data,
                             B_num,
                             test_stat,
                             alpha_residual,
                             beta_residual,
                             z_alpha,
                             z_beta,
                             lambda_alpha,
                             lambda_beta,
                             g_alpha_fit_values,
                             s = 1, s_star = 0, covariates_new = rep(0, length(covariates)), ncpus = 7) {
  ab_single <- function (my_data, boot_index) {
    n <- nrow(my_data)
    #1. bootstrap indexes and then data
    boot_data <- my_data[boot_index,]
    boot_res_alpha <- alpha_residual[boot_index]
    boot_res_beta <- beta_residual[boot_index]
    boot_g_alpha_fit_values <- g_alpha_fit_values[boot_index]

    #2. compute classical bootstrap statistic
    tmp_boot_res <- compute_stat_Binary_LM(boot_data, s = s, s_star = s_star, covariates_new = covariates_new)
    mediation_est_boot <-
      tmp_boot_res$me #classical bootstrap estimate
    z_alpha_boot <- tmp_boot_res$z_alpha
    z_beta_boot <- tmp_boot_res$z_beta
    M_proj_boot <- tmp_boot_res$M_proj_boot
    # g_alpha_fit_values <- tmp_boot_res$g_alpha_fit_values
    alpha_hat_boot <- tmp_boot_res$alpha_hat
    alpha_int_hat_boot <- tmp_boot_res$alpha_int_hat
    alpha_vec_hat_boot <- tmp_boot_res$alpha_vec_hat

    #3. depends on passing threshold or not
    t_alpha <-
      (abs(z_alpha) <= lambda_alpha) &
      (abs(z_alpha_boot) <= lambda_alpha)
    t_beta <-
      (abs(z_beta) <= lambda_beta) & (abs(z_beta_boot) <= lambda_beta)

    #4. compute local expansion bootstrap
    if (t_alpha & t_beta) {
      S_boot <- boot_data$exposure
      M_boot <- boot_data$mediator
      ab_stat <-
        compute_local_stat_Binary_LM(
          boot_data,
          M_proj_boot,
          boot_res_alpha,
          boot_res_beta,
          alpha_hat_boot,
          alpha_int_hat_boot,
          boot_g_alpha_fit_values,
          alpha_vec_hat = alpha_vec_hat_boot, s = s, s_star = s_star, covariates_new = covariates_new
        ) - test_stat
    } else {
      ab_stat <- mediation_est_boot
    }

    return(c(mediation_est_boot, ab_stat))
  }

  one_boot_res <- boot::boot(my_data, ab_single,
                       # weights = weights,
                       # strata = strata,
                       R = B_num, parallel = "multicore", ncpus=ncpus)

  return(one_boot_res)
}



### Main function


PoC_AB_Binary_LM <- function(S, M, Y, X, s = 1, s_star = 0, covariates_new = rep(0, ncol(X) - 1), B = 500, lambda = 2) {
  data <- data.frame(
    S = S,
    M = M,
    Y = Y,
    X = X
  )
  p <- ncol(X)
  n <- nrow(X)

  colnames(data) <- c("exposure", "mediator", "outcome", paste(paste("X", 0:(p-1), sep = "")))
  covariates <<- paste("X", 1:(p-1), sep = "") # make it global


  ## AB TEST
  # get indirect effects
  tmp_res <- compute_stat_Binary_LM(data, s = s,
                          s_star = s_star,
                          covariates_new = covariates_new)

  test_stat <- tmp_res$me
  mediation_estimate <- test_stat/sqrt(n)




  alpha_residual <- tmp_res$alpha_residual
  beta_residual <- tmp_res$beta_residual
  z_alpha <- tmp_res$z_alpha
  z_beta <- tmp_res$z_beta
  g_alpha_fit_values <- tmp_res$g_alpha_fit_values


  lambda1 <- lambda2 <- lambda
  lambda_alpha <-
    lambda1 * sqrt(n) / log(n)
  lambda_beta <- lambda2 * sqrt(n) / log(n)
  # set.seed(47)
  one_boot_res <- one_ab_bootstrap_Binary_LM(
    data,
    B,
    test_stat,
    alpha_residual,
    beta_residual,
    z_alpha,
    z_beta,
    lambda_alpha,
    lambda_beta,
    g_alpha_fit_values,
    s = s,
    s_star = s_star,
    covariates_new = covariates_new,
    ncpus = parallel::detectCores() - 1
  )

  tmp_p_ab <- mean(one_boot_res$t[, 2] > 0)
  p_value <- 2 * min(tmp_p_ab, 1 - tmp_p_ab)
  return(list(
    mediation_effect = mediation_estimate,
    p_value=p_value
  ))
}
