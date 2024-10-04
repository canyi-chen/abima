################################################################
#fit alpha logistic regression, beta linear regression, and ME
################################################################
compute_stat_GLM_LM <- function(boot_data, M.family, s = 1, s_star = 0, covariates_new = rep(0, length(covariates))) {
  ########################################
  #1. fit logistic: alpha coefficient
  ########################################
  ########################################
  #equivalent closed-form
  # S_boot <- boot_data$exposure
  # M_boot <- boot_data$mediator
  logit_transform_fun <- M.family$linkinv
  logit_deriv <- function(x) M.family$mu.eta(x)
  # logit_deriv <- function(x) deriv_inverse.link(M.family$linkinv(x),M.family$link)

  outcome_reg_formula <- as.formula(paste0("outcome", "~",
                                           paste(c("exposure", "mediator", covariates), collapse = "+")))
  mediator_reg_formula <- as.formula(paste0("mediator", "~",
                                            paste(c("exposure", covariates), collapse = "+")))

  lm_alpha_tmp <-
    glm(mediator_reg_formula,
        family = M.family,
        data = boot_data)
  lm_alpha <- summary(lm_alpha_tmp)
  result_alpha <- lm_alpha$coefficients
  alpha_hat <- result_alpha["exposure", "Estimate"]
  alpha_int_hat <- result_alpha["(Intercept)", "Estimate"]
  # z_alpha <- result_alpha["exposure", "t value"]
  z_alpha <- result_alpha["exposure", 3]
  alpha_vec_hat <- result_alpha[-c(1, 2), "Estimate"]
  dispersion_hat <- summary(lm_alpha_tmp)$dispersion

  # logit_transform_fun <- M.family$linkinv
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
           alpha_int_hat,
           g_alpha_fit_values,
           dispersion_hat,
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

    ## < ccy
    # dispersion_hat <- sqrt(sum(alpha_residual^2)/(length(alpha_residual) - 1 - length(covariates) - 1))
    # dispersion_hat <- 1

    # V_alpha <-
    #   t(D_boot_mat) %*% sweep(D_boot_mat,
    #                           MARGIN = 1,
    #                           #g_alpha_fit_values * (1 - g_alpha_fit_values),
    #                           # 1/dispersion_hat*M.family$variance(g_alpha_fit_values),
    #                           # 1/(dispersion_hat*M.family$variance(g_alpha_fit_values)),
    #                           M.family$variance(g_alpha_fit_values),
    #                           '*')
    #
    # alpha_stat_res <-
    #   W_alpha_t %*%
    #   solve(V_alpha) %*%
    #   t(D_boot_mat) %*%
    #   alpha_residual

    # OPT2
    V_alpha <-
      t(D_boot_mat) %*% sweep(D_boot_mat,
                              MARGIN = 1,
                              #g_alpha_fit_values * (1 - g_alpha_fit_values),
                              # 1/dispersion_hat*M.family$variance(g_alpha_fit_values),
                              1/(dispersion_hat*M.family$variance(g_alpha_fit_values))*M.family$mu.eta(M.family$linkfun(g_alpha_fit_values))^2,
                              # M.family$variance(g_alpha_fit_values),
                              '*')

    alpha_stat_res <-
      W_alpha_t %*%
      solve(V_alpha,
      t(D_boot_mat) %*% diag(as.numeric(1/(dispersion_hat*M.family$variance(g_alpha_fit_values))*M.family$mu.eta(M.family$linkfun(g_alpha_fit_values))))%*%
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
                             s = 1, s_star = 0, covariates_new = rep(0, length(covariates)), ncpus = 7) {
  ab_single <- function (my_data, boot_index) {
    logit_transform_fun <- M.family$linkinv
    logit_deriv <- function(x) M.family$mu.eta(x)
    # logit_deriv <- function(x) deriv_inverse.link(M.family$linkinv(x),M.family$link)

    n <- nrow(my_data)
    #1. bootstrap indexes and then data
    boot_data <- my_data[boot_index,]
    boot_res_alpha <- alpha_residual[boot_index]
    boot_res_beta <- beta_residual[boot_index]
    boot_g_alpha_fit_values <- g_alpha_fit_values[boot_index]

    #2. compute classical bootstrap statistic
    tmp_boot_res <- compute_stat_GLM_LM(boot_data, M.family = M.family, s = s, s_star = s_star, covariates_new = covariates_new)
    mediation_est_boot <-
      tmp_boot_res$me #classical bootstrap estimate
    z_alpha_boot <- tmp_boot_res$z_alpha
    z_beta_boot <- tmp_boot_res$z_beta
    M_proj_boot <- tmp_boot_res$M_proj_boot
    # g_alpha_fit_values <- tmp_boot_res$g_alpha_fit_values
    alpha_hat_boot <- tmp_boot_res$alpha_hat
    alpha_int_hat_boot <- tmp_boot_res$alpha_int_hat
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
      S_boot <- boot_data$exposure
      M_boot <- boot_data$mediator
      ab_stat <-
        compute_local_stat_GLM_LM(
          boot_data,
          M.family = M.family,
          M_proj_boot,
          boot_res_alpha,
          boot_res_beta,
          alpha_hat_boot,
          alpha_int_hat_boot,
          boot_g_alpha_fit_values,
          # dispersion_hat = dispersion_hat,
          dispersion_hat = dispersion_hat_boot,
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


PoC_AB_GLM_LM <- function(S, M, Y, X, M.family = binomial(link = 'logit'), s = 1, s_star = 0, covariates_new = rep(0, ncol(X) - 1), B = 500, lambda = 2) {
  data <- data.frame(
    S = S,
    M = M,
    Y = Y,
    X = X
  )
  p <- ncol(X)
  n <- nrow(X)

  colnames(data) <- c("exposure", "mediator", "outcome", paste(paste("X", 0:(p-1), sep = "")))
  covariates <<- paste("X", 1:(p-1), sep = "")


  ## AB TEST
  # get indirect effects
  tmp_res <- compute_stat_GLM_LM(data,
                          M.family = M.family,
                          s = s,
                          s_star = s_star,
                          covariates_new = covariates_new)

  test_stat <- tmp_res$me
  mediation_estimate <- test_stat/sqrt(n)




  alpha_residual <- tmp_res$alpha_residual
  beta_residual <- tmp_res$beta_residual
  z_alpha <- tmp_res$z_alpha
  z_beta <- tmp_res$z_beta
  g_alpha_fit_values <- tmp_res$g_alpha_fit_values
  dispersion_hat <- tmp_res$dispersion_hat

  # B_num <- 3e4
  # lambda1 <- 2
  # lambda2 <- 2
  B_num <- B
  lambda1 <- lambda2 <- lambda
  lambda_alpha <-
    lambda1 * sqrt(n) / log(n)
  lambda_beta <- lambda2 * sqrt(n) / log(n)
  # set.seed(47)
  one_boot_res <- one_ab_bootstrap_GLM_LM(
    data,
    M.family = M.family,
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



################################################################
#generate data
################################################################
data_generate_model_GLM_LM<-function(n, alpha_coef, beta_coef, M.family = binomial("logit")){

  v0 <- 1
  # alpha_int <- -1; beta_int <- -1
  alpha_int <- 1; beta_int <- -1
  alpha_cov <- 1; beta_cov <- 1
  tauS <- 1

  #1. generate exposure
  S <- rbinom(n, 1, 0.5)
  # X <- rbinom(n, 1, 0.5)
  X <- rnorm(n, 0, 2)

  #2. generate mediator
  mu_vec <- S * alpha_coef + alpha_int + alpha_cov * X
  mu_logit <- M.family$linkinv(mu_vec) #logit_transform_fun(mu_vec)
  if(M.family$family == "binomial"| M.family$family=="quasibinomial")
    M <- rbinom(n, size = 1, prob = mu_logit)
  if(M.family$family == "Gamma")
   M <- rgamma(n, shape = 1, scale = mu_logit)
  if(M.family$family == "poisson")
    M <- rpois(n, lambda = mu_logit)
  if(M.family$family == "gaussian")
    M <- rnorm(n, mean = mu_logit)
  # generate_random_by_family <- function(family, n = 100) {
  #   if (M.family$family == "gaussian") {
  #     return(rnorm(n, mu = mu_logit))
  #   } else if (M.family$family == "binomial") {
  #     return(rbinom(n, size = 1, prob = mu_logit))
  #   } else if (M.family$family == "poisson") {
  #     return(rpois(n, lambda = mu_logit))
  #   } else if (M.family$family == "Gamma") {
  #     return(rgamma(n, shape = 2, rate = 1))
  #   } else if (M.family$family == "inverse.gaussian") {
  #     return(statmod::rinvgauss(n, mean = 1, shape = 2))
  #   } else {
  #     stop("Unsupported family")
  #   }
  # }


  #3. generate outcome
  sigma_sd = 0.5
  mean_y <- M * beta_coef + beta_int + beta_cov * X + tauS * S
  Y <- mean_y + rnorm(n, 0, sigma_sd)

  my_data = data.frame(exposure=S, mediator=M, outcome=Y, covariate1 = X, covariate0 = rep(1,n))

  return(my_data)
}

#
# library(parallel)
# library(abind)
#
# #bootstrap of logistic regression
#
# #data generation
# #0. parameters
# n <- 500
#
# alpha_coef <- 0
# beta_coef <- 0
#
# #repeat simulation
# B_num <- 100 #bootstrap repetition
# n_repeat <- 200
# num_cores <- 8 #number of cores for parallel computing
# # set.seed(123)
# all_seeds = floor(runif(n_repeat) * 10^5)
#
#
# lambda1 <- 1.9; lambda2 <- 1.9
# lambda_alpha <- lambda1*sqrt(n)/log(n); lambda_beta <- lambda2*sqrt(n)/log(n)
#
#
# set.seed(4)
#
# library(future.apply)
# plan(multisession, workers = 12)
#
# # M.family <- gaussian()
# # M.family <- poisson()
# # M.family <- Gamma("log") #Gamma("log")
# # M.family <- quasibinomial(link = "logit")
# M.family <- binomial("logit")
# out <- future_replicate(n_repeat, {
#   my_data <- data_generate_model_GLM_LM(n, alpha_coef, beta_coef, M.family = M.family)
#
#   S <- my_data$exposure
#   M <- my_data$mediator
#   Y <- my_data$outcome
#   X <- my_data[,c("covariate0", "covariate1")]
#
#   PoC_AB_GLM_LM(my_data$exposure, my_data$mediator,my_data$outcome, my_data[,c("covariate0", "covariate1")], M.family = M.family, #binomial("logit"),
#                        lambda = 2, s = 2, s_star = 0, covariates_new = c(2),
#                 B = 200)$p_value
# })
#
# hist(out)
# plot((1:n_repeat)/n_repeat, sort(out), xlim = c(0,1), ylim=c(0,1))
# abline(0,1, col = "red")
# mean(out<=0.05)
# mean(out)
# var(out)
# 1/12


