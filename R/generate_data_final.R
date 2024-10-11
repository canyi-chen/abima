#' Generate simulated data
#'
#' @param n the sample size
#' @param alpha_S the parameter in M~S
#' @param beta_M the parameter in Y~S+M
#' @param M.family a description of the error distribution and link function to be used in the mediator model. This is the
#' result of a call to a family function. (See \code{\link{family}} for detials of family functions.) The default
#' family is \code{\link{gaussian}}.
#' @param Y.family a description of the error distribution and link function to be used in the outcome model. This is the
#' result of a call to a family function. (See \code{\link{family}} for detials of family functions.) The default
#' family is \code{\link{gaussian}}.
#' @param sigma_M the noise level for M~S
#' @param sigma_Y the noise level for Y~S+M
#'
#' @returns A list with the following components:
#' \item{S}{exposure}
#' \item{M}{mediator}
#' \item{Y}{outcome}
#' \item{X}{confounder}
#' @export
generate_all_data <- function(n = 200,
                              alpha_S = 0,
                              beta_M = 0,
                              M.family = stats::gaussian(),
                              Y.family = stats::gaussian(),
                              sigma_M = 0.5,
                              sigma_Y = 0.5) {
  if ((M.family$family == "gaussian") &
      (Y.family$family == "gaussian")) {
    S <- stats::rbinom(n, 1, 1 / 2)
    X1 <- stats::rnorm(n)
    X2 <- stats::rbinom(n, 1, 1 / 2)
    eps_M <- stats::rnorm(n, 0, sd = sigma_M)
    eps_Y <- stats::rnorm(n, 0, sd = sigma_Y)

    alpha_vec <- c(1, 1, 1)
    beta_vec <- c(1, 1, 1)
    tau_S <- 1
    # generate M
    M <- alpha_S * S +
      alpha_vec[1] +
      alpha_vec[2] * X1 +
      alpha_vec[3] * X2 +
      eps_M
    Y <- beta_M * M +
      beta_vec[1] +
      beta_vec[2] * X1 +
      beta_vec[2] * X2 +
      tau_S * S +
      eps_Y
    X <- cbind(X1, X2)
    return(list(
      S = S,
      M = M,
      Y = Y,
      X = X
    ))
  }
  logit_transform_fun <- M.family$linkinv
  # logit_deriv <- function(x) deriv_inverse.link(M.family$linkinv(x),M.family$link)
  logit_deriv <- function(x)
    M.family$mu.eta(x)

  if (M.family$family == "binomial" &
      Y.family$family == "binomial") {
    v0 <- 1
    alpha_int <- -1
    beta_int <- -1
    alpha_cov <- v0
    beta_cov <- 1
    tauS <- 1

    #1. generate exposure
    S <- stats::rbinom(n, 1, 0.5)
    # S <- rnorm(n)
    X <- stats::rbinom(n, 1, 0.5)

    #2. generate mediator
    mu_vec <- S * alpha_S + alpha_int + alpha_cov * X
    mu_logit <- logit_transform_fun(mu_vec)
    M <- stats::rbinom(n, size = 1, prob = mu_logit)

    #3. generate outcome
    mean_y <- M * beta_M +
      beta_int +
      beta_cov * X +
      tauS * S
    mean_y_logit <- logit_transform_fun(mean_y)
    Y <- stats::rbinom(n, size = 1, prob = mean_y_logit)

    return(list(
      S = S,
      M = M,
      Y = Y,
      X = X
    ))
  }
  if (M.family$family == "binomial" &
      Y.family$family == "gaussian") {
    alpha_int <- -1
    beta_int <- -1
    alpha_cov <- 1
    beta_cov <- 1
    tauS <- 1

    #1. generate exposure
    S <- stats::rbinom(n, 1, 0.5)
    X <- stats::rbinom(n, 1, 0.5)

    #2. generate mediator
    mu_vec <- S * alpha_S + alpha_int + alpha_cov * X
    mu_logit <- logit_transform_fun(mu_vec)
    M <- stats::rbinom(n, size = 1, prob = mu_logit)

    #3. generate outcome
    sigma_sd = 0.5
    mean_y <- M * beta_M + beta_int + beta_cov * X + tauS * S
    Y <- mean_y + stats::rnorm(n, 0, sigma_sd)

    return(list(
      S = S,
      M = M,
      Y = Y,
      X = X
    ))
  }
  if (M.family$family != "binomial" &
      Y.family$family == "gaussian") {
    # alpha_int <- -1; beta_int <- -1
    alpha_int <- 1
    beta_int <- -1
    alpha_cov <- 1
    beta_cov <- 1
    tauS <- 1

    #1. generate exposure
    S <- stats::rbinom(n, 1, 0.5)
    # X <- rbinom(n, 1, 0.5)
    X <- stats::rnorm(n, 0, 2)

    #2. generate mediator
    mu_vec <- S * alpha_S + alpha_int + alpha_cov * X
    mu_logit <- M.family$linkinv(mu_vec) #logit_transform_fun(mu_vec)
    if (M.family$family == "binomial" |
        M.family$family == "quasibinomial") {
      M <- stats::rbinom(n, size = 1, prob =
                    mu_logit)
    } else if (M.family$family == "Gamma") {
      X <- stats::rbinom(n, 1, 0.5)
      mu_vec <- S * alpha_S + alpha_int + alpha_cov * X
      mu_logit <- M.family$linkinv(mu_vec) #logit_transform_fun(mu_vec)
      M <- stats::rgamma(n, shape = 1, scale = mu_logit)
    } else if (M.family$family == "poisson") {
      M <- stats::rpois(n, lambda = mu_logit)
    } else if (M.family$family == "gaussian") {
      M <- stats::rnorm(n, mean = mu_logit)
    }  else {
      stop("Unsupported families.")
    }


    #3. generate outcome
    sigma_sd = 0.5
    mean_y <- M * beta_M + beta_int + beta_cov * X + tauS * S
    Y <- mean_y + stats::rnorm(n, 0, sigma_sd)

    return(list(
      S = S,
      M = M,
      Y = Y,
      X = X
    ))
  }
  stop("Unsupported families.")
}
