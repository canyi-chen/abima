#' Generate simulated data
#'
#'
#' @importFrom stats rnorm
#'
#' @param n the sample size
#' @param alpha_S the parameter in M~S+X for S
#' @param beta_M the parameter in Y~S+M+X for M
#' @param tau_S the parameter in Y~S+M+X for S
#' @param alpha_vec the parameter in M~S+X for X
#' @param beta_vec the parameter in Y~S+M+X for X
#' @param M.family a description of the error distribution and link function to be used in the mediator model. This is the
#' result of a call to a family function. (See \code{\link{family}} for details of family functions.) The default
#' family is \code{\link{gaussian}}.
#' @param Y.family a description of the error distribution and link function to be used in the outcome model. This is the
#' result of a call to a family function. (See \code{\link{family}} for details of family functions.) The default
#' family is \code{\link{gaussian}}.
#' @param sigma_M the noise level for M~S+X
#' @param sigma_Y the noise level for Y~S+M+X
#'
#' @examples
#' out <- generate_all_data(n = 200)
#' summary(out$S)
#' summary(out$M)
#' summary(out$Y)
#' summary(out$X)
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
                              tau_S = 1,
                              alpha_vec = rep(1, 3),
                              beta_vec = rep(1, 3),
                              M.family = stats::gaussian(),
                              Y.family = stats::gaussian(),
                              sigma_M = 0.5,
                              sigma_Y = 0.5) {
  stopifnot(
    "The length of alpha_vec should equal that of beta_vec" = length(alpha_vec) ==
      length(beta_vec)
  )
  stopifnot("The length of alpha_vec should equal that of beta_vec" = length(alpha_S) ==
              length(beta_M))
  p <- length(alpha_vec)
  J <- length(alpha_S)
  linkinv <- M.family$linkinv

  if (Y.family$family == "gaussian") {
    alpha_vec <- as.matrix(alpha_vec)
    beta_vec <- as.matrix(beta_vec)
    alpha_S <- as.matrix(alpha_S)
    beta_M <- as.matrix(beta_M)

    S <- matrix(stats::rbinom(n, 1, 0.5), n, 1)
    X <- cbind(1, matrix(rnorm((p - 1) * n, 0, 1), n, p - 1))


    eps_M <- matrix(stats::rnorm(n * J, 0, sigma_M), n, J)
    eps_Y <- stats::rnorm(n, 0, sigma_Y)

    mu_glm <- M.family$linkinv(S %*% t(alpha_S) + X %*% alpha_vec %*% matrix(1, 1, J))

    if (M.family$family == "gaussian") {
      M <- S %*% t(alpha_S) + X %*% alpha_vec %*% matrix(1, 1, J) + eps_M
    } else if (M.family$family == "Gamma") {
      M <- stats::rgamma(n * J, shape = 1, scale = mu_glm)
      M <- matrix(M, n, J)
    } else if (M.family$family == "poisson") {
      M <- stats::rpois(n * J, lambda = mu_glm)
      M <- matrix(M, n, J)
    } else if (M.family$family == "binomial") {
      M <- stats::rbinom(n * J, 1, prob = mu_glm)
      M <- matrix(M, n, J)
    } else {
      stop("We only support Gaussian, Gamma, Poisson and Binomial family for the mediator M.")
    }


    Y <-   M %*% beta_M +
      X %*% beta_vec +
      tau_S * S +
      eps_Y


    return(structure(list(
      S = S,
      M = M,
      Y = Y,
      X = X[, -1, drop = F]
    ), class = "sim_data"))
  } else {
    stop("We only support Gaussian family for the outcome Y.")
  }

}
