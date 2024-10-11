#' Adaptive bootstrap for mediation analysis
#'
#' @description
#' \code{ab.mediation} is used for testing no mediation effect via adaptive bootstrap under several (generalized)
#' linear model.
#'
#' @param S a n-by-1 matrix for exposure.
#' @param M a n-by-m matrix for mediator, each row corresponding to an observation. The dimension m could be 1 or
#' larger than 1.
#' @param Y a n-by-1 matrix for outcome.
#' @param X a n-by-p matrix for confounder.
#' @param intercept logical, default at FALSE; if the intercept has been included in the confounder. If the intercept
#' is included, please put it in the first column of X.
#' @param s exposure level, default at 1
#' @param s_star another exposure level, default at 0
#' @param B the number of bootstrap, default at 199
#' @param lambda the constant used in the pretest when conducing adaptive bootstrap, default at 2.
#'
#' @returns mediation_effect the estimated mediation effect between s and s_star (conditioning on the covariates_new).
#' @returns p_value the p value.
#' @references \url{https://academic.oup.com/jrsssb/advance-article/doi/10.1093/jrsssb/qkad129/7420463}
#' @example man/examples/example_abYlm.Mlm.R
#'
#' @export
abYlm.Mlm <- function(S, M, Y, X, intercept = FALSE, s = 1, s_star = 0, B = 199, lambda = 2) {

  # ============================================================ #
  # Parameters checking and cleaning
  # ============================================================ #
  S <- as.matrix(S)
  M <- as.matrix(M)
  Y <- as.matrix(Y)
  X <- as.matrix(X)

  # only support scalar exposure and outcome
  stopifnot((ncol(S) == 1) & ncol(Y) == 1)

  # check intercept
  if (!intercept) {
    X <- cbind(1, X)
  }

  if (ncol(M) > 1) {
      out <- PoC_AB_multi(S = S, M = M, Y = Y, X = X, B = B, lambda = lambda)
      return(
        list(
          mediation_effect = as.numeric(out$mediation_effect) * (s - s_star),
          p_value = out$p_value
        )
      )
  } else {
      out <- PoC_AB(S = S, M = M, Y = Y, X = X, B = B, lambda = lambda)
      return(
        list(
          mediation_effect = as.numeric(out$mediation_effect) * (s - s_star),
          p_value = out$p_value
        )
      )
  }
}