#' Adaptive Bootstrap for Mediation Analysis with Linear Models
#'
#' @description
#' \code{abYlm.Mlm} conducts adaptive bootstrap for mediation analysis with linear models on both the mediator M and the outcome Y.
#'
#' @param S an n-by-1 vector for exposure.
#' @param M an n-by-m matrix for mediator, each row corresponding to an observation. The dimension m could be 1 or
#' larger than 1.
#' @param Y an n-by-1 vector for outcome.
#' @param X an n-by-p matrix for confounder. Please do not include intercept in X. If you do not have confounder, you do not need specify this parameter.
#' @param s exposure level, default is 1
#' @param s_star another exposure level, default is 0
#' @param B the number of bootstrap samples, default is 199
#' @param lambda the constant used in the pretest when conducting adaptive bootstrap, default is 2.
#'
#' @returns NIE the estimated mediation effect between s and s_star.
#' @returns p_value_NIE the p value for NIE.
#' @returns NDE the estimated natural direct effect between s and s_star.
#' @returns p_value_NDE the p value for NDE.
#' @returns NTE the estimated natural total (treatment) effect between s and s_star.
#' @returns p_value_NTE the p value for NTE.
#' @references He, Y., Song, P. X. K., and Xu, G. (2023), “Adaptive bootstrap tests for composite null hypotheses in the mediation pathway analysis,” Journal of the Royal Statistical Society Series B: Statistical Methodology, qkad129. <doi:10.1093/jrsssb/qkad129>.
#' @example man/examples/example_abYlm.Mlm.R
#'
#' @export
abYlm.Mlm <- function(S, M, Y, X = NULL, s = 1, s_star = 0, B = 199, lambda = 2) {

  # ============================================================ #
  # Parameters checking and cleaning
  # ============================================================ #
  S <- as.matrix(S)
  M <- as.matrix(M)
  Y <- as.matrix(Y)
  # X <- as.matrix(X)

  # Check dimensions to ensure matrices have the same number of rows
  stopifnot(nrow(S) == nrow(M), nrow(M) == nrow(Y))

  # Ensure scalar exposure and outcome
  stopifnot(ncol(S) == 1, ncol(Y) == 1)

  if(is.null(X)) {
    X <- matrix(1, nrow(S), 1)
  } else {
    # Check dimensions to ensure matrices have the same number of rows
    stopifnot(nrow(Y) == nrow(X))
    X <- cbind(1, as.matrix(X))  # Add intercept column if it is not present
  }

  # ============================================================ #
  # Main computation based on the number of mediators
  # ============================================================ #
  if (ncol(M) > 1) {
    out <- PoC_AB_multi(S = S, M = M, Y = Y, X = X, B = B, lambda = lambda)
  } else {
    out <- PoC_AB(S = S, M = M, Y = Y, X = X, B = B, lambda = lambda)
  }

  # ============================================================ #
  # Return structured output
  # ============================================================ #
  return(structure(list(
    NIE = as.numeric(out$NIE) * (s - s_star),
    p_value_NIE = out$p_value_NIE,
    NDE = as.numeric(out$NDE) * (s - s_star),
    p_value_NDE = out$p_value_NDE,
    NTE = as.numeric(out$NTE) * (s - s_star),
    p_value_NTE = out$p_value_NTE
  ), class = "abYlmMlmResult"))
}
