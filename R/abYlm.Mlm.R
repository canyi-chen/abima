#' Adaptive Bootstrap Test in Mediation Analysis with Linear Models
#'
#' @description
#' \code{abYlm.Mlm} performs an adaptive bootstrap test in mediation analysis in which both mediator M and outcome Y are modeled by linear models. The mediators can be multi-dimensional.
#'
#' @param S an n-by-1 vector for exposure.
#' @param M an n-by-m matrix for m mediators, where each row corresponds to a m-element vector of mediators for one subject. The dimension m can be 1 or an integer larger than 1.
#' @param Y an n-by-1 vector for outcome.
#' @param X an n-by-p matrix for p confounders. Do not include intercept in matrix X. In the absence of confounders, matrix is not required to be specified.
#' @param s exposure level with default set at 1.
#' @param s_star baseline exposure level with default set at 0.
#' @param B the number of bootstrap samples with default set at 199.
#' @param lambda the tuning parameter used in the adaptive bootstrap procedure with default set at 2.
#'
#' @returns NIE refers to natural indirect effect, which is the estimated when a change of exposure S occurs from s_star to s.
#' @returns p_value_NIE refers to the AB p value for NIE.
#' @returns NDE refers to natural direct effect, which is estimated when exposure S changes from s_star to s.
#' @returns p_value_NDE is the p value for NDE.
#' @returns NTE is the estimated natural total effect when exposure S changes from s_star to s.
#' @returns p_value_NTE is the p value for NTE.
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
