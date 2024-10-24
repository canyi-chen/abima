#' Adaptive Bootstrap Test in Mediation Analysis with Linear Models or Generalized Linear Models
#'
#' @description
#' \code{abYlm.Mglm} conducts an adaptive bootstrap test in a mediation analysis in which mediator M is modeled by a generalized linear model (GLM) and outcome Y is modeled by a linear model (LM). The family of the generalized linear model is specified according to \code{M.family}.
#'
#' @param S an n-by-1 vector for exposure.
#' @param M an n-by-1 vector for mediator.
#' @param Y an n-by-1 vector for outcome.
#' @param X an n-by-p matrix for p confounders. Do not include intercept in matrix X. In the absence of confounders, matrix is not required to be specified.
#' @param covariates_cfder a vector of p confounders to adjust mediation effects, and default is zero (or void). The length of this vector is p, the same as the column dimension of X.
#' @param M.family the error distribution and link function for the mediator model. The default family is \code{gaussian()}.
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
#' @example man/examples/example_abYlm.Mglm.R
#'
#' @export
abYlm.Mglm <- function(S, M, Y, X = NULL, covariates_cfder = NULL,
                       M.family = stats::gaussian(), s = 1, s_star = 0, B = 199,
                       lambda = 2) {

  # ============================================================ #
  # Parameters checking and cleaning
  # ============================================================ #

  # Check M.family
  if (is.character(M.family))
    M.family <- get(M.family, mode = "function", envir = parent.frame())
  if (is.function(M.family))
    M.family <- M.family()
  if (is.null(M.family$family)) {
    # message(deparse(M.family))
    stop("'M.family' not recognized")
  }

  S <- as.matrix(S)
  M <- as.matrix(M)
  Y <- as.matrix(Y)
  # Ensure scalar exposure and outcome
  stopifnot("We currently only support a single exposure S and a single outcome Y." =
              (ncol(S) == 1) & ncol(Y) == 1)

  if(is.null(X)) {
    X <- matrix(1, nrow(S), 1) # Add intercept column
    if(is.null(covariates_cfder)) {
      covariates_cfder <- 1
    } else {
      warning("X is null while covariates_new is not null. We will default covariates_new at zero.")
      covariates_cfder <- 1
    }
  } else {
    X <- cbind(1, as.matrix(X))
    if(is.null(covariates_cfder)) {
      covariates_cfder <- c(1, rep(0, ncol(X) - 1))
    } else {
      stopifnot("The length of covariates_new must match the number of confounders, i.e., the column dimension of X." = length(covariates_cfder) ==  (ncol(X) - 1))
      covariates_cfder <- c(1, covariates_cfder)
    }
  }



  # Ensure support for a single mediator
  if (ncol(M) > 1) {
    stop("We currently only support a single mediator M.")
  }

  # ============================================================ #
  # Main computation
  # ============================================================ #
  out <- PoC_AB_GLM_LM(S = S, M = M, Y = Y, X = X, M.family = M.family,
                       lambda = lambda, s = s, s_star = s_star,
                       covariates_new = covariates_cfder, B = B)

  # ============================================================ #
  # Return structured output
  # ============================================================ #
  return(structure(list(
    NIE = as.numeric(out$NIE),
    p_value_NIE = out$p_value_NIE,
    NDE = as.numeric(out$NDE),
    p_value_NDE = out$p_value_NDE,
    NTE = as.numeric(out$NTE),
    p_value_NTE = out$p_value_NTE
  ), class = "abYlmMglmResult"))
}
