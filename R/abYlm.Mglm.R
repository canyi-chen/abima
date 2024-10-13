#' Adaptive Bootstrap for Mediation Analysis with Linear Models and Generalized Linear Models
#'
#' @description
#' \code{abYlm.Mglm} conducts adaptive bootstrap for mediation analysis with the generalized linear model on the mediator M and the linear model on the outcome Y. The family of the generalized linear model is specified with \code{M.family}.
#'
#' @param S a n-by-1 matrix for exposure.
#' @param M a n-by-m matrix for mediator, each row corresponding to an observation. The dimension m could be 1 or
#' larger than 1.
#' @param Y a n-by-1 matrix for outcome.
#' @param X a n-by-p matrix for confounder.
#' @param M.family a description of the error distribution and link function to be used in the mediator model. This is the
#' result of a call to a family function. (See \code{\link{family}} for detials of family functions.) The default
#' family is gaussian().
#' @param intercept logical, default at FALSE; if the intercept has been included in the confounder. If the intercept
#' is included, please put it in the first column of X.
#' @param s exposure level, default at 1
#' @param s_star another exposure level, default at 0
#' @param B the number of bootstrap, default at 199
#' @param lambda the constant used in the pretest when conducing adaptive bootstrap, default at 2.
#' @param covariates_new the level of confounders you want to condition on, default at zero.
#'
#' @returns mediation_effect the estimated mediation effect between s and s_star (conditioning on the covariates_new).
#' @returns p_value the p value.
#' @references He, Y., Song, P. X. K., and Xu, G. (2023), “Adaptive bootstrap tests for composite null hypotheses in the mediation pathway analysis,” Journal of the Royal Statistical Society Series B: Statistical Methodology, qkad129. \url{https://doi.org/10.1093/jrsssb/qkad129}.
#' @example man/examples/example_abYlm.Mglm.R
#'
#' @export
abYlm.Mglm <- function(S, M, Y, X, M.family = stats::gaussian(),
                       intercept = FALSE, s = 1, s_star = 0, B = 199,
                       lambda = 2, covariates_new = NULL) {

  # ============================================================ #
  # Parameters checking and cleaning
  # ============================================================ #
  S <- as.matrix(S)
  M <- as.matrix(M)
  Y <- as.matrix(Y)
  X <- as.matrix(X)

  # only support scalar exposure and outcome
  stopifnot("We currently only support a single exposure S and a single outcome Y." =
              (ncol(S) == 1) & ncol(Y) == 1)

  if (is.null(covariates_new)) {
    if (intercept) {
      covariates_new <- rep(0, ncol(X) - 1)
    } else {
      covariates_new <- rep(0, ncol(X))
    }
  }
  # check intercept
  if (!intercept) {
    X <- cbind(1, X)
  }

  if (ncol(M) > 1) {
    stop("We currently only support a single mediator M.")
  } else {
    out <- PoC_AB_GLM_LM(S = S, M = M, Y = Y, X = X, M.family =
      M.family, lambda = lambda, s = s, s_star = s_star, covariates_new = covariates_new,
                         B = B)
    return(
      list(
        mediation_effect = as.numeric(out$mediation_effect),
        p_value = out$p_value
      )
    )
  }

}
