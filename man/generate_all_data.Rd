% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_all_data.R
\name{generate_all_data}
\alias{generate_all_data}
\title{Data simulated from a directed acyclic graph (DAG) model}
\usage{
generate_all_data(
  n = 200,
  alpha_S = 0,
  beta_M = 0,
  tau_S = 1,
  alpha_vec = rep(1, 3),
  beta_vec = rep(1, 3),
  M.family = stats::gaussian(),
  Y.family = stats::gaussian(),
  sigma_M = 0.5,
  sigma_Y = 0.5
)
}
\arguments{
\item{n}{a sample size}

\item{alpha_S}{a parameter in mediator model M~S+X for S}

\item{beta_M}{a parameter in outcome model Y~S+M+X for M}

\item{tau_S}{a parameter in outcome model Y~S+M+X for S}

\item{alpha_vec}{a parameter vector in mediator model M~S+X for X}

\item{beta_vec}{a parameter vector in outcome model Y~S+M+X for X}

\item{M.family}{an error distribution and link function to be used in the mediator model. See \code{\link{family}} for details of family functions. The default family is \code{\link{gaussian}}.}

\item{Y.family}{an error distribution and link function to be used in the outcome model. See \code{\link{family}} for details of family functions. The default family is \code{\link{gaussian}}.}

\item{sigma_M}{the variance of the error term in the mediator model M~S+X}

\item{sigma_Y}{the variance of the error term in the outcome model Y~S+M+X}
}
\value{
A list with the following components:
\item{S}{exposure}
\item{M}{mediator}
\item{Y}{outcome}
\item{X}{confounder}
}
\description{
Simulate data from a directed acyclic graph (DAG) model.
}
\examples{
out <- generate_all_data(n = 200)
summary(out$S)
summary(out$M)
summary(out$Y)
summary(out$X)

}
