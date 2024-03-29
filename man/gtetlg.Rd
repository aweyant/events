% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtetlg.R
\name{gtetlg}
\alias{gtetlg}
\alias{rgtetlg}
\alias{pgtetlg}
\alias{cdf_gtetlg}
\title{Title}
\usage{
rgtetlg(n, beta, prob_p, prob_q, rounding = FALSE)

pgtetlg(
  x_lower = 0,
  x_upper = Inf,
  y_lower = 0,
  y_upper = Inf,
  n_lower = 1,
  n_upper = Inf,
  beta,
  prob_p,
  prob_q
)

cdf_gtetlg(x, y, n, alpha, beta, prob_p, prob_q)
}
\arguments{
\item{n}{when provided to rsmp, this is the number of observations. Otherwise
it represents the duration of an event.}

\item{prob_p}{probability of success in each trial after the first trial, if
the first trial was not a success. \eqn{0 <} prob_p \eqn{\le 1}.}

\item{prob_q}{probability of success in the first trial. \eqn{0 \le} prob_q \eqn{\le 1}}

\item{rounding}{logical; if TRUE, generated Pareto II values are rounded to
the nearest whole number. This feature is present in case this is important
for a comparison of actual observations to some theoretical output.}

\item{x_lower, x_upper, y_lower, y_upper}{each are individual numbers; limits of
integration of the SMP CDF}

\item{n_lower, n_upper}{each are individual numbers; limits of summation of
for the duration of the a trivariate distribution, where n_lower is the
lowest duration included in the summation and n_upper is the largest included
in the summation; for example, n_lower = 3, n_upper = 5 sums over {3,4,5},
whereas n_lower = 3, n_upper = 3 sums over {3}.}

\item{x, y}{numbers; magnitude and maxima of a Pareto II vector}

\item{alpha, beta}{tail/shape parameter and scale parameters of a Multivariate
Pareto II distribution (see details for the parameterization used).}
}
\description{
Title
}
