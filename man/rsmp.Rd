% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smp.R
\name{rsmp}
\alias{rsmp}
\title{Sums & Maxima of Multivariate Pareto II Vectors}
\usage{
rsmp(n, alpha, beta, m, rounding = FALSE)
}
\arguments{
\item{n}{number of observations.}

\item{alpha, beta}{tail/shape parameter and scale parameters of a Multivariate
Pareto II distribution (see details for the parameterization used).}

\item{m}{the number of elements in the vectors. This is expected to either be
a single number (all generated vectors have the same number of elements) or a
vector of length n (in which case the number of elements in each vector is
specified)}

\item{rounding}{logical; if TRUE, generated Pareto II values are rounded to
the nearest whole number. This feature is present in case this is important
for a comparison of actual observations to some theoretical output.}
}
\value{
a data.frame() with three columns:
\itemize{
\item n the number of elements in the Pareto II vector
\item x the sum of the elements in the Pareto II vector
\item y the maximal value in the Pareto II vector
}
}
\description{
Sums & Maxima of Multivariate Pareto II Vectors
}
\examples{
rsmp(n = 3, m = c(3,4,2), alpha = 0.03, beta = 1/100)
\dontrun{
  n        x        y
1 3 241.6733 157.2649
2 4 357.7778 187.8859
3 2 219.1081 140.8779
}
}
