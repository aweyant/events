% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ted_utils.R
\name{ted}
\alias{ted}
\alias{pteddt}
\alias{psdeddt}
\alias{construct_cdf_ted}
\alias{construct_pted}
\title{\strong{P}robability of \strong{T}rivariate \strong{E}vent \strong{D}istribution \strong{D}erived
\strong{T}otal (P TED DT) and \strong{P}robability of \strong{S}et \strong{D}uration \strong{E}vent
\strong{D}istribution \strong{D}erived \strong{T}otal}
\usage{
pteddt(
  q = NULL,
  threshold,
  event_duration_marginal_pmf,
  event_duration_marginal_pmf_args = NULL,
  event_length_arg_name = "n",
  event_magnitude_conditional_cdf,
  event_magnitude_conditional_cdf_args,
  lower.tail = TRUE,
  tol = 1e-06,
  max_N = 10 * (2^9)
)

psdeddt(
  q = NULL,
  length = 1,
  threshold,
  event_magnitude_conditional_cdf,
  event_magnitude_conditional_cdf_args,
  event_length_arg_name = "n",
  lower.tail = TRUE,
  return_NA_outside_support = FALSE
)

construct_cdf_ted(
  args,
  event_bivariate_conditional_cdf,
  event_duration_marginal_pmf,
  conditional_cdf_argument_transformations,
  env = parent.frame()
)

construct_pted(
  args,
  pbed,
  event_duration_marginal_pmf,
  pbed_argument_transformations,
  env = parent.frame(),
  tol = 1e-06,
  max_N = 10 * (2^9)
)
}
\arguments{
\item{q}{a number; the event total, including the below-threshold portion}

\item{threshold}{a number; the minimum threshold used to define events}

\item{event_duration_marginal_pmf}{a function; the marginal CDF for event durations}

\item{event_duration_marginal_pmf_args}{a list; arguments for the marginal duration CDF}

\item{event_magnitude_conditional_cdf}{a function; the CDF of event magnitudes
conditioned on durations, e.g. psmp_marginal_x()}

\item{lower.tail}{logical; if TRUE (default), probabilities are
\eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}

\item{tol}{a very small number; minimal exceedance probability desired}

\item{max_N}{a number; highest N to search}

\item{length}{integer; the duration of the event}

\item{return_NA_outside_support}{logical; should an unsupported event cause
a return of NA or 0? Decide here.}
}
\value{
a number in the interval [0,1], the probability of not exceeding or
exceeding q
}
\description{
This evaluates the CDF of the total accumulation of an amount over an event.
}
\details{
pteddt() by way of example:
a 3 day precipitation event might have a magnitude of 430m. If
the threshold for defining the event was 75mm/day, then the actual total
precipitation over this event was 655mm [430 + 3*75].

Note that durations and daily accumulations are, in general, both random
variables, so if we are interested in probabilties (return intervals) of
uninterrupted precipitation accumulations, we must construct an event total
CDF from (at least) two components: the (discrete) PMF describing event
durations and the CDF for magnitudes, conditioned upon a duration. The Law of
Total Probability lets us create a CDF for the event total, which is quite
meaningful if we are studying precipitation.

psdeddt() evaluates the CDF for event totals of a given duration (hence, "set
duration").
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
