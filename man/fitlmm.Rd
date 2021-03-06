% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/mainfunc.R
\name{fitlmm}
\alias{fitlmm}
\title{Learn RLL or Log-posterior for two variance MM.}
\usage{
fitlmm(lines, startbox, eps = 0, delE = 0, delS = 0, M = Inf,
  maxit = 10, ratio = FALSE, lognote = "summary")
}
\arguments{
\item{lines}{a dataframe that contains the constants that define the
line represented by each term in the sum. Created as output from
\code{\link{findlines}}.}

\item{startbox}{a list of boxes and their bounds, possibly the output
from a previous call to \code{fitlmm}. If left empty, will create
a new startbox from a call to \code{makebox}.}

\item{eps}{a non-negative numeric indicating the tolerance within which you
would like to learn the function. Default is 0.}

\item{delE}{a non-negative numeric indicating the width of \eqn{\sigma^2_e}
beyond which the algorithm will stop branching boxes. Default is 0.}

\item{delS}{a non-negative numeric indicating the width of \eqn{\sigma^2_s}
beyond which the algorithm will stop branching boxes. Default is 0.}

\item{M}{a non-negative numeric indicating the size of the buffer between
the highest global lower bound and the upper bound of a given box past
that box will now be branched further.  Default is \code{Inf}.}

\item{maxit}{a positive integer indicating the maximum number of iterations
of the algorithm.  Default is 10.}

\item{ratio}{a logical indicating if \code{delE} and \code{delS} are specified as
ZZQ.  Default is \code{FALSE}.}

\item{lognote}{a string to append to the "log.out" log file to annotate
  the purpose of each run of the algorithm. For use in benchmarking computation
  time. Default is \code{"summary"}.}
}
\value{
A list of boxes including their limits in \eqn{\sigma^2_e} and
  \eqn{\sigma^2_e} as well as their bounds. Running this function also results
  in the creation of a log file called "log.out" containing box counts at every
  iteration.
}
\description{
Learns the shape of the objective function within epsilon.
}
\details{
This is the primary function that implements the branch-bound-kill
algorithm described in "Approximately Exact Calculations for Linear
Mixed Models" by Lavine & Hodges (2015). This an interative algorithm
that can require substantial computation time. It is recommended that
the user start with a conservative \code{maxit} and do additional
computations as necessary.

The arguments of this function include five settings to control when
to stop branching (subdividing) a box:
\enumerate{
  \item when the upper and lower bounds are within \code{eps} of one another.
  \item when the width of \eqn{\sigma^2_e} is less than \code{delE}.
  \item when the width of \eqn{\sigma^2_s} is less than \code{delS}.
  \item when the upper bound is more than \code{M} below the highest global
  lower bound.
  \item when the number of iterations reaches \code{maxit}.
  }
}

