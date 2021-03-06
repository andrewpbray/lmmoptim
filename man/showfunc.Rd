% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/vizfuncs.R
\name{showfunc}
\alias{showfunc}
\title{Plot of function}
\usage{
showfunc(boxes)
}
\arguments{
\item{box}{a list of boxes and their bounds, usually output from
\code{\link{fitlmm}}.}
}
\description{
Generates a plot of the restricted log-likelihood or log-posterior as
estimated by \code{\link{fitlmm}}.
}
\details{
The plot produced is an image plot of the lower bound within each box shown
in grayscale. The red dot indicates the box with the highest lower bound.
Boxes outlined in blue have upper bounds that are higher than the red dot
(i.e. it is possible that the blue boxes contain the point at which the
true function achieves its maximum value).
}

