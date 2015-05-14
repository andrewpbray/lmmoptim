#' HMO premiums
#'
#' A dataset containing has one row for each of 341 healthcare plans offered by
#' health maintenance organizations (HMOs). For state-level data, see
#' \code{\link{hmo_states}}.
#'
#' This dataset was assembled in the mid-1990s to estimate the cost of moving families
#' from a Department of Defense healthcare plan to private plans. See James Hodges,
#' "Richly Parameterized Linear Models" (2013) for details and analysis.
#'
#' @format A data frame with 340 rows and 6 variables:
#' \describe{
#'   \item{state}{the state (two-letter abbreviation).}
#'   \item{plan}{a two-character code identifying the plan.}
#'   \item{indPrem}{total premium for an individual (dollars).}
#'   \item{famPrem}{total premium for a family (dollars).}
#'   \item{indEnr}{total enrollment of federal employees as individuals (number
#'   of persons).}
#'   \item{famEnr}{total enrollment of federal employees as famlilies (number of
#'   families).}
#'   }
"hmo"

#' HMO states
#'
#' A dataset containing state- or jurisdiction-level information on the states
#' in which health maintainance organization (HMO) plans are offered. See
#' \code{\link{hmo}} for plan-level information.
#'
#' Regional Abbreviations: MA = Middle Atlantic, MT = Mountain,
#' NC = North Central, NE = New England, PA = Pacific, SA = South Atlantic,
#' SC = South Central. See James Hodges, "Richly Parameterized Linear Models"
#' (2013) for details and analysis of this dataset.
#'
#' @format A data frame with 45 rows and 4 variables:
#' \describe{
#'   \item{state}{the state (two-letter abbreviation).}
#'   \item{expPerAdm}{state average expenses per admission
#'   (dollars;  American Medical Association 1991 Annual Survey of Hospitals).}
#'   \item{pop}{population (1990 US Census).}
#'   \item{region}{region, from the Marion Merrill Dow Managed Care Digest
#'   1991. See details for abbreviations.}
#'   }
"hmo_states"

#' Global mean surface temperature
#'
#' A dataset containing the global mean surface temperature deviations for the
#' years 1881 through 2005 inclusive.
#'
#' This version was downloaded in 2006 and has been superseded by updated
#' series based on more temperature measurements. See James Hodges, "Richly
#' Parameterized Linear Models" (2013) for details and analysis.
#'
#' @format A data frame with 125 rows and 2 variables:
#' \describe{
#'   \item{Year}{integer year.}
#'   \item{temp.dev}{mean temperature deviation for that year (in units 0.01C).}
#' }
#' @source \url{http://data.giss.nasa.gov/gistemp/tabledata/GLB.Ts.txt}
"gmst"
