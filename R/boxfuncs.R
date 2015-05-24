#' Make a box with bounds
#'
#' Given the limits of the box, calculate its status relative to each line and
#' calculate the bounds on the objective function within that box.
#'
#' @param lims.sigsqe numeric vector of length 2 of the form \code{c(lower, upper)}
#'   containing the x-limits of the box.
#' @param lims.sigsqs numeric vector of length 2 of the form \code{c(lower, upper)}.
#'   containing the y-limits of the box.
#' @param status factor of length of \code{lines} denoting the location of the box
#'   relative to each line. If the box has a parent, these correspond to the status
#'   of that parent. The levels are \{above, below, straddle\}.
#' @param lines dataframe containing the constants that define the shape of each
#'   term in the sum. Output from \code{\link{findlines}}.
#'
#' @return A list containing:
#' \itemize{
#'   \item limits of the box in \eqn{\sigma^2_e} (numeric vector of length 2)
#'   \item limits of the box in \eqn{\sigma^2_s} (numeric vector of length 2)
#'   \item the status of the box relative to each line (factor of length of the
#'   number of terms in the sum)
#'   \item the bounds on the objective function (numeric vector of length 2 of
#'   the form \code{c(lower, upper)})
#'   }

makebox <- function(lims.sigsqe = NA, lims.sigsqs = NA, status = NA, lines) {

    # sanity checks
    if (missing(lims.sigsqs) || missing(lims.sigsqe))
        print("please supply lims.sigsqs and lims.sigsqe")
    if (missing(status))
        print("please supply status")
    if (missing(lines))
        print("please supply lines")

    # If the box has a parent then it inherits the parent's status.  But if the
    # parent straddles a line, we must check whether the child also straddles the
    # line.
    strad <- which(status == "straddle")
    status[strad] <- getstatus(lims.sigsqe = lims.sigsqe, lims.sigsqs = lims.sigsqs,
        lines = lines[strad, ])

    # we could get some of the bounds from the parent, but it's just as easy to
    # recalculate them
    bounds <- getbounds(lims.sigsqe = lims.sigsqe, lims.sigsqs = lims.sigsqs,
        status = status, lines = lines)

    return(list(lims.sigsqe = lims.sigsqe, lims.sigsqs = lims.sigsqs, status = status,
        bounds = colSums(bounds)))
}

#' Split a parent box into 4 children.
#'
#' Split a parent box into four equally sized boxes as part of the branching part
#' of the algorithm.
#'
#' @param box list containing the properties of a parent box. Created as output
#'   from \code{\link{makebox}}.
#' @param lines dataframe containing the constants that define the shape of each
#'   term in the sum. Output from \code{\link{findlines}}.
#'
#' @return a list of the four child boxes.

splitbox <- function(box, lines) {
    # split a box into four children
    NW <- with(box, makebox(lims.sigsqe = c(lims.sigsqe[1], mean(lims.sigsqe)),
        lims.sigsqs = c(mean(lims.sigsqs), lims.sigsqs[2]), status = status,
        lines = lines))
    NE <- with(box, makebox(lims.sigsqe = c(mean(lims.sigsqe), lims.sigsqe[2]),
        lims.sigsqs = c(mean(lims.sigsqs), lims.sigsqs[2]), status = status,
        lines = lines))
    SW <- with(box, makebox(lims.sigsqe = c(lims.sigsqe[1], mean(lims.sigsqe)),
        lims.sigsqs = c(lims.sigsqs[1], mean(lims.sigsqs)), status = status,
        lines = lines))
    SE <- with(box, makebox(lims.sigsqe = c(mean(lims.sigsqe), lims.sigsqe[2]),
        lims.sigsqs = c(lims.sigsqs[1], mean(lims.sigsqs)), status = status,
        lines = lines))
    return(list(NW, NE, SW, SE))
}

#' Get the location of a box relative to lines.
#'
#' Computes the location of a box with given limits relative to each of
#' the lines.
#'
#' @param lims.sigsqe numeric vector of length 2 of the form \code{c(lower, upper)}
#'   containing the x-limits of the box.
#' @param lims.sigsqs numeric vector of length 2 of the form \code{c(lower, upper)}
#'   containing the y-limits of the box.
#' @param lines dataframe containing the constants that define the shape of each
#'   term in the sum. Output from \code{\link{findlines}}.
#'
#' @return A factor of length of \code{lines} denoting the location of the box
#'   relative to each line.

getstatus <- function(lims.sigsqe, lims.sigsqs, lines) {
    # Is a box above, below, or straddling the lines with these slopes and
    # intercepts?

    # value of the lines at the left side of the box
    tmp1 <- lines$int.sigsqs + lines$slope * lims.sigsqe[1]
    tmp1[is.infinite(tmp1)] <- NA

    # value of the lines at the right side of the box
    tmp2 <- lines$int.sigsqs + lines$slope * lims.sigsqe[2]
    tmp2[is.infinite(tmp2)] <- NA

    # where is the box relative to the lines?
    above <- with(lines, ifelse(slope > -Inf, lims.sigsqs[1] > tmp1, lims.sigsqe[1] > int.sigsqe))
    below <- with(lines, ifelse(slope > -Inf, lims.sigsqs[2] < tmp2, lims.sigsqe[2] < int.sigsqe))
    status <- rep("straddle", nrow(lines))
    status[above] <- "above"
    status[below] <- "below"

    return(status)
}

#' Get bounds for a given box
#'
#' Computes the lower and upper bounds on every term in the sum of
#' the objective function within a box.
#'
#' @param lims.sigsqe numeric vector of length 2 of the form \code{c(lower, upper)}
#'   containing the x-limits of the box.
#' @param lims.sigsqs numeric vector of length 2 of the form \code{c(lower, upper)}
#'   containing the y-limits of the box.
#' @param lines dataframe containing the constants that define the shape of each
#'   term in the sum. Output from \code{\link{findlines}}.
#' @param status factor of length of \code{lines} denoting the location of the box
#'   relative to each line.
#'
#' @return A matrix with a row corresponding each term in the sum. The first column
#' contains the lower bound of that term within the box and the second column
#' contains the upper bound. Note that the \code{colSums} of this matrix yields the
#' \code{c(lower, upper)} bounds on the full objective function.

getbounds <- function(lims.sigsqe, lims.sigsqs, lines, status) {
    # small sanity check
    if (missing(lims.sigsqe) || missing(lims.sigsqs))
        print("Please supply lims.sigsqe and lims.sigsqs")
    if (missing(lines))
        print("Please supply lines")
    if (missing(status))
        print("Please supply status")
    if (length(status) != nrow(lines))
        print("length(status) != nrow(lines)")

    # evaluate each line at the upper-right corner of the box
    ur <- with(lines, a * lims.sigsqs[2] + b * lims.sigsqe[2])
    eval.ur <- with(lines, -0.5 * (multiplier.log * log(ur) + multiplier.inv/ur))
    # evaluate each line at the lower-left corner of the box
    ll <- with(lines, a * lims.sigsqs[1] + b * lims.sigsqe[1])
    eval.ll <- ifelse(ll == 0, -Inf, with(lines, -0.5 * (multiplier.log * log(ll) +
        multiplier.inv/ll)))

    bounds <- matrix(rep(eval.ur, 2), ncol = 2)
    # The next two lines of code are for lines that are not straddled.  'above'
    # means the box is above the line
    bounds[status != "above", 1] <- eval.ll[status != "above"] # lower
    bounds[status == "above", 2] <- eval.ll[status == "above"] # upper
    # now we'll take care of straddled lines
    strad <- status == "straddle"
    bounds[strad, 1] <- pmin(eval.ur[strad], eval.ll[strad]) # lower
    # for the upper bound, we can evaluate anywhere on the line, so we might as
    # well evaluate at (int.sigsqe,0)
    bounds[strad, 2] <- with(lines[strad, ], ifelse(is.na(int.sigsqe), -0.5 *
       (multiplier.log * log(int.sigsqs) + multiplier.inv/int.sigsqs), -0.5 *
       (multiplier.log * log(int.sigsqe) + multiplier.inv/int.sigsqe)))

    return(bounds)
}
