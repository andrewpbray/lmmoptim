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
#' @param lines matrix containing the constants that define the shape of each
#'   term in the sum. Output from \code{\link{findlines}}.
#'
#' @return A matrix containing:
#' \itemize{
#'   \item limits of the box in \eqn{\sigma^2_e} (numeric vector of length 2)
#'   \item limits of the box in \eqn{\sigma^2_s} (numeric vector of length 2)
#'   \item the status of the box relative to each line (factor of length of the
#'   number of terms in the sum)
#'   \item the bounds on the objective function (numeric vector of length 2 of
#'   the form \code{c(lower, upper)})
#'   }

makebox <- function(lims.sigsqe = NA, lims.sigsqs = NA, status = NA, lines) {
    # If the box has a parent then it inherits the parent's status.  But if the
    # parent straddles a line, we must check whether the child also straddles the
    # line.
    strad <- which(status == 2)
    status[strad] <- getstatus(lims.sigsqe = lims.sigsqe, lims.sigsqs = lims.sigsqs,
        lines = lines[strad, ])

    # we could get some of the bounds from the parent, but it's just as easy to
    # recalculate them
    bounds <- getbounds(lims.sigsqe = lims.sigsqe, lims.sigsqs = lims.sigsqs,
        status = status, lines = lines)

    return(c(lims.sigsqe, lims.sigsqs, status, bounds))
}

#' Split a parent box into 4 children.
#'
#' Split a parent box into four equally sized boxes as part of the branching part
#' of the algorithm.
#'
#' @param box numeric vector containing the properties of a parent box. Created as output
#'   from \code{\link{makebox}}.
#' @param lines matrix containing the constants that define the shape of each
#'   term in the sum. Output from \code{\link{findlines}}.
#'
#' @return a matrix of the four child boxes.

splitbox <- function(box, lines) {
  # split a box into four children
  boxes <- matrix(NA, nrow = 4, ncol = length(box))
  boxes[1, ] <- makebox(lims.sigsqe = c(box$lims.sigsqe.l, mean(box$lims.sigsqe.l, box$lims.sigsqe.u)), # NW
                lims.sigsqs = c(mean(box$lims.sigsqs.l, box$lims.sigsqs.u), box$lims.sigsqs.u),
                status = status, lines = lines)
  boxes[2, ] <- makebox(lims.sigsqe = c(mean(box$lims.sigsqe.l, box$lims.sigsqe.u), box$lims.sigsqe.u), # NE
                lims.sigsqs = c(mean(box$lims.sigsqs.l, box$lims.sigsqs.u), box$lims.sigsqs.u),
                status = status, lines = lines)
  boxes[3, ] <- makebox(lims.sigsqe = c(box$lims.sigsqe.l, mean(box$lims.sigsqe.l, box$lims.sigsqe.u)), # SW
                lims.sigsqs = c(box$lims.sigsqs.l, mean(box$lims.sigsqs.l, box$lims.sigsqs.u)),
                status = status, lines = lines)
  boxes[4, ] <- makebox(lims.sigsqe = c(mean(box$lims.sigsqe.l, box$lims.sigsqe.u), box$lims.sigsqe.u), # SE
                lims.sigsqs = c(box$lims.sigsqs.l, mean(box$lims.sigsqs.l, box$lims.sigsqs.u)),
                status = status, lines = lines)
  return(boxes)
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
#' @param lines matrix containing the constants that define the shape of each
#'   term in the sum. Output from \code{\link{findlines}}.
#'
#' @return A factor of length of \code{lines} denoting the location of the box
#'   relative to each line.

getstatus <- function(lims.sigsqe, lims.sigsqs, lines) {
    # Is a box above, below, or straddling the lines with these slopes and
    # intercepts?

    # value of the lines at the left side of the box
    tmp1 <- lines$int.sigsqs + lines$slope * lims.sigsqe.l
    tmp1[is.infinite(tmp1)] <- NA

    # value of the lines at the right side of the box
    tmp2 <- lines$int.sigsqs + lines$slope * lims.sigsqe.u
    tmp2[is.infinite(tmp2)] <- NA

    # where is the box relative to the lines?
    status <- rep(2, nrow(lines)) # straddle
    status[lims.sigsqs.l > tmp1] <- 1 # above
    status[lims.sigsqs.u < tmp2] <- 0 # below

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
#' @param lines matrix containing the constants that define the shape of each
#'   term in the sum. Output from \code{\link{findlines}}.
#' @param status factor of length of \code{lines} denoting the location of the box
#'   relative to each line.
#'
#' @return A matrix with a row corresponding each term in the sum. The first column
#' contains the lower bound of that term within the box and the second column
#' contains the upper bound. Note that the \code{colSums} of this matrix yields the
#' \code{c(lower, upper)} bounds on the full objective function.

getbounds <- function(lims.sigsqe, lims.sigsqs, lines, status) {
    # evaluate each line at the upper-right corner of the box
    ur <- lines[, "a"] * lims.sigsqs[2] + lines[, "b"] * lims.sigsqe[2]
    eval.ur <- -0.5 * (lines[, "multiplier.log"] * log(ur) + lines[, "multiplier.inv"]/ur)
    # evaluate each line at the lower-left corner of the box
    ll <- lines[, "a"] * lims.sigsqs[1] + lines[, "b"] * lims.sigsqe[1]
    eval.ll <- ifelse(ll == 0, -Inf, -0.5 * (lines[, "multiplier.log"] * log(ll) +
        lines[, "multiplier.inv"]/ll))

    bounds <- matrix(rep(eval.ur, 2), ncol = 2)
    # The next two lines of code are for lines that are not straddled.  'above'
    # means the box is above the line
    bounds[status != 1, 1] <- eval.ll[status != 1]
    bounds[status == 1, 2] <- eval.ll[status == 1]
    # now we'll take care of straddled lines
    strad <- status == 2
    bounds[strad, 1] <- pmin(eval.ur[strad], eval.ll[strad])
    # for the upper bound, we can evaluate anywhere on the line, so we might as
    # well evaluate at (int.sigsqe,0)
    bounds[strad, 2] <- ifelse(is.na(lines[strad, "int.sigsqe"]),
                               -0.5 * (lines[strad, "multiplier.log"] * log(lines[strad, "int.sigsqs"]) +
                                       lines[strad, "multiplier.inv"]/lines[strad, "int.sigsqs"]),
                               -0.5 * (lines[strad, "multiplier.log"] * log(lines[strad, "int.sigsqe"]) +
                                       lines[strad, "multiplier.inv"]/lines[strad, "int.sigsqe"]))

    return(colSums(bounds))
}
