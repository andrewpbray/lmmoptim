#' Learn RLL or Log-posterior for two variance MM.
#'
#' Learns the shape of the objective function within epsilon.
#'
#' This is the primary function that implements the branch-bound-kill
#' algorithm described in "Approximately Exact Calculations for Linear
#' Mixed Models" by Lavine & Hodges (2015). This an interative algorithm
#' that can require substantial computation time. It is recommended that
#' the user start with a conservative \code{maxit} and do additional
#' computations as necessary.
#'
#' The arguments of this function include five settings to control when
#' to stop branching (subdividing) a box:
#' \enumerate{
#'   \item when the upper and lower bounds are within \code{eps} of one another.
#'   \item when the width of \eqn{\sigma^2_e} is less than \code{delE}.
#'   \item when the width of \eqn{\sigma^2_s} is less than \code{delS}.
#'   \item when the upper bound is more than \code{M} below the highest global
#'   lower bound.
#'   \item when the number of iterations reaches \code{maxit}.
#'   }
#'
#' @param lines a dataframe that contains the constants that define the
#'   line represented by each term in the sum. Created as output from
#'   \code{\link{findlines}}.
#' @param startbox a list of boxes and their bounds, possibly the output
#'   from a previous call to \code{fitlmm}. If left empty, will create
#'   a new startbox from a call to \code{makebox}.
#' @param eps a non-negative numeric indicating the tolerance within which you
#'   would like to learn the function. Default is 0.
#' @param delE a non-negative numeric indicating the width of \eqn{\sigma^2_e}
#'   beyond which the algorithm will stop branching boxes. Default is 0.
#' @param delS a non-negative numeric indicating the width of \eqn{\sigma^2_s}
#'   beyond which the algorithm will stop branching boxes. Default is 0.
#' @param M a non-negative numeric indicating the size of the buffer between
#'   the highest global lower bound and the upper bound of a given box past
#'   that box will now be branched further.  Default is \code{Inf}.
#' @param maxit a positive integer indicating the maximum number of iterations
#'   of the algorithm.  Default is 10.
#' @param ratio a logical indicating if \code{delE} and \code{delS} are specified as
#'   ZZQ.  Default is \code{FALSE}.
#' @param lognote a string to append to the "log.out" log file to annotate
#'   the purpose of each run of the algorithm. For use in benchmarking computation
#'   time. Default is \code{"summary"}.
#'
#' @return A list of boxes including their limits in \eqn{\sigma^2_e} and
#'   \eqn{\sigma^2_e} as well as their bounds. Running this function also results
#'   in the creation of a log file called "log.out" containing box counts at every
#'   iteration.

fitlmm <- function(lines, startbox, eps = 0, delE = 0, delS = 0, M = Inf, maxit = 10,
                   ratio = FALSE, lognote = "summary") {

    start_time <- Sys.time()

    if (missing(lines)) {
        print("please supply lines")
        return
    }
    if (missing(startbox)) {
        startbox <- makebox(lines = lines, lims.sigsqs = c(0, max(lines$int.sigsqs[is.finite(lines$int.sigsqs)])),
            lims.sigsqe = c(0, max(lines$int.sigsqe[is.finite(lines$int.sigsqe)])),
            status = rep("straddle", nrow(lines)))
    }
    inactive <- list()
    ninact <- 0
    # check whether startbox is a box or a list of boxes
    active <- ifelse(length(startbox) == 4 && identical(names(startbox), c("lims.sigsqe",
        "lims.sigsqs", "status", "bounds")), list(startbox), startbox)
    nact <- length(active)

    lowbound <- -Inf
    iter <- 0

    # conditions under which a box becomes inactive
    killfunc <- function(box, lb, M, eps, delE, delS, ratio) {
        # lb is a lower bound on max f; it changes at each iteration. M, eps,
        # delE, delS stay constant throughout the iterations.
        cond.low <- box$bounds[2] < lb - M
        cond.eps <- diff(box$bounds) < eps
        cond.E <- ifelse(ratio, diff(log(box$lims.sigsqe)) < delE, diff(box$lims.sigsqe) <
            delE)
        cond.S <- ifelse(ratio, diff(log(box$lims.sigsqs)) < delS, diff(box$lims.sigsqs) <
            delS)
        return(cond.low || cond.eps || cond.E || cond.S)
    }

    while (nact > 0 && iter < maxit) {
        # get new global highest lower bound
        low.act <- max(vapply(X = active, FUN = function(box) {
            box$bounds[1]
        }, FUN.VALUE = 0.1))
        lowbound <- max(lowbound, low.act)
        # look at active list for boxes to be killed. return logical vector to kill
        kill <- vapply(X = active, FUN = killfunc, FUN.VALUE = TRUE, lb = lowbound,
            M = M, eps = eps, delE - delE, delS = delS, ratio = ratio)
        nkill <- sum(kill)
        if (nkill > 0) {
            inactive[(ninact + 1):(ninact + nkill)] <- active[kill]
        }
        ninact <- length(inactive)
        kids <- list()
        length(kids) <- sum(!kill) * 4
        nkids <- 0
        for (i in which(!kill)) {
            kids[(nkids + 1):(nkids + 4)] <- splitbox(active[[i]], lines)  # boxes are split into 4 parts
            nkids <- nkids + 4
        }
        active <- kids
        nact <- length(active)

        iter <- iter + 1
        write(c("iteration", iter, "nact", nact, "ninact", ninact, "lowbound",
            lowbound), file = "log.out", ncolumns = 8, append = TRUE)
    }


    tmp <- t(vapply(X = c(active, inactive), FUN = function(box) {
        c(box$lims.sigsqs, box$lims.sigsqe, box$bounds)
    }, FUN.VALUE = c(sigsqs.lo = 0.1, sigsqs.hi = 0.1, sigsqe.lo = 0.1, sigsqe.hi = 0.1,
        rll.lower = 0.1, rll.upper = 0.1)))
    end_time <- Sys.time()
    write(paste(lognote, "runtime: ", round(end_time - start_time, digits = 3)),
          file = "log.out", ncolumns = 1, append = TRUE)
    return(data.frame(tmp))
}
