fitlmm <- function(lines, startbox, eps = 0, delE = 0, delS = 0, M = Inf, maxit = 10,
                   ratio = FALSE, lognote = "summary") {
    # five control settings don't subdivide boxes whose upper and lower bounds
    # are within eps don't subdivide boxes whose sigsqE spans less than delE
    # don't subdivide boxes whose sigsqS spans less than delS don't subdivide
    # boxes whose max RLL is more than M below the MRLE don't do more than maxit
    # iterations ratio determines whether the sigsqE and sigsqS conditions are
    # for ratios or absolute differences
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
    # active <- list() active[[1]] <- startbox
    nact <- length(active)
    lowbound <- -Inf
    iter <- 0

    # conditions under which a box becomes inactive
    killfunc <- function(box, lb, M, eps, delE, delS, ratio) {
        # lb is a lower bound on max logRL; it changes at each iteration M, eps,
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
        # Find the lower bound of each box and the maximum of the lower bounds.  For
        # each active box, either make it inactive or divide it.
        low.act <- max(vapply(X = active, FUN = function(box) {
            box$bounds[1]
        }, FUN.VALUE = 0.1))
        lowbound <- max(lowbound, low.act)
        kill <- vapply(X = active, FUN = killfunc, FUN.VALUE = TRUE, lb = lowbound,
            M = M, eps = eps, delE - delE, delS = delS, ratio = ratio)
        nkill <- sum(kill)
        if (nkill > 0) {
            inactive[(ninact + 1):(ninact + nkill)] <- active[kill]
        }
        ninact <- length(inactive)
        kids <- list()
        nkids <- 0
        for (i in which(!kill)) {
            kids[(nkids + 1):(nkids + 4)] <- splitbox(active[[i]], lines)  # boxes are split into 4 parts
            nkids <- nkids + 4
        }
        active <- kids
        nact <- length(active)

        iter <- iter + 1
        write(c("iteration", iter, "nact", nact, "ninact", ninact, "lowbound",
            lowbound), file = "bigcode.out", ncolumns = 8, append = TRUE)
    }


    tmp <- t(vapply(X = c(active, inactive), FUN = function(box) {
        c(box$lims.sigsqs, box$lims.sigsqe, box$bounds)
    }, FUN.VALUE = c(sigsqs.lo = 0.1, sigsqs.hi = 0.1, sigsqe.lo = 0.1, sigsqe.hi = 0.1,
        rll.lower = 0.1, rll.upper = 0.1)))
    end_time <- Sys.time()
    write(paste(lognote, "runtime: ", round(end_time - start_time, digits = 3)),
          file = "bigcode.out", ncolumns = 1, append = TRUE)
    return(data.frame(tmp))
}
