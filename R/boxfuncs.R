# When dividing the plane into boxes, a box should have the following
# attributes: xlims=sigsqelims, ylims=sigsqslims for each line, an indicator
# of whether the box is above, below, or straddling the line for each line,
# upper and lower bounds on RLL

makebox <- function(lims.sigsqs = NA, lims.sigsqe = NA, status = NA, lines) {

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
}  # end makebox

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
    return(list(NW, NE, SW, SE))  # list of four boxes
}

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
    status <- rep("straddle", nrow(lines))
    status[lims.sigsqs[1] > tmp1] <- "above"
    status[lims.sigsqs[2] < tmp2] <- "below"

    return(status)
}

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
    bounds[status != "above", 1] <- eval.ll[status != "above"]
    bounds[status == "above", 2] <- eval.ll[status == "above"]
    # now we'll take care of straddled lines
    strad <- status == "straddle"
    bounds[strad, 1] <- pmin(eval.ur[strad], eval.ll[strad])
    # for the upper bound, we can evaluate anywhere on the line, so we might as
    # well evaluate at (int.sigsqe,0)
    #bounds[strad, 2] <- with(lines-0.5 * (multiplier.log * log(int.sigsqe) + multiplier.inv/int.sigsqe)
    #bounds[is.na(lines$int.sigsqe)] <- with(lines[is.na(lines$int.sigsqe)],
     #                                       -0.5 * (multiplier.log * log(int.sigsqs) +
      #                                                multiplier.inv/int.sigsqs))
    bounds[strad, 2] <- with(lines[strad, ], ifelse(is.na(int.sigsqe), -0.5 *
       (multiplier.log * log(int.sigsqs) + multiplier.inv/int.sigsqs), -0.5 *
       (multiplier.log * log(int.sigsqe) + multiplier.inv/int.sigsqe)))

    return(bounds)
}
