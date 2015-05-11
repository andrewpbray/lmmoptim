

library(ggplot2)

showlines <- function(lines) {
    p <- ggplot(subset(lines, is.finite(slope) & slope < 0), aes(x = int.sigsqe,
        y = 0, xend = 0, yend = int.sigsqs)) + xlab(expression(sigma[e]^2)) +
        ylab(expression(sigma[s]^2)) + geom_segment() + geom_vline(data = subset(lines,
        slope == -Inf), aes(xintercept = int.sigsqe), lty = 2) + theme_bw()

    # Is there a horizontal line?
    if (any(lines$slope == 0)) {
        p <- p + geom_hline(data = subset(lines, slope == 0), aes(yintercept = int.sigsqs),
            lty = 2)
    }
    return(p)
}

showboxes <- function(boxes) {
    p <- ggplot(boxes, aes(xmin = sigsqe.lo, xmax = sigsqe.hi, ymin = sigsqs.lo,
        ymax = sigsqs.hi))
    p + geom_rect(color = "black", fill = "white") + xlab("Sigma_e^2") + ylab("Sigma_s^2")

}

showfunc <- function(boxes) {
    mle <- which.max(boxes$rll.lower)
    mle.val <- boxes$rll.lower[mle]
    hi <- which(boxes$rll.upper > mle.val)
    p <- ggplot(boxes, aes(xmin = sigsqe.lo, xmax = sigsqe.hi, ymin = sigsqs.lo,
        ymax = sigsqs.hi))
    return(p + geom_rect(aes(fill = rll.lower)) + scale_fill_continuous(limits = c(mle.val -
        20, mle.val), low = "black", high = "white", na.value = "black", name = expression(L^b)) +
        geom_rect(data = boxes[hi, ], aes(fill = rll.lower), color = "blue") +
        geom_point(data = boxes[mle, ], aes(x = sigsqe.lo, y = sigsqs.lo), color = "red",
            size = 5))
}
