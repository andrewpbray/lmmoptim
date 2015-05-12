#' Plot of linear terms
#'
#' Generate a plot of the lines that correspond to where each term
#' in the sum of the objective function is maximized.
#'
#' @param lines dataframe containing the constants that define the shape of each
#'   term in the sum. Output from \code{\link{findlines}}.

showlines <- function(lines) {
  p <- ggplot(subset(lines, is.finite(slope) & slope < 0),
              aes(x = int.sigsqe, y = 0, xend = 0, yend = int.sigsqs)) +
    xlab(expression(sigma[e]^2)) +
        ylab(expression(sigma[s]^2)) +
    geom_segment() +
    geom_vline(data = subset(lines, slope == -Inf), aes(xintercept = int.sigsqe), lty = 2) +
    theme_bw()

  # Is there a horizontal line?
  if (any(lines$slope == 0)) {
    p <- p + geom_hline(data = subset(lines, slope == 0),
                        aes(yintercept = int.sigsqs), lty = 2)
    }
  return(p)
}

#' Plot of box partition
#'
#' Generates a plot showing the partitioning of the parameter space
#' into boxes.
#'
#' @param box a list of boxes, usually output from \code{\link{fitlmm}}.

showboxes <- function(boxes) {
  p <- ggplot(boxes, aes(xmin = sigsqe.lo, xmax = sigsqe.hi,
                         ymin = sigsqs.lo, ymax = sigsqs.hi))
  p + geom_rect(color = "black", fill = "white") + xlab("Sigma_e^2") + ylab("Sigma_s^2")
}

#' Plot of function
#'
#' Generates a plot of the restricted log-likelihood or log-posterior as
#' estimated by \code{\link{fitlmm}}.
#'
#' The plot produced is an image plot of the lower bound within each box shown
#' in grayscale. The red dot indicates the box with the highest lower bound.
#' Boxes outlined in blue have upper bounds that are higher than the red dot
#' (i.e. it is possible that the blue boxes contain the point at which the
#' true function achieves its maximum value).
#'
#' @param box a list of boxes and their bounds, usually output from
#'   \code{\link{fitlmm}}.

showfunc <- function(boxes) {
  mle <- which.max(boxes$rll.lower)
  mle.val <- boxes$rll.lower[mle]
  hi <- which(boxes$rll.upper > mle.val)
  p <- ggplot(boxes, aes(xmin = sigsqe.lo, xmax = sigsqe.hi,
                         ymin = sigsqs.lo, ymax = sigsqs.hi))
  return(p + geom_rect(aes(fill = rll.lower)) +
           scale_fill_continuous(limits = c(mle.val - 20, mle.val), low = "black",
                                 high = "white", na.value = "black", name = expression(L^b)) +
           geom_rect(data = boxes[hi, ], aes(fill = rll.lower), color = "blue") +
           geom_point(data = boxes[mle, ], aes(x = sigsqe.lo, y = sigsqs.lo), color = "red", size = 5))
}
