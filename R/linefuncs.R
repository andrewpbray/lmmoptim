#' Find the constants that define each line.
#'
#' Given X, Y, Z, find the constants {a_j, b_j, c_j, d_j} that define the shape
#' of each term in the sum.
#'
#' @param x a matrix with n rows corresponding to the fixed effects.
#' @param z a matrix with n rows corresponding to the random effects.
#' @param y a numeric vector of length n.
#' @param SigE an n x n covariance matrix for the random error.
#' @param SigS an n x n covariance matrix for the random effects.
#'
#' @return A dataframe of containing the constants that define the shape of each
#'   term in the sum.

findlines <- function(x, z, y, SigE, SigS) {
    if (!is.matrix(x))
        print("x should be a matrix")
    if (!is.matrix(z))
        print("z should be a matrix")
    if (!is.vector(y))
        print("y should be a vector")

    n <- length(y)
    nx <- ncol(x)
    nz <- ncol(z)

    if (nrow(x) != n)
        print("nrow(x) != length(y)")
    if (nrow(z) != n)
        print("nrow(z) != length(y)")

    if (!is.matrix(SigE))
        print("SigE should be a matrix")
    if (!is.matrix(SigS))
        print("SigS should be a matrix")
    if (nrow(SigE) != n || ncol(SigE) != n)
        print("SigE should be a square matrix to match length(y)")
    if (nrow(SigS) != nz || ncol(SigS) != nz)
        print("SigS should be a square matrix to match ncol(z)")

    # Is SigE the identity?  If not, make it so.
    if (!identical(SigE, diag(n))) {
        tmp <- inv.sqrt(SigE)
        y <- tmp %*% y
        x <- tmp %*% x
        z <- tmp %*% z
        SigE <- diag(n)
    }

    # Is SigS the identity?  If not, make it so.
    if (!identical(SigS, diag(nz))) {
        tmp <- sqrt.m(SigS)
        z <- z %*% tmp
        SigS <- diag(nz)
    }

    # sx, sz, Gamma_x, Gamma_z, Gamma_c # Is Gamma_c really needed?
    qrx <- qr(x, LAPACK = FALSE)  # use LINPACK to get rank(x)
    sx <- qrx$rank
    Gamx <- qr.Q(qrx)[, 1:sx]
    if (sx == 1)
        Gamx <- matrix(Gamx, ncol = 1)

    tmp <- qr.resid(qrx, z)
    qrz <- qr(tmp, LAPACK = FALSE)  # use LINPACK to get rank(x)
    sz <- qrz$rank
    Gamz <- qr.Q(qrz)[, 1:sz]
    if (sz == 1)
        Gamz <- matrix(Gamz, ncol = 1)

    M <- qr.solve(cbind(Gamx, Gamz), cbind(x, z))
    M.zz <- M[-(1:sx), -(1:ncol(x)), drop = FALSE]
    tmp <- svd(M.zz)
    a <- tmp$d^2  # follows Eq (15.1)
    v <- t(tmp$u) %*% t(Gamz) %*% y  # follows Eq (15.4)
    if (length(a) != sz)
        print("length(a) != sz")
    if (length(v) != sz)
        print("length(v) != sz")
    rss <- sum(resid(lm(y ~ cbind(Gamx, Gamz)))^2)

    lines <- data.frame(a = c(a, 0), v = c(v, sqrt(rss)), int.sigsqs = c(v^2/a,
        NA), int.sigsqe = c(v^2, rss/(n - (sx + sz))), slope = c(-1/a, -Inf),
        multiplier.log = c(rep(1, sz), n - (sx + sz)),
        multiplier.inv = c(v^2, rss), b = 1)
    return(lines)
}

#' Add a prior to the RLL
#'
#' Modifies the \code{lines} dataframe that is output from the
#' \code{\link{findlines}} function to add a term corresponding to the
#' prior for use in a Bayesian analysis.
#'
#' Working with the posterior density of \eqn{\sigma^2_e, \sigma^2_e} instead
#' of the RLL requires this incorporation of a prior distribution on those same
#' parameters. We work here with the inverse Gamma, which is the conjugate prior,
#' \enumerate{
#'  \item \eqn{\sigma^2_e ~ IG(\alpha_e, \beta_e)}
#'  \item \eqn{\sigma^2_s ~ IG(\alpha_s, \beta_s)}
#'  }
#' Given a dataframe output from \code{findlines},
#'
#' @param alpha_e a positive numeric. The shape parameter for the IG prior on
#'   \eqn{\sigma^2_e}.
#' @param beta_e a positive numeric. The scale parameter for the IG prior on
#'   \eqn{\sigma^2_e}.
#' @param alpha_s a positive numeric. The shape parameter for the IG prior on
#'   \eqn{\sigma^2_s}.
#' @param beta_s a positive numeric. The scale parameter for the IG prior on
#'   \eqn{\sigma^2_s}.
#'
#' @return Returns a new \code{lines} dataframe with the constants associated
#' with the prior term appended to the bottom.

addprior <- function(lines, alpha_e = 0, beta_e = 0, alpha_s = 0, beta_s = 0) {
    if ((alpha_e != 0) || (beta_e != 0)) {
        vert <- which(lines$slope == -Inf)
        lines$multiplier.log[vert] <- lines$multiplier.log[vert] + 2 * (alpha_e +
            1)
        lines$multiplier.inv[vert] <- lines$multiplier.inv[vert] + beta_e
        lines$v[vert] <- sqrt(lines$multiplier.inv[vert])
        lines$int.sigsqe[vert] <- lines$multiplier.inv[vert]/lines$multiplier.log[vert]
    }
    if ((alpha_s != 0) || (beta_s != 0)) {
        horiz <- data.frame(a = 1, multiplier.inv = 2 * beta_s, multiplier.log = 2 *
            (alpha_s + 1), v = sqrt(2 * beta_s), int.sigsqs = beta_s/(alpha_s + 1),
            int.sigsqe = NA, slope = 0, b = 0)
        lines <- rbind(lines, horiz)
    }
    return(lines)
}
