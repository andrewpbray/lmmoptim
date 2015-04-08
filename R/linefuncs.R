# Given X, Z, Y, find the constants in Hodges' reexpression

findlines <- function(x, z, y, SigE, SigS) {
    # Find the constants in Hodges' reexpression.  Just those constants related
    # to the model, not the prior.  sanity checks
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


addprior <- function(lines, a.E = 0, b.E = 0, a.S = 0, b.S = 0) {
    if ((a.E != 0) || (b.E != 0)) {
        # if sigesq ~ IG(a.E,b.E)
        vert <- which(lines$slope == -Inf)
        lines$multiplier.log[vert] <- lines$multiplier.log[vert] + 2 * (a.E +
            1)
        lines$multiplier.inv[vert] <- lines$multiplier.inv[vert] + b.E
        lines$v[vert] <- sqrt(lines$multiplier.inv[vert])
        lines$int.sigsqe[vert] <- lines$multiplier.inv[vert]/lines$multiplier.log[vert]
    }
    if ((a.S != 0) || (b.S != 0)) {
        # if sigssq ~ IG(a.S,b.S)
        horiz <- data.frame(a = 1, multiplier.inv = 2 * b.S, multiplier.log = 2 *
            (a.S + 1), v = sqrt(2 * b.S), int.sigsqs = b.S/(a.S + 1), int.sigsqe = NA,
            slope = 0, b = 0)
        lines <- rbind(lines, horiz)
    }
    return(lines)
}
