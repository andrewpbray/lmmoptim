sqrt.m <- function(m) {
    # m should be a square-rootable matrix
    m.eig <- eigen(m)
    tmp <- m.eig$vectors %*% diag(sqrt(m.eig$values)) %*% t(m.eig$vectors)
    return(tmp)
}
inv.sqrt <- function(m) {
    # m should be an invertible matrix
    m.eig <- eigen(m)
    tmp <- m.eig$vectors %*% diag(1/sqrt(m.eig$values)) %*% t(m.eig$vectors)
    return(tmp)
}

# orthogonal projection operator onto the column space of matrix x
proj <- function(x) {
    x %*% solve(crossprod(x)) %*% t(x)
}

# calculate RLL from the reexpression
RLL <- function(sigsqe, sigsqs) {
    terms.lin <- with(lines, a * sigsqs + sigsqe)
    return(sum(with(lines, -0.5 * (multiplier.log * log(terms.lin) + multiplier.inv/terms.lin))))
}

# calculate RLL from Hodges' (1.16)
RLL2 <- function(x, y, z, SigE, SigS, sigsqe, sigsqs) {
    V <- z %*% (sigsqs * SigS) %*% t(z) + sigsqe * SigE
    Vinv <- solve(V)
    xVinv <- crossprod(x, Vinv)
    -0.5 * (det(V, log = TRUE) + det(xVinv %*% x, log = TRUE) + t(y) %*% (Vinv - 
        t(xVinv) %*% solve(xVinv %*% x) %*% xVinv) %*% y)
}

whichbox <- function(sigsqe, sigsqs, boxes) {
    # Which box contains ( sigsqe, sigsqs )
    tmp <- t(sapply(X = boxes, FUN = function(box) {
        c(box$lims.sigsqe, box$lims.sigsqs)
    }))
    return(which(sigsqe > tmp[, 1] & sigsqe < tmp[, 2] & sigsqs > tmp[, 3] & 
        sigsqs < tmp[, 4]))
} 
