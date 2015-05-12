# global mean surface temperature

library(nlme) # for the standard analysis, by REML
library(devtools)
devtools::install_github("andrewpbray/lmmoptim")
library("lmmoptim")

data(gmst)
y <- gmst$temp.dev

# plot the data
pdf("gmst-scatter.pdf")
p <- ggplot(data = gmst, aes(x = Year, y = temp.dev)) + geom_point()
p + ylab("temperature deviation") + theme_bw()
dev.off()

# Hodges' code for centered X, Z, GamX, GamZ
x <- 1:125
xs <- (x - mean(x))/sqrt(var(x))
X <- matrix(c(rep(1, 125), xs, xs^2), 125, 3)
knots <- seq(4, 121, 4)
ks <- (knots - mean(x))/sqrt(var(x))
kp <- matrix(0 == 1, 125, 30)
for (i in 1:30) {
  for (j in 1:125) {
    kp[j, i] <- xs[j] >= ks[i]
    }
  }
Z <- kp * (matrix(xs, 125, 30) - matrix(ks, 125, 30, byrow = T))^2
GX <- svd(X)$u
CX <- t(svd(cbind(X,Z))$u) %*% GX
QCX <- diag(rep(1, 33)) - CX %*% solve(t(CX) %*% CX) %*% t(CX)
CZ <- eigen(QCX, sym = T)$vec[, 1:30]
GZ <- svd(cbind(X, Z))$u %*% CZ
x <- X
z <- Z

# Lavine's code for uncentered  X, Z, GamX, GamZ
yr <- gmst$year
nyr <- length(yr)
knots <- 1880 + seq(4, 120, by = 4)
x <- cbind(rep(1, nyr), yr, yr^2)
# x <- cbind ( rep(1,nyr), scale ( x[,-1] ) ) # either use either this line or not.  it shouldn't matter
z <- pmax(outer(yr, knots, "-" ), 0)^2

# Compute H_ZZ and PDP';  Compute D and canon regressors GZ P

# HZZ <- t(GZ) %*% Z
# P <- eigen(HZZ %*% t(HZZ),sym=T)$vec
# D <- eigen(HZZ %*% t(HZZ),sym=T)$val

# D matches previous work!

# > D
#  [1] 3.596518e+01 3.151748e+00 5.615607e-01 1.471589e-01 4.926579e-02 1.953047e-02 8.762126e-03
#  [8] 4.320602e-03 2.295369e-03 1.295259e-03 7.682428e-04 4.751401e-04 3.045424e-04 2.013100e-04
# [15] 1.367069e-04 9.507613e-05 6.755040e-05 4.893385e-05 3.608908e-05 2.707027e-05 2.064154e-05
# [22] 1.600112e-05 1.261906e-05 1.013951e-05 8.320570e-06 6.996742e-06 6.055456e-06 5.422789e-06
# [29] 5.055363e-06 3.749042e-06

# CR <- GZ %*% P
# y <- as.matrix(c( -13, -1, -5, -42, -23, -26, -46, -23, 6, -20, -56, -40, -39, -32, -32, -27, -15, -20, -25, -6, -5, -30, -36, -42, -25, -15, -40, -30, -30, -20, -25, -33, -28, -2, 5, -20, -46, -34, -9, -18, -4, -9, -15, -11, -15, 5, -4, 1, -22, -3, 2, 4, -11, 5, -8, 1, 12, 15, -2, 14, 11, 10, 6, 10, -1, 1, 12, -3, -9, -17, -2, 3, 12, -9, -8, -18, 8, 9, 5, -2, 10, 5, 3, -25, -15, -8, -2, -9, 0, 4, -10, -5, 18, -6, -2, -21, 16, 7, 14, 28, 40, 9, 34, 15, 13, 19, 35, 39, 26, 48, 44, 15, 19, 32, 47, 39, 41, 72, 46, 42, 58, 69, 68, 61, 77))

# tp <- cbind(GX,GZ)
# GC <- eigen(diag(rep(1,125)) - tp %*% solve(t(tp) %*% tp) %*% t(tp),sym=T)$vec[,1:92]
# RSS <- sum((t(GC) %*% y)^2)
# hv <- t(P) %*% t(GZ) %*% y

###################################################



# y <- c( -13, -1, -5, -42, -23, -26, -46, -23, 6, -20, -56, -40, -39, -32, -32, -27, -15, -20, -25, -6, -5, -30, -36, -42, -25, -15, -40, -30, -30, -20, -25, -33, -28, -2, 5, -20, -46, -34, -9, -18, -4, -9, -15, -11, -15, 5, -4, 1, -22, -3, 2, 4, -11, 5, -8, 1, 12, 15, -2, 14, 11, 10, 6, 10, -1, 1, 12, -3, -9, -17, -2, 3, 12, -9, -8, -18, 8, 9, 5, -2, 10, 5, 3, -25, -15, -8, -2, -9, 0, 4, -10, -5, 18, -6, -2, -21, 16, 7, 14, 28, 40, 9, 34, 15, 13, 19, 35, 39, 26, 48, 44, 15, 19, 32, 47, 39, 41, 72, 46, 42, 58, 69, 68, 61, 77)

# the standard reml analysis
# df <- data.frame ( cbind(y,x,z) )
# names(df) <- c ( "y", paste("x",1:3,sep=""), paste("z",1:30,sep="") )
# mod <- lme ( y ~ x1 + x2 + x3 + z1, data=df, random = ~1|state )
# summary(mod)

# find lines
y <- gmst$temp.dev
SigE <- diag(length(y))
SigS <- diag(length(knots))
lines <- findlines(x, z, y, SigE, SigS)

# plot the lines
pdf("gmst-lines.pdf")
showlines(lines)
dev.off()

# To find the MRLE the initial box is determined by the maxima of the x- and y-intercepts

# fit the model
boxes <- fitlmm(lines = lines, eps = 10, delE = 1, delS = 1, M = 10, maxit = 15) # first try

# view the boxes
pdf("gmst-boxes1.pdf")
p <- showboxes(boxes)
p + scale_x_continuous(name = expression(sigma[e]^2)) +
    scale_y_continuous(name = expression(sigma[s]^2)) +
    theme_bw()
dev.off()

p <- showfunc(boxes)
pdf("gmst-rll1.pdf")
p + scale_x_log10(name = expression(sigma[e]^2),
                  breaks = c(30, 300, 3000)) +
    scale_y_log10(name = expression(sigma[s]^2),
                  breaks = c(10, 1000, 100000, 10000000)) +
    theme_bw()
dev.off()

# after the first try we can narrow down the search region
startbox <- makebox(lines = lines,
	                  lims.sigsqe = c(50, 500),
                    lims.sigsqs = c(0, 1000000),
	                  status = rep("straddle", nrow(lines)))

boxes <- fitlmm(lines = lines, startbox = startbox, eps = 2, M = 10, maxit = 20,
                delE = 0, delS = 0)

# view the boxes
jpeg("gmst-boxes2.jpg")
p <- showboxes(boxes)
p + scale_x_continuous(name = expression(sigma[e]^2)) +
    scale_y_continuous(name = expression(sigma[s]^2)) +
    theme_bw()
dev.off()

jpeg("gmst-rll2.jpg")
p <- showfunc(boxes)
p + scale_x_log10(name = expression(sigma[e]^2),
                  breaks = c(30, 100, 300)) +
    scale_y_log10(name = expression(sigma[s]^2),
                  breaks = c(1, 10, 100, 1000, 10000, 100000)) +
    theme_bw()
dev.off()
