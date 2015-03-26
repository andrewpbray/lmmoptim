library(ggplot2)
library(dplyr)
library(nlme)
library(microbenchmark)
library(devtools)
#devtools::install_github("andrewpbray/lmmoptim", ref = "a773c8dfbf54d9f388801c9fa9ffb4e49b92774b")
devtools::install_github("andrewpbray/lmmoptim")
library("lmmoptim")

data(hmo)
data(hmo_states)
hmobig <- left_join(hmo, hmo_states, by="state")

y <- hmobig$indPrem
n <- length(y)
mod <- lm(indPrem ~ expPerAdm + (region=="NE") + state, data=hmobig, x=TRUE )
x <- mod$x[,1:3]
z <- mod$x[,-(1:3),drop=FALSE]
SigE <- diag(n)
SigS <- diag(44)

# the standard reml analysis
lines <- findlines(x, z, y, SigE, SigS)

mlrebox <- with(lines,
                  makebox(lines = lines,
                          lims.sigsqs = c(0, max(int.sigsqs[is.finite(int.sigsqs)])),
                          lims.sigsqe = c(0, max(int.sigsqe[is.finite(int.sigsqe)])),
                          status = rep("straddle", nrow(lines))
                  )
)

####################
boxes.HH11 <- fitlmm ( lines=lines, mlrebox, eps=5, delE=log(10), delS=log(10),
                       ratio=TRUE, M=5, maxit=10, lognote = "First ME model v2")
####################

# p <- showboxes(boxes.HH11)
# pdf ( "development/hmo-HH11-boxes.pdf" )
# p + scale_x_continuous ( name = expression(sigma[e]^2) ) +
#   scale_y_continuous ( name = expression(sigma[s]^2) ) +
#   theme_bw()
# dev.off()
#
# p <- showfunc(boxes.HH11)
# pdf ( "development/hmo-HH11-rll.pdf" )
# p + scale_x_log10 ( name = expression(sigma[e]^2),
#                     breaks = c (100,300,1000,3000,10000)
# ) +
#   scale_y_log10 ( name = expression(sigma[s]^2),
#                   breaks = c (30,100,300,1000,3000,10000)
#   ) +
#   theme_bw()
# dev.off()

box2 <- with(lines,
             makebox(lines = lines,
                     lims.sigsqe = c(300, 1000 ),
                     lims.sigsqs = c(0, 1000 ),
                     status = rep("straddle", nrow(lines))))

####################
boxes.HH11 <- fitlmm(lines=lines, box2, eps=1, M=10, maxit=7, lognote = "Refined ME model v2")
####################

# p <- showboxes (boxes.HH11)
# pdf ( "development/hmo-HH11-boxes2.pdf" )
# p + scale_x_continuous ( name = expression(sigma[e]^2), limits=c(300,1000) ) +
#   scale_y_continuous ( name = expression(sigma[s]^2) ) +
#   theme_bw()
# dev.off()
#
# p <- showfunc(boxes.HH11)
# pdf ( "development/hmo-HH11-rll2.pdf" )
# p + scale_x_log10 ( name = expression(sigma[e]^2),
#                     breaks = c (300,500,1000)
# ) +
#   scale_y_log10 ( name = expression(sigma[s]^2),
#                   breaks = c (1, 3, 10, 30,100,300,1000)
#   ) +
#   theme_bw()
# dev.off()

# Bayesian analysis using Hodges 1998 prior
lines <- addprior ( lines, a.E=1, b.E=0, a.S=1.1, b.=.1 )
# pdf("hmolines.HH11.Bayes.pdf")
# showlines(lines)
# dev.off()

box2 <- with(lines,
             makebox(lines = lines,
                     lims.sigsqe = c(300, 1000 ),
                     lims.sigsqs = c(0, 1000 ),
                     status = rep("straddle", nrow(lines))))

boxes.HH11Bayes <- fitlmm ( lines=lines, box2, eps=1, M=10, maxit=6, lognote = "bayes v2")

# p <- showboxes (boxes.HH11Bayes)
# jpeg ( "hmo_HH11Bayes_boxes.jpg" )
# p + scale_x_continuous ( name = expression(sigma[e]^2), limits=c(300,1000) ) +
#   scale_y_continuous ( name = expression(sigma[s]^2) ) +
#   theme_bw()
# dev.off()
#
# p <- showfunc(boxes.HH11Bayes)
# jpeg ( "hmo_HH11Bayes_rll.jpg" )
# p + scale_x_log10 ( name = expression(sigma[e]^2),
#                     breaks = c (300,500,1000)
# ) +
#   scale_y_log10 ( name = expression(sigma[s]^2),
#                   breaks = 10^seq(-3,3,by=1)
#   ) +
#   theme_bw()
# dev.off()

# The figure reveals a ridge.  Let's investigate
ridge <- with ( boxes.HH11Bayes,
                rll.upper > max(rll.lower) & (rll.upper-rll.lower) > 5
)
with ( boxes.HH11Bayes, range(sigsqs.lo[good]) )
with ( boxes.HH11Bayes, range(sigsqs.hi[good]) )
with ( boxes.HH11Bayes, range(rll.upper[good]) )
with ( boxes.HH11Bayes, range(rll.lower[good]) )
with ( boxes.HH11Bayes, range(rll.upper[good]-rll.lower[good]) )

tmp <- seq(0,.1,length=1000)
tmp2 <- dinvgamma(tmp,shape=1.1,rate=0.1)
plot(tmp,tmp2,type="l")
abline (v=0.04761905)




# testing a new way to initialize Lb

linemaxes <- with(lines, ifelse(is.na(int.sigsqe),
                                -0.5 * (multiplier.log * log(int.sigsqs) +
                                          multiplier.inv/int.sigsqs),
                                -0.5 * (multiplier.log * log(int.sigsqe) +
                                          multiplier.inv/int.sigsqe)))

# try visualizing the contribution of a single term
fterm <- function(x, y, lines, which.line) {
  linearterm <- lines[which.line, "a"] * y + x
  ifelse(linearterm == 0, NA,
         -0.5 * (lines[which.line, "multiplier.log"] * log(linearterm) +
                   lines[which.line, "multiplier.inv"]/linearterm))
}

plotTerm <- function(lines, which.line = 1, npix = 100) {
  x <- seq(0, lines[which.line, "int.sigsqe"], length.out = npix)
  y <- seq(0, lines[which.line, "int.sigsqs"], length.out = npix)
  z <- outer(x, y, FUN = fterm, lines = lines, which.line = which.line)

  par(mfrow = c(1, 2))
  image(x, y, z, zlim = c(.05 * min(as.numeric(z), na.rm = T),
                          max(as.numeric(z), na.rm = T)))
  lines(c(lines[which.line, "int.sigsqe"], 0),
        c(0, lines[which.line, "int.sigsqs"]), col = "blue", lwd = 2)
  hist(z, main = "Histogram of f")
}

plotTerm(lines, which.line = 1, npix = 200)

