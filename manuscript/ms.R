# dinvgamma: density function of the inverse-gamma distribution
dinvgamma = function(x, shape = 1, rate = 1, scale = 1/rate, log = FALSE) {
    # return( rate^shape / gamma(shape) * x^(-shape-1) * exp(-rate/x) )
    logval = shape*log(rate) - lgamma(shape) - (shape+1)*log(x) - rate/x
    if (log)
	return(logval)
    else
	return(exp(logval))
}

# pinvgamma: cumulative distribution function of the inverse-gamma distribution
pinvgamma = function(q, shape = 1, rate = 1, scale = 1/rate,
			lower.tail = TRUE, log.p = FALSE) {
    return( pgamma(1/q,shape,rate,scale,!lower.tail,log.p) )
}

# qinvgamma: quantile function of the inverse-gamma distribution
qinvgamma = function(p, shape = 1, rate = 1, scale = 1/rate,
			lower.tail = TRUE, log.p = FALSE) {
    return( 1 / qgamma(p,shape,rate,scale,!lower.tail,log.p) )
}

# rinvgamma: sample a vector of inverse-gamma random variables
rinvgamma = function(n, shape = 1, rate = 1, scale = 1/rate) {
    return( 1 / rgamma(n,shape,rate,scale) )
}


baseplot <- function(...){
  plot ( x=c(0,1.3), y=c(0,1.3), type="n", bty="l",
         xlab = expression(sigma[e]^2), ylab = expression(sigma[s]^2),
         line=2,
         xaxt="n", yaxt="n",
         ...
       )
}

line1 <- function ( label=TRUE, ... ){
  lines ( x=c(0,.4), y=c(.7,0), ... )
  if(label) text ( x=.23, y=.3, expression(j==1) )
}
line2 <- function ( label=TRUE, ... ){
  lines ( x=c(0,.1), y=c(1,0), ... )
  if(label) text ( x=.1, y=.1, expression(j==2) )
}
line3 <- function ( label=TRUE, ... ){
  lines ( x=c(0,1), y=c(.6,0), ... )
  if(label) text ( x=.5, y=.3, expression(j==3) )
}
line0 <- function ( label=TRUE, ... ){
  lines ( x=c(.6,.6), y=c(1,0), ... )
  if(label) text ( x=.6, y=.7, expression(j==0) )
}
linev <- function ( label=TRUE, ... ){
  lines ( x=c(.6,.6), y=c(1,0), ... )
  if(label) text ( x=.61, y=.4, expression(j==v) )
}
lineh <- function ( label=TRUE, ... ){
  lines ( x=c(0,1), y=c(.6,.6), ... )
  if(label) text ( x=.3, y=.62, expression(j==h) )
}


boxA <- function(){
  lines ( x=c(.1, .15, .15, .1, .1), y=c(.2, .2, .25, .25, .2) )
  text ( x=.125, y=.225, "A", cex=1.5)
}
boxB <- function(){
  lines ( x=c(.2, .25, .25, .2, .2), y=c(.5, .5, .55, .55, .5) )
  text ( x=.225, y=.525, "B", cex=1.5)
}
boxC <- function(){
  lines ( x=c(.3,.35,.35,.3,.3), y=c(.1,.1,.15,.15,.1) )
  text ( x=.325, y=.125, "C", cex=1.5)
}
boxb <- function(){
  lines ( x=c(.3,.35,.35,.3,.3), y=c(.1,.1,.15,.15,.1) )
  text ( x=.325, y=.125, expression(italic(b)), cex=1.5)
}

# Figure "oneline"
pdf("figs/oneline.pdf",width=5,height=5)
baseplot()
line1()
text ( x=0.1, y=0.2, "+", cex=2 )
text ( x=0.3, y=0.4, "-", cex=2 )
dev.off()

# Figure "multilines"
pdf("figs/multilines.pdf")
baseplot()
line1()
line2()
line3()
line0()
boxA()
boxB()
boxC()
dev.off()

#Figure "boundingbox"
pdf("figs/boundingbox.pdf")
baseplot()
line1(lty=3)
line2(lty=3)
line3(lty=3)
#line0(lty=3)
linev(lty=3)
lineh(lty=3)
# lines ( x=c(1,.6,.6,0,0,2/47,.1,1), y=c(0,.24,1,1,.6,27/47,0,0), lty=2 )
lines ( x=c(0,1,1,0,0), y=c(0,0,1,1,0), lwd=2 )
text(x=.2, y=.98, expression(B[1]))
text(x=.98, y=.2, expression(B[1]))
points ( x=0.5, y=1.2, pch=19 ); text ( x=.55, y=1.2, expression(bolditalic(p)[1]))
points ( x=0.5, y=1.0, pch=19 ); text ( x=.55, y=1.05, expression(bolditalic(p)[1]^"*"))
points ( x=1.05, y=0.5, pch=19 ); text ( x=1.1, y=0.5, expression(bolditalic(p)[2]))
points ( x=1.0, y=0.5, pch=19 ); text ( x=0.95, y=0.5, expression(bolditalic(p)[2]^"*"))
points ( x=1.05, y=1.1, pch=19 ); text ( x=1.05, y=1.15, expression(bolditalic(p)[3]))
points ( x=1, y=1, pch=19 ); text ( x=1.05, y=1, expression(bolditalic(p)[3]^"*"))
dev.off()

#Figure "smallboundingregion"
pdf("figs/smallboundingregion.pdf")
baseplot()
line1(lty=3)
line2(lty=3)
line3(lty=3)
linev(lty=3)
lineh(lty=3)
#lines ( x=c(1,.6,.6,0,0,2/47,.1,1), y=c(0,.24,1,1,.6,27/47,0,0), lwd=2 )
lines ( x=c(1,1,.6,.6,0,0,2/47,.1,1), y=c(0,.6,.6,1,1,.6,27/47,0,0), lwd=2 )
lines ( x=c(0,1,1,0,0), y=c(0,0,1,1,0), lty=2 )
dev.off()


#Figure "boundingbox2"
pdf("figs/boundingbox2.pdf")
baseplot()
line1(lty=3,label=FALSE)
line2(lty=3,label=FALSE)
line3(lty=3,label=FALSE)
line0(lty=3,label=FALSE)
lines ( x=c(0,1,1,0,0), y=c(0,0,1,1,0), lwd=1 )
text (x=.97, y=.5, expression(B[1]))
text (x=.3, y=.97, expression(B[1]))
points ( x=0.25, y=.7, pch=19 ); text ( x=.25, y=.75, expression(bolditalic(p)[4]))
lines (x=c(0,1.2,1.2,0,0), y=c(1.2,1.2,0,0,1.2), lwd=2)
text (x=1.17, y=.5, expression(B))
text (x=.3, y=1.17, expression(B))
axis(side=1,at=1.2,labels=expression(sigma[e]^"2*"))
axis(side=2,at=1.2,labels=expression(sigma[s]^"2*"))
lines(x=c(1.2,1.2),y=c(1.2,1.3),lty=2)
lines ( x=c(.6,.6), y=c(1,1.3), lty=2)
lines(x=c(0,0),y=c(1.2,1.3),lty=2)
points(x=1.25,y=1.25,pch=19); text(x=1.25,y=1.25,expression(bolditalic(q)[1]),pos=4)
points(x=1.2,y=0,pch=19); text(x=1.2,y=0,expression(bolditalic(q)[1]^"*"),pos=4)
#points(x=1,y=0,pch=19); text(x=1,y=0.05,expression(bolditalic(q)[1]^"**"),pos=4)
points (x=.8, y=1.25, pch=19); text (x=.8, y=1.25, expression(bolditalic(q)[2]),pos=4)
points (x=.6, y=1.2, pch=19); text(x=.6,y=1.225, expression(bolditalic(q)[2]^"*"),pos=4)
#points (x=.6, y=1, pch=19); text(x=.6,y=1.025, expression(bolditalic(q)[2]^"**"),pos=4)
points(x=.2,y=1.25,pch=19); text(x=.2, y=1.25, expression(bolditalic(q)[3]),pos=4)
points(x=.0,y=1.2,pch=19); text(x=0, y=1.225, expression(bolditalic(q)[3]^"*"),pos=4)
#points(x=.0,y=1,pch=19); text(x=0, y=1.025, expression(bolditalic(q)[3]^"**"),pos=4)
dev.off()

# Figure "abcboxes"
pdf("figs/abcboxes.pdf")
baseplot()
line1()
boxA()
boxB()
boxC()
dev.off()

# Figure "boxb"
pdf("figs/boxb.pdf")
baseplot()
lines ( x=c(0,1,1,0,0), y=c(0,0,1,1,0), lwd=2 )
text(x=.2, y=.98, expression(B[1]))
text(x=.98, y=.2, expression(B[1]))
line1(lty=3)
line2(lty=3)
line3(lty=3)
boxb()
dev.off()


#Figure "contour box"
pdf("figs/contourbox.pdf")
baseplot(usr=c(0,1.2,0,1.2))
line1(label=F,lty=2)
line2(label=F,lty=2)
line3(label=F,lty=2)
line0(label=F,lty=2)
arrows ( 1, 0, 1, 1.04 ); text ( x=.98, y=.5, expression(bolditalic(l)[1]))
points ( x=1, y=0, pch=19 ); text ( x=.97, y=0, expression(bolditalic(p)[1]))
arrows ( 0, 1, 1, 1 ); text ( x=.5, y=.98, expression(bolditalic(l)[2]))
points ( x=0, y=1, pch=19 ); text ( x=.022, y=.98, expression(bolditalic(p)[2]))
points ( .6, 1, pch=19 ); text ( x=.62, y=.98, expression(bolditalic(p)[3]))
dev.off()
