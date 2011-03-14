# TODO: Add comment
# 
# Author: twutz
###############################################################################

windows(width=3.4,height=2.8, pointsize=10, record=TRUE)
par( las=1 )					#also y axis labels horizontal
par(mar=c(2.2,2.3,0,0)+0.3 )  #margins
par(tck=0.02 )				#axe-tick length inside plots             
par(mgp=c(1.4,0.2,0) )  #positioning of axis title, axis labels, axis


a=0.5
b=0.05
x=0:4
xClose <- seq(min(x),max(x),length.out=30)
fModel <- function(x,a,b){a+b*exp(x)}
yTrue = fModel(x,a,b)
sd1 = 0.2

yObs <- yTrue + rnorm(length(x), sd=sd1)
plot( 1, xlim=c(0,4.6), ylim=c(0,max(yTrue)+0.5), ylab="", xlab="x", type="n" )
mtext("y",2,1,las=1)
lines( fModel(xClose,a,b) ~ xClose, col="gray" )
points( yObs ~ x)

yres0 <- seq(-3*sd1,+3*sd1, length.out=30)
xres0d <- dnorm(yres0,sd=sd1)	#scaled so that maximum is 0.5
xres0 <- xres0d*0.5/max(xres0d)
#plot(xres0~yres0,ylim=c(0,0.6), type="l")
#plot(yres0~xres0,xlim=c(0,0.6), type="l")

i=2
for( i in seq_along(x) ){
	yi <- yTrue[i]
	xi <- x[i] 
	yres = yi+yres0
	xres = xi+xres0
	lines(yres~xres)
	lines(range(yres)~rep(xi,2),col="lightgray")
}


#--------- logNormal multiplicative error
yObs <- yTrue * rlnorm(length(x), meanlog=sd1^2/2, sdlog=sd1)
plot( 1, xlim=c(0,4.6), ylim=c(0,max(yTrue)+2.5), ylab="", xlab="x", type="n" )
mtext("y",2,1,las=1)
lines( fModel(xClose,a,b) ~ xClose, col="gray" )
points( yObs ~ x)

yres0 <- seq(1/2,1*2, length.out=30)
xres0d <- dlnorm(yres0,meanlog=sd1^2/2, sdlog=sd1)	#scaled so that maximum is 0.5
xres0 <- xres0d*0.5/max(xres0d)
#plot(xres0~yres0, type="l")
#plot(yres0~xres0,xlim=c(0,0.6), type="l")

i=2
for( i in seq_along(x) ){
	yi <- yTrue[i]
	xi <- x[i] 
	yres = yi*yres0
	xres = xi+xres0
	lines(yres~xres)
	lines(range(yres)~rep(xi,2),col="lightgray")
}

#--------- logNormal additive error
# Y ~ lognormal( f(x), sigma^2 )

yObs <- rlnorm(length(x), meanlog=log(yTrue)-sd1^2/2, sdlog=sd1)
plot( 1, xlim=c(0,4.6), ylim=c(0,max(yTrue)+2.5), ylab="", xlab="x", type="n" )
mtext("y",2,1,las=1)
lines( fModel(xClose,a,b) ~ xClose, col="gray" )
points( yObs ~ x)

i=2
for( i in seq_along(x) ){
	yi <- yTrue[i]
	xi <- x[i]
	#yr <- qlnorm(c(0.025,0.975), meanlog=log(yi)-sd1^2/2, sdlog=sd1 )
	yr <- qlnorm(c(0.025/10,1-(0.025/10)), meanlog=log(yi)-sd1^2/2, sdlog=sd1 )
	yres <- seq(yr[1],yr[2],length.out=30)
	xres0d <- dlnorm(yres, meanlog=log(yi)-sd1^2/2, sdlog=sd1 )
	xres0 <- xres0d*0.5/max(xres0d)
	xres <- x[i] + xres0
	lines(yres~xres)
	lines(range(yres)~rep(xi,2),col="lightgray")
}

#acutally the same as multiplicative error



