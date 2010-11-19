library(sos)
(fKD <- findFn('kernal density estimation'))

#----------- 
kde.contour <- function(x,y=NULL, conf=0.95, xdiv=100, ydiv=100, kernel.sd ) {
	xy <- xy.coords(x,y)
	xr <- range(xy$x)
	yr <- range(xy$y)
	
	xr <- xr + c(-1,1)*0.1*diff(xr)
	yr <- yr + c(-1,1)*0.1*diff(yr)
	
	if(missing(kernel.sd)) {
		kernel.sd <- c( diff(xr)/6, diff(yr)/6 )
	} else if (length(kernel.sd)==1) {
		kernel.sd <- rep(kernel.sd, 2)
	}
	
	xs <- seq(xr[1], xr[2], length.out=xdiv)
	ys <- seq(yr[1], yr[2], length.out=ydiv)
	mydf <- expand.grid( xx=xs, yy=ys )
	
	tmpfun <- function(xx,yy) {
		sum( dnorm(xx, xy$x, kernel.sd[1]) * dnorm(yy, xy$y, kernel.sd[2]) )
	}
	
	z <- mapply(tmpfun, xx=mydf$xx, yy=mydf$yy)
	
	sz <- sort(z, decreasing=TRUE)
	cz <- cumsum(sz)
	cz <- cz/cz[length(cz)]
	
	cutoff <- sz[ which( cz > conf )[1] ]
	
	plot(xy, xlab='x', ylab='y', xlim=xr, ylim=yr)
	#contour( xs, ys, matrix(z, nrow=xdiv), add=TRUE, col='blue')
	contour( xs, ys, matrix(z, nrow=xdiv), add=TRUE, col='red',
		levels=cutoff, labels='')
	
	invisible(NULL)
}

# test
kde.contour( rnorm(100), rnorm(100) )

# correlated data
my.xy <- MASS::mvrnorm(100, c(3,10), matrix( c(1,.8,.8,1), 2) )

kde.contour( my.xy, kernel.sd=.5 )

# compare to theoritical
lines(ellipse::ellipse( 0.8, scale=c(1,1), centre=c(3,10)), col='green')


# bimodal

new.xy <- rbind( MASS::mvrnorm(65, c(3,10), matrix( c(1,.6,.6,1),2) ),
	MASS::mvrnorm(35, c(6, 7), matrix( c(1,.6,.6,1), 2) ) )

kde.contour( new.xy, kernel.sd=.75 )
