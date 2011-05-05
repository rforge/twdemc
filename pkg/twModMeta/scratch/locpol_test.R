install.packages("locpol")

library("locpol")
?locpol

example(locpol)

# seems that only the f = f(x) with x one-dimensional is supported.

#------- two dimensional example
# bimodal
new.xy <- rbind( MASS::mvrnorm(65, c(3,10), matrix( c(1,.6,.6,1),2) ),
	MASS::mvrnorm(35, c(6, 7), matrix( c(1,.6,.6,1), 2) ) )
#colnames(new.xy) <- c("x","y")
locpol( )

#--------- seven dimensions
library(mvtnorm)
help(package="mvtnorm", help_type = "html")
?dmvnorm
n=7
new.dat <- rbind( rmvnorm(400, rep(0,n), diag(1,n) ),
	MASS::mvrnorm(200, rep(1,n), diag(0.5,n)) )
bw2 <- npudensbw(dat=new.dat)
npplot(bws=bw)

function (formula, data, weig = rep(1, nrow(data)), bw = NULL, 
	kernel = EpaK, deg = 1, xeval = NULL, xevalLen = 100) 
{
	stopifnot(nrow(data) == length(weig))
	res <- list()
	res$bw <- bw
	res$KName <- match.call()
	res$kernel <- kernel
	res$deg <- deg
	res$xeval <- xeval
	res$mf <- model.frame(formula, data)
	datCla <- attr(attr(res$mf, "terms"), "dataClasses")
	varNames <- names(datCla)[datCla == "numeric"]
	stopifnot(length(varNames) == 2)
	res$Y <- varNames[1]
	res$X <- varNames[2]
	xo <- order(res$mf[, res$X])
	res$mf <- res$mf[xo, ]
	res$weig <- weig[xo]
	if (is.null(xeval)) 
		res$xeval <- seq(min(res$mf[, res$X]), max(res$mf[, res$X]), 
			len = xevalLen)
	else res$xeval <- sort(xeval)
	if (is.null(res$bw)) 
		res$bw <- regCVBwSelC(data[, res$X], data[, res$Y], res$deg, 
			res$kernel, res$weig)
	res$lpFit <- locPolSmootherC(res$mf[, res$X], res$mf[, res$Y], 
		res$xeval, res$bw, res$deg, res$kernel, DET = TRUE, res$weig)
	names(res$lpFit)[] <- c(res$X, res$Y, paste(res$Y, 1:deg, 
			sep = ""), "xDen")
	res$lpFit$xDen <- res$lpFit$xDen^(1/(deg + 1))/(nrow(data) * 
			res$bw)
	nu <- 0
	res$CIwidth <- computeRK(equivKernel(kernel, nu, deg), lower = dom(res$kernel)[[1]], 
		upper = dom(res$kernel)[[2]], subdivisions = 25) * factorial(nu)^2
	res$CIwidth <- res$CIwidth/(nrow(data) * res$bw)
	res$residuals <- res$mf[, res$Y] - locLinSmootherC(res$mf[, 
			res$X], res$mf[, res$Y], res$mf[, res$X], res$bw, res$kernel, 
		res$weig)$beta0
	res$lpFit$var <- locCteSmootherC(res$mf[, res$X], res$residuals^2, 
		res$xeval, 1.2 * res$bw, res$kernel, res$weig)$beta0
	class(res) <- "locpol"
	return(res)
}
B

