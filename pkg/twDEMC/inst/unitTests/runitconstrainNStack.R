.setUp <- function(){
	data(twdemcEx1)
	data(twLinreg1)
	suppressWarnings( rm( list=names(twLinreg1) ) )
	suppressWarnings( rm( list=names(twdemcEx1) ) )
	attach( twLinreg1 )
	attach( twdemcEx1 )
}

.tearDown <- function(){
	detach( twdemcEx1 )
	detach( twLinreg1 )
}

test.allVars <- function(){
	#str(twdemcEx1 )
	.pss1 <- stackChains(twdemcEx1)
	.thetaPrior <- thetaTrue
	.invCovarTheta <- invCovarTheta
	.n <- nrow(.pss1)%/%4
	.res1 <- constrainNStack( .pss1, .thetaPrior, n=.n )
	checkTrue( is.matrix(.res1) )
	#checkEquals( length(.thetaPrior), ncol(.res1))	#.thetaPrior only to constrain output
	checkEquals( ncol(.pss1), ncol(.res1))
	checkEquals( .n, nrow(.res1))
	#
	.res2 <- constrainNStack( .pss1, .thetaPrior, invCovarTheta=.invCovarTheta, n=.n )
	checkTrue( is.matrix(.res2) )
	#checkEquals( length(.thetaPrior), ncol(.res1))	#.thetaPrior only to constrain output
	checkEquals( ncol(.pss1), ncol(.res2))
	checkEquals( .n, nrow(.res2))
}
#mtrace(test.allVars); {.setUp(); test.allVars(); .tearDown() }


test.constrainCfStack <- function(){
	#str(twdemcEx1 )
	.pss1 <- stackChains(twdemcEx1)
	.thetaPrior <- thetaTrue
	.invCovarTheta <- invCovarTheta
	.alpha <- 0.7
	.res1 <- constrainCfStack( .pss1, .thetaPrior, alpha=.alpha )
	checkTrue( is.matrix(.res1) )
	#checkEquals( length(.thetaPrior), ncol(.res1))	#.thetaPrior only to constrain output
	checkEquals( ncol(.pss1), ncol(.res1))
	
	.res2 <- constrainCfStack( .pss1, .thetaPrior, invCovarTheta=.invCovarTheta, alpha=.alpha )
	checkTrue( is.matrix(.res2) )
	#checkEquals( length(.thetaPrior), ncol(.res1))	#.thetaPrior only to constrain output
	checkEquals( ncol(.pss1), ncol(.res2))
}
#mtrace(constrainCfStack)
#mtrace(test.allVars); {.setUp(); test.allVars(); .tearDown() }


.tmp.f <- function(){
	library("KernSmooth")
	.d <- 	bkde2D(.pss1[,-1], bandwidth=0.2)
	names(.d) <- c("x","y","z")
	contour(.d)
	points( thetaTrue[1], thetaTrue[2] )
	
	.dres <- 	bkde2D(.res1[,-1], bandwidth=0.2)
	names(.dres) <- c("x","y","z")
	contour(.dres, add=TRUE, col="red")
	
	.dres2 <- 	bkde2D(.res2[,-1], bandwidth=0.2)
	names(.dres2) <- c("x","y","z")
	contour(.dres2, add=TRUE, col="blue")
	
	plot(density(.pss1[,"a"]))
	lines(density(.res1[,"a"]), col="red")
}


