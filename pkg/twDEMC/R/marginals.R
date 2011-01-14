.cut2dim2 <- function(x,y,c,xn=round(sqrt(length(c))),yn=round(sqrt(length(c))), xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), fun=median ){
	# cut2dim2
	#
	#matrix mean c over cells of unregular grid over x and y
	#grid is defined by classes of equal number of observations in x and y
	#xn and yn number of classes in x and y dimension
	xcuts=cutQuantiles(x,g=xn, levels.mean=TRUE)
	ycuts=cutQuantiles(y,g=yn, levels.mean=TRUE)
	tmp <- list(xcuts,ycuts); names(tmp) <- c(xlab,ylab)
	cagg <- tapply(c, tmp, fun)
	list( cagg=cagg, xmeans=as.numeric(levels(xcuts)), ymeans=as.numeric(levels(ycuts)) )
}
.cut2dim3 <- function(x,y,z,c,xn=round(length(c)^(1/3)),yn=xn, zn=xn, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), fun=mean ){
	#matrix mean c over cells of unregular grid over x,y, and z
	#grid is defined by classes of equal number of observations in x and y
	#xn and yn number of classes in x and y dimension
	xcuts=cutQuantiles(x,g=xn, levels.mean=TRUE)
	ycuts=cutQuantiles(y,g=yn, levels.mean=TRUE)
	zcuts=cutQuantiles(z,g=zn, levels.mean=TRUE)
	tmp <- list(xcuts,ycuts,zcuts); names(tmp) <- c(xlab,ylab,zlab)
	cmean <- tapply(c, tmp, fun)
	list( cagg=cmean
		, xmeans=as.numeric(levels(xcuts))
		, ymeans=as.numeric(levels(ycuts)) 
		, zmeans=as.numeric(levels(zcuts)) 
	)
}

.colNames <- function(x,index){
	if( is.numeric(index) ) colnames(x)[index] else index
}
	
marginals1d <- function(
	### integrate over all others but 1 column
	x					##<<  numeric matrix
	,vars=2				##<< index (integer or string): column to define integration boxes
	,intCol=1			##<< index (integer or string): column to average
	,n=round(nrow(x)/12)##<< number of groups, defaults to about 10 observations per group
	,dimnames= .colNames(x,vars)	##<< dimension names of the result array
	,fAgg=mean			##<< function applied over the subsets of x[,intCol]
	,...				##<< further arguments to fAgg
){
	##seealso<<  
	## \code{\link{twDEMCInt}}
	
	##details<< 
	## There are several methods to get calculate marginal densities for multivariate sample. \itemize{
	## \item{ integrate over all others but 1 column: this method  } 
	## \item{ integrate over all others but 2 column: \code{\link{marginals2d}}  } 
	## \item{ integrate over all others but 3 column: \code{\link{marginals3d}}  } 
	## \item{ plotting a 2D marginal distribution: \code{\link{plotMarginal2D}}  } 
	## \item{ plotting a 2D kernel density estimate: \code{\link{plot2DKDMarginals}}  } 
	##}
	
	xcuts <- cutQuantiles(x[,vars],g=n, levels.mean=TRUE)
	tmp <- list(xcuts); names(tmp) <- dimnames	
	#cagg <- tapply( x[,intCol], tmp, fAgg )
	cagg <- tapply( x[,intCol], tmp, fAgg, ...)
	### numeric vector with each item representing application of fAgg to x[,intCol]
	### and dimnames equal to means of the classes
}
attr(marginals1d,"ex") <- function(){
	res <- marginals1d(sample,vars="kO",n=20)
	plot( res ~ as.numeric(names(res)) )
}

marginals2d <- function(
	### integrate over all others but 2 columns
	x					##<<  numeric matrix
	,vars=c(2,3)		##<< index (integer or string): column to define integration boxes
	,intCol=1			##<< index (integer or string): column to average
	,n=round(sqrt(nrow(x))/12) ##<< numeric vector: number of groups per variable
	,dimnames= .colNames(x,vars)	##<< dimension names of the result array
	,fAgg=mean			##<< function applied over the subsets of x[,intCol]
	,...				##<< further arguments to fAgg
){
	##seealso<<  
	## \code{\link{marginals1d}}
	## \code{\link{twDEMCInt}}
	
	if( length(n)==1 ) n=rep(n,2)
	xcuts <- cutQuantiles(x[,vars[1] ],g=n[1], levels.mean=TRUE)
	ycuts <- cutQuantiles(x[,vars[2] ],g=n[2], levels.mean=TRUE)
	tmp <- list(xcuts,ycuts); names(tmp) <- dimnames	
	#cagg <- tapply( x[,intCol], tmp, fAgg )
	cagg <- tapply( x[,intCol], tmp, fAgg, ...)
	### numeric matrix with each item representing application of fAgg to x[,intCol]
	### and dimnames equal to means of the classes
}
attr(marginals2d,"ex") <- function(){
	res <- marginals2d(sample,vars=c("kY","kO"),n=10)
	tmp.f <- function(){
		library(lattice)
		levelplot( value~kY*kO, data=melt(res), col.regions=rev(heat.colors(100)), xlab=names(dimnames(res))[1], ylab=names(dimnames(res))[2] )
	}
}

marginals3d <- function(
	### integrate over all others but 3 columns
	x					##<<  numeric matrix
	,vars=c(2,3,4)		##<< index (integer or string): column to define integration boxes
	,intCol=1			##<< index (integer or string): column to average
	,n=round(nrow(x)^(1/3)/12) ##<< numeric vector: number of groups per variable
	,dimnames= .colNames(x,vars)	##<< dimension names of the result array
	,fAgg=mean			##<< function applied over the subsets of x[,intCol]
	,...				##<< further arguments to fAgg
){
	##seealso<<  
	## \code{\link{marginals1d}}
	## \code{\link{twDEMCInt}}
	
	if( length(n)==1 ) n=rep(n,3)
	xcuts <- cutQuantiles(x[,vars[1] ],g=n[1], levels.mean=TRUE)
	ycuts <- cutQuantiles(x[,vars[2] ],g=n[2], levels.mean=TRUE)
	zcuts <- cutQuantiles(x[,vars[3] ],g=n[3], levels.mean=TRUE)
	tmp <- list(xcuts,ycuts,zcuts); names(tmp) <- dimnames	
	#cagg <- tapply( x[,intCol], tmp, fAgg )
	cagg <- tapply( x[,intCol], tmp, fAgg, ...)
	### numeric array with each item representing application of fAgg to x[,intCol]
	### and dimnames equal to means of the classes
}
attr(marginals3d,"ex") <- function(){
	res <- marginals3d(sample,vars=c("kY","kO","tLagLeaf"),n=8)
	tmp.f <- function(){
		library(Rcmdr)
		ds <- melt(res)
		scatter3d(ds$kY, ds$tLagLeaf, ds$kO
			, surface=FALSE
			,bg="white", axis.scales=TRUE, grid=TRUE, ellipsoid=FALSE
			, xlab="kY",ylab="tLagLeaf", zlab="kO"
			, point.col=rev(heat.colors(100))[round(rescale(ds$value,to=c(1,100)))]
			,sphere.size=1.5
		)
	}
	
}

plot2DKDMarginals <- function(
	### Evaluates the Kernel-density regression on a grid and displays a figure.
	form	##<< formula of regression e.g. rLogLik~dO+tvrO
	,ds		##<< dataframe or matrix with column names corresponding to variables in the formula
	,dims=40	##<< vector of spaces across predictor space where 
	,bandwidth = NULL	 ##<< integer vector: bandwidth for each dimension
		### see adaptive bandwidth of npreg.
		### about Records to include for estimating bandwidth.
		### If only one value is supplied, then it is repeated for each predictor
	,nSample=200 ##<< the sample size to estimate the bandwidth and the marginal
		### The larger the sample, the longer it takes
	,bwtype="adaptive_nn" ,tol=.1, ftol=.1	##<< parameters for \code{\link{npregbw}}
	,bws=NULL	##<< may directly supply a bandwidht object
	,...		##<< further arguments to \code{\link{plot.twApply2DMesh}}
){
	##details<<
	## Uses package snowfall to split calculation to subprocesses.
	## Requires package np loaded in cluster (sfLibrary(np))
	
	if( !is.data.frame(ds) ) ds <- as.data.frame(ds)
	dsPred <- model.frame(form,ds)
	nVars <- ncol(dsPred)-1
	if( nVars != 2) error("plot2DKDMarginals: number of predictor variables must be 2.")
	if( length(bandwidth) == 1) bandwidth <- rep(bandwidth, nVars)
	#done in twApply3DMesh if( length(dims) == 1) dims <- rep(dims, nVars)
	
	if( nSample < nrow(ds))
		dss <- ds[sample.int(nrow(ds),nSample),]	# tradeoff 200: larger sample sizes become really slow
	else
		dss <- ds
	dsPred <- model.frame(form,dss)	# constrain to sample
	
	#bw0 <- npregbw(formula=form, ds100, bwmethod="normal-reference" )
	#bw3 <- bw30 <- npregbw(formula=form, ds100, bwtype="generalized_nn" ,tol=.1, ftol=.1 )
	dssbw <- if( 0<length(bandwidth) ){
		#do not need to estimate bandwidht, just take 3 samples
		dss[sample.int(nrow(dss),3),]
	}else
		dss
	sfLibrary(np)
	tmpf <- function(i, dssbw, ...){
		#generate different starting values
		set.seed(i*1000)
		bwsi <- pmax(3,round( nrow(dssbw)/10+rnorm(2,sd=nrow(dssbw)/20)) )
		npregbw( nmulti=1, bws=bwsi, data=dssbw, ... )
	}
	#mtrace(tmpf)
	#tmpf(1, formula=form, dssbw=dssbw, bwtype=bwtype ,tol=tol, ftol=ftol)
	bw3L <- sfLapply( 1:4, tmpf, formula=form, dssbw=dssbw, bwtype=bwtype ,tol=tol, ftol=ftol)			
	fval <- sapply(bw3L,"[","fval")
	#sapply(bw3L,"[","bw")
	bw3 <- bw30 <- 	bw3L[[which.min(fval)]]

	#bw3 <- bw30 <- npregbw(formula=form, dssbw, bwtype=bwtype ,tol=tol, ftol=ftol )
	# refine (but does not change)
	#bw3 <- npregbw(bws=bw30, formula=form, ds100, bwtype="adaptive_nn", tol=0.01, ftol=0.01 )
	#bw3$bw <- bw30$bw*2
	#plot(bw32)
	#est <- npreg(bws=bw3, tdat=dsPred, edat=dsPred)		
	lMin <- min(dsPred[,1])
	tmpf <- function(x,y, bws, minLogLik=lMin){
		ind <- splitIndices(length(x),4)
		exdat <- data.frame(x=as.numeric(x),y=as.numeric(y)); names(exdat) <- bws$xnames
		fList <- list(
			function(bws,dsPred,exdat,ind){ npreg(bws=bws,txdat=dsPred[,-1], tydat=dsPred[,1], exdat=exdat[ind[[1]],]) }
			,function(bws,dsPred,exdat,ind){ npreg(bws=bws,txdat=dsPred[,-1], tydat=dsPred[,1], exdat=exdat[ind[[2]],]) }
			,function(bws,dsPred,exdat,ind){ npreg(bws=bws,txdat=dsPred[,-1], tydat=dsPred[,1],exdat=exdat[ind[[3]],]) }
			,function(bws,dsPred,exdat,ind){ npreg(bws=bws,txdat=dsPred[,-1], tydat=dsPred[,1], exdat=exdat[ind[[4]],]) }
		)
		argsList <- list( bws=bws,dsPred=dsPred,exdat=exdat,ind=ind)
		#tmpf2<-fList[[2]]
		#mtrace(tmpf2)
		#resL <- tmpf2(bws=bws,dsPred=dsPred,exdat=exdat,ind=ind)
		resL <- sfPar(fList,sfParArgsList=argsList)
		tmpMean <- do.call(c, lapply( resL, "[[", "mean" ))
		ifelse( tmpMean>minLogLik, tmpMean, as.numeric(NA))
	}
	#tmp2 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=bw3$bw*8;tmp}), dims=60, col=rev(heat.colors(20)) )
	#tmp1 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=c(200,200);tmp}), dims=60, col=rev(heat.colors(20)) )
	#mtrace(tmpf)
	if( 0==length(bandwidth) ) bandwidth<-bw3$bw
	tmp2 <- twApply2DMesh(dsPred[,-1],FUN=tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=bandwidth;tmp}), dims=dims, label=colnames(dsPred)[1] )
	#tmp2 <- twPlot2DFun( dsPred[,2], dsPred[,3], tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=bandwidth;tmp}), dims=dim[1], ydiv=dim[2], col=rev(heat.colors(20)), xlab=bw3$xnames[1], ylab=bw3$xnames[2],key.title=title(sub=paste(colnames(dsPred)[1],"\n",sep="")) )
	names(dimnames(tmp2)) <- bw3$xnames
	attr(tmp2,"bws") <- bw3
	plot(tmp2, ...)
	return(invisible(tmp2))
	### Matrix of predictions for the grid.
}
attr(plot2DKDMarginals,"ex") <- function(){
	.tmp.f <- function(){
		#wrapped inside function because it takes long
		data(twdemcEx1)
		sample <- stackChains(twdemcEx1)
		sample0 <- sample[ sample[,1] >= quantile(sample[,1],0.05), ]
		# just plot the sample
		plot( sample0[,-1], col=rev(heat.colors(20))[round(twRescale(sample0[,1],c(1,20)))] )
		
		#simple plot
		#mtrace(plotMarginal2D)
		plotMarginal2D( sample0, "a", "b", minN=1 )	# actually not a marginal
		
		# kernel density regression estimate
		#mtrace(plot2DKDMarginals)
		rl <- sample0[,"rLogLik"]
		set.seed(0815)
		sample0b <- cbind(sample0, rLogLik2 = rl+rnorm(length(rl),sd=diff(range(rl))/4 ))
		plot( rLogLik2 ~ rLogLik, sample0b)
		plot( sample0b[,c("a","b")], col=rev(heat.colors(20))[round(twRescale(sample0b[,"rLogLik2"],c(1,20)))] )
		#refit disturbed sample
		tmp2 <- tmp3 <- plot2DKDMarginals( rLogLik2 ~ a+b, sample0b, dim=20  )
		attributes(tmp2)$bws$bandwidth
		# the bandwidth seems to be underestimated for the second dimension, repeat with higher bandwidth
		#tmp3 <- plot2DKDMarginals( rLogLik ~ a+b, sample0b, dim=20, bandwidth=30  )
		#refine plot
		plot(tmp3, xlab="Intercept a", ylab="Slope b", zlab="Log-Like-\nlihood", contour=TRUE )
	}
}

plot3DKDMarginals <- function(
	### Evaluates the Kernel-density regression on a grid and displays a figure.
	form	##<< formula of regression e.g. rLogLik~tvrY+tvrO+biasLitterLeaf
	,ds		##<< dataframe or matrix with column names corresponding to variables in the formula
	,probs=c(0.2,0.5,0.75,0.95)		##<< the percentiles at which to plot contour surfaces for the response variable
		## the levels are calculated from the subsample 
	,...		##<< further argument to \code{\link{plot.twApply3DMesh}}
	,dims=8	##<< vector of points across predictor space within each dimension: scales bad: n=prod(dims)~dims[1]^3  
	,bandwidth = NULL	 ##<< integer vector: bandwidth for each dimension
	,knotSpacing = "quantile"	##<< see argument in \code{\link{twApply3DMesh}}
	### see adaptive bandwidth of npreg.
	### about Records to include for estimating bandwidth.
	### If only one value is supplied, then it is repeated for each predictor
	,nSampleBandwidth=500 ##<< the sample size to estimate the bandwidth
	,nSample=1000		##<< sample used for interpolation
	### The larger the sample, the longer it takes
	,bwtype="adaptive_nn" ,tol=.1, ftol=.1	##<< parameters for \code{\link{npregbw}}
){
	##details<<
	## Uses package snowfall to split calculation to subprocesses.
	## Requires package np loaded in cluster (sfLibrary(np))
	
	if( !is.data.frame(ds) ) ds <- as.data.frame(ds)
	dsPredAll <- model.frame(form,ds)
	nVars <- ncol(dsPredAll)-1
	if( nVars != 3) error("plot3DKDMarginals: number of predictor variables must be 3.")
	if( length(bandwidth) == 1) bandwidth <- rep(bandwidth, nVars)
	if( length(dims) == 1) dims <- rep(dims, nVars)
	
	# the sample to draw (and to predict new values)
	if( nSample < nrow(ds))
		dss <- ds[sample.int(nrow(ds),nSample),]	# tradeoff 200: larger sample sizes become really slow
	else
		dss <- ds
	dsPred <- model.frame(form,dss)	# constrain to sample
	
	#bw0 <- npregbw(formula=form, ds100, bwmethod="normal-reference" )
	#bw3 <- bw30 <- npregbw(formula=form, ds100, bwtype="generalized_nn" ,tol=.1, ftol=.1 )
	# the sample to estimate the bandwidth
	dssbw <- if( 0<length(bandwidth) ){
			#do not need to estimate bandwidht, just take 3 samples
			dss[sample.int(nrow(dss),3),]
		}else if( nSampleBandwidth<nrow(dss) ){
			#do not need to estimate bandwidht, just take 3 samples
			dss[sample.int(nrow(dss),nSampleBandwidth),]
		}else 
			dss
	
	#bw0 <- npregbw(formula=form, ds100, bwmethod="normal-reference" )
	#bw3 <- bw30 <- npregbw(formula=form, ds100, bwtype="generalized_nn" ,tol=.1, ftol=.1 )
	#sfLibrary(np) must be done before
	tmpf <- function(i, dssbw, ...){
		#generate different starting values
		set.seed(i*1000)
		bwsi <- pmax(3,round( nrow(dssbw)/10+rnorm(3,sd=nrow(dssbw)/20)) )
		npregbw( nmulti=1, bws=bwsi, data=dssbw, ... )
	}
	#mtrace(tmpf)
	#tmpf(1, formula=form, dssbw=dssbw, bwtype=bwtype ,tol=tol, ftol=ftol)
	bw3L <- sfLapply( 1:4, tmpf, formula=form, dssbw=dssbw, bwtype=bwtype ,tol=tol, ftol=ftol)			
	fval <- sapply(bw3L,"[","fval")
	#sapply(bw3L,"[","bw")
	bw3 <- bw30 <- 	bw3L[[which.min(fval)]]
	# refine (but does not change)
	#bw3 <- npregbw(bws=bw30, formula=form, ds100, bwtype="adaptive_nn", tol=0.01, ftol=0.01 )
	#bw3$bw <- bw30$bw*2
	#plot(bw32)
	#est <- npreg(bws=bw3, tdat=dsPred, edat=dsPred)		
	lMin <- min(dsPred[,1])
	tmpf <- function(x,y,z, bws, minLogLik=lMin){
		ind <- splitIndices(length(x),4)
		exdat <- data.frame(x=as.numeric(x),y=as.numeric(y),z=as.numeric(z)); names(exdat) <- bws$xnames
		fList <- list(
			function(bws,dsPred,exdat,ind){ npreg(bws=bws,txdat=dsPred[,-1], tydat=dsPred[,1], exdat=exdat[ind[[1]],]) }
			,function(bws,dsPred,exdat,ind){ npreg(bws=bws,txdat=dsPred[,-1], tydat=dsPred[,1], exdat=exdat[ind[[2]],]) }
			,function(bws,dsPred,exdat,ind){ npreg(bws=bws,txdat=dsPred[,-1], tydat=dsPred[,1],exdat=exdat[ind[[3]],]) }
			,function(bws,dsPred,exdat,ind){ npreg(bws=bws,txdat=dsPred[,-1], tydat=dsPred[,1], exdat=exdat[ind[[4]],]) }
		)
		resL <- sfPar(fList,sfParArgsList=list( bws=bws,dsPred=dsPred,exdat=exdat,ind=ind))
		tmpMean <- do.call(c, lapply( resL, "[[", "mean" ))
		ifelse( tmpMean>minLogLik, tmpMean, as.numeric(NA))
	}
	#tmp2 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=bw3$bw*8;tmp}), dims=60, col=rev(heat.colors(20)) )
	#tmp1 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=c(200,200);tmp}), dims=60, col=rev(heat.colors(20)) )
	#mtrace(tmpf)
	if( 0==length(bandwidth) ) bandwidth<-bw3$bw
	#tmp2 <- twPlot3DContourFun( dsPred[,-1], FUN=tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=bandwidth;tmp}), dims=dim[1], ydiv=dim[2], col=rev(heat.colors(20)), xlab=bw3$xnames[1], ylab=bw3$xnames[2],key.title=title(sub=paste(colnames(dsPred)[1],"\n",sep="")) )
	tmp2 <- twApply3DMesh(dsPred[,-1],FUN=tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=bandwidth;tmp}), dims=dims, nSample=nSample, label=colnames(dsPred)[1], knotSpacing=knotSpacing )
	#names(dimnames(tmp2)) <- bw3$xnames
	tmp2$bandwidht <- bandwidth
	#levels <- quantile( dsPredAll[,1], probs=probs)
	levels <- quantile( tmp2$sample[,4], probs=probs)
	plot(tmp2, levels=levels, ...)
	return(invisible(tmp2))
	### 3D Array of Kernel density regression predictions over the mesh. See \code{\link{twApply3DMesh}}.
}
attr(plot3DKDMarginals,"ex") <- function(){
	#Example: Nested contours of mixture of three tri-variate normal densities
	# wrapped inside a function, because it takes so long.
	.tmp.f <- function(){
		nmix3 <- function(x, y, z, m, s) {
			0.4 * dnorm(x, m, s) * dnorm(y, m, s) * dnorm(z, m, s) +
				0.3 * dnorm(x, -m, s) * dnorm(y, -m, s) * dnorm(z, -m, s) +
				0.3 * dnorm(x, m, s) * dnorm(y, -1.5 * m, s) * dnorm(z, m, s)
		}
		f <- function(x,y,z) nmix3(x,y,z,.5,.5)
		
		n <- 250
		x <- rnorm(n,.5,.5)
		y <- c(rnorm(n/2,.5,.5), rnorm(n/2,-.5,.5)) 
		zz <- rnorm(n,.5,.5)
		
		plot(tmp <- twApply3DMesh(x,y,zz,f, nSample=1000, dims=10))	# just the points
		plot(tmp, probs=seq(0.5, 0.95, len=4))
		tmpObs <- attributes(tmp)$sample
		tmpObs$f2 <- tmpObs$f + rnorm(nrow(tmpObs),sd=0.025)
		plot(tmpObs$f2~tmpObs$f)
		plot3d( tmpObs, col=rev(heat.colors(20))[round(twRescale(tmpObs$f2,c(1,20)))]  )
		
		#reestimate density from disturbed sample
		# kernel density regression estimate
		#mtrace(plot3DKDMarginals)
		tmp2 <- plot3DKDMarginals( f2 ~ x+y+zz, tmpObs, dims=10, nSampleBandwidth=500  )
		tmp2b <- tmp2; attr(tmp2b,"sample")<-NULL; plot(tmp2b)	# show the mesh instead of the example
		plot(tmp2, probs=c(0.25,0.5,0.75,0.9), nDrawPoints=0)		#add contours
	}
}





