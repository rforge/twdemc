.cut2dim2 <- function(x,y,c,xn=round(sqrt(length(c))),yn=round(sqrt(length(c))), xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), fun=median ){
	# cut2dim2
	#
	#matrix mean c over cells of unregular grid over x and y
	#grid is defined by classes of equal number of observations in x and y
	#xn and yn number of classes in x and y dimension
	xcuts=cut2(x,g=xn, levels.mean=TRUE)
	ycuts=cut2(y,g=yn, levels.mean=TRUE)
	tmp <- list(xcuts,ycuts); names(tmp) <- c(xlab,ylab)
	cagg <- tapply(c, tmp, fun)
	list( cagg=cagg, xmeans=as.numeric(levels(xcuts)), ymeans=as.numeric(levels(ycuts)) )
}
.cut2dim3 <- function(x,y,z,c,xn=round(length(c)^(1/3)),yn=xn, zn=xn, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), fun=mean ){
	#matrix mean c over cells of unregular grid over x,y, and z
	#grid is defined by classes of equal number of observations in x and y
	#xn and yn number of classes in x and y dimension
	xcuts=cut2(x,g=xn, levels.mean=TRUE)
	ycuts=cut2(y,g=yn, levels.mean=TRUE)
	zcuts=cut2(z,g=zn, levels.mean=TRUE)
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
	
	xcuts <- cut2(x[,vars],g=n, levels.mean=TRUE)
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
	xcuts <- cut2(x[,vars[1] ],g=n[1], levels.mean=TRUE)
	ycuts <- cut2(x[,vars[2] ],g=n[2], levels.mean=TRUE)
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
	xcuts <- cut2(x[,vars[1] ],g=n[1], levels.mean=TRUE)
	ycuts <- cut2(x[,vars[2] ],g=n[2], levels.mean=TRUE)
	zcuts <- cut2(x[,vars[3] ],g=n[3], levels.mean=TRUE)
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
	,dim=40	##<< vector of spaces across predictor space where 
	,bandwidth = NULL	 ##<< integer vector: bandwidth for each dimension
		### see adaptive bandwidth of npreg.
		### about Records to include for estimating bandwidth.
		### If only one value is supplied, then it is repeated for each predictor
	,nSampleBandwidth=200 ##<< the sample size to estimate the bandwidth
		### The larger the sample, the longer it takes
	,bwtype="adaptive_nn" ,tol=.1, ftol=.1	##<< parameters for \code{\link{npregbw}}
){
	##details<<
	## Uses package snowfall to split calculation to subprocesses.
	## Requires package np loaded in cluster (sfLibrary(np))
	
	if( !is.data.frame(ds) ) ds <- as.data.frame(ds)
	dsPred <- model.frame(form,ds)
	nVars <- ncol(dsPred)-1
	if( nVars != 2) error("plot2DKDMarginals: number of predictor variables must be 2.")
	if( length(bandwidth) == 1) bandwidth <- rep(bandwidth, nVars)
	if( length(dim) == 1) dim <- rep(dim, nVars)
	
	if( nSampleBandwidth < nrow(ds))
		ds100 <- ds[sample.int(nrow(ds),nSampleBandwidth),]	# tradeoff 200: larger sample sizes become really slow
	else
		ds100 <- ds
	#bw0 <- npregbw(formula=form, ds100, bwmethod="normal-reference" )
	#bw3 <- bw30 <- npregbw(formula=form, ds100, bwtype="generalized_nn" ,tol=.1, ftol=.1 )
	bw3 <- bw30 <- npregbw(formula=form, ds100, bwtype=bwtype ,tol=tol, ftol=ftol )
	# refine (but does not change)
	#bw3 <- npregbw(bws=bw30, formula=form, ds100, bwtype="adaptive_nn", tol=0.01, ftol=0.01 )
	#bw3$bw <- bw30$bw*2
	#plot(bw32)
	#est <- npreg(bws=bw3, tdat=dsPred, edat=dsPred)		
	lMin <- min(dsPred[,1])
	#sfLibrary(np)
	tmpf <- function(x,y, bws, minLogLik=lMin){
		ind <- splitIndices(length(x),4)
		exdat <- data.frame(x=x,y=y); names(exdat) <- bws$xnames
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
	#tmp2 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=bw3$bw*8;tmp}), xdiv=60, col=rev(heat.colors(20)) )
	#tmp1 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=c(200,200);tmp}), xdiv=60, col=rev(heat.colors(20)) )
	#mtrace(tmpf)
	if( 0==length(bandwidth) ) bandwidth<-bw3$bw
	tmp2 <- twPlot2DFun( dsPred[,2], dsPred[,3], tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=bandwidth;tmp}), xdiv=dim[1], ydiv=dim[2], col=rev(heat.colors(20)), xlab=bw3$xnames[1], ylab=bw3$xnames[2],key.title=title(sub=paste(colnames(dsPred)[1],"\n",sep="")) )
	names(dimnames(tmp2)) <- bw3$xnames
	tmp2
	### Matrix of predictions for the grid.
}
attr(plot2DKDMarginals,"ex") <- function(){
	data(twdemcEx1)
	sample <- stackChains(twdemcEx1)
	sample0 <- sample[ sample[,1] >= quantile(sample[,1],0.05), ]
	#simple plot
	plotMarginal2D( sample0, "a", "b", minN=1 )	# actually not a marginal
	
	# kernel density regression estimate
	#mtrace(plot2DKDMarginals)
	tmp2 <- plot2DKDMarginals( rLogLik ~ a+b, sample0, dim=20  )
	
	#refine plot
	xy <- apply(do.call(cbind,dimnames(tmp2)),2,as.numeric  )
	twPlot2DFun.contour( xy[,1], xy[,2], tmp2, xdiv=nrow(tmp2)
		, xlab="Intercept a"
		, ylab="Slope b"
		, key.title=title(sub="Log-Like-\nlihood\n\n")
		, color.palette=function(n){rev(heat.colors(n))
	})
}




