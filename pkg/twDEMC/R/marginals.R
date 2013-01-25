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
.cut2dim3 <- function(x,y,z,c,xn=round(length(c)^(1/3)),yn=xn, zn=xn, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), zlab=deparse(substitute(z)), fun=mean ){
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
	## \code{\link{twDEMCBlockInt}}
	
	##details<< 
	## There are several methods to get calculate marginal densities for multivariate sample. \itemize{
	## \item{ integrate over all others but 1 column: this method  } 
	## \item{ integrate over all others but 2 column: \code{\link{marginals2d}}  } 
	## \item{ integrate over all others but 3 column: \code{\link{marginals3d}}  } 
	##}
    # moved to plot package
	#? \item{ plotting a 2D marginal distribution: \code{\link{plotMarginal2D}}  } 
    # \item{ plotting a 2D kernel density estimate: \code{\link{plot2DKDMarginals}}  } 
    
	xcuts <- cutQuantiles(x[,vars],g=n, levels.mean=TRUE)
	tmp <- list(xcuts); names(tmp) <- dimnames	
	#cagg <- tapply( x[,intCol], tmp, fAgg )
	cagg <- tapply( x[,intCol], tmp, fAgg, ...)
	### numeric vector with each item representing application of fAgg to x[,intCol]
	### and dimnames equal to means of the classes
}
attr(marginals1d,"ex") <- function(){
	data(twdemcEx1)
	sample <- stackChains(twdemcEx1)
	#res <- marginals1d(sample,vars="kO",n=20)
	res <- marginals1d(sample,vars="b",n=20)
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
	## \code{\link{twDEMCBlockInt}}
	
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
	data(twdemcEx1)
	sample <- stackChains(twdemcEx1)
	#res <- marginals2d(sample,vars=c("kY","kO"),n=10)
	res <- marginals2d(sample,vars=c("a","b"),n=10)
	if( require(lattice) ){
		levelplot( value~a*b, data=melt(res), col.regions=rev(heat.colors(100)), xlab=names(dimnames(res))[1], ylab=names(dimnames(res))[2] )
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
	## \code{\link{twDEMCBlockInt}}
	
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
	.tmp.f <- function(){
		res <- marginals3d(sample,vars=c("kY","kO","tLagLeaf"),n=8)
		#library(Rcmdr)
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




