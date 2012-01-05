plotMarginal2D <- function(
	### Plot the 2D marginal Density of a sample
	smp			##<< numeric matrix: first column logDensity, other columns free parameters, see \code{\link{stackChains.twDEMC}}
	,xCol=2		##<< index (column number or column name) of column for x ordinate
	,yCol=3		##<< index (column number or column name) of column for y ordinate
	, intCol=1  ##<< index (column number or column name) of column of LogDensity values
	, grains=18 ##<< vector of length 1 or 2 giving the number of groups for x and y classes respectively
	, FUN=mean	 ##<< the function applied over logDensitys in the classes
	, argsFUN=list(na.rm=TRUE)	##<< additional arguments to FUN 
	, col=rev(heat.colors(20))	##<< vector of colors
	, minN=7	##<< minimum number of items in classpixel to be plotted
	, ...		##<< additional arguments to twPlot2D
){
	##details<< 
	## The entire sample is split into bins of about equal number of observations and all observations regarding x and y.
	## Within the bin the values are aggregated. By default the mean is calculated.
	
	##seealso<<  
	## \code{\link{twDEMCInt}}
	
	##details<< 
	## There are several plotting methods related to twDEMC run. \itemize{
	## \item{ TODO: link methods  } 
	## \item{ the Gelman criterion: this method  } 
	## \item{ the theorectical minimum logDen-Value for significant model difference : \code{\link{getRLogDenQuantile}}  } 
	##}
	
	if( length(grains)==1) grains=c(grains,grains)	
	smpGrained <- cbind(smp[,c(xCol,yCol)],smp[,intCol,drop=FALSE])
	smpGrained[,1] <- {tmp<-cutQuantiles( smpGrained[,1],g=grains[1],levels.mean=TRUE); as.numeric(levels(tmp))[tmp]}
	smpGrained[,2] <- {tmp<-cutQuantiles( smpGrained[,2],g=grains[2],levels.mean=TRUE); as.numeric(levels(tmp))[tmp]}
	tmpGrp <- (as.factor(smpGrained[,1])):(as.factor(smpGrained[,2]))
	xs <- as.numeric(levels(as.factor(smpGrained[,1])))
	ys <- as.numeric(levels(as.factor(smpGrained[,2])))
	smp2 <- data.frame(
		x=rep(xs, each=grains[1])
		,y=rep(ys, grains[2])
		,marginalLogDen = do.call( tapply, c(list(smpGrained[,3], tmpGrp, FUN),argsFUN)  )
		,n=as.vector(table(tmpGrp))
	)
	names(smp2)[1:2] <- colnames(smpGrained)[1:2]
	smp3 <- smp2
	smp3[ smp2[,"n"]<minN, 3] <- NA
	#plot( tmp <- twApply2DMesh(xs, ys, FUN=smp3$marginalLogDen, knotSpacing="all"), xlab=colnames(smpGrained)[1], ylab=colnames(smpGrained)[2], col=col )
	#plot.twApply2DMesh(smp3$marginalLogDen, xlab=colnames(smpGrained)[1], ylab=colnames(smpGrained)[2], col=col ) 
	twPlot2D(xs, ys, z=matrix(smp3$marginalLogDen,nrow=length(xs)), xlab=colnames(smpGrained)[1], ylab=colnames(smpGrained)[2], zlab=colnames(smpGrained)[3], col=col, ... ) 
	return(invisible(smp2))
}
attr(plotMarginal2D,"ex") <- function(){
	data(twdemcEx1)
	sample <- stackChains(twdemcEx1)
	sample0 <- sample[ sample[,1] >= quantile(sample[,1],0.05), ]
	#simple plot
	#mtrace(plotMarginal2D)
	plotMarginal2D( sample0, "a", "b", minN=1 )	# actually not a marginal
}

.plotMarginal2D.levelplot <- function(
	### Plot the 2D marginal Density of a sample
	smp			##<< numeric matrix: first column logDensity, other columns free parameters, see \code{\link{stackChains.twDEMC}}
	,xCol=2		##<< index (column number or column name) of column for x ordinate
	,yCol=3		##<< index (column number or column name) of column for y ordinate
	, intCol=1  ##<< index (column number or column name) of column of LogDensity values
	, ...		##<< additional arguments to FUN 
	, grains=20	 ##<< vector of length 1 or 2 giving the number of groups for x and y classes respectively
	, FUN=mean	 ##<< the function applied over logDensitys in the classes
	, col=rev(heat.colors(100))	##<< vector of colors
){
	##details<< 
	## The entire sample is split into bins of about equal number of observations and all observations regarding x and y.
	## Within the bin the values are aggregated. By default the mean is calculated.
	
	##seealso<<  
	## \code{\link{twDEMCInt}}
	
	##details<< 
	## There are several plotting methods related to twDEMC run. \itemize{
	## \item{ TODO: link methods  } 
	## \item{ the Gelman criterion: this method  } 
	## \item{ the theorectical minimum logDen-Value for significant model difference : \code{\link{getRLogDenQuantile}}  } 
	##}
	
	if( length(grains)==1) grains=c(grains,grains)	
	smpGrained <- cbind(smp[,c(xCol,yCol)],smp[,intCol])
	smpGrained[,1] <- {tmp<-cutQuantiles( smpGrained[,1],g=grains[1],levels.mean=TRUE); as.numeric(levels(tmp))[tmp]}
	smpGrained[,2] <- {tmp<-cutQuantiles( smpGrained[,2],g=grains[2],levels.mean=TRUE); as.numeric(levels(tmp))[tmp]}
	#smpMarginal <- aggregate(smpGrained[,3], as.data.frame(smpGrained[,1:2]), FUN=FUN )
	smpMarginal <- aggregate(smpGrained[,3], as.data.frame(smpGrained[,1:2]), FUN=FUN, ... )
	names(smpMarginal) <- c("x","y","z")
	levelplot(z~x*y, data=smpMarginal, col.regions=col, xlab=colnames(smpGrained)[1], ylab=colnames(smpGrained)[2])
}


plotConditional2D <- function(
	### Plot the 2D conditional profile-Density of a sample: calling \code{\link{plotMarginal2D}} with aggregating function max.
	smp			
	,xCol=2
	,yCol=3
	,intCol=1
	,...
){
	##seealso<<  
	## \code{\link{plotMarginal2D}}
	## \code{\link{twDEMCInt}}
	
	plotMarginal2D(smp,xCol,yCol,intCol,FUN=max,...)
}

