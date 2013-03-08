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
	## \code{\link{twDEMCBlockInt}}
	
	##details<< 
	## There are several plotting methods related to twDEMC run. \itemize{
	## \item{ TODO: link methods  } 
	## \item{ the Gelman criterion: this method  } 
	##}

    # \item{ the theorectical minimum logDen-Value for significant model difference : \code{\link{getRLogDenQuantile}}  } 
    
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
	## \code{\link{twDEMCBlockInt}}
	
	##details<< 
	## There are several plotting methods related to twDEMC run. \itemize{
	## \item{ TODO: link methods  } 
	## \item{ the Gelman criterion: this method  } 
	##}
    
    ## \item{ the theorectical minimum logDen-Value for significant model difference : \code{\link{getRLogDenQuantile}}  } 
    
	
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
	## \code{\link{twDEMCBlockInt}}
	
	plotMarginal2D(smp,xCol,yCol,intCol,FUN=max,...)
}




plot2DKDMarginals <- function(
        ### Evaluates the Kernel-density regression on a grid and displays a figure.
        form	##<< formula of regression e.g. rLogDen~dO+tvrO
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
    #
    if( !is.data.frame(ds) ) ds <- as.data.frame(ds)
    dsPred <- model.frame(form,ds)
    nVars <- ncol(dsPred)-1
    if( nVars != 2) stop("plot2DKDMarginals: number of predictor variables must be 2.")
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
    tmpf <- function(x,y, bws, minLogDen=lMin){
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
        ifelse( tmpMean>minLogDen, tmpMean, as.numeric(NA))
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
        rl <- sample0[,"rLogDen"]
        set.seed(0815)
        sample0b <- cbind(sample0, rLogDen2 = rl+rnorm(length(rl),sd=diff(range(rl))/4 ))
        plot( rLogDen2 ~ rLogDen, sample0b)
        plot( sample0b[,c("a","b")], col=rev(heat.colors(20))[round(twRescale(sample0b[,"rLogDen2"],c(1,20)))] )
        #refit disturbed sample
        tmp2 <- tmp3 <- plot2DKDMarginals( rLogDen2 ~ a+b, sample0b, dim=20  )
        attributes(tmp2)$bws$bandwidth
        # the bandwidth seems to be underestimated for the second dimension, repeat with higher bandwidth
        #tmp3 <- plot2DKDMarginals( rLogDen ~ a+b, sample0b, dim=20, bandwidth=30  )
        #refine plot
        plot(tmp3, xlab="Intercept a", ylab="Slope b", zlab="LogDensity", contour=TRUE )
    }
}

plot3DKDMarginals <- function(
        ### Evaluates the Kernel-density regression on a grid and displays a figure.
        form	##<< formula of regression e.g. rLogDen~tvrY+tvrO+biasLitterLeaf
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
    tmpf <- function(x,y,z, bws, minLogDen=lMin){
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
        ifelse( tmpMean>minLogDen, tmpMean, as.numeric(NA))
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

