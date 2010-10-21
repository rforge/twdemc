plotMarginal2D <- function(
	### Plot the 2D marginal Likelihood of a sample
	smp			##<< numeric matrix: first column log-Likelihood, other columns free parameters, see \code{\link{stackChains.twDEMC}}
	,xCol=2		##<< index (column number or column name) of column for x ordinate
	,yCol=3		##<< index (column number or column name) of column for y ordinate
	, intCol=1  ##<< index (column number or column name) of column of LogLikelihood values
	, ...		##<< additional arguments to FUN 
	, grains=20	 ##<< vector of length 1 or 2 giving the number of groups for x and y classes respectively
	, FUN=mean	 ##<< the function applied over log-Likelihoods in the classes
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
	## \item{ the theorectical minimum logLik-Value for significant model difference : \code{\link{getRLogLikQuantile}}  } 
	##}
	
	if( length(grains)==1) grains=c(grains,grains)	
	smpGrained <- cbind(smp[,c(xCol,yCol)],smp[,intCol])
	smpGrained[,1] <- {tmp<-cut2( smpGrained[,1],g=grains[1],levels.mean=TRUE); as.numeric(levels(tmp))[tmp]}
	smpGrained[,2] <- {tmp<-cut2( smpGrained[,2],g=grains[2],levels.mean=TRUE); as.numeric(levels(tmp))[tmp]}
	#smpMarginal <- aggregate(smpGrained[,3], as.data.frame(smpGrained[,1:2]), FUN=FUN )
	smpMarginal <- aggregate(smpGrained[,3], as.data.frame(smpGrained[,1:2]), FUN=FUN, ... )
	names(smpMarginal) <- c("x","y","z")
	levelplot(z~x*y, data=smpMarginal, col.regions=col, xlab=colnames(smpGrained)[1], ylab=colnames(smpGrained)[2])
}

plotConditional2D <- function(
	### Plot the 2D conditional profile-Likelihood of a sample: calling \code{\link{plotMarginal2D}} with aggregating function max.
	smp			
	,xCol
	,yCol
	,intCol=1
	,...
){
	##seealso<<  
	## \code{\link{plotMarginal2D}}
	## \code{\link{twDEMCInt}}
	
	plotMarginal2D(smp,xCol,yCol,intCol,FUN=max,...)
}


ggplotChainPopMoves <- function(
	### Plot boxplots of the distribution of first and last fifth of steps for each population 
	resB			##<< the twDEMC to examine
	, iChains = rep(1:ncol(resB$temp), each=ncol(resB$rLogLik)%/%ncol(resB$temp))
		### mapping from chain to population
	, doSort=TRUE	##<< if TRUE result $rLogLik is sorted
){
	# ggplotChainPopMoves
	##seealso<<  
	## \code{\link{plotMarginal2D}}
	## \code{\link{twDEMCInt}}
	
	iGen <- cbind( floor(c(0,1/5)*nrow(resB$rLogLik))+1, floor(c(4/5,1)*nrow(resB$rLogLik)) )
	iLabel <- apply(iGen,2,function(iGeni){paste(range(iGeni), collapse=" to ")}) 
	tmp1 <- melt(resB$rLogLik[iGen[1,1]:iGen[2,1],]); tmp1$pos="from"
	tmp2 <- melt(resB$rLogLik[iGen[1,2]:iGen[2,2],]); tmp2$pos="to"
	
	dfLogLik2 <- rbind(tmp1,tmp2)
	dfLogLik2$Pop <- iChains[dfLogLik2$X2]  #(dfLogLik2$X2-1) %% ncol(resB$temp) +1
	dfLogLik2 <- dfLogLik2[order(dfLogLik2$Pop,dfLogLik2$pos),]
	dfLogLik2$PopPos <- paste(dfLogLik2$Pop,dfLogLik2$pos, sep="_")
	dfLogLik2$PopPos <- factor( dfLogLik2$PopPos, levels=unique(dfLogLik2$PopPos) )
	
	dfMed <- aggregate(dfLogLik2$value,list(Pop=dfLogLik2$Pop,pos=dfLogLik2$pos), median)
	#dfMed$q25 <- aggregate(dfLogLik2$value,list(Pop=dfLogLik2$Pop,pos=dfLogLik2$pos), function(values){ quantile(values, probs=0.25)})
	dfMed$q25 <- aggregate(dfLogLik2$value,list(Pop=dfLogLik2$Pop,pos=dfLogLik2$pos), quantile, probs=0.25)$x
	dfMed <- dfMed[order(dfMed$Pop,dfMed$pos),]
	dfMed$PopPos <- paste(dfMed$Pop,dfMed$pos, sep="_")
	dfMed$PopPos <- factor( dfMed$PopPos, levels=unique(dfMed$PopPos) )
	dfMed$medFrom <- dfMed$x 
	dfMed$medTo <- c(dfMed$x[-1],NA)
	dfMedSingle <- subset(dfMed, pos=="from")
	
	p5 <- ggplot( dfLogLik2, aes(x=pos,y=value) )+ xlab("Population") +ylab("Log-Likelihood") +	facet_grid(.~Pop)
	p6 <- p5+
	geom_boxplot(aes(fill=pos))+
	geom_point( aes(x=1,y=medFrom), data=dfMedSingle, shape=1 )+
	geom_point( aes(x=2,y=medTo), data=dfMedSingle, shape=1, colour="blue" )+
	geom_segment(aes(x=1, xend=2, y=medFrom, yend=medTo ), data=dfMedSingle,
		colour="gray30", arrow=arrow(length=unit(0.18,"cm")))+
	scale_x_discrete(breaks=as.character(dfMed$PopPos), labels=as.vector(rbind(dfMedSingle$Pop,"")))+	
	scale_fill_manual("Records",value=c("white","aliceblue"), breaks=c("from","to"),labels=iLabel )+
	#scale_fill_brewer("Generations", breaks=c("from","to"),labels=iLabel,pal="Paired" )+
	coord_cartesian(ylim = {tmp<-c(min(dfMed$q25), max(dfLogLik2$value));tmp+diff(tmp)*0.05*c(-1,-1)},  )
	
	### List with components 
	list( ##describe<< 
		plot=p6,				##<< the ggplot object 
		rLogLik=if(doSort) sort(dfMedSingle$medTo,decreasing=TRUE) else dfMedSingle$medTo	##<< Median of Log-Likelihoods of last period
	)##end<< 
}

plotChainPopMoves <-function(	
	### Plot boxplots of the distribution of first and last fifth of steps for each population 
	resB			##<< the twDEMC to examine
	, iChains = rep(1:ncol(resB$temp), each=ncol(resB$rLogLik)%/%ncol(resB$temp))
	### mapping from chain to population
	, doSort=TRUE	##<< if TRUE result $rLogLik is sorted
	,...			##<< further arguements passed to matplot
	, xlab="Chains"
	, ylab="Log-Likelihood"
){
	##details<< 
	## \code{\link{ggplotChainPopMoves}} gives nicer results, but this functions i faster.
	
	##seealso<<  
	## \code{\link{plotMarginal2D}}
	## \code{\link{twDEMCInt}}
	
	iGen <- cbind( floor(c(0,1/5)*nrow(resB$rLogLik))+1, floor(c(4/5,1)*nrow(resB$rLogLik)) )
	#iLabel <- apply(iGen,2,function(iGeni){paste(range(iGeni), collapse=" to ")}) 
	rLogLik1 <- colMeans(twDEMCPopMeans(resB$rLogLik[iGen[,1],],ncol(resB$temp)))
	rLogLik2 <- colMeans(twDEMCPopMeans(resB$rLogLik[iGen[,2],],ncol(resB$temp)))
	
	matplot( cbind(rLogLik1,rLogLik2), pch=c(as.character(1:9),LETTERS), type="n", xlab=xlab, ylab=ylab, ... )
	points(rLogLik1, pch=c(as.character(1:9),LETTERS))
	points(rLogLik2, pch=c(as.character(1:9),LETTERS), col="red")
	arrows( 1:length(rLogLik1),rLogLik1,1:length(rLogLik1),rLogLik2, length=0.1)
	if( doSort) sort(rLogLik2, decr=TRUE) else rLogLik2
	
}

.tmp.f <- function(){
	matplot( cbind(rLogLik1,rLogLik2), pch=c(as.character(1:9),LETTERS), type="n", ... )
	points(rLogLik1, pch=c(as.character(1:9),LETTERS))
	points(rLogLik2, pch=c(as.character(1:9),LETTERS), col="red")
	arrows( 1:length(rLogLik1),rLogLik1,1:length(rLogLik1),rLogLik2, length=0.1)
	if( doSort) sort(rLogLik2, decr=TRUE) else rLogLik2
	
	iLabel <- apply(iGen,2,function(iGeni){paste(range(iGeni), collapse=" to ")}) 
	rLogLik1 <- colMeans(twDEMCPopMeans(resB$rLogLik[iGen[,1],],ncol(resB$temp)))
	rLogLik2 <- colMeans(twDEMCPopMeans(resB$rLogLik[iGen[,2],],ncol(resB$temp)))
	#tmp1 <- tapply(res$rLogLik[nrow(res$rLogLik)-nBack,], iChains, mean)
	#tmp2 <- tapply(res$rLogLik[nrow(res$rLogLik),], iChains, mean)
	#windows(); 
	dfLogLik <- data.frame(Chain=seq_along(rLogLik1), from=rLogLik1, to=rLogLik2) 
	p1 <- ggplot( dfLogLik, aes(x = Chain, y = from))+ xlab("Chain")+ ylab("Log-Likelihood")+
		geom_segment(aes(xend=Chain, yend = to),
			colour="gray40",	arrow=arrow(length=unit(0.15,"cm")))+
		geom_point(aes(colour=iLabel[1]))+
		geom_point(aes(y=to, colour=iLabel[2]))+
		scale_colour_manual("Generations",structure(c("black","red"),names=iLabel) )
	p1
	
	p2 <- ggplot( dfLogLik2, aes(x=PopPos,y=value) )+ xlab("Population") +ylab("Log-Likelihood") + 
		geom_boxplot(aes(fill=pos))+
		scale_colour_manual("Generations",breaks=c("from","to"),labels=iLabel )
	#scale_x_discrete(breaks=as.character(dfMed$PopPos), labels=as.vector(rbind(dfMedSingle$Pop,"")))+	
	#geom_segment(aes(x=(Pop-1)*2+1, xend=(Pop-1)*2+2, y=medFrom, yend=medTo ), data=dfMedSingle,
	#	colour="gray40", arrow=arrow(length=unit(0.15,"cm")))+
	#p3 <- p2 + ylim( min(dfMed$q25), max(dfLogLik2$value) )
	p3 <- p2 + coord_cartesian(ylim = {tmp<-c(min(dfMed$q25), max(dfLogLik2$value));tmp+diff(tmp)*0.05*c(-1,-1)},  )+
		p3	
	
	
}

plotThinned <- function (
	### Thin dataset before plotting.
	x, ...
	) UseMethod("plotThinned")

plotThinned.mcmc.list <- function(
	### Plot of thinned coda's \code{mcmc.list}.
	x,
	maxVars=3,	##<< maximum number of variable to put in one plot window
	vars=1:min(maxVars,ncol(x[[1]])), ##<< index of variables to plot
	maxN=500, 	##<< maximum number of cases to plot, (thin before)
	suppressPlotWarnings=TRUE,	##<< if true, plot command is executed with suppressWarnings
	...
){
	# plotThinned.mcmc.list
	
	##seealso<<  
	## \code{\link{plotMarginal2D}}
	## \code{\link{twDEMCInt}}
	
	##details<< 
	## default coda plot for mcmc.list squeezes all variables in one window
	
	# before plotting, take only up to maxVars variables and thin so that each chain has no more than maxN entries
	# may specify vars: the columns of x to plot
	thinFac <- max(1, ceiling(nrow(x[[1]]) / maxN )) 
	tmp.post <- window( x[,vars], thin=thin(x)*thinFac )
	if( suppressPlotWarnings )
		suppressWarnings(plot(tmp.post))
	else
		plot(tmp.post)
}





### ggplot2 statistics: adds columns prior and priorScaled to the output
StatPrior <- StatDensity$proto( calculate<- function(.,data, scales, poptDistr, doTransOrig=FALSE, ...){
	parname <- as.character(data$parName[1])
	range <- scales$x$output_set()
	xGrid <- seq(range[1], range[2], length = 100)
	#mtrace(dDistr)
	if( (parname %in% names(poptDistr$mu)) ){
		if( doTransOrig )
			prior <- dDistr(xGrid, mu=poptDistr$mu[parname], sigma=poptDistr$sigmaDiag[parname], trans=poptDistr$trans[parname])
		else
			prior <- dnorm(xGrid, mean=poptDistr$mu[parname], sd=poptDistr$sigmaDiag[parname])
		densdf <- data.frame(x=xGrid,prior=prior)
		densdf$priorScaled <- densdf$prior/max(densdf$prior, na.rm = TRUE)
		#xgrid <- seq( poptDistr$mu[parname]-sqrt(poptDistr$sigmaDiag[parname]), poptDistr$mu[parname]+sqrt(poptDistr$sigmaDiag[parname]), length.out=100)	
		#plot( dnorm(xgrid, mean=poptDistr$mu[parname], sd=sqrt(poptDistr$sigmaDiag[parname])) )
		#densdf$priorScaledOpt <- densdf$prior/ dnorm(poptDistr$mu[parname], mean=poptDistr$mu[parname], sd=sqrt(poptDistr$sigmaDiag[parname]))
		maxDen <- if( doTransOrig ) switch(as.character(poptDistr$trans[parname])
			,logitnorm = dlogitnorm(modeLogitnorm(mu=poptDistr$mu[parname], sigma=poptDistr$sigmaDiag[parname]),mu=poptDistr$mu[parname], sigma=poptDistr$sigmaDiag[parname])
			,lognorm = dlnorm(exp(poptDistr$mu[parname]-poptDistr$sigmaDiag[parname]^2), meanlog=poptDistr$mu[parname], sdlog=poptDistr$sigmaDiag[parname])
			,norm = dnorm(poptDistr$mu[parname], mean=poptDistr$mu[parname], sd=poptDistr$sigmaDiag[parname])
		) else dnorm(poptDistr$mu[parname], mean=poptDistr$mu[parname], sd=poptDistr$sigmaDiag[parname]) 
		densdf$priorScaledOpt <- densdf$prior/ maxDen
	}else{
		densdf <- data.frame(x=xGrid,prior=0,priorScaled=0,priorScaledOpt=0)
	}
	#qqplot( x, priorScaledOpt, data=densdf )
	#plot(densdf$x, densdf$priorScaledOpt )
	#if( doTransOrig ) 
	#	densdf$x <- transOrigPopt(xgrid,poptDistr$trans[parname])
	densdf
}, required_aes=c("x","parName"))
#with(StatPrior, mtrace(calculate)) 
#with(StatPrior, mtrace(calculate,F)) 

stat_prior <- function (
	### constructs a new StatPrior statistics based on aesthetics x and parName
	mapping = NULL, data = NULL, geom = "line", position = "stack", 
	adjust = 1,	poptDistr, doTransOrig=FALSE, ...
){ 
	StatPrior$new(mapping = mapping, data = data, geom = geom, 
		position = position, adjust = adjust, poptDistr=poptDistr, doTransOrig=doTransOrig, ...)
}
#stat_prior(

ggplotDensity.twDEMC <- function(
	### Plotting the densities for each parameter.
	res				##<< the twDEMC whose densities to plot
	,poptDistr=NULL	##<< parameter Distributions for the prior, usually \code{poptDistr <- \link{twConstrainPoptDistr}(poptNames,HamerParameterPriors$parDistr )}
	,pMin=0.05		##<< if > 0, the results are constrained to quantiles of rLogLik>percMin. Can avoid extremes
	,doTransOrig=FALSE	##<< if TRUE, parameters are translated to original scale
	,doDispLogLik=TRUE	##<< include density of LogLikelihoods
){
	##seealso<<  
	## \code{\link{plotMarginal2D}}
	## \code{\link{twDEMCInt}}
	
	nPop = 	ncol(res$temp)
	nChainsPop = ncol(res$rLogLik)%/%nPop
	#thin result to about 500 cases per constrainted population, to save calculation time
	resT <- resT0 <- thin(res, newThin=max(1,floor((nrow(res$rLogLik)*res$thin*nChainsPop*(1-pMin))/500.0)%/%res$thin)*res$thin )
	if( 0 < length(poptDistr) )
		poptDistr2 <- twConstrainPoptDistr(rownames(resT0$parms),poptDistr )	#also used for prior
	if( doTransOrig ){
		if(0==length(poptDistr)) stop("for doTransOrig one must provide argument poptDistr")
		resT <- transOrigPopt.twDEMC(resT0,poptDistr2$trans)
	}
	# stack populations
	pParms <- twDEMCPopApply( resT$parms, nPop, function(x){ abind(twListArrDim(x),along=2) })
	pLogLik <- twDEMCPopApply( resT$rLogLik, nPop, as.vector)
	tmp <- abind(pLogLik, pParms, along=1); 
	pTmp3 <- if( pMin > 0){
			# remove cases with lowest rLogLik
			nDrop <- floor(ncol(tmp)*pMin)
			abind( lapply( 1:nPop, function(iPop){
						iKeep <- order(pLogLik[,iPop,drop=TRUE])[-(1:nDrop)]
						twExtractFromLastDims(adrop(tmp[,,iPop,drop=FALSE],3),iKeep) 
					}), rev.along=0 )
		}else tmp
	dimnames(pTmp3)<-list( parms=c("rLogLik",rownames(res$parms)),steps=NULL, pops=NULL)
	#pTmp3 <- pTmp3[c("rLogLik","epsA","kS"),,][,,1:2]
	#pTmp3 <- pTmp3[c("rLogLik","tvr"),,][,,1:2]
	#poptDistr2 <- twConstrainPoptDistr(rownames(pTmp3)[-1],poptDistr )
	if( !doDispLogLik ){
		pTmp3 <- pTmp3[rowNames(pTmp3)!="rLogLik",,]		
	}		
	tmpDs4 <- melt(pTmp3)
	tmpDs4$pops <- as.factor(tmpDs4$pops)
	
	p1 <- p2 <- ggplot(tmpDs4,aes(x=value))+ #, colour=pops, fill=pops
		opts(axis.title.x = theme_blank())
	#print(p1 + geom_density(aes(y=..scaled..,colour=pops)))
	if( 0 < length(poptDistr) ){
		.nPop <- length(levels(tmpDs4$pops))
		cols <- c("gray40",scale_colour_hue(limits=1:.nPop)$output_set() )
		names(cols) <- c("prior",as.character(1:.nPop))
		p2 <- p1+ 
			stat_prior(aes(y=..priorScaledOpt..,parName=parms,colour="prior")
				#,data=tmpDs3	#will confuse the x-ranges
				,poptDistr=poptDistr2,doTransOrig = doTransOrig
				,size=0.8,linetype="twodash")+
			scale_colour_manual("Populations",cols )
	}
	p3 <- p2 + stat_density(aes(y=..scaled..,colour=pops), geom="line")+
		facet_wrap(~parms, scales="free_x")
	#p3 <- p2 + stat_density(aes(y=..density..,colour=pops), geom="line")+
	#		facet_wrap(~parms, scales="free")
	#p4 <- p3 + stat_prior(aes(y=..priorScaled..,parName=parms),poptDistr=poptDistr,col="blue",)
	p5 <- p3 + xlab("Parameter") + ylab("Scaled posterior density")
	p5
}

ggplotDensity.poptDistr <- function(
	### Plotting the densities for each parameter.
	poptDistr		##<< parameter Distributions for the prior, usually \code{poptDistr <- twConstrainPoptDistr(poptNames,HamerParameterPriors$parDistr )}
	,pMin=0.005		##<< range of the distribution from pMin to 1-pMin
	,parmsBounds=NULL	##<< list parName <- c(mode, upperQuantile)
	#,upperBoundProb = 0.99	##<< percentile of the upper bound
	,plotUpperQuantile=TRUE	##<< wheter to include upper quantile (set to FALSE if this inflates the displayed domain)
	,doTransOrig=TRUE		##<< set to FALSE to display transform to normal scale
){
	##seealso<<  
	## \code{\link{plotMarginal2D}}
	## \code{\link{twDEMCInt}}
	
	pRange = c(pMin,1-pMin)
	pNames <- names(poptDistr$mu) 
	#iPar=1
	qRange <- t(sapply( seq_along(poptDistr$mu), function(iPar){
			qNorm <- qnorm(pRange,mean=poptDistr$mu[iPar],sd=poptDistr$sigmaDiag[iPar])
			if(doTransOrig) transOrigPopt(qNorm, poptDistr$trans[iPar]) else qNorm
		}))
	dimnames(qRange)<-list(parms=pNames, iRec=NULL)
	tmpDs4 <- melt(qRange)
	p1 <- p2 <- ggplot(tmpDs4,aes(x=value),geom="blank")+ #, colour=pops, fill=pops
		opts(axis.title.x = theme_blank())
	p2 <- p1 + facet_wrap(~parms, scales="free_x") 
	p4 <- p2 + stat_prior(aes(y=..priorScaledOpt..,parName=parms),poptDistr=poptDistr,doTransOrig=doTransOrig)
	p5 <- p6 <- p4 + xlab("Parameter") + ylab("Scaled density")
	if( is.list(parmsBounds)){
		tmpDs5 <- if(plotUpperQuantile)
			melt(parmsBounds[names(poptDistr$mu)])
		else
			melt(lapply(parmsBounds[names(poptDistr$mu)],"[",1,drop=FALSE))
		colnames(tmpDs5) <- c("value","parms")
		if( !doTransOrig ) tmpDs5$value <- transNormPopt(tmpDs5$value, poptDistr$trans[tmpDs5$parms])
		p6 <- p5 + geom_vline(aes(xintercept=value), tmpDs5, color="blue") 
	}
	p6
}



