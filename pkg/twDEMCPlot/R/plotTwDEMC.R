ggplotChainPopMoves <- function(
	### Plot boxplots of the distribution of first and last fifth of steps for each population 
	resB			##<< the twDEMC to examine
	, iChains = rep(1:ncol(resB$temp), each=ncol(resB$rLogDen)%/%ncol(resB$temp))
		### mapping from chain to population
	, doSort=TRUE	##<< if TRUE result $rLogDen is sorted
){
	# ggplotChainPopMoves
	##seealso<<  
	## \code{\link{plotMarginal2D}}
	## \code{\link{twDEMCBlockInt}}
	
	iGen <- cbind( floor(c(0,1/5)*nrow(resB$rLogDen))+1, floor(c(4/5,1)*nrow(resB$rLogDen)) )
	iLabel <- apply(iGen,2,function(iGeni){paste(range(iGeni), collapse=" to ")}) 
	tmp1 <- reshape2::melt(resB$rLogDen[iGen[1,1]:iGen[2,1],]); tmp1$pos="from"
	tmp2 <- reshape2::melt(resB$rLogDen[iGen[1,2]:iGen[2,2],]); tmp2$pos="to"
	
	dfLogDen2 <- rbind(tmp1,tmp2)
	dfLogDen2$Pop <- iChains[dfLogDen2$X2]  #(dfLogDen2$X2-1) %% ncol(resB$temp) +1
	dfLogDen2 <- dfLogDen2[order(dfLogDen2$Pop,dfLogDen2$pos),]
	dfLogDen2$PopPos <- paste(dfLogDen2$Pop,dfLogDen2$pos, sep="_")
	dfLogDen2$PopPos <- factor( dfLogDen2$PopPos, levels=unique(dfLogDen2$PopPos) )
	
	dfMed <- aggregate(dfLogDen2$value,list(Pop=dfLogDen2$Pop,pos=dfLogDen2$pos), median)
	#dfMed$q25 <- aggregate(dfLogDen2$value,list(Pop=dfLogDen2$Pop,pos=dfLogDen2$pos), function(values){ quantile(values, probs=0.25)})
	dfMed$q25 <- aggregate(dfLogDen2$value,list(Pop=dfLogDen2$Pop,pos=dfLogDen2$pos), quantile, probs=0.25)$x
	dfMed <- dfMed[order(dfMed$Pop,dfMed$pos),]
	dfMed$PopPos <- paste(dfMed$Pop,dfMed$pos, sep="_")
	dfMed$PopPos <- factor( dfMed$PopPos, levels=unique(dfMed$PopPos) )
	dfMed$medFrom <- dfMed$x 
	dfMed$medTo <- c(dfMed$x[-1],NA)
	dfMedSingle <- subset(dfMed, pos=="from")
	
	p5 <- ggplot( dfLogDen2, aes(x=pos,y=value) )+ xlab("Population") +ylab("LogDensity") +	facet_grid(.~Pop)
	p6 <- p5+
	geom_boxplot(aes(fill=pos))+
	geom_point( aes(x=1,y=medFrom), data=dfMedSingle, shape=1 )+
	geom_point( aes(x=2,y=medTo), data=dfMedSingle, shape=1, colour="blue" )+
	geom_segment(aes(x=1, xend=2, y=medFrom, yend=medTo ), data=dfMedSingle,
		colour="gray30", arrow=arrow(length=unit(0.18,"cm")))+
	scale_x_discrete(breaks=as.character(dfMed$PopPos), labels=as.vector(rbind(dfMedSingle$Pop,"")))+	
	scale_fill_manual("Records",value=c("white","aliceblue"), breaks=c("from","to"),labels=iLabel )+
	#scale_fill_brewer("Generations", breaks=c("from","to"),labels=iLabel,pal="Paired" )+
	coord_cartesian(ylim = {tmp<-c(min(dfMed$q25), max(dfLogDen2$value));tmp+diff(tmp)*0.05*c(-1,-1)},  )
	
	### List with components 
	list( ##describe<< 
		plot=p6,				##<< the ggplot object 
		rLogDen=if(doSort) sort(dfMedSingle$medTo,decreasing=TRUE) else dfMedSingle$medTo	##<< Median of LogDensitys of last period
	)##end<< 
}


.tmp.f <- function(){
	matplot( cbind(rLogDen1,rLogDen2), pch=c(as.character(1:9),LETTERS), type="n", ... )
	points(rLogDen1, pch=c(as.character(1:9),LETTERS))
	points(rLogDen2, pch=c(as.character(1:9),LETTERS), col="red")
	arrows( 1:length(rLogDen1),rLogDen1,1:length(rLogDen1),rLogDen2, length=0.1)
	if( doSort) sort(rLogDen2, decr=TRUE) else rLogDen2
	
	iLabel <- apply(iGen,2,function(iGeni){paste(range(iGeni), collapse=" to ")}) 
	rLogDen1 <- colMeans(popMeansTwDEMC(resB$rLogDen[iGen[,1],],ncol(resB$temp)))
	rLogDen2 <- colMeans(popMeansTwDEMC(resB$rLogDen[iGen[,2],],ncol(resB$temp)))
	#tmp1 <- tapply(res$rLogDen[nrow(res$rLogDen)-nBack,], iChains, mean)
	#tmp2 <- tapply(res$rLogDen[nrow(res$rLogDen),], iChains, mean)
	#windows(); 
	dfLogDen <- data.frame(Chain=seq_along(rLogDen1), from=rLogDen1, to=rLogDen2) 
	p1 <- ggplot( dfLogDen, aes(x = Chain, y = from))+ xlab("Chain")+ ylab("LogDensity")+
		geom_segment(aes(xend=Chain, yend = to),
			colour="gray40",	arrow=arrow(length=unit(0.15,"cm")))+
		geom_point(aes(colour=iLabel[1]))+
		geom_point(aes(y=to, colour=iLabel[2]))+
		scale_colour_manual("Generations",structure(c("black","red"),names=iLabel) )
	p1
	
	p2 <- ggplot( dfLogDen2, aes(x=PopPos,y=value) )+ xlab("Population") +ylab("LogDensity") + 
		geom_boxplot(aes(fill=pos))+
		scale_colour_manual("Generations",breaks=c("from","to"),labels=iLabel )
	#scale_x_discrete(breaks=as.character(dfMed$PopPos), labels=as.vector(rbind(dfMedSingle$Pop,"")))+	
	#geom_segment(aes(x=(Pop-1)*2+1, xend=(Pop-1)*2+2, y=medFrom, yend=medTo ), data=dfMedSingle,
	#	colour="gray40", arrow=arrow(length=unit(0.15,"cm")))+
	#p3 <- p2 + ylim( min(dfMed$q25), max(dfLogDen2$value) )
	p3 <- p2 + coord_cartesian(ylim = {tmp<-c(min(dfMed$q25), max(dfLogDen2$value));tmp+diff(tmp)*0.05*c(-1,-1)},  )+
		p3	
}

### ggplot2 statistics: adds columns prior and priorScaled to the output
StatPrior <- {
	require(ggplot2)
	proto::proto( ggplot2:::StatDensity, {
	objname <- "prior"
	calculate<- function(.,data, scales, poptDistr, doTransOrig=FALSE, nGrid=30, ...){
		parname <- as.character(data$parName[1])
		#range <- scales$x$output_set()
		range <- scales$x$range$range
		xGrid <- seq(range[1], range[2], length = nGrid)
        if( is.data.frame(poptDistr) ){ 
            orig <- poptDistr
            poptDistr <- as.list(poptDistr)
            poptDistr <- lapply( poptDistr, function(item){ structure(item, names=rownames(orig))})
        }
		#mtrace(dDistr)
		if( (parname %in% names(poptDistr$mu)) ){			if( doTransOrig )
				prior <- dDistr(xGrid, mu=poptDistr$mu[parname], sigma=poptDistr$sigmaDiag[parname], trans=poptDistr$trans[parname])
			else
				prior <- dnorm(xGrid, mean=poptDistr$mu[parname], sd=poptDistr$sigmaDiag[parname])
			densdf <- data.frame(x=xGrid,prior=prior)
			densdf$priorScaled <- densdf$prior/max(densdf$prior, na.rm = TRUE)
			#xgrid <- seq( poptDistr$mu[parname]-sqrt(poptDistr$sigmaDiag[parname]), poptDistr$mu[parname]+sqrt(poptDistr$sigmaDiag[parname]), length.out=100)	
			#plot( dnorm(xgrid, mean=poptDistr$mu[parname], sd=sqrt(poptDistr$sigmaDiag[parname])) )
			#densdf$priorScaledOpt <- densdf$prior/ dnorm(poptDistr$mu[parname], mean=poptDistr$mu[parname], sd=sqrt(poptDistr$sigmaDiag[parname]))
			maxDen <- if( doTransOrig ) switch( as.character(poptDistr$trans[parname])
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
	}
	required_aes=c("x","parName")
})}
#with(StatPrior, mtrace(calculate)) 
#with(StatPrior, mtrace(calculate,F)) 


stat_prior <- function (
	### constructs a new StatPrior statistics based on aesthetics x and parName
	mapping = NULL, data = NULL, geom = "line", position = "stack", 
	adjust = 1,	poptDistr, doTransOrig=FALSE, nGrid=100, ...
){ 
	StatPrior$new(mapping = mapping, data = data, geom = geom, 
		position = position, adjust = adjust, poptDistr=poptDistr, doTransOrig=doTransOrig, nGrid=nGrid, ...)
}
#stat_prior(

ggplotDensity.twDEMC <- function(
	### Plotting the densities for each parameter.
	res				##<< the twDEMC whose densities to plot
	,poptDistr=NULL	##<< parameter Distributions for the prior, usually \code{poptDistr <- \link{twConstrainPoptDistr}(poptNames,HamerParameterPriors$parDistr )}
	,pMin=0.05		##<< if > 0, the results are constrained to quantiles of rLogDen>percMin. Can avoid extremes
	,doTransOrig=FALSE	##<< if TRUE, parameters are translated to original scale
	,doDispLogDen=TRUE	##<< include density of LogDensitys
    ,idInfo=names(res$dInfos)[1]
){
	##seealso<<  
	## \code{\link{plotMarginal2D}}
	## \code{\link{twDEMCBlockInt}}
	#
	nPop = 	ncol(res$temp)
	nChainsPop = ncol(res$rLogDen)%/%nPop
	#thin result to about 500 cases per constrainted population, to save calculation time
	resT <- resT0 <- thin(res, newThin=max(1,floor((nrow(res$rLogDen)*res$thin*nChainsPop*(1-pMin))/500.0)%/%res$thin)*res$thin )
	if( 0 < length(poptDistr) )
		poptDistr2 <- twConstrainPoptDistr(rownames(resT0$parms),poptDistr )	#also used for prior
	if( doTransOrig ){
		if(0==length(poptDistr)) stop("for doTransOrig one must provide argument poptDistr")
		resT <- transOrigPopt.twDEMC(resT0,poptDistr2$trans)
	}
	# each population one single chain
    m <- stackChainsPop( resT )
    mA <- stackChains(m)    # all samples
    # remove cases with 
    m <- m[ m[,idInfo] >=  quantile(m[,idInfo],pMin) ,   ]
    # pick columns
    .colNames <- c( {if(doDispLogDen) idInfo else NULL}, colnames(res$parms) )
	tmpDs4 <- structure( reshape2::melt(m[,.colNames]), names=c("pops","parms","value"))
	tmpDs4$pops <- as.factor(tmpDs4$pops)
	#
	p1 <- p2 <- ggplot(tmpDs4,aes(x=value))+ ylim(0,1) +#, colour=pops, fill=pops
		theme(axis.title.x = element_blank())
	#print(p1 + geom_density(aes(y=..scaled..,colour=pops)))
	if( 0 < length(poptDistr) ){
		.nPop <- length(levels(tmpDs4$pops))
		cols <- c("gray40",scale_colour_hue(limits=1:.nPop)$output_set() )
		names(cols) <- c("prior",as.character(1:.nPop))
		p2 <- p1+ 
			stat_prior(aes(y=..priorScaledOpt..,parName=parms,nGrid=30, colour="prior")
				#,data=tmpDs3	#will confuse the x-ranges
				,poptDistr=poptDistr2,doTransOrig = doTransOrig
				,size=0.8,linetype="twodash")+
			scale_colour_manual("Populations",cols )
	}
	p3 <- p2 + stat_density(aes(y=..scaled..,colour=pops,ymax=1), geom="line")+
		facet_wrap(~parms, scales="free_x") +
        theme()
	#p3 <- p2 + stat_density(aes(y=..density..,colour=pops), geom="line")+
	#		facet_wrap(~parms, scales="free")
	#p4 <- p3 + stat_prior(aes(y=..priorScaled..,parName=parms),poptDistr=poptDistr,col="blue",)
	p5 <- p3 + xlab("Parameter") + ylab("Scaled posterior density")
	p5
}
attr( ggplotDensity.twDEMC, "ex") <- function(){
    data(twdemcEx1)  # from package twDEMC
    res <- concatPops(twdemcEx1)
    ggplotDensity.twDEMC( res )
}

ggplotDensity.twDEMCPops <- function(
        ### Plotting the densities for each parameter.
        res				##<< the twDEMCPops object whose densities to plot
        ,poptDistr=NULL	##<< parameter Distributions for the prior, usually \code{poptDistr <- \link{twConstrainPoptDistr}(poptNames,HamerParameterPriors$parDistr )}
        ,pMin=0.05		##<< if > 0, the results are constrained to quantiles of rLogDen>percMin. Can avoid extremes
        ,doTransOrig=FALSE	##<< if TRUE, parameters are translated to original scale
        ,doDispLogDen=TRUE	##<< include density of LogDensitys
        ,idInfo=names(res$dInfos)[1]
        ,nSamplesPop=500      ##<< thin to about these number of samples within each population
        ,popNames=as.character(seq_along(res$pops)) ##<< character vector (nPop): names of the populations displayed in colour legend
        ,popCols=scale_colour_hue(1:.nPop)$palette(.nPop)   ##<< colors of the populations
        ,legendTitle="Populations"
){
    ##seealso<<  
    ## \code{\link{plotMarginal2D}}
    ## \code{\link{twDEMCBlockInt}}
    #
    nPop = 	length(res$pops)
    # one chain per population
    resPops <- lapply( 1:nPop, function(iPop){ stackChainsPop(subPops(res,iPop=iPop)) })
    resPopsI <- resPops[[1]]
    parNames <- colnames(resPopsI$pops[[1]]$parms)
    # thin result to about 500 cases per constrainted population, to save calculation time and retrieve matrix
    resMPops <- lapply(resPops, function(mPopsI){
                stackChains(
                thin(mPopsI, newThin=max(1,
                    floor( (nrow(mPopsI$pops[[1]]$parms)*mPopsI$thin*(1-pMin))/nSamplesPop)  %/% mPopsI$thin
                            ) * mPopsI$thin
                    ))
            })
    # pick columns
    .colNames <- c( {if(doDispLogDen) idInfo else NULL}, parNames )
    # omit the cases of very low density
    #minDen <- quantile( do.call( c, lapply(resMPops,"[",,idInfo)), pMin)
    resM2Pops <- lapply( resMPops, function(m){
                m[ m[,idInfo] >=  quantile(m[,idInfo],pMin) ,  .colNames ]      # better for each pop separately
            })
    if( 0 < length(poptDistr) )
        poptDistr2 <- twConstrainPoptDistr(parNames,poptDistr )	#also used for prior
    if( doTransOrig ){
        if(0==length(poptDistr)) stop("for doTransOrig one must provide argument poptDistr")
        for( iPop in seq_along(resM2Pops) ){
            resM2Pops[[iPop]][,-1] <- transOrigPopt(resM2Pops[[iPop]][,-1], poptDistr2$trans)
        }
    }
    tmpDs4 <- structure( reshape2::melt(resM2Pops) , names=c("cases","parms","value","pops"))
    tmpDs4$pops <- as.factor(tmpDs4$pops)
    .nPop <- length(levels(tmpDs4$pops))
    #
    p1 <- p2 <- ggplot(subset(tmpDs4, TRUE),aes(x=value,colour=pops))#+ ylim(0,1) +#, colour=pops, fill=pops
            theme(axis.title.x = element_blank())
    #print(p1 + geom_density(aes(y=..scaled..,colour=pops)))
    if( 0 < length(poptDistr) ){
        #1:.nPop
        p2 <- p1+ 
                stat_prior(aes(y=..priorScaledOpt..,parName=parms,linetype="prior",ymax=1) #, position="identity"
                        #,data=tmpDs3	#will confuse the x-ranges
                        ,poptDistr=poptDistr2,doTransOrig = doTransOrig
                        ,size=I(0.8)
                        #,linetype=I("twodash")
                        , colour="gray40"
                        ) 
    }
    #cols <- structure( c("gray40",popCols  ), names=c("prior",seq_along(res$pops)) )
    cols <- structure( popCols  , names=seq_along(res$pops) )
    lts <- structure( c("twodash", rep("solid",nPop) ), names=c("prior",seq_along(res$pops)) )
    p3 <- p2 + stat_density(aes(y=..scaled..,ymax=1), geom="line", position = "identity")+
            facet_wrap(~parms, scales="free_x") +
            #scale_colour_manual(legendTitle,values=cols, labels=c(popNames,"prior") )+
            scale_colour_manual(legendTitle,values=cols, labels=popNames )+
            scale_linetype_manual("",values=lts, labels=c("prior",popNames) )+
            theme()
    #p3 <- p2 + stat_density(aes(y=..density..,colour=pops), geom="line")+
    #		facet_wrap(~parms, scales="free")
    #p4 <- p3 + stat_prior(aes(y=..priorScaled..,parName=parms),poptDistr=poptDistr,col="blue",)
    p5 <- p3 + xlab("Parameter") + ylab("Scaled posterior density")
    p5
}
attr( ggplotDensity.twDEMCPops, "ex") <- function(){
    data(twdemcEx1)  # from package twDEMC
    res <- twdemcEx1
    ggplotDensity.twDEMCPops( res )
    ggplotDensity.twDEMCPops( res, doDispLogDen=FALSE )
}


ggplotDensity.poptDistr <- function(
	### Plotting the densities for each parameter.
	parDistr		##<< data.frame (nParm x c(mu, sigmaDiag)): parameter Distributions for the prior, usually \code{poptDistr <- twConstrainPoptDistr(poptNames,HamerParameterPriors$parDistr )}
	,pMin=0.005		##<< range of the distribution from pMin to 1-pMin
	,parmsBounds=NULL	##<< numeric vector (nParm x nLines) to be plottes as lines
	,doTransOrig=TRUE		##<< set to FALSE to display transform to normal scale
){
	##seealso<<  
	## \code{\link{plotMarginal2D}}
	## \code{\link{twDEMCBlockInt}}
	
	pRange = c(pMin,1-pMin)
	pNames <- rownames(parDistr) 
	#iPar=1
	qRange <- t(sapply( seq_along(parDistr$mu), function(iPar){
			qNorm <- qnorm(pRange,mean=parDistr$mu[iPar],sd=parDistr$sigmaDiag[iPar])
			if(doTransOrig) transOrigPopt(qNorm, parDistr$trans[iPar]) else qNorm
		}))
	dimnames(qRange)<-list(parms=pNames, iRec=NULL)
	tmpDs4 <- reshape2::melt(qRange)
	p1 <- p2 <- ggplot(tmpDs4,aes(x=value),geom="blank")+ #, colour=pops, fill=pops
		theme(axis.title.x = element_blank())
	p2 <- p1 + facet_wrap(~parms, scales="free_x") 
	#p4 <- p2 + stat_prior(aes(y=..priorScaledOpt..,parName=parms),poptDistr=poptDistr,doTransOrig=doTransOrig)
	p4 <- p2 + stat_prior(aes(y=..priorScaled..,parName=parms),poptDistr=parDistr,doTransOrig=doTransOrig)
	p5 <- p6 <- p4 + xlab("Parameter") + ylab("Scaled density")
	if( 0 < length(parmsBounds) ){
		tmpDs5 <- reshape2::melt(parmsBounds[pNames,])
		colnames(tmpDs5) <- c("parms","quantile","value")
		if( !doTransOrig ) tmpDs5$value <- transNormPopt(tmpDs5$value, parDistr[as.character(tmpDs5$parms),"trans"])
		p6 <- p5 + geom_vline(aes(xintercept=value), tmpDs5, color="blue") 
	}
	p6
}
attr(ggplotDensity.poptDistr,"ex") <- function(){
	parmsBounds = do.call( rbind,list(		# mode and upper bound
		A0 = c(10,15)		
		,D0 = c(10, 100)
		,C0 = c(0.6,0.8)
	))
	colnames(parmsBounds) <- c("mode","upper")
	varDistr <- twVarDistrVec( rownames(parmsBounds) )	# by default assumed normal
	varDistr["D0"] <- "lognorm"
	varDistr["C0"] <- "logitnorm"
	parDistr <- twQuantiles2Coef( parmsBounds, varDistr, upperBoundProb=0.975 )

	#mtrace(ggplotDensity.poptDistr)
	ggplotDensity.poptDistr( parDistr, parmsBounds=parmsBounds )	
	ggplotDensity.poptDistr( parDistr, parmsBounds=parmsBounds, doTransOrig=FALSE )	
	
}



