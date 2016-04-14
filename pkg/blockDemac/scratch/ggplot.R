.tmp.f <- function(){
	subvector<-subset(vector, Devel==0,
		select = c(Code, GDP50, GDP80,
			GDP2000,PRIM50,PRIM80,PRIM2000))
	
	p1 <- ggplot(subvector, aes(x = GDP50, y = PRIM50))
	
	p1 <- p1 + xlab("GDP")
	p1 <- p1 + ylab("Primacy")
	
	p1 + geom_segment(aes(xend = GDP80, yend = PRIM80),
		colour="blue",linetype=5,
		arrow=arrow(length=unit(0.15,"cm")))
	
	last_plot() + geom_segment(aes(xend = GDP2000, yend = PRIM2000),
		colour="dark red", size=0.72,
		arrow=arrow(length=unit(0.15,"cm")))
	
	last_plot() + geom_segment(aes(x = GDP80, y = PRIM80,
			xend = GDP2000, yend = PRIM2000),
		colour="blue",linetype=4, size=0.75,
		arrow=arrow(length=unit(0.15,"cm")))
	
	last_plot() + geom_text(aes(x=GDP2000, y=PRIM2000, label=Code),
		data=subvector, vjust=1, hjust=0,
		size=3) 
}

ggplotDensity.twDEMC <- function(
	### Plotting the densities for each parameter.
	res				##<< the twDEMC whose densities to plot
	,poptDistr=NULL	##<< parameter Distributions for the prior, usually \code{poptDistr <- twConstrainPoptDistr(poptNames,HamerParameterPriors$parDistr )}
	,pMin=0.05		##<< if > 0, the results are constrained to quantiles of rLogDen>percMin. Can avoid extremes
	,doTransOrig=FALSE	##<< if TRUE, parameters are translated to original scale
){
	nPop = 	ncol(res$temp)
	nChainPop = ncol(res$rLogDen)%/%nPop
	#thin result to about 200 cases per constrainted population, to save calculation time
	resT <- resTStart <- thin(res, newThin=(floor((nrow(res$rLogDen)*res$thin*nChainPop*(1-pMin))/200.0)%/%res$thin)*res$thin )
	if( doTransOrig )
		resT <- transOrigPopt.twDEMC(resTStart) 
	pLogDen <- popApplyTwDEMC( resT$rLogDen, nPop, as.vector)
	# include rLogDen as variable
	tmp <- abind(resT$rLogDen, resT$parms, along=1); dimnames(tmp)[[1]][1] <- "rLogDen"; names(dimnames(tmp))<-names(dimnames(resT$parms))
	#rownames(tmp)[1] <- c("rLogDen",rownames(resT$parms))
	# stack populations
	#mtrace(popApplyTwDEMC)
	pTmp <- popApplyTwDEMC( tmp, nPop, function(x){ abind(twListArrDim(x),along=2) })
	dimnames(pTmp) <- c(dimnames(tmp)[1:2], list(pops=NULL)) 
	pTmp3 <- if( pMin > 0){
			# remove cases with lowest rLogDen
			#popQuantiles <- as.vector(popApplyTwDEMC( resT$rLogDen, nPop, quantile, probs=pMin ))
			#tmpDs2 <- ddply(tmpDs, .(pop), function(df,pop)subset, rLogDen>popQuantiles[.(pop)] )
			nDrop <- round(ncol(pTmp)*pMin)
			#aaply takes very long: pTmp2 <- aaply( pTmp, 3, function(A){A[, -(order(A["rLogDen",])[1:nDrop]) ]})
			abind( lapply( 1:nPop, function(iPop){pTmp[,-(order(pTmp["rLogDen",,iPop])[1:nDrop]),iPop] }), rev.along=0 )
		}else pTmp
	dimnames(pTmp3)<-dimnames(pTmp)
	#pTmp3 <- pTmp3[c("rLogDen","epsA","kS"),,][,,1:2]
	#poptDistr2 <- twConstrainPoptDistr(rownames(pTmp3)[-which(rownames(pTmp3)=="rLogDen")],poptDistr )
	poptDistr2 <- twConstrainPoptDistr(rownames(pTmp3)[-1],poptDistr )
	tmpDs3 <- tmpDs4 <- melt(pTmp3)
	tmpDs3$pops <- as.factor(tmpDs3$pops)
	if( doTransOrig ){
		origMat <- function(pop,A){ t(transOrigPopt(t(adrop(A[,,pop,drop=FALSE],3)),poptDistr2$trans))}
		#mtrace(transOrigPopt.matrix)
		#origMat(1,pTmp3[-which(rownames(pTmp3)=="rLogDen"),,])
		pTmp4 <- abind( lapply(1:dim(pTmp3)[3], origMat, A=pTmp3[-which(rownames(pTmp3)=="rLogDen"),,,drop=FALSE]), rev.along=0)
		pTmp4 <- abind( adrop(pTmp3["rLogDen",,,drop=FALSE],1), pTmp4, along=1) 
		dimnames(pTmp4)<-dimnames(pTmp3)
		tmpDs4 <- melt(pTmp4)
	} 
	tmpDs4$pops <- as.factor(tmpDs4$pops)
	
	#popTrans <- matrix( 1:ncol(resT$rLogDen), ncol=ncol(resT$temp) )
	p1 <- p2 <- ggplot(tmpDs4,aes(x=value))+ #, colour=pops, fill=pops
		scale_colour_hue("Populations")
	#p1 + geom_density(aes(y=..scaled..,colour=pops))
	if( 0 < length(poptDistr2) ){
		.nPop <- length(levels(tmpDs4$pops))
		cols <- c("gray40",scale_colour_hue(limits=1:.nPop)$output_set() )
		names(cols) <- c("prior",as.character(1:.nPop))
		p2 <- p1+ 
			stat_prior(aes(y=..priorScaledOpt..,parName=parms,colour="prior")
				#,data=tmpDs3	#will confuse the x-ranges
				,poptDistr=poptDistr2,doTransNorm = doTransOrig
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

