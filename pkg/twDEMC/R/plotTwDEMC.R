
plotThinned <- function (
	### Thin dataset before plotting.
	x		##<< object to be thinned
	, ...
	) UseMethod("plotThinned")

plotThinned.mcmc.list <- function(
	### Plot of thinned coda's \code{mcmc.list}.
	x
	,...
	,maxVars=3	##<< maximum number of variable to put in one plot window
	,vars=1:min(maxVars,ncol(x[[1]])) ##<< index of variables to plot
	,maxN=500 	##<< maximum number of cases to plot, (thin before)
	,suppressPlotWarnings=TRUE	##<< if true, plot command is executed with suppressWarnings
){
	# plotThinned.mcmc.list
	
	##seealso<<  
	## \code{\link{plotChainPopMoves}}
	## \code{\link{twDEMC}}
	
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


plotChainPopMoves <-function(	
	### Plot boxplots of the distribution of first and last fifth of steps for each population 
	resB			##<< the twDEMC to examine
	, iChains = rep(1:getNPops(resB), each=getNChainsPop(resB))
	### mapping from chain to population
	, doSort=TRUE	##<< if TRUE result $rLogDen is sorted
	,...			##<< further arguements passed to matplot
	, xlab="Chains"
	, ylab="LogDensity"
){
	##details<< 
	## \code{ggplotChainPopMoves} in twDEMCPlot gives nicer results, but this functions i faster.
	
	##seealso<<  
	## \code{\link{plotThinned.mcmc.list}}
	## \code{\link{twDEMC}}
	nPop <- getNPops(resB)
	nSamp <- getNSamples(resB) 
	iGen <- cbind( floor(c(0,1/5)*nSamp)+1, floor(c(4/5,1)*nSamp) )
	#iLabel <- apply(iGen,2,function(iGeni){paste(range(iGeni), collapse=" to ")}) 
	rLogDen1 <- colMeans(popMeansTwDEMC(resB$logDen[iGen[,1],1,],nPop))
	rLogDen2 <- colMeans(popMeansTwDEMC(resB$logDen[iGen[,2],1,],nPop))
	
	matplot( cbind(rLogDen1,rLogDen2), pch=c(as.character(1:9),LETTERS), type="n", xlab=xlab, ylab=ylab, ... )
	points(rLogDen1, pch=c(as.character(1:9),LETTERS))
	points(rLogDen2, pch=c(as.character(1:9),LETTERS), col="red")
	arrows( 1:length(rLogDen1),rLogDen1,1:length(rLogDen1),rLogDen2, length=0.1)
	### the logDensity of the last fifth of all the populations
	if( doSort) sort(rLogDen2, decreasing=TRUE) else rLogDen2
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


