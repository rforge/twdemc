# slowly decreasing Temperature (cost reduction factor)

twDEMCSA <- function(
	### simulated annealing DEMC 
	thetaPrior			##<< vector of parameters, point estimate
	,covarTheta			##<< the a prior covariance of parameters 
	,nChainPop=4		##<< number of chains within population
	,nPop=2				##<< number of populations
	,doIncludePrior=TRUE	##<< should the prior be part of initial population
	#
	,dInfos				##<< argument to \code{\link{twDEMCBlockInt}}
	,qTempInit=0.75		##<< quantile of logDensities used to calculate initial beginning and end temperature, with default 0.75 - 1/4 of the space is accepted 
	#
	,nObs				##<< integer vector (nResComp) specifying the number of observations for each result component
	,nGen=512			##<< number of generations in the initial batch
	,TFix=c(1,NA,NA)	##<< numeric vector (nResComp) specifying a finite value for components with fixed Temperatue, and a non-finite for others
	,nBatch=4			##<< number of batches with recalculated Temperature
	,maxRelTChange=0.025	##<< if Temperature of the components changes less than specified value, the algorithm can finish
	,... 				##<< further argument to \code{\link{twDEMCBlockInt}}
){
	#mtrace(initZtwDEMCNormal)
	Zinit0 <- initZtwDEMCNormal( thetaPrior, covarTheta, nChainPop=nChainPop, nPop=nPop, doIncludePrior=doIncludePrior)
	ss <- stackChains(Zinit0)
	
	logDenL <- lapply( dInfos, function(dInfo){
		resLogDen <- twCalcLogDenPar( dInfo$fLogDen, ss, argsFLogDen=dInfo$argsFLogDen)$logDenComp
	})
	# replace missing cases
	logDenDS <- abind( logDenL )
	boFinite <- apply(is.finite(logDenDS), 1, all )
	m0FiniteFac <- sum(boFinite) / nrow(logDenDS)
	if( m0FiniteFac < 0.9 ){
		## XXTODO extend Zinit0
	}
	Zinit <- Zinit0
	
	#------ intial temperatures: 
	#sapply( seq_along(expLogDenBest), function(i){ max(logDenDS[,i]) })
	#iResComp <- 
	temp0 <- sapply( seq_along(nObs), function(iResComp){
			pmax(1, -2/nObs[iResComp]* quantile( logDenDS[,iResComp],c(qTempInit) ))
	})
	names(temp0) <- colnames(logDenDS)
	iFixTemp <- which( is.finite(TFix) )
	temp0[iFixTemp] <- TFix[iFixTemp]
	print(paste("initial T=",paste(signif(temp0,2),collapse=","),"    ", date(), sep="") )
	
	res0 <-  res <- twDEMCBlock( Zinit
		, nGen=nGen0
		#,nGen=16, debugSequential=TRUE
		, dInfos=dInfos
		, TSpec=cbind( T0=temp0, TEnd=temp0 )
		, nPop=.nPop
		# XXTODO: replace lines below later on by ...
		#, blocks = blocks
		,...
	)
	.tmp.f <- function(){
		mc0 <- concatPops(res)
		matplot( mc0$temp, type="l" )
		matplot( mc0$pAccept[,1,], type="l" )
		matplot( mc0$pAccept[,2,], type="l" )
		plot( as.mcmc.list(mc0), smooth=FALSE )
	}
	TCurr <- getCurrentTemp(res)
	print(paste("finished ",1*nGen," out of ",nBatch*nGen," gens. T=",paste(signif(TCurr,2),collapse=" "),"    ", date(), sep="") )
	
	if( nBatch == 1 ) return(res0)
	#iBatch=2
	for(iBatch in (2:nBatch)){
		resEnd <- thin(res, start=getNGen(res)%/%2 )	# neglect the first half part
		ssc <- stackChainsPop(resEnd)	# combine all chains of one population
		mcl <- as.mcmc.list(ssc)
		#plot( as.mcmc.list(mcl), smooth=FALSE )
		#plot( res$pops[[1]]$parms[,"b",1] ~ res$pops[[1]]$parms[,"a",1], col=rainbow(255)[round(twRescale(-res$pops[[1]]$resLogDen[,"logDen1",1],c(10,200)))] )
		TEnd <- if( gelman.diag(mcl)$mpsrf <= 1.2 ){
				logDenT <- calcTemperatedLogDen(resEnd, TCurr)
				iBest <- getBestModelIndex( logDenT, resEnd$dInfos )
				maxLogDenT <- logDenT[iBest, ]
				TEnd <- pmax(1, -2*maxLogDenT/nObs )
				TEnd[iFixTemp] <- TFix[iFixTemp]
				relTChange <- abs(TEnd - TCurr)/TEnd
				if( max(relTChange) <= maxRelTChange ){
					print(paste("twDEMCSA: Maximum Temperture change only ",signif(max(relTChange)*100,2),"%. Finishing early.",sep=""))
					break
				}
				TEnd
			}else{
				TEnd <-TCurr
			}
		res1 <- res <- twDEMCBlock( resEnd
			, nGen=nGen0
			#, debugSequential=TRUE
			, dInfos=dInfos
			, TEnd=TEnd
			# replace lines below later on by ...
			#, blocks = blocks
			,...
		)
		TCurr <- getCurrentTemp(res)
		print(paste("finished ",iBatch*nGen," out of ",nBatch*nGen," gens. T=",paste(signif(TCurr,2),collapse=","),"    ", date(), sep="") )
	}
	res
}
attr(twDEMCSA,"ex") <- function(){
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	dInfos=list(
		dSparce=list(fLogDen=denSparcePrior, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
		,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
	)
	blocks = list(
		a=list(dInfoPos="dSparce", compPos="a")
		,b=list(dInfoPos="dRich", compPos="b")
	)
	
	do.call( dInfos$dSparce$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparce$argsFLogDen))
	do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen))
	
	#str(twTwoDenEx1)
	nObs <- c( parmsSparce=1, obsSparce=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )
	
	
	resPops <- res <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs, nBatch=8 )
	
	(TCurr <- getCurrentTemp(resPops))
	mc0 <- concatPops(res)
	logDenT <- calcTemperatedLogDen(stackChains(mc0$resLogDen), TCurr)
	iBest <- getBestModelIndex( logDenT, res$dInfos )
	maxLogDenT <- logDenT[iBest, ]
	ss <- stackChains(mc0$parms)
	(thetaBest <- ss[iBest, ])
	(.qq <- apply(ss,2,quantile, probs=c(0.025,0.5,0.975) ))
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(rowSums(logDenT),c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"obsSparce"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"logDen1"],c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rgb(
			twRescale(logDenT[,"obsSparce"]),0, twRescale(logDenT[,"logDen1"]) ))
	apply( apply( logDenT, 2, quantile, probs=c(0.1,0.9) ),2, diff )

	# density of parameters
	plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
	plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")
	
	# predictive posterior (best model only)
	pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
	plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
	plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 
	
}

.tmp.ggplotResults <- function(){
	mc0 <- concatPops(res)
	.nSample <- 128
	dfDen <- rbind(
		cbind( data.frame( scenario="S1", {tmp <- stackChains(res1)[,-(1:2)]; tmp[ seq(1,nrow(tmp),length.out=.nSample),] }))
	)
	dfDenM <- melt(dfDen)
	str(dfDenM)
	
	require(ggplot2)
	StatDensityNonZero <- StatDensity$proto( calculate<- function(.,data, scales, dMin=0.001, ...){
			res <- StatDensity$calculate(data,scales,...)
			# append a zero density at the edges
			if( res$density[1] > 0) res <- rbind( data.frame(x=res$x[1]-diff(res$x[1:2]),density=0,scaled=0,count=0),res)
			if( res$density[nrow(res)] > 0) res <- rbind( res, data.frame(x=res$x[nrow(res)]+diff(res$x[nrow(res)-c(1,0)]),density=0,scaled=0,count=0))
			bo <- (res$density < dMin) & (c(res$density[-1],0) < dMin) & (c(0,res$density[-length(res$density)]) < dMin)
			res$density[ bo ] <- NA
			res$count[ bo ] <- NA
			res$scaled[ bo ] <- NA
			res
		}, required_aes=c("x"))
#with(StatDensityNonZero, mtrace(calculate)) 
#with(StatDensityNonZero, mtrace(calculate,F)) 
	
	stat_densityNonZero <- function (
		### constructs a new StatPrior statistics based on aesthetics x and parName
		mapping = NULL, data = NULL, geom = "line", position = "stack", 
		adjust = 1,	...
	){ 
		StatDensityNonZero$new(mapping = mapping, data = data, geom = geom, 
			position = position, adjust = adjust, ...)
	}
	
#pgm <- geom_ribbon( alpha=0.8, aes(ymax = ..density.., ymin = -..density..), stat = "density")
#pgm <- geom_ribbon(  aes(ymax = ..density.., ymin = -..density..), stat = "density")
#pgm <- stat_densityNonZero(  aes(ymax = -..density..),  size=1, geom="line")
	pgm <- geom_line(aes(y = +..scaled..), stat="densityNonZero", size=1)
	pgm2 <- geom_line(aes(y = -..scaled..), stat="densityNonZero", size=1) 
	optsm <- opts(axis.title.x = theme_blank(), axis.text.y = theme_blank(), axis.ticks.y = theme_blank() ) 
	pa <- ggplot(dfDen, aes(x = a, colour=scenario, linetype=scenario)) + pgm + pgm2 + optsm + scale_y_continuous('a')
	pb <- ggplot(dfDen, aes(x = b, colour=scenario, linetype=scenario)) + pgm + pgm2 + optsm + scale_y_continuous('b')
#pb
	windows(width=7, height=3)
	grid.newpage()
	pushViewport( viewport(layout=grid.layout(2,1)))
	print(pa , vp = viewport(layout.pos.row=1,layout.pos.col=1))	
	print(pb + opts(legend.position = "none") , vp = viewport(layout.pos.row=2,layout.pos.col=1))	
	
	#----------- ggplot predictive posterior
	#scenarios <- c("R","RS","RSw","DG","DM")
	scenarios <- c("R")
	nScen <- length(scenarios)
	# infer quantiles of predictions
	.nSample=128
	y1M <- array( NA_real_, dim=c(.nSample, length(twTwoDenEx1$obs$y1),nScen), dimnames=list(sample=NULL,iObs=NULL,scenario=scenarios)  )
	y2M <- array( NA_real_, dim=c(.nSample, length(twTwoDenEx1$obs$y2),nScen), dimnames=list(sample=NULL,iObs=NULL,scenario=scenarios) )
	resScen <- #list( res1, res2, res2b, res3a, res3 ); names(resScen) <- scenarios
	resScen <- list( res1 ); names(resScen) <- scenarios
	#scen <- "R"
	for( scen in scenarios ){
		ss <- stackChains(resScen[[scen]])[,-(1:getNBlocks(resScen[[scen]]))]
		ssThin <- ss[round(seq(1,nrow(ss),length.out=.nSample)),]
		#i <- .nSample
		for( i in 1:.nSample){
			pred <-  twTwoDenEx1$fModel(ssThin[i,], xSparce=twTwoDenEx1$xSparce, xRich=twTwoDenEx1$xRich, thresholdCovar=thresholdCovar) 
			y1M[i,,scen] <- pred$y1
			y2M[i,,scen] <- pred$y2
		}
	}
	predQuantilesY1 <- lapply( scenarios, function(scen){
			t(apply( y1M[,,scen],2, quantile, probs=c(0.025,0.5,0.975) ))
		}); names(predQuantilesY1) <- scenarios
	predQuantilesY2 <- lapply( scenarios, function(scen){
			t(apply( y2M[,,scen],2, quantile, probs=c(0.025,0.5,0.975) ))
		}); names(predQuantilesY2) <- scenarios
	str(predQuantilesY1)
	rm( y1M, y2M)	# free space, we only need the summary
	
	#str(pred1)
	dfPred <- rbind(
		cbind( data.frame( scenario="R", variable="y1", value = pred1$y1 ), predQuantilesY1[["R"]] )
		,cbind( data.frame( scenario="R", variable="y2", value = pred1$y2 ), predQuantilesY2[["R"]] )
	)
	dfPred$observations <- NA_real_
	dfPred$observations[dfPred$variable=="y1"] <- twTwoDenEx1$obs$y1
	dfPred$observations[dfPred$variable=="y2"] <- twTwoDenEx1$obs$y2
	colnames(dfPred)[ match(c("2.5%","50%","97.5%"),colnames(dfPred)) ] <- c("lower","median","upper")
	str(dfPred)
	
	p1 <- ggplot(dfPred, aes(x=value, y=observations, colour=scenario) ) +
		geom_errorbarh(aes(xmax = upper, xmin = lower)) +
		geom_point() +
		facet_wrap( ~ variable, scales="free") +
		opts(axis.title.x=theme_blank() ) +
		geom_abline(colour="black") +
		c()
	p1
	
	
}

getBestModelIndex <- function(
	### select the best model based on (temperated) logDensity components
	logDenT		##<< numeric matrix (nStep x nResComp): logDensity (highest are best)
	, dInfos 	##<< list of lists with entry resCompPos (integer vector) specifying the position of result components for each density 
){
	iBest <- if( length(dInfos) > 1){
		# with several densities, each parameter vector is ranked differently
		# select the case where the maximum of all the ranks across densities is lowest
		logDenTDens <- do.call( cbind, 	lapply( seq_along(dInfos), function(iDen){
					dInfo <- dInfos[[iDen]]
					logDenInfoT <- rowSums(logDenT[,dInfo$resCompPos ,drop=FALSE])	
				}))
		rankDen <- apply( -logDenTDens,2, rank) # starting with the highest density
		rankDenMax <- apply( rankDen, 1, max )   
		iBest <- which.min(rankDenMax)
	}else{
		iBest <- which.max( rowSums(logDenT) )
	}
	### the index within logDenT with the best rank
	iBest
}
attr(getBestModelIndex,"ex") <- function(){
	logDenT <- cbind( -sample(5)/2, -sample(5), -sample(5) )
	#dInfos <- list( d1=list(resCompPos=1:2), d2=list(resCompPos=3) )
	dInfos <- list( d1=list(resCompPos=2), d2=list(resCompPos=3) )
	getBestModelIndex(logDenT, dInfos)
	-logDenT
}

