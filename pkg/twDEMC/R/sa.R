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
	,TFix=numeric(0)	##<< numeric vector (nResComp) specifying a finite value for components with fixed Temperatue, and a non-finite for others
	,nBatch=4			##<< number of batches with recalculated Temperature
	,maxRelTChange=0.025	##<< if Temperature of the components changes less than specified value, we do not need further batches because of Temperature
	,maxLogDenDrift=0.3	##<< if difference between mean logDensity of first and fourth quartile of the sample is less than this value, we do not need further batches because of drift in logDensity
	,restartFilename	##<< filename to write intermediate results to 
	,... 				##<< further argument to \code{\link{twDEMCSA}}
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
	# now the legnth of resComp is known (colnames(logDenDS)), so check the TFix and TMax parameters
	nResComp <- ncol(logDenDS)
	if( 0 == length( TFix) ) TFix <- structure( rep( NA_real_, nResComp), names=colnames(logDenDS) )
	if( nResComp != length( TFix) ) stop("twDEMCSA: TFix must be of the same length as number of result Components.")
	iFixTemp <- which( is.finite(TFix) )
	#if( 0 == length( TMaxInc) ) TMaxInc <- structure( rep( NA_real_, nResComp), names=colnames(logDenDS) )
	#if( nResComp != length( TMaxInc) ) stop("twDEMCSA: TMaxInc must be of the same length as number of result Components.")
	#iMaxIncTemp <- which( is.finite(TMaxInc) )
	#
	#------ intial temperatures: 
	#sapply( seq_along(expLogDenBest), function(i){ max(logDenDS[,i]) })
	#iResComp <- 
	temp0 <- sapply( seq_along(nObs), function(iResComp){
			pmax(1, -2/nObs[iResComp]* quantile( logDenDS[,iResComp],c(qTempInit) ))
	})
	names(temp0) <- colnames(logDenDS)
	temp0[iFixTemp] <- TFix[iFixTemp]
	print(paste("initial T=",paste(signif(temp0,2),collapse=","),"    ", date(), sep="") )
	#
	res0 <-  res <- twDEMCBlock( Zinit
		, nGen=nGen
		#,nGen=16, debugSequential=TRUE
		, dInfos=dInfos
		, TSpec=cbind( T0=temp0, TEnd=temp0 )
		, nPop=nPop
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
		tmp <- calcTemperatedLogDen(res$pops[[1]]$resLogDen[,,1], getCurrentTemp(res) )
		matplot( tmp, type="l" )
		plot(rowSums(tmp))
		#
		matplot( mc0$resLogDen[,3,], type="l" )
		matplot( mc0$parms[,"a",], type="l" )
		matplot( mc0$parms[1:20,"b",], type="l" )
		bo <- 1:10; iPop=1
		plot( mc0$parms[bo,"a",iPop], mc0$parms[bo,"b",iPop], col=rainbow(100)[twRescale(mc0$resLogDen[bo,"parmsSparce",iPop],c(10,100))] )
	}
	print(paste("finished initial ",1*nGen," out of ",nBatch*nGen," gens. T=",paste(signif(temp0,2),collapse=" "),"    ", date(), sep="") )
	#
	if( nBatch == 1 ) return(res0)
	twDEMCSACont( mc=res0, nObs=nObs, nGen=nGen, TFix=TFix, nBatch=nBatch-1, maxRelTChange = maxRelTChange, maxLogDenDrift=maxLogDenDrift, restartFilename=restartFilename, ...)
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
	
	#trace(twDEMCSACont, recover )
	resPops <- res <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs
		, TFix=c(1,NA,NA)
		, nGen=256
		, nBatch=5 
		, restartFilename=file.path("tmp","example_twDEMCSA.RData")
	)

	(TCurr <- getCurrentTemp(resPops))
	mc0 <- concatPops(res)
	plot( as.mcmc.list(mc0) , smooth=FALSE )
	matplot( mc0$temp, type="l" )
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

twDEMCSACont <- function(
	### continuing simulated annealing DEMC based on previous reslt
	mc					##<< result of twDEMCBlock
	,nObs				##<< integer vector (nResComp) specifying the number of observations for each result component
	,nGen=512			##<< number of generations in the initial batch, default 512
	,TFix=numeric(0)	##<< numeric vector (nResComp) specifying a finite value for components with fixed Temperatue, and a non-finite for others
	,nBatch=4			##<< number of batches with recalculated Temperature
	,maxRelTChange=0.025 ##<< if Temperature of the components changes less than specified value, the algorithm can finish
	,maxLogDenDrift=0.3		##<< if difference between mean logDensity of first and fourth quartile of the sample is less than this value, we do not need further batches because of drift in logDensity
	,restartFilename=NULL	##<< filename to write intermediate results to 
	,... 				##<< further argument to \code{\link{twDEMCBlockInt}}
){
	# save calling arguments to allow an continuing an interrupted run
	args <- c( list( nObs=nObs, nGen=nGen, TFix=TFix, maxRelTChange=maxRelTChange, maxLogDenDrift=maxLogDenDrift,  restartFilename=restartFilename), list(...) )
	#
	nResComp <- ncol(mc$pops[[1]]$resLogDen)
	if( 0 == length( TFix) ) TFix <- structure( rep( NA_real_, nResComp), names=colnames(logDenDS) )
	if( nResComp != length( TFix) ) stop("twDEMCSA: TFix must be of the same length as number of result Components.")
	iFixTemp <- which( is.finite(TFix) )
	#
	res <- mc
	.tmp.f <- function(){
		mc0 <- concatPops(res)
		matplot( mc0$temp, type="l" )
		matplot( mc0$pAccept[,1,], type="l" )
		matplot( mc0$pAccept[,2,], type="l" )
		plot( as.mcmc.list(mc0), smooth=FALSE )
		tmp <- calcTemperatedLogDen(res$pops[[1]]$resLogDen[,,1], getCurrentTemp(res) )
		matplot( tmp, type="l" )
		plot(rowSums(tmp))
		#
		matplot( mc0$resLogDen[,3,], type="l" )
		matplot( mc0$parms[,"a",], type="l" )
		matplot( mc0$parms[1:20,"b",], type="l" )
		bo <- 1:10; iPop=1
		plot( mc0$parms[bo,"a",iPop], mc0$parms[bo,"b",iPop], col=rainbow(100)[twRescale(mc0$resLogDen[bo,"parmsSparce",iPop],c(10,100))] )
	}
	# sum nObs within density
	iDens <- seq_along(mc$dInfos)
	#iDen=1
	iCompsNonFixDen <- lapply( iDens, function(iDen){ 
			irc <-  mc$dInfos[[iDen]]$resCompPos
			irc <- irc[ !(irc %in% iFixTemp) ]
		})
	nObsDen <- sapply( iDens, function(iDen){ sum( nObs[iCompsNonFixDen[[iDen]] ]) })
	TCurr <- getCurrentTemp(mc)
	iBatch=1
	for(iBatch in (1:nBatch)){
		if((0 < length(restartFilename)) && is.character(restartFilename) && restartFilename!=""){
			resRestart.twDEMCSA = res #avoid variable confusion on load by giving a longer name
			resRestart.twDEMCSA$iBatch <- iBatch	# also store updated calculation of burnin time
			resRestart.twDEMCSA$args <- args		# also store updated calculation of burnin time
			save(resRestart.twDEMCSA, file=restartFilename)
			cat(paste("Saved resRestart.twDEMCSA to ",restartFilename,"\n",sep=""))
		}
		resEnd <- thin(res, start=getNGen(res)%/%2 )	# neglect the first half part
		ssc <- stackChainsPop(resEnd)	# combine all chains of one population
		mcl <- as.mcmc.list(ssc)
		#plot( as.mcmc.list(mcl), smooth=FALSE )
		#plot( res$pops[[1]]$parms[,"b",1] ~ res$pops[[1]]$parms[,"a",1], col=rainbow(255)[round(twRescale(-res$pops[[1]]$resLogDen[,"logDen1",1],c(10,200)))] )
		TEnd <- if( gelman.diag(mcl)$mpsrf <= 1.2 ){
				logDenT <- calcTemperatedLogDen(resEnd, TCurr)
				#mtrace(getBestModelIndex)
				iBest <- getBestModelIndex( logDenT, resEnd$dInfos )
				maxLogDenT <- logDenT[iBest, ]
				#thetaBest <- stackChains(concatPops(resEnd)$parms)[iBest, ]
				#calcTemperatedLogDen( do.call( dInfos[[1]]$fLogDen, c(list(thetaBest), dInfos[[1]]$argsFLogDen) ), TCurr )
				#sum logDenT within density
				#logDenTDen <- lapply( iDens, function(iDen){
				#			rcpos <- iCompsNonFixDen[[iDen]]
				#			if( length(rcpos)==1) logDenT[,rcpos ,drop=TRUE] else						
				#				rowSums( logDenT[ ,rcpos ,drop=FALSE])
				#		})
				maxLogDenTDen <- sapply(iDens, function(iDen){ 
						sum(maxLogDenT[iCompsNonFixDen[[iDen]] ])
					})
				#TEnd0 <- pmax(1, -2*maxLogDenT/nObs )
				TEndDen <- pmax(1, -2*maxLogDenTDen/nObsDen )
				TEnd <- TCurr
				for( iDen in iDens){
					TEnd[ resEnd$dInfos[[iDen]]$resCompPos ] <- TEndDen[iDen]
				} 
				TEnd[iFixTemp] <- TFix[iFixTemp]
				relTChange <- abs(TEnd - TCurr)/TEnd
				#if( (max(relTChange) <= maxRelTChange) ) recover()
				#trace(isLogDenDrift, recover )
				if( (max(relTChange) <= maxRelTChange) && !isLogDenDrift(logDenT, resEnd$dInfos, maxDrift=maxLogDenDrift) ){
					res <- resEnd
					print(paste("twDEMCSA: Maximum Temperture change only ",signif(max(relTChange)*100,2),"% and no drift in logDensity. Finishing early.",sep=""))
					break
				}
				TEnd
			}else{
				TEnd <-TCurr
			}
		#TMin <- pmin(TMin, TEnd)
		res1 <- res <- twDEMCBlock( resEnd
			, nGen=nGen
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

.tmp.ggplotResults <- function(){
	mc0 <- concatPops(res)
	.nSample <- 128
	dfDen <- rbind(
		cbind( data.frame( scenario="S1", {tmp <- stackChains(mc0)[,-(1:getNBlocks(mc0))]; tmp[ seq(1,nrow(tmp),length.out=.nSample),] }))
	)
	dfDenM <- melt(dfDen)
	#str(dfDenM)
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
	resScen <- list( mc0 ); names(resScen) <- scenarios
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

.tmp.AllDensities <- function(){ # does not work
	# same as example but with each parameter updated against both densities
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	# for updating each in turn against both densities
	dInfos=list(
		dSparce=list(fLogDen=denSparcePrior, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
		,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
	)
	blocks = list(
		bSparce=list(dInfoPos="dSparce", compPos=c("a","b") )
		,bRich=list(dInfoPos="dRich", compPos=c("a","b") )
	)
	
	do.call( dInfos$dSparce$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparce$argsFLogDen))
	do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen))
	
	#str(twTwoDenEx1)
	nObs <- c( parmsSparce=1, obsSparce=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )
	
	#trace(twDEMCSA, recover )
	resPops <- res <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs, nBatch=5, debugSequential=TRUE )
	
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

.tmp.AllSparce <- function(){ # 
	# same as example but with b also updated against sparce, a only against sparce 
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	# for updating each in turn against both densities
	dInfos=list(
		dSparce=list(fLogDen=denSparcePrior, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
		,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
	)
	blocks = list(
		bSparce=list(dInfoPos="dSparce", compPos=c("a","b") )
		,bRich=list(dInfoPos="dRich", compPos=c("a") )
	)
	
	do.call( dInfos$dSparce$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparce$argsFLogDen))
	do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen))
	#str(twTwoDenEx1)
	nObs <- c( parmsSparce=1, obsSparce=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )
	
	#untrace(twDEMCSA )
	#trace(twDEMCSA, recover )
	resPops <- res <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs, nBatch=5, debugSequential=TRUE, TFix=c(1,NA,NA) )
	
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


.tmp.2DenSwitched <- function(){ # does not work
	# same as example but with each parameter b updated against sparce
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	# for updating each in turn against both densities
	dInfos=list(
		dSparce=list(fLogDen=denSparcePrior, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta))
		,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar,twTwoDenEx=twTwoDenEx1, theta0=thetaPrior))
	)
	blocks = list(
		bSparce=list(dInfoPos="dSparce", compPos=c("b") )
		,bRich=list(dInfoPos="dRich", compPos=c("a") )
	)
	
	do.call( dInfos$dSparce$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dSparce$argsFLogDen))
	do.call( dInfos$dRich$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos$dRich$argsFLogDen))
	
	#str(twTwoDenEx1)
	nObs <- c( parmsSparce=1, obsSparce=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )
	
	#untrace(twDEMCSA )
	#trace(twDEMCSA, recover )
	resPops <- res <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, blocks=blocks, nObs=nObs, nBatch=5, debugSequential=TRUE )
	
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


.tmp.oneDensity <- function(){
	# same as example but with only one combined density for both parameters
	data(twTwoDenEx1)
	
	thetaPrior <- twTwoDenEx1$thetaTrue
	covarTheta <- diag((thetaPrior*0.3)^2)
	invCovarTheta <- (thetaPrior*0.3)^2		# given as independent variances for faster calculation
	
	thresholdCovar = 0.3	# the true value
	thresholdCovar = 0		# the effective model that glosses over this threshold
	
	# for only one logDensity - Temperature of the strongest component goes to zero and other increase according to mismatch
	 dInfos=list( list(fLogDen=denBoth, argsFLogDen=list(thresholdCovar=thresholdCovar, twTwoDenEx=twTwoDenEx1, theta0=thetaPrior, thetaPrior=thetaPrior, invCovarTheta=invCovarTheta)) )
	
	do.call( dInfos[[1]]$fLogDen, c(list(theta=twTwoDenEx1$theta0),dInfos[[1]]$argsFLogDen))
	#str(twTwoDenEx1)
	nObs <- c( parmsSparce=1, y1=length(twTwoDenEx1$obs$y1), y2=length(twTwoDenEx1$obs$y2) )
	
	#trace(twDEMCSA, recover)
	resPops <- res <- twDEMCSA( thetaPrior, covarTheta, dInfos=dInfos, nObs=nObs
		,TFix = c(1,NA,NA)
		,nGen=256
		#,TMax = c(NA,1.2,NA)	# do not increase sparce observations again too much
		, nBatch=5 
		, debugSequential=TRUE
	)
	
	(TCurr <- getCurrentTemp(res))
	mc0 <- concatPops(res)
	logDenT <- calcTemperatedLogDen(stackChains(mc0$resLogDen), TCurr)
	iBest <- getBestModelIndex( logDenT, res$dInfos )
	maxLogDenT <- logDenT[iBest, ]
	ss <- stackChains(mc0$parms)
	(thetaBest <- ss[iBest, ])
	(.qq <- apply(ss,2,quantile, probs=c(0.025,0.5,0.975) ))
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(rowSums(logDenT),c(1,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"y1"],c(10,100))] )
	plot( ss[,"a"], ss[,"b"], col=rainbow(100)[twRescale(logDenT[,"y2"],c(10,100))] )
	plot( ss[,"a"], ss[,"b"], col=rgb(
			twRescale(logDenT[,"y1"]),0, twRescale(logDenT[,"y2"]) ))
	apply( apply( logDenT, 2, quantile, probs=c(0.2,0.8) ),2, diff )
	
	# density of parameters
	plot( density(ss[,"a"])); abline(v=thetaPrior["a"]); abline(v=thetaBest["a"], col="blue")
	plot( density(ss[,"b"])); abline(v=thetaPrior["b"]); abline(v=thetaBest["b"], col="blue")
	
	# predictive posterior (best model only)
	pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
	plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
	plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 
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

isLogDenDrift <- function(
	### check whether first quartile all the logDensities is significantly smaller than last quartile 
	logDenT		##<< numeric matrix (nStep x nResComp): logDensity (highest are best)
	, dInfos 	##<< list of lists with entry resCompPos (integer vector) specifying the position of result components for each density
	, alpha=0.05	##<< the significance level for a difference
	, maxDrift=0.3	##<< difference in LogDensity, below which no drift is signalled
){
	##details<<
	## Because of large sample sizes, very small differences may be significantly different.
	## Use argument minDiff to specify below which difference a significant difference is not regarded as drift.
	nr <- nrow(logDenT)
	nr4 <- nr %/% 4
	#iDen <- 1
	boL <- sapply( seq_along(dInfos), function(iDen){
			dInfo <- dInfos[[iDen]]
			rs1 <- rowSums(logDenT[1:nr4,  dInfo$resCompPos ,drop=FALSE])
			rs4 <- rowSums(logDenT[(nr+1-nr4):nr, dInfo$resCompPos ,drop=FALSE])
			resTTest <- t.test(rs1,rs4,"less")
			bo <- (diff(resTTest$estimate) >= maxDrift ) && (resTTest$p.value <= alpha) 
		})
	any( boL )
	### TRUE if any of the logDensities are a significantly greater in the fourth quantile compared to the first quantile of the samples
}
attr(isLogDenDrift,"ex") <- function(){
	data(twdemcEx1)
	logDenT <- calcTemperatedLogDen( twdemcEx1, getCurrentTemp(twdemcEx1))
	isLogDenDrift(logDenT, twdemcEx1$dInfos )
}


twRunDEMCSA <- function(
	### wrapper to twDEMCSA compliant to runCluster.R
	argsTwDEMCSA	 
	### Arguments passed to twDEMCSA -> twDEMCBlockInt
	### It is updated by \dots.
	,...				 	##<< further arguments passed to twDEMCSA -> twDEMCBlockInt
	,prevResRunCluster=NULL	##<< results of call to twRunDEMCSA, argument required to be called from runCluster.R
	,restartFilename=NULL	##<< name of the file to store restart information, argument required to be called from runCluster.R 
){
	#update argsDEMC to args given 
	argsDEMC <- argsTwDEMCSA
	# update the restartFilename argument
	.dots <- list(...)
	argsDEMC[ names(.dots) ] <- .dots
	#for( argName in names(.dots) ) argsDEMC[argName] <- .dots[argName]
	if( 0 != length(restartFilename)) argsDEMC$restartFilename <- restartFilename
	
	# do the actual call
	res <- do.call( twDEMCSA, argsDEMC )
}
tmp.testRestart <- function(){
	rFileName <- file.path("tmp","example_twDEMCSA.RData")
	load(rFileName)	# resRestart.twDEMCSA
	#untrace(twDEMCSACont)
	#trace(twDEMCSACont, recover )
	res1 <- do.call( twDEMCSACont, c( list(mc=resRestart.twDEMCSA), resRestart.twDEMCSA$args ))
}





