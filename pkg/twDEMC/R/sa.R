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
	,qTempInit=c(0.4,0.05)	##<< quantile of logDensities used to calculate initial beginning and end temperature
	#
	,expLogDenBest		##<< numeric vector (nResComp) specifying the expected logDen of the best model
	,nGen0=1024			##<< number of generations in the initial batch
	
	
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
	
	#------ intial temperatures: 0.4 percentile of rLogDen
	#sapply( seq_along(expLogDenBest), function(i){ max(logDenDS[,i]) })
	temp0 <- sapply( seq_along(expLogDenBest), function(iResComp){
		pmax(1, quantile( (expLogDenBest[iResComp] - logDenDS[,iResComp]),c(qTempInit) ))
	})
	colnames(temp0) <- colnames(logDenDS)
	temp0[,"parmsSparce"] <- 1
	
	res0 <-  res <- twDEMCBlock( Zinit
		, nGen=nGen0
		#,nGen=16, debugSequential=TRUE
		, dInfos=dInfos
		, TSpec=cbind( T0=temp0[1,], TEnd=temp0[2,] )
		, nPop=.nPop
		# XXTODO: replace lines below later on by ...
		, blocks = blocks
	)
	.tmp.f <- function(){
		mc0 <- concatPops(res)
		matplot( mc0$temp, type="l" )
		matplot( mc0$pAccept[,1,], type="l" )
		matplot( mc0$pAccept[,2,], type="l" )
		plot( as.mcmc.list(mc0), smooth=FALSE )
	}
	resEnd <- thin(res, start=getNGen(res)%/%2 )	# neglect the first half part
	ssc <- stackChainsPop(resEnd)	# combine all chains of one population
	mcl <- as.mcmc.list(ssc)
	#plot( mcl, smooth=FALSE )
	#plot( mc0$parms[,"b",1] ~ mc0$parms[,"a",1], col=rainbow(255)[round(twRescale(-mc0$resLogDen[,"logDen1",1],c(10,200)))] )
	#plot( mc0$parms[,"b",1] ~ mc0$parms[,"a",1], col=rainbow(255)[round(twRescale(-mc0$resLogDen[,"obsSparce",1],c(10,200)))] )
	#plot( res$pops[[1]]$parms[,"b",1] ~ res$pops[[1]]$parms[,"a",1], col=rainbow(255)[round(twRescale(-res$pops[[1]]$resLogDen[,"logDen1",1],c(10,200)))] )
	#plot( res$pops[[1]]$parms[,"b",1] ~ res$pops[[1]]$parms[,"a",1], col=rainbow(255)[round(twRescale(-res$pops[[1]]$resLogDen[,"obsSparce",1],c(10,200)))] )
	TEnd <- if( gelman.diag(mcl)$mpsrf <= 1.2 ){
			.tmp.estimateLMax.Raftery <- function(){
				# estimate lmax according to Raftery07: l_{max} = \bar{l} + Var{l}
				# does not work out, varL is much too high for strongly skewed distr.
				# need to thin because of autocorrlation before estimating the variance
				sscT <- squeeze( ssc, length.out=effectiveSize(mcl) )
				logDen <- stackChains(concatPops(sscT)$logDen)
				varL <- apply(logDen,2,var)
				meanL <- apply(logDen,2,mean)
				lMax1 <- meanL + varL
			}
			#logDen0 <- stackChains(concatPops(res0)$resLogDen); sapply( seq_along(expLogDenBest), function(i){ max(logDen0[,i]) })
			logDen <- stackChains(concatPops(ssc)$resLogDen)
			#iResComp <- 3
			#sapply( seq_along(expLogDenBest), function(i){ max(logDen[,i]) })
			#sapply( seq_along(expLogDenBest), function(i){ max(logDen[,i]) })
			#mc0$temp[ nrow(mc0$temp), ]
			TEnd <- sapply( seq_along(expLogDenBest), function(iResComp){
					pmax(1, (expLogDenBest[iResComp] - max(logDen[,iResComp])) )
				})
			names(TEnd) <- colnames(logDen)
			TEnd
		}else{
			# stay at given Temperature
			TEnd <-getCurrentTemp(ssc)
		}
	
	res1 <- res <- twDEMCBlock( resEnd
		, nGen=nGen0
		#, debugSequential=TRUE
		, dInfos=dInfos
		, TEnd=TEnd
		# XXTODO: replace lines below later on by ...
		, blocks = blocks
	)
	

	
	
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
	expLogDenBest <- -1/2 * c( parmsSparce=0, obsSparce=length(twTwoDenEx1$obs$y1), logDen1=length(twTwoDenEx1$obs$y2) )
	
	
	#SADEMC(
	

	res <- res1
	logDen <- stackChains(concatPops(res)$resLogDen)
	TCurr <- getCurrentTemp(res)
	logDenT <- apply( logDen, 1, function(logDenRow){ sum(logDenRow / TCurr)} )
	iDen <- 1:getNBlocks(res1)
	
	ss <- stackChains(res1)
	(thetaBest <- thetaBest1 <- ss[ which.max(sLogDen), -iDen])		# sum of logDen
	(thetaBest <- thetaBest1 <- ss[ which.max(ss[,2]), -iDen])		# par b
	(thetaBest <- thetaBest1 <- ss[ which.max(logDenT), -iDen])		# Temperated sum of lodDen
	(.qq <- apply(ss[,-iDen],2,quantile, probs=c(0.025,0.5,0.975) ))
	
	pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
	plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) 
	plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1)
	
	plot( density(ss[,"b"])); abline(v=thetaPrior["b"])
	
}

.tmp.f <- function(){
	mc0 <- concatPops(res1)
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

