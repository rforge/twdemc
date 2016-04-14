# TODO: Add comment
# 
# Author: twutz
###############################################################################


twDEMCBatchInt <- function(
	### Calls \code{\link{twDEMCBlockInt}} successively with \code{nBatch} generations.
	Zinit
	### the initial population
	, nGen=10
	### number of generations in total
	, nBatch=512
	### number of generations between saving results
	, ...
	### further arguments passed to twDEMC
	, restartFilename=NULL
	### name of the file to save intermediate results resRestart.twDEMC, NULL for not saving
	, fCheckConvergence=function(res, addArgs){ FALSE }
	### checking convergence of a DEMC and interrupting 
	, fCheckConvergenceArgs=list()
	### additional arguments to the DEMC convergence checking function
	, TStart=1
	### initial temperature of burnin , defaults to 1 or if Zinit is twDEMC to the temperature of the last row
	, nGenBurnin=0
	### number of generations of burnin (Temperature decrease to 0, and probUpDirBurnin)
	, doResetOutlierN=256
	### if > 0, outlier chains are reset to best chain
	, nPops = 1
	### number of independent populations
	, probUpDirBurnin=0.75
	### probUbDir during burnin (see twDEMC argument propuUpDir)
	, controlTwDEMC=list()
	### controls to twDEMC, some items are overwritten
	#,minPCompAcceptTempDecr=0.2	
	### if maximum di=Lpi-Lai drops below this rate 
	### then Temperature is not decreased within the next batch.
	, doRepeatLowAcceptanceChains=TRUE
	, maxNGenBurnin=50000	##<< maximum burnin beyond which can not be extendend on too low acceptance rate
){
	# twDEMCBatchInt
	##seealso<<   
	## \code{\link{twDEMCBlockInt}}
	## \code{\link{twDEMCBatch}}
	##details<< 
	## Usually invoked by \code{\link{twDEMCBatch}}.
	if( !hasArg(nBatch))nBatch=512
	if( !hasArg(restartFilename))restartFilename=NULL
	cl <- match.call()
	for( i in (2:length(cl)) )
		try( cl[i] <- list(eval.parent(cl[[i]])) )
	#cl[ names(cl)[-1] ] <- lapply(as.list(cl)[-1],eval.parent) 	#substitute all variables by their values, -1 needed for not substituting the function name
	iRun <- 0    #ialready completed generations
	nRun <- min(nBatch, (if(iRun<nGenBurnin) min(nGenBurnin,nGen) else nGen) -iRun)
	
	ctrl <- controlTwDEMC
	ctrl$TStart=TStart
	ctrl$Tend=TStart 	#no Temp decrease in first batch (if Zinit is not twDEMC see below) 
	ctrl$probUpDir=(if(nRun<=nGenBurnin) probUpDirBurnin else NULL)
	minAccepRateTempDecrease <- minPCompAcceptTempDecr <- if(is.numeric(ctrl$minPCompAcceptTempDecr)) ctrl$minPCompAcceptTempDecr else 0.16
	TFix <- if(is.numeric(ctrl$TFix)) ctrl$TFix else numeric(0) 
	thin <- if(is.numeric(ctrl$thin)) 	ctrl$thin else 1
	nGen <- (nGen %/% thin)*thin
	
	##details<< 
	## If Zinit is of class twDEMC, initial temperature is set to the temperature of the last row
	## and the number of generations already in Zinit are skipped.
	if( is(Zinit,"twDEMC") ){
		res <- Zinit
		iRun <- (nrow(res$rLogDen)-1)*res$thin
		if( iRun >= nGen ) return(res)
		ctrl$initialAcceptanceRate <- popMeansTwDEMC( res$pAccept[nrow(res$pAccept),],nPops )
		ctrl$TStart <- TStartc <- res$temp[ nrow(res$temp), ,drop=FALSE ]
		if( length(TStartc) != nPops) stop(paste("twDEMCBlockInt: encoutered temperature recored with",length(TStartc),"columns but argument nPops=",nPops))
		
		#calculate optimal end temperature
		resCols <- match( rownames(res$logDenComp), rownames(res$Y))
		#nPops <- ncol(res$temp)
		nChains <- dim(res$parms)[3]
		nChainPop <- nChains %/% nPops
		chain2Pop <- rep(1:nPops, each=nChainPop )	#mapping of chain to population
		
		diffLogDen <- getDiffLogDen.twDEMCPops(res$Y, resCols, nLastSteps=ceiling(128/nChainPop)) 	#in twDEMC S3twDEMC.R
		diffLogDenPops <- popApplyTwDEMC( diffLogDen, nPops=nPops, function(x){ abind(twListArrDim(x),along=2,new.names=dimnames(x)) })	#stack param columns by population
		
		pAcceptChains <- res$pAccept[ nrow(res$pAccept), ]
		pAcceptPops <- tapply( pAcceptChains, chain2Pop, mean) 
		#boPopCoolingTooFast <- (pAcceptPops < minAccepRateTempDecrease)
		boPopCoolingTooFast <- (pAcceptPops < minPCompAcceptTempDecr)
		
		.tmp.f <- function(){
			#mtrace(.calcTemperatedDiffLogDen)
			diffLogDenPopsT <- .calcTemperatedDiffLogDen( diffLogDenPops, TFix, TStartc)
			# acceptance rate per component		
			tmpPercAcc <- 1-apply(diffLogDenPopsT,c(1,3),function(d){ ecdf(d)(log(0.5))} ) #comp x pops
			boPopCoolingTooFast <- (apply(tmpPercAcc,2,min) < minPCompAcceptTempDecr)
		}
		
		if( any(boPopCoolingTooFast) ){
			nGenBurnin=min(maxNGenBurnin,nGenBurnin+nRun)
			#recalculate nRun with changed nGenBurnin
			nRun <- min(nBatch, (if(iRun<nGenBurnin) min(nGenBurnin,nGen) else nGen) -iRun)		
		}
		pTarget=minPCompAcceptTempDecr+0.02
		temp <- TStartc	# keep Temperature, only for others cool down further
		if( any(!boPopCoolingTooFast)) temp[!boPopCoolingTooFast] <- {
				tempExp <- calcDEMCTemp( TStartc[!boPopCoolingTooFast], 1, nGenBurnin-iRun, nRun) #recalculate with initial temperature
				#mtrace(calcDEMCTempDiffLogDenConst)
				#tempEmp <- sapply( seq_along(TStartc)[!boPopCoolingTooFast], function(iPop){ calcDEMCTempDiffLogDenConst(diffLogDenPops[,,iPop , drop=FALSE], TFix=ctrl$TFix, Tmax=TStartc[iPop], pTarget=pTarget)})
				#cool faster than tempExp but give more than 1/3 weight to the empirical estimate to avoid much too fast cooling 
				#pmax(1, pmin( tempExp,(2*tempExp+tempEmp)/3 ))
				tempExp
			}
		ctrl$Tend <- temp
		
		##details<< \describe{\item{Temperature estimate from proposal distribution}{
		## The distribution of differences between Density of proposals Lp and of accepted state La
		## can be used to estimate an optimal temperature per data stream, so that each
		## datastream contributes to rejections in about the same magnitude and the overall
		## acceptance rate is aobut a specified value.
		## The proportions of the so calculated datastream specific temperature are multiplied 
		## with the global temperature on Metropolis decisions.
		## Further, if the temperatures of the datasteams are all below the goal of the global
		## temperature Tend, Tend is also lowered.
		## }}
		if( (0<length(ctrl$useMultiT)) ) if( ctrl$useMultiT ){
				ctrl$Tprop <- TDiffLogDen <- sapply( seq_along(TStartc), function(i){calcDEMCTempDiffLogDen2(diffLogDenPops[,,i], pTarget=pTarget, TFix=TFix, Tmax=TStartc[i])})  # will be scaled in twDEMCBlockInt
				#ctrl$Tend <- pmin( ctrl$Tend, apply(TDiffLogDen,2,max))
			}
	}
	
	#--------- do the twDEMC ------
	cat(paste(iRun," out of ",nGen," generations completed. T=",paste({T<-ctrl$TStart;round(T,digits=ifelse(T>20,0,1))},collapse="  "),"     ",date(),"\n",sep=""))
	res <- twDEMC( Zinit=Zinit, nGen=nRun, nPops=nPops, controlTwDEMC=ctrl, ... )
	attr(res,"batchCall") <- cl
	boConverged=FALSE
	if( hasArg(fCheckConvergence))
		boConverged = (all(res$temp[nrow(res$temp),]<1.1)) & fCheckConvergence(res, fCheckConvergenceArgs)
	
	iRun <- iRun + nRun		#current number of runs
	nChains <- dim(res$parms)[3]
	nChainPop <- nChains %/% nPops
	chain2Pop <- rep(1:nPops, each=nChainPop )	#mapping of chain to population
	resCols <- match( rownames(res$logDenComp), rownames(res$Y))	#index of columns of results components in Y
	while( !boConverged & (iRun < nGen) ){
		cat(paste(iRun," out of ",nGen," generations completed. T=",paste({T<-res$temp[nrow(res$temp),];round(T,digits=ifelse(T>20,0,1))},collapse="  "),"     ",date(),"\n",sep=""))
		##details<<
		## Saves the result as variable \code{resRestart.twDEMC} to file \code{restartFilename}.
		if(is.character(restartFilename)){
			resRestart.twDEMC = res #avoid variable confusion on load by giving a longer name
			save(resRestart.twDEMC, file=restartFilename)
			cat(paste("Saved resRestart.twDEMC to ",restartFilename,"\n",sep=""))
		}
		# for populations with loo low acceptance rate repeat with extended nGenBurnin and recalculated temperature
		nRunPrev <- nRun
		if((iRun+nRun)<nGenBurnin) {
			##details<< \describe{\item{cooling and acceptance rate: repeating}{ 
			## If acceptance rate of some population drops below 0.1 then cooling was much too fast.
			## In this case, repeat the populations with low acceptance 
			## with extended nGenBurnin and recalculated temperature
			## }}
			pAcceptChains <- res$pAccept[ nrow(res$pAccept), ]
			pAcceptPops <- tapply( pAcceptChains, chain2Pop, mean) 
			#boPopCoolingTooFast <- (pAcceptPops < minAccepRateTempDecrease)
			boPopCoolingMuchTooFast <- (pAcceptPops < 0.1)
			#boPopCoolingMuchTooFast[2] <- TRUE	#for debugging
			i=0
			
			cl <-attr(res,"batchCall")
			relTempDecr <- (cl$TStart-res$temp[nrow(res$temp),])/cl$TStart
			while( doRepeatLowAcceptanceChains & any(boPopCoolingMuchTooFast) & any(relTempDecr[boPopCoolingMuchTooFast] > 0.05) & (nGenBurnin < maxNGenBurnin) & i<4){
				#XX TODO: repeat chains that cooled much too fast
				nGenBurnin <- min(maxNGenBurnin, nGenBurnin*1.5)	#extend the entire burnin
				cat(paste("repeating period for chains with too low acceptance rate: ",paste(round(pAcceptPops[boPopCoolingMuchTooFast],2),collapse=","),"; nGenBurnin=",nGenBurnin,"\n"))
				resF <- subChains(res,iPops=which(boPopCoolingMuchTooFast))
				iRunPrev <- iRun-nRunPrev
				resFt <- thin(resF, end=res$thin+iRunPrev)
				cl[["Zinit"]] <- resFt	
				cl$nPops <- ncol(resFt$temp)
				cl$nGen <- iRun
				cl$nGenBurnin <- nGenBurnin	#important to update for recalculating temperature
				#will be done in batch Callcl$controlTwDEMC$Tend <- calcDEMCTemp(resFt$temp[nrow(resFt$temp),,drop=FALSE], 1, nGenBurnin-iRunPrev, nRunPrev) #recalculate with initial temperature
				cl$doRepeatLowAcceptanceChains=FALSE	#else may be recursive
				#cl$initialAcceptanceRate <- popMeansTwDEMC( resFt$pAccept[nrow(resFt$pAccept),],nPops )
				cl$restartFilename<-NULL
				cl$X <- numeric(0)
				attr(resFt,"batchCall") <- cl
				resFtnew <- eval(cl)
				# matplot(resFtnew$pAccept,type="l")
				#mtrace(replacePops.twDEMC)
				res <- replacePops( res, resFtnew, iPops=which(boPopCoolingMuchTooFast))
				#recalculate acceptance 
				pAcceptChains <- res$pAccept[ nrow(res$pAccept), ]
				pAcceptPops <- tapply( pAcceptChains, chain2Pop, mean) 
				boPopCoolingMuchTooFast <- (pAcceptPops < 0.1)
			}
			if( any(boPopCoolingMuchTooFast & doRepeatLowAcceptanceChains &  (nGenBurnin < maxNGenBurnin)) ) stop("too low acceptance rate. Increasing nGenBurnin did not help.")
		}			
		
		zGen <- dim(res$parms)[2]
		if( (doResetOutlierN>0) & (iRun <= nGenBurnin) ){
			iGenOmega <- max(1,zGen-doResetOutlierN+1):zGen #(zGen%/%2):zGen
			# according to Vrugt09
			omega <- sapply( 1:nChains, function(iChain){mean(res$rLogDen[iGenOmega,iChain], na.rm=TRUE)}) #mean logDen across last halv of chain
			for( iPop in 1:nPops ){
				iChains <- ((iPop-1)*nChainPop+1):(iPop*nChainPop)
				q13 <- quantile( omega[iChains], c(1/4,3/4) )	#lower and upper quartile
				bo <- omega[iChains] < q13[1] -2*diff(q13)		#outside 2 interquartile ranges, see Vrugt09
				if( any(bo) ){
					#reset state of outliers to the best sampled parameter state
					tmp.best <- which( res$rLogDen[,iChains] == max(res$rLogDen[,iChains]), arr.ind = TRUE )[1,]	
					res$parms[,zGen,iChains[bo]] <- res$parms[,tmp.best[1],iChains[ tmp.best[2] ] ]
				}
			}
		}
		nRun <- min(nBatch, (if(iRun<nGenBurnin) min(nGenBurnin,nGen) else nGen) -iRun)		#iRun: Generation after batch run
		.dots <- list(...)
		.dots[c("logDenX","logDenCompX")] <- NULL;	#those will be inferred from res
		clArgs <- c(list(Zinit=res), .dots)	#Zinit must be first argument 
		clArgs$nGen<-nRun
		clArgs$nPops<-nPops
		#clArgs$TStart=max(1,b*exp(-a*iRun))		# if temp did not decrease start from this temperature
		clArgs$controlTwDEMC <- controlTwDEMC
		clArgs$controlTwDEMC$initialAcceptanceRate <- popMeansTwDEMC( res$pAccept[nrow(res$pAccept),],nPops )
		clArgs$controlTwDEMC$TStart<-TStartc<-res$temp[ nrow(res$temp), ,drop=FALSE]
		#clArgs$Tend=max(1,b*exp(-a*(iRun+nRun)))
		##--calculating end temperature
		clArgs$controlTwDEMC$Tend <- 1
		if((iRun+nRun)<nGenBurnin) {
			##details<< \describe{\item{cooling and acceptance rate}{ 
			## If acceptance rate of some population drops below rate=minAccepRateTempDecrease then cooling is too fast.
			## In this moderate case do not repeat the rund but keep the current temperature for the next period for this population.
			## and extend the burnin phase by the length of this period.
			## }}
			pAcceptChains <- res$pAccept[ nrow(res$pAccept), ]
			pAcceptPops <- tapply( pAcceptChains, chain2Pop, mean) 
			boPopCoolingTooFast <- (pAcceptPops < minAccepRateTempDecrease) 
			
			##details<< \describe{\item{cooling and expected difference in LogDensity}{ 
			## If the difference of temperated LogDensitys between Proposed stepds and accepted steps
			## of the component with highest difference (which is negative)
			## drops below rate=minPCompAcceptTempDecr then cooling is too fast.
			## The median of the last 128 steps is used
			## In this case keep the current temperature for the next period for this population.
			## and extend the burnin phase by the length of this period.
			## }}
			#construct Temp for results and populations
			tempResPops <- matrix( rep(TStartc, length(resCols)), ncol=length(TStartc), byrow=TRUE, dimnames=list(rownames(res$Y)[resCols],NULL))
			tempResPops[names(ctrl$TFix),] <- matrix( rep(TFix, length(TStartc)), ncol=length(TStartc) )
			#mtrace(getDiffLogDen.twDEMCPops)
			#diffLogDenT <- getDiffLogDen.twDEMCPops(res$Y, resCols, temp=tempResPops, nLastSteps=ceiling(128/nChainPop)) 	#in twDEMC S3twDEMC.R
			diffLogDen <- getDiffLogDen.twDEMCPops(res$Y, resCols, nLastSteps=ceiling(128/nChainPop)) 	#in twDEMC S3twDEMC.R
			diffLogDenPops <- popApplyTwDEMC( diffLogDen, nPops=nPops, function(x){ abind(twListArrDim(x),along=2,new.names=dimnames(x)) })	#stack param columns by population
			#diffLogDenPops[!is.finite(diffLogDenPops)] <- NA
			#XXTODO: think about criterion for too fast cooling
			#tmp <- ecdf(diffLogDenPops["amendm",,])
			# calculate temperated diffLogDen, i.e. divided by the component and population specific temperature
			.tmp.f <- function(){ 
				diffLogDenPopsT <- .calcTemperatedDiffLogDen( diffLogDenPops, TFix, TStartc)
				#acceptance rate for each single parameter (percentil > log(0.5)
				tmpPercAcc <- 1-apply(diffLogDenPopsT,c(1,3),function(d){ ecdf(d)(log(0.5))} ) #comp x pops
				boPopCoolingTooFast <- (apply(tmpPercAcc,2,min) < minPCompAcceptTempDecr)
			}
			
			if( any(boPopCoolingTooFast) )
				nGenBurnin=min(maxNGenBurnin, nGenBurnin+nRun)
			temp <- TStartc	# keep Temperature, only for others cool down further
			if( any(!boPopCoolingTooFast)) temp[!boPopCoolingTooFast] <- {
					tempExp <- calcDEMCTemp( TStartc[!boPopCoolingTooFast], 1, nGenBurnin-iRun, nRun) #recalculate with initial temperature
					#mtrace(calcDEMCTempDiffLogDenConst)
					#tempEmp <- sapply( seq_along(TStartc)[!boPopCoolingTooFast], function(iPop){ calcDEMCTempDiffLogDenConst(diffLogDenPops[,,iPop , drop=FALSE], TFix=clArgs$controlTwDEMC$TFix, Tmax=TStartc[iPop], pTarget=minPCompAcceptTempDecr+0.02)})
					#pmax(1, pmin( tempExp,(2*tempExp+tempEmp)/3 ))
					tempExp
				}
			clArgs$controlTwDEMC$Tend <- temp
			
			pTarget=minPCompAcceptTempDecr+0.02
			if( (0<length(clArgs$controlTwDEMC$useMultiT)) ) if( clArgs$controlTwDEMC$useMultiT ){
					#mtrace(calcDEMCTempDiffLogDen2)
					clArgs$controlTwDEMC$Tprop <- TDiffLogDen <- sapply( seq_along(TStartc), function(i){calcDEMCTempDiffLogDen2(diffLogDenPops[,,i], pTarget=pTarget, TFix=TFix, Tmax=TStartc[i])})  # will be scaled in twDEMCBlockInt
					#clArgs$controlTwDEMC$Tprop <- TDiffLogDen <- apply( diffLogDenPops, 3, calcDEMCTempDiffLogDen, pTarget=pTarget)  # will be scaled in twDEMCBlockInt
					#clArgs$controlTwDEMC$Tend <- pmin( tempExp, apply(TDiffLogDen,2,max))
				}
			
		}
		clArgs$controlTwDEMC$probUpDir <- (if((iRun+nRun)<=nGenBurnin) probUpDirBurnin else NULL)	#set to NULL after burnin
		res <- do.call( twDEMC, clArgs, quote=TRUE )
		attr(res,"batchCall") <- cl
		#res <- twDEMC( Zinit=res, nGen=nRun, ... ) problems with double Zinit 
		if( hasArg(fCheckConvergence))
			boConverged = (all(res$temp[nrow(res$temp),]<1.1)) & fCheckConvergence(res, fCheckConvergenceArgs)
		iRun = iRun + nRun
	}
	cat(paste(iRun," out of ",nGen," generations completed. T=",paste({T<-res$temp[nrow(res$temp),];round(T,digits=ifelse(T>20,0,1))},collapse="  "),"     ",date(),"\n",sep=""))
	res$nGenBurnin <- nGenBurnin
	res
	### List of class twDEMC (see \code{\link{twDEMCBlockInt}}) with additional entry nGenBurnin
}

