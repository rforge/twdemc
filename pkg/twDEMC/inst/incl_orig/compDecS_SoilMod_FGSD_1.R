.tmp.f <- function(){
	#sfLibrary(mnormt)
	#sfLibrary(MASS) #better random number generator
	#sfStop()
	#source("../mcmc/tw.DEMCzs/R/sensitivity.R")  #  corOrdered
}



.tmp.f <- function(){
	library(R.utils)
	source(file.path("inst","cluster","setupCluster.R"))
	data(HamerLigninGlucoseRespUnc)
	data(HamerParameterPriors)

	#state 1.7:
	#FGSD1: long run with only 6 population by 8 chains
	#FGSD2 (iproc21): short run with 16 populationas by 4 chains to 3*1024 
	paramFilenameAsom <- file.path("parms","MDIHamer_FGSD2.RData")
	paramFilenameAsom <- file.path("parms","MDIHamer_FGSD1.RData")
	paramFilename <- file.path("..",paramFilenameAsom)
	bsubOptions="-q SLES10"
	
	poptNamesFGSD <- c("A0", "pA01toA0", "F0", "kF", "mFF", "epsF", "sA", "epsA"
			,"tvr", "tvrP", "tvrA"
			,"kG", "mGG", "decGPotAStar", "epsG", "betaA", "betaP", "pStorageOnActivation"
			,"D0"
			,"biasRespSum","biasRespC14Sum"
		)
	
	parms <- within( HamerParameterPriors$parms0, {
			fSDec="schimel"
			f.sDec=decompFList[[fSDec]]
		})
	poptNames <- c( poptNamesFGSD, decompParList[[parms$fSDec]])	#schimel and bias in observations
	
	poptDistr <- twConstrainPoptDistr(poptNames,HamerParameterPriors$parDistr )
	normpopt <- poptDistr$mu
	transOrigPopt(normpopt, poptDistr$trans)
	
	.tmp.f <- function(){
		mtrace(initStateHamer.SoilMod_FGSD)
		mtrace(initState.SoilMod_FGSD)
		mtrace(initState.SoilMod)
	}
	
	model <- list(
		modMeta = modMeta.SoilMod_FGSD()
		,fInitState = initStateHamer.SoilMod_FGSD
		,fSolve = solve.SoilMod_FGSD
	) 
	
	argsFLogLik <- argsFLogLik0 <- list(
		model=model		 
		,poptDistr=poptDistr
		,obs=HamerLigninGlucoseRespUnc
		,parms=HamerParameterPriors$parms0
		,fCalcBiasedObs=calcBiasedObs.hamer	#include creation of biased observations
		,argsFCalcBiasedObs=list(namesDs=calcBiasedObs.hamer.namesDs)	#not sourced on remote hosts
		#,includeStreams=c("parms",includeStreams.hamer.resp)
		,includeStreams=c("parms",includeStreams.hamer.resp,includeStreams.hamer.accresp)
	#,useRk4=TRUE		#need to specify here, becuase also used in replaceNonFiniteZinit
	#,useRImpl=.useRImpl
	)
	#mtrace(ofb.hamer)
	nObs <-max( sapply( argsFLogLik$obs, function(dstr) length(dstr$obs) ) ) #expected misfit
	Tstd <- nObs / length(normpopt)		#Temperature for which totally free model of n parameters expected to yield Log-Likelihood L=-1
		
	# tried to increase number of chains per population because acceptance rate dropped to low, but this did not work
	.nPops=16	
	#.nPops=6	
	#.nPops=3
	.internalComponents<- c("parms","amendm","c14obs","amendm_small","c14obs_small","amendmSum","c14obsSum","amendmSum_small","c14obsSum_small")
	#Zinit <- initZtwDEMCNormal( normpopt, poptDistr$sigma, nChains=4*.nPops, nPops=.nPops, doIncludePrior=FALSE)
	Zinit <- initZtwDEMCNormal( normpopt, poptDistr$sigma, nChains=4*.nPops, nPops=.nPops, doIncludePrior=FALSE)
	
	#------- calculate the LogLik of the proposals 
	str(tmp <- t(abind(twListArrDim(Zinit),along=2,new.names=dimnames(Zinit)[1:2])))
	argsTwDEMCBatchInit <- within( list(),{
			fLogLik<-ofb.hamer
			argsFLogLik<-argsFLogLik
			xProp<-tmp
		})
	runClusterParms <- runClusterParmsInit <- within( list(),{
		fSetupCluster <- setupClusterMDIHamerDev
		fRun <- twRunFLogLikPar
		#specify intResCompNames here, so that intResCompNames from previous results is not overwritten
		argsFRun <- argsTwDEMCBatchInit
		#, intResCompNames=intersect(.internalComponents, argsFLogLik$includeStreams)
	})
	.tmp.f <- function(){
		argsFRun2 <- runClusterParmsInit$argsFRun
		argsFRun2$xProp <- argsFRun2$xProp[1:500,]
		tmp <- do.call(twRunFLogLikPar, argsFRun2)
		resLogLik <- tmp
	}
	save( runClusterParms, file=paramFilename )
	iproc=0;outFilename0 <- twResultFilename(paramFilename,iproc=iproc); outFilenameAsom1 <- twResultFilename(paramFilenameAsom,iproc=iproc)
	unlink( outFilename0 )
	nProc=min(7,dim(Zinit)[3])
	nProc=min(14,dim(Zinit)[3])
	nProc=min(16,dim(Zinit)[3])
	#	nP = min(8,dim(demcDrivers$Zinit[3])); copy2clip(as.character(GString("bsub -n ${nP} ./bsubr_i.sh ${`getHamerDemcFilename('FGSD')`} iproc=41 nprocSingle=${nP} 'paramFile=\"parms/demcDrivers_1tFGSD.RData\"'")))
	#copy2clip(as.character(GString("~/bsubr_iLocal.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\\\"${paramFilenameAsom}\\\"' 'argsFRun=list(nGen=16)'" )))	#testing with 16 generations on pc026
	#copy2clip(as.character(GString("~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=16)'" )))	#testing with 16 generations on pc026
	copy2clip(as.character(GString("./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=4 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=16)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=16)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"'" )))
	
	load(outFilename0)
	#load(file.path("..","parms","res_res_MDIHamer_FGSD1_1_1.RData"))
	resLogLik <- runClusterRes$res
	
	boFin <- is.finite(resLogLik$logLik)
	rL <- resLogLik$logLik[ boFin ]
	p <- c(0.05, 0.1, 0.2, 100/length(rL) )
	q <- quantile(rL, probs=1-p)
	bo <- sapply(q, function(qi) rL > qi) #is.finite(resLogLik$logLik); 
	apply(bo,2,sum)	#145 cases is ok for 5plevel is ok
	i=4
	resLogLikQ <- resLogLik
	resLogLikQ$logLik <- resLogLik$logLik[ boFin ][bo[,i]]
	resLogLikQ$resFLogLik <- resLogLik$resFLogLik[ boFin, ][bo[,i],]
	tmp <- melt(resLogLikQ$resFLogLik)
	p1 <- qplot( X2, value, geom="boxplot", data=tmp)+
		opts(axis.text.x=theme_text(angle=30, hjust=1, size=8))
	p1

	#apply( resLogLikQ$resFLogLik, 2, range )
	#sample from the 5% best as accepted LogLik
	Lp <- resLogLikQ$resFLogLik
	ord <- order(resLogLikQ$logLik, decreasing = TRUE)
	#resLogLikQ$logLik[ord]
	La <- Lp[sample( ord[1:round(length(ord)*0.05)], nrow(Lp), replace=TRUE ), ]
	
	d <- Lp-La
	tmp <- melt(d)
	p2 <- qplot( X2, value, geom="boxplot", data=tmp)+
		opts(axis.text.x=theme_text(angle=30, hjust=1, size=8))
	#p2
	

	q20 <- apply( d, 2, quantile, probs=1-c(0.2))
	T <- pmax(1,q20/log(0.2))
	
	ps <- 0.8
	pps <- function(ps,pTarget=0.2){
		qps <- apply( d, 2, quantile, probs=1-ps)
		sum(apply( d, 1, function(dj){ all(dj>=qps)} ))/nrow(d) - pTarget
	}
	#pps(0.82,0.2)
	ps <- uniroot( pps, c(0.2, 1), tol=0.01 )$root
	qps <- apply( d, 2, quantile, probs=1-ps)
	T <- structure(pmax(1,qps/log(ps), names=names(qps)))
	
	ZinitR <- replaceZinitNonFiniteLogLiks( Zinit, resLogLik$logLik)
	
	res$Y[names(T),1:10,]
	
	
	
	#-------  first 3*1024 with high temperature and multistep Metropolis
	#experiment with T0, so that chains are not parallel wit 4 pops, when found repeat with 8
	minNBatch <- ceiling(calcM0twDEMC( length(normpopt),.nPops,ncol(Zinit))/0.06)
	argsTwDEMCBatch0 <- within( list(),{
			#Zinit<-Zinit
			Zinit<-ZinitR
			nPops<-.nPops
			fLogLik<-ofb.hamer
			argsFLogLik<-argsFLogLik
			#loosing variability with > 20 parameters nBatch<-512
			#nBatch=max(minNBatch,1024)
			#nGen<-3*1024
			nGen <- 512
			nBatch<-nGen
			controlTwDEMC<-list(thin=8,Tprop=T)
			#T0<-nObs*16
			T0<-max(T)
			nGenBurnin<-nGen*1000	#essentially not decreasing temp
			doRecordProposals=TRUE
	})
	runClusterParms <- runClusterParms0 <- within( list(),{
			fSetupCluster <- setupClusterMDIHamerDev
			fRun <- twRunDEMC
			#specify intResCompNames here, so that intResCompNames from previous results is not overwritten
			argsFRun <- list(argsTwDEMCBatch=argsTwDEMCBatch0
				#, intResCompNames=intersect(.internalComponents, argsFLogLik$includeStreams)
			)
		})
	#runClusterParms$argsFRun$fLogLik
	.tmp.f <- function(){
		runClusterParms$argsFRun$argsTwDEMCBatch$Zinit <- ZinitR[,,1:8]
		runClusterParms$argsFRun$argsTwDEMCBatch$nPops <- 2
		
		.prevDir <- setwd("..")
		#do.call( fSetupCluster, list()  )
		argsFRunDebug <- runClusterParms$argsFRun 
		argsFRunDebug$argsTwDEMCBatch <- within(argsFRunDebug$argsTwDEMCBatch, {
				Zinit=Zinit[,,1:8,drop=FALSE]
				nPops=2
				nGen=16
				controlTwDEMC$thin=4
				debugSequential=TRUE
				#nBatch=8
				#doStopOnError=TRUE
				#argsFLogLik$useRImpl=TRUE
			})
		#options(error=dump.frames)
		#mtrace(twRunDEMC)
		#mtrace(replaceZinitNonFiniteLogLiksLastStep)
		#mtrace(twCalcLogLikPar)
		#tmp<-argsFRunDebug$argsTwDEMCBatch$argsFLogLik$model$fSolve; mtrace(tmp); argsFRunDebug$argsTwDEMCBatch$argsFLogLik$model$fSolve<-tmp
		#tmp<-ofb.hamer; mtrace(tmp); argsFRunDebug$argsTwDEMCBatch$fLogLik<-tmp
		#mtrace(deriv.SoilMod_FGSD)
		#mtrace(twDEMCInt)
		#mtrace(.doDEMCStep)
		tmp <- NULL; tmp <-do.call( twRunDEMC, argsFRunDebug )
		tmp$rLogLik[ nrow(tmp$rLogLik), ]
		setwd(.prevDir)
	}
	
	iproc=21;outFilename1 <- twResultFilename(paramFilename,iproc=iproc); outFilenameAsom1 <- twResultFilename(paramFilenameAsom,iproc=iproc)
	iproc=1;outFilename1 <- twResultFilename(paramFilename,iproc=iproc); outFilenameAsom1 <- twResultFilename(paramFilenameAsom,iproc=iproc)
		
	save( runClusterParms, file=paramFilename )
	unlink( outFilename1 )
	nmax <- dim(runClusterParms$argsFRun$argsTwDEMCBatch$Zinit)[3]
	nProc=min(7,nmax)
	nProc=min(14,nmax)
	nProc=min(16,nmax)
	#	nP = min(8,dim(demcDrivers$Zinit[3])); copy2clip(as.character(GString("bsub -n ${nP} ./bsubr_i.sh ${`getHamerDemcFilename('FGSD')`} iproc=41 nprocSingle=${nP} 'paramFile=\"parms/demcDrivers_1tFGSD.RData\"'")))
	#copy2clip(as.character(GString("~/bsubr_iLocal.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\\\"${paramFilenameAsom}\\\"' 'argsFRun=list(nGen=16)'" )))	#testing with 16 generations on pc026
	#copy2clip(as.character(GString("~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=16)'" )))	#testing with 16 generations on pc026
	copy2clip(as.character(GString("./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=4 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=16)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=16)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"'" )))

	
	load(outFilename1)
	#load(file.path("..","parms","res_res_MDIHamer_FGSD1_1_1.RData"))
	res <- runClusterRes$res
	nRun1 <- (nrow(res$rLogLik)-1)*res$thin
	#to extend the run, copy the results file to param file (res will be used) and change nGen
	#file.copy(outFilename1,paramFilename,overwrite=TRUE);	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilename}\"' 'argsFRun=list(nGen=${`nRun1+1*1024`})'" )))
	
	.start=16
	#.start=145
	#windows(record=TRUE)
	matplot( res$temp[-(1:.start),], type="l")
	matplot( popMeansTwDEMC(res$pAccept[-(1:.start),], ncol(res$temp), 10), type="l" )
	abline(h=0.1,col="gray")
	matplot( res$rLogLik[-(1:.start),], type="l")
	tmp <- res$rLogLik[nrow(res$rLogLik),]
	
	rest <- thin(res, start=nRun1-512)
	matplot( rest$rLogLik[,], type="l")
	matplot( popMeansTwDEMC(rest$pAccept, ncol(rest$temp), 5), type="l" )
	abline(h=0.1,col="gray")
	#rescoda <- as.mcmc.list(rest)
	#windows(record=TRUE)
	#plot(thinN(rescoda), smooth=FALSE)
	#rest <- thin(res, start=(1024+512))
	p1 <- ggplotDensity.twDEMC( rest, poptDistr, pMin=0.1)
	#p1
	#rescoda <- as.mcmc.list(rest)
	#windows(record=TRUE)
	#plot(thinN(rescoda), smooth=FALSE)

	.tmp.f <- function(){
		#-------  second remove outlier chains
		# now do not use multistep Metropolis any more
		# tried with going down in 4*1024 steps with holding temp if acceptance rate drops to low 
		
		iproc=1;outFilename1 <- twResultFilename(paramFilename,iproc=iproc)
		#iproc=21;outFilename1 <- twResultFilename(paramFilename,iproc=iproc)
		load(outFilename1)
		#load(file.path("..","parms","res_res_MDIHamer_FGSD1_1_1.RData"))
		res <- runClusterRes$res
		#res0 <- runClusterRes$argsFRun$argsTwDEMCBatch$Zinit
		runClusterRes$res <- runClusterRes$argsFRun$prevResRunCluster <- NULL	#else interferes with continued run, also save space
		
		#nRun0 <- (nrow(res0$rLogLik)-1)*res0$thin
		nRun1 <- (nrow(res$rLogLik)-1)*res$thin
		runClusterParms <- runClusterParms0
		runClusterParms$argsFRun$argsTwDEMCBatch <- within( runClusterRes$argsFRun$argsTwDEMCBatch,{
				Zinit<-res
				nBatch<-256	#	check on outlier chains quite often 
				nGen<-nRun1+2*1024	
				nGenBurnin<-nGen*1000	
				nPops<-ncol(res$temp)
			})
		runClusterParms$argsFRun$resFLogLikX<-character(0)
		runClusterParms$fRun <- twRunDEMC		# update functions to represent current dev state 
		runClusterParms$argsFRun$fLogLik <- ofb.hamer
		iproc=2; outFilename2 <- twResultFilename(paramFilename,iproc=iproc)
		#iproc=22; outFilename2 <- twResultFilename(paramFilename,iproc=iproc)
		save( runClusterParms, file=paramFilename )
		unlink( outFilename2 )
		
		#nProc=min(15,dim(Zinit)[3])
		nProc=min(16,dim(runClusterParms$argsFRun$Zinit)[3])
		copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${`nRun1+16`},doRecordProposals=TRUE)'" )))
		copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${nRun1}+1024*1)'" )))
		copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"'" )))
		
		load(outFilename2)
		res <- runClusterRes$res
		res0 <- runClusterRes$argsFRun$argsTwDEMCBatch$Zinit
		runClusterRes$res <- runClusterRes$argsFRun$prevResRunCluster <- NULL	#else interferes with continued run, also save space
		nRun0 <- (nrow(res0$rLogLik)-1)*res0$thin
		nRun1 <- (nrow(res$rLogLik)-1)*res$thin
		
		.start=nRun0; .starti <- .start%/% res$thin
		#tmp<-apply(res$temp,1,max); .starti= min(which(tmp<50)); .start <- .starti*res$thin
		#windows(record=TRUE)
		matplot( res$temp[-(1:.starti),], type="l")
		matplot( popMeansTwDEMC(res$pAccept[-(1:.starti),], ncol(res$temp), 2), type="l" )
		abline(h=0.1,col="gray")
		matplot( res$rLogLik[-(1:.starti),], type="l")
		
		
		copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${`nRun1+5*1024`},nGenBurnin=${`res$nGenBurnin`})' doRestart=TRUE" )))
		
		rest0 <- thin(res, start=nRun1-1024)
		p0 <- ggplotDensity.twDEMC( rest0, poptDistr)
	}

	#-------  second decrease temperature quite fast and observe retardation by acceptance rate
	# now do not use multistep Metropolis any more
	# tried with going down in 4*1024 steps with holding temp if acceptance rate drops to low 
	
	iproc=1;outFilename1 <- twResultFilename(paramFilename,iproc=iproc)
	#iproc=21;outFilename1 <- twResultFilename(paramFilename,iproc=iproc)
	load(outFilename1)
	#load(file.path("..","parms","res_res_MDIHamer_FGSD1_1_1.RData"))
	res <- runClusterRes$res
	res0 <- runClusterRes$argsFRun$argsTwDEMCBatch$Zinit
	runClusterRes$res <- runClusterRes$argsFRun$prevResRunCluster <- NULL	#else interferes with continued run, also save space
	
	nRun0 <- (nrow(res0$rLogLik)-1)*res0$thin
	nRun1 <- (nrow(res$rLogLik)-1)*res$thin
	runClusterParms <- runClusterParms0
	runClusterParms$argsFRun$argsTwDEMCBatch <- within( runClusterRes$argsFRun$argsTwDEMCBatch,{
			Zinit<-res
			nBatch<-512	#	want to see output of progress and remove outlier chains 
			nGen<-nRun1+32*1024	
			nGenBurnin<-nRun1+16*1024	
			nPops<-ncol(res$temp)
	})
	runClusterParms$argsFRun$resFLogLikX<-character(0)
	runClusterParms$fRun <- twRunDEMC		# update functions to represent current dev state 
	runClusterParms$argsFRun$fLogLik <- ofb.hamer
	iproc=2; outFilename2 <- twResultFilename(paramFilename,iproc=iproc)
	#iproc=22; outFilename2 <- twResultFilename(paramFilename,iproc=iproc)
	save( runClusterParms, file=paramFilename )
	unlink( outFilename2 )

	#nProc=min(15,dim(Zinit)[3])
	nProc=min(16,dim(runClusterParms$argsFRun$Zinit)[3])
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${`nRun1+16`},doRecordProposals=TRUE)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${nRun1}+1024*1)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"'" )))

	load(outFilename2)
	res <- runClusterRes$res
	res0 <- runClusterRes$argsFRun$argsTwDEMCBatch$Zinit
	runClusterRes$res <- runClusterRes$argsFRun$prevResRunCluster <- NULL	#else interferes with continued run, also save space
	nRun0 <- (nrow(res0$rLogLik)-1)*res0$thin
	nRun1 <- (nrow(res$rLogLik)-1)*res$thin
	
	.start=nRun0; .starti <- .start%/% res$thin
	#tmp<-apply(res$temp,1,max); .starti= min(which(tmp<50)); .start <- .starti*res$thin
	#windows(record=TRUE)
	matplot( res$temp[-(1:.starti),], type="l")
	matplot( popMeansTwDEMC(res$pAccept[-(1:.starti),], ncol(res$temp), 2), type="l" )
	abline(h=0.1,col="gray")
	matplot( res$rLogLik[-(1:.starti),], type="l")
	
	
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${`nRun1+5*1024`},nGenBurnin=${`res$nGenBurnin`})' doRestart=TRUE" )))
	
	rest0 <- thin(res, start=nRun1-512)
	p0 <- ggplotDensity.twDEMC( rest0, poptDistr)
	#p0
	tmp <- plotChainPopMoves(rest0)
	rest <- subChains( rest0, iPops=(1:8)[-c(3,7)])
	rest <- subChains( rest0, iPops=(1:8)[-c(1,2,3,7)])
	tmp <- plotChainPopMoves(rest)
	rescoda <- as.mcmc.list(rest)
	#windows(record=TRUE)
	plot(thinN(rescoda), smooth=FALSE)
	p1 <- ggplotDensity.twDEMC( rest, poptDistr)
	#p1
	
	
	#-------- third decrease Temp further, speed based on observed temp decrease
	iproc=2; outFilename2 <- twResultFilename(paramFilename,iproc=iproc)
	load(outFilename2)
	#load(outFilename3)
	resOrig <- runClusterRes$res
	res <- subChains( resOrig, iPops=(1:8)[-c(3,7)])		#will cancel attribute batchCall
	res0 <- runClusterRes$argsFRun$argsTwDEMCBatch$Zinit
	runClusterRes$res <- runClusterRes$argsFRun$prevResRunCluster <- NULL	#else interferes with continued run, also save space
	nRun1 <- (nrow(res$rLogLik)-1)*res$thin
	nRun0 <- (nrow(res0$rLogLik)-1)*res0$thin

	iStep <- nRun1-nRun0	#nRun0 from outFilename1 above
	.burnin <- max( res$nGenBurnin,ceiling(nRun0+calcDEMCnGenBurnin(mean(res$temp[nRun0/res$thin,]), mean(res$temp[nRun1/res$thin,]), iStep)))
		
	runClusterParms <- runClusterRes
	runClusterParms$argsFRun$argsTwDEMCBatch <- within( runClusterParms$argsFRun$argsTwDEMCBatch,{
			Zinit<-res
			nBatch<-512	#	want to see output of progress and remove outlier chains 
			nGen<-.burnin+5*1024	
			nGenBurnin<-.burnin
			nPops <- ncol(res$temp)
		})
	runClusterParms$argsFRun$resFLogLikX<-character(0)
	runClusterParms$fRun <- twRunDEMC		# update functions to represent current dev state 
	#runClusterParms$argsFRun$fLogLik <- ofb.hamer
	iproc=3; outFilename3 <- twResultFilename(paramFilename,iproc=iproc)
	save( runClusterParms, file=paramFilename )
	unlink( outFilename3 )
	
	#nProc=min(15,dim(Zinit)[3])
	nProc=min(16,dim(runClusterParms$argsFRun$Zinit)[3])
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${`nRun1+16`},doRecordProposals=TRUE)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${nRun1}+1024*1)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"'" )))
	
	load(outFilename3)
	res <- runClusterRes$res
	res0 <- runClusterRes$argsFRun$argsTwDEMCBatch$Zinit
	runClusterRes$res <- runClusterRes$argsFRun$prevResRunCluster <- NULL	#else interferes with continued run, also save space
	nRun0 <- (nrow(res0$rLogLik)-1)*res0$thin
	nRun1 <- (nrow(res$rLogLik)-1)*res$thin
	
	.start=nRun0; .starti <- .start%/% res$thin
	#windows(record=TRUE)
	matplot( res$temp[-(1:.starti),], type="l")
	matplot( popMeansTwDEMC(res$pAccept[-(1:.starti),], ncol(res$temp), 10), type="l" )
	matplot( res$rLogLik[-(1:.starti),], type="l")
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ./bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${`nRun1+5*1024`},nGenBurnin=${`res$nGenBurnin`})' doRestart=TRUE" )))
	
	rest0 <- thin(res, start=nRun1-512)
	p0 <- ggplotDensity.twDEMC( rest0, poptDistr)
	#p0
	tmp <- plotChainPopMoves(rest0)
	rest <- subChains( rest0, iPops=(1:8)[-c(1,2,6)])
	rest <- subChains( rest0, iPops=(1:8)[c(1,2,6)])
	tmp <- plotChainPopMoves(rest)
	rescoda <- as.mcmc.list(rest)
	#windows(record=TRUE)
	plot(thinN(rescoda), smooth=FALSE)
	p1 <- ggplotDensity.twDEMC( rest, poptDistr)
	#p1
	p2 <- ggplotDensity.twDEMC( rest, poptDistr, doTransOrig=TRUE)
	#p2

	normpoptBest2 <- twExtractFromLastDims(rest$parms, which.max( rest$rLogLik) )[,1]
	tmp <- attr(res,"batchCall")	
	#tmp2 <- hamerPlot3.of( normpoptBest2, tmp$fLogLik, argsFLogLik=tmp$argsFLogLik)
	tmp2 <- hamerPlot3.of( normpoptBest2, ofb.hamer, argsFLogLik=tmp$argsFLogLik)
	barplot(unlist(as.list(tmp2)))
	barplot(attributes(tmp2)$logLikParms)
	
	#recovery
	obsB <- attributes(tmp2)$obs
	parmsAdj <- attributes(tmp2)$parms
	unlist(parmsAdj[c("biasRespSum","biasRespC14Sum")])
	obsB$c14obsSum$obs / parmsAdj$Fa0 
	obsB$c14obsSum_small$obs / parmsAdj$Fa0s 
	

	
	normp <- twExtractFromLastDims( rest$parms, which.max(rest$rLogLik) )[,1] 
	transOrigPopt(normpopt, poptDistr$trans)
	(popt <- transOrigPopt(normp, poptDistr$trans))
	
	#plot best results
	#resOf <- hamerPlot3.of( normp ) 
	tmp <- attr(res,"batchCall")	
	resOf <- hamerPlot3.of( normp, argsFLogLik=tmp$argsFLogLik) 
	unlist(as.list(resOf))
	sort(attr(resOf,"logLikParms"), decreasing=TRUE)

	#mtrace(hamerPlotResCols.resOf)
	hamerPlotResCols.resOf( resOf, paste("csums",model$modMeta$rowNames,sep="_"), "control" )
	hamerPlotResCols.resOf( resOf, paste("csums",c("F","A"),sep="_"), "control" )
	
	
}










