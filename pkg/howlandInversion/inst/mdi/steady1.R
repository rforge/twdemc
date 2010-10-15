.tmp.f <- function(){
	#sfLibrary(mnormt)
	#sfLibrary(MASS) #better random number generator
	#sfStop()
	#source("../mcmc/tw.DEMCzs/R/sensitivity.R")  #  corOrdered
}



.tmp.f <- function(){
	library(R.utils)
	source(file.path("inst","cluster","setupCluster.R"))
	data(HamerLigninGlucoseRespAggrGrowth)
	data(HamerParameterPriors)

	#state 1.7:
	#FGSD1: long run with only 6 population by 8 chains
	#FGSD2 (iproc21): short run with 16 populationas by 4 chains to 3*1024 
	paramFilenameAsom <- file.path("parms","MDIHamer_FGSD2.RData")
	paramFilenameAsomB <- file.path("parms","MDIHamer_FGSD2_B.RData")
	paramFilename <- file.path("..",paramFilenameAsom)
	paramFilenameB <- file.path("..",paramFilenameAsomB)
	#bsubOptions="-q SLES10"
	bsubOptions="-q mpi_large"	# will get several nodes instead one with 16 procs
	bsubOptions="-q BTM"		
	copy2clip("bqueues -l mpi_large | grep HOSTS | cut -d ' ' -f 3- | xargs bhosts")
	
	poptNamesFGSD <- c("A0", "pA01toA0", "kF", "mFF", "epsF", "sA"
			,"tvr", "tvrP", "tvrA"
			#, "F0"
			#,"kG", "mGG", "decGPotAStar", "epsG", "betaA", "betaP", "pStorageOnActivation", "epsA"
			,"kG", "mGG", "decGPotAStar", "epsG"
			#,"D0"
			,"biasRespSum","biasRespC14Sum"
		)
	
	parms <- within( HamerParameterPriors$parms0, {
			fSDec="schimel"
			f.sDec=decompFList[[fSDec]]
			betaA=0; D0=0			# no activation/deactivation
		})
	poptNames <- c( poptNamesFGSD, decompParList[[parms$fSDec]])	#schimel and bias in observations
	
	poptDistr <- twConstrainPoptDistr(poptNames,HamerParameterPriors$parDistr )
	normpopt <- poptDistr$mu
	transOrigPopt(normpopt, poptDistr$trans)
	
	TFix <- c(parms=1)	#do not use increased temperate for parms data stream
	
	model <- list(
		modMeta = modMeta.SoilMod_FGSD()
		,fInitState = initStateHamer.SoilMod_FGSD
		,fSolve = solve.SoilMod_FGSD
	) 
	
	.includeStreams <- c("parms",includeStreams.hamer.resp,includeStreams.hamer.accresp)
	.obs <- HamerLigninGlucoseRespAggrGrowth
	#.rescompList <- getOfbRescompList.hamer(.includeStreams,.obs) 
	#.internalComponents<- c("parms","amendm","c14obs","amendm_small","c14obs_small","amendmSum","c14obsSum","amendmSum_small","c14obsSum_small")
	.internalComponents<- c("parms")
	argsFLogLik <- argsFLogLik0 <- list(
		model=model		 
		,poptDistr=poptDistr
		,obs=.obs
		,parms=HamerParameterPriors$parms0
		,fCalcBiasedObs=calcBiasedObs.hamer	#include creation of biased observations
		,argsFCalcBiasedObs=list(namesDs=calcBiasedObs.hamer.namesDs)	#not sourced on remote hosts
		,includeStreams=.includeStreams
		#,useRk4=TRUE		#need to specify here, becuase also used in replaceNonFiniteZinit
		#,useRImpl=.useRImpl
		#,rescompList=.rescompList
		#,logLikAcceptPos=sapply(names(.rescompList),function(resComp){match(.rescompList[[resComp]], .internalComponents )})
		#,ofb.subF=ofb.subF.hamer
	)
	
	#mtrace(ofb.hamer)
	nObsSum <-sum( sapply( argsFLogLik$obs, function(dstr) length(dstr$obs) ) ) #expected misfit
	nObs <-max( sapply( argsFLogLik$obs, function(dstr) length(dstr$obs) ) ) #expected misfit
	Tstd <- nObs / length(normpopt)		#Temperature for which totally free model of n parameters expected to yield Log-Likelihood L=-1
		
	# tried to increase number of chains per population because acceptance rate dropped to low, but this did not work
	.nPops=16	
	#.nPops=6	
	#.nPops=3
	Zinit <- initZtwDEMCNormal( normpopt, poptDistr$sigma, nChains=4*.nPops, nPops=.nPops, doIncludePrior=FALSE)
	# try to not spread parameters too wide so that parms datastream takes not half of the proposals away in the beginning
	# does now work: parameters spread during steps and temperature has been calculated too low and acceptance rate goes down
	#Zinit <- initZtwDEMCNormal( normpopt, poptDistr$sigma/2, nChains=4*.nPops, nPops=.nPops, doIncludePrior=FALSE)
	
	#------- calculate the LogLik of the proposals 
	str(.xProp <- t(abind(twListArrDim(Zinit),along=2,new.names=dimnames(Zinit)[1:2])))
	#mtrace(ofb.hamer)
	#
	argsTwDEMCBatchInit <- within( list(),{
			fLogLik<-ofb.hamer
			argsFLogLik<-argsFLogLik
			xProp<-.xProp
			debugSequential<-TRUE
		})
	runClusterParms <- runClusterParmsInit <- within( list(),{
		fSetupCluster <- setupClusterMDIHamerDev
		fRun <- twRunFLogLikPar
		#specify intResCompNames here, so that intResCompNames from previous results is not overwritten
		argsFRun <- argsTwDEMCBatchInit
		intResCompNames<-.internalComponents#intersect(.internalComponents, argsFLogLik$includeStreams)
	})
	.tmp.f <- function(){
		argsFRun2 <- runClusterParmsInit$argsFRun
		argsFRun2$xProp <- argsFRun2$xProp[1:500,]
		#tmp <- argsFRun2$fLogLik; mtrace(tmp); argsFRun2$fLogLik <- tmp
		#mtrace(deriv.SoilMod_FGSD); 
		#argsFRun2$argsFLogLik$useRImpl=TRUE
		#tmp <- argsFRun2$argsFLogLik$model$fSolve; mtrace(tmp); argsFRun2$argsFLogLik$model$fSolve <- tmp
		tmp <- do.call(twRunFLogLikPar, argsFRun2)
		resLogLik <- tmp
	}
	save( runClusterParms, file=paramFilename )
	iproc=0;outFilename0 <- twResultFilename(paramFilename,iproc=iproc); outFilenameAsom1 <- twResultFilename(paramFilenameAsom,iproc=iproc)
	#unlink( outFilename0 )
	nProc=min(7,nrow(.xProp))
	nProc=min(15,nrow(.xProp))
	nProc=min(16,nrow(.xProp))
	#	nP = min(8,dim(demcDrivers$Zinit[3])); copy2clip(as.character(GString("bsub -n ${nP} ~/bsubr_i.sh ${`getHamerDemcFilename('FGSD')`} iproc=41 nprocSingle=${nP} 'paramFile=\"parms/demcDrivers_1tFGSD.RData\"'")))
	#copy2clip(as.character(GString("~/bsubr_iLocal.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\\\"${paramFilenameAsom}\\\"' 'argsFRun=list(nGen=16)'" )))	#testing with 16 generations on pc026
	#copy2clip(as.character(GString("~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=16)'" )))	#testing with 16 generations on pc026
	copy2clip(as.character(GString("~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=4 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=16)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=16)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"'" )))
	
	
	load(outFilename0)
	#load(file.path("..","parms","res_res_MDIHamer_FGSD1_1_1.RData"))
	resLogLik <- runClusterRes$res
	ZinitR <- replaceZinitNonFiniteLogLiks( Zinit, resLogLik$logLik)
	Zinit2 <- ZinitR[,,1:8]	#only two populations
	Zinit3 <- ZinitR[,,1:12]	#three populations
	
	TInit <-  calcDEMCTempDiffLogLik3Init( resLogLik, TFix=TFix, pTarget=0.16 )	# in twDEMC twDEMC.R
	sort(TInit, decreasing=TRUE)
	
	#-------  first 3*1024 with high temperature and multistep Metropolis
	#experiment with T0, so that chains are not parallel wit 4 pops, when found repeat with 8
	.thin=8
	argsTwDEMCBatch0 <- within( list(),{
			#Zinit<-ZinitR;nPops<-.nPops
			#Zinit<-Zinit2; nPops<-2			
			Zinit<-Zinit3; nPops<-3			
			fLogLik<-ofb.hamer
			argsFLogLik<-argsFLogLik
			#loosing variability with > 20 parameters nBatch<-512
			#nBatch=max(minNBatch,1024)
			#nGen<-3*1024
			#enough recored steps after first batch, and enough cases per population to estimate Temperature from DiffLogLik
			#nBatch<- max(512,.minNBatchT <- max( ncol(Zinit)*.thin, ceiling(128/(.nChainsPop<-dim(Zinit)[3]/nPops)) ))
			nBatch<-512*2	#with less nBatch there is a deep drop of acceptance rate at second batch, when the initial states have been discarded
			nGen <- 512*2.2
			#controlTwDEMC<-list(thin=.thin,Tprop=TInit,TFix=TFix)
			controlTwDEMC<-list(thin=.thin,TFix=TFix,useMultiT=TRUE,Tprop=TInit,minPCompAcceptTempDecr=0.16)
			#T0<-nObs*16
			T0<-max(TInit)/2
			#T0=40000
			nGenBurnin<-nGen*1000	#essentially not decreasing temp
			#doRecordProposals=TRUE
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
		
		#mtrace(thin.twDEMC)
		#nGen0 <- res$thin+2*1024
		#resTmp3 <- thin(res, end=nGen0)	#from out1
		#argsFRunDebug$argsTwDEMCBatch
		.prevDir <- setwd("..")
		#do.call( fSetupCluster, list()  )
		#mtrace(ofb.hamer)
		argsFRunDebug <- runClusterParms$argsFRun 
		argsFRunDebug$argsTwDEMCBatch <- within(argsFRunDebug$argsTwDEMCBatch, {
				Zinit=Zinit[,,1:8,drop=FALSE]
				#Zinit=subChains(resTmp3,iPops=1:2)
				nPops=2
				nBatch=64
				nGen=2*64+3*controlTwDEMC$thin
				#nGen=nGen0+64
				#controlTwDEMC$thin=4
				#debugSequential=TRUE
				#intResCompNames=.internalComponents
				#model$fLogLik<-ofb.hamer
				#doStopOnError=TRUE
				#argsFLogLik$useRImpl=TRUE
				doRecordProposals=TRUE
			})
		#mtrace(ofb.hamer)
		argsFRunDebug$fLogLik<-ofb.hamer
		#options(warn=2)
		#options(error=dump.frames)
		#mtrace(twRunDEMC)
		#mtrace(replaceZinitNonFiniteLogLiksLastStep)
		#mtrace(twCalcLogLikPar)
		#tmp<-argsFRunDebug$argsTwDEMCBatch$argsFLogLik$model$fSolve; mtrace(tmp); argsFRunDebug$argsTwDEMCBatch$argsFLogLik$model$fSolve<-tmp
		#tmp<-ofb.hamer; mtrace(tmp); argsFRunDebug$argsTwDEMCBatch$fLogLik<-tmp
		#mtrace(deriv.SoilMod_FGSD)
		#mtrace(.doDEMCStep)
		#mtrace(.doDEMCSteps)
		#mtrace(twDEMCInt)
		#mtrace(calcDEMCTempDiffLogLikConst)
		#mtrace(twDEMCBatchInt)
		#mtrace(calcDEMCTempDiffLogLik2)
		tmpRes <- NULL; tmpRes <-do.call( twRunDEMC, argsFRunDebug )
		
		
		tmpRes$Y["parms",,]
		
		
		matplot( tmpRes$temp , type="l")
		matplot( tmpRes$pAccept , type="l"); abline(h=0.1,col="gray"); abline(v=nGen0/8)
		#matplot( tmpRes3$pAccept , type="l"); abline(h=0.1,col="gray")
		argsFRunDebug2 <- argsFRunDebug; argsFRunDebug2$argsTwDEMCBatch$nGen = (nrow(tmpRes$rLogLik)-1)*tmpRes$thin+2*32; argsFRunDebug2$argsTwDEMCBatch$Zinit <- tmpRes
		#mtrace(twDEMCBatchInt)
		tmpRes2 <- NULL; tmpRes2 <- do.call( twRunDEMC, argsFRunDebug2)
		matplot( tmpRes2$temp , type="l")
		matplot( tmpRes2$pAccept , type="l"); abline(h=0.1,col="gray")
		setwd(.prevDir)
	}
	
	iproc=21;outFilename1 <- twResultFilename(paramFilename,iproc=iproc); outFilenameAsom1 <- twResultFilename(paramFilenameAsom,iproc=iproc)
	iproc=1;outFilename1 <- twResultFilename(paramFilename,iproc=iproc); outFilenameAsom1 <- twResultFilename(paramFilenameAsom,iproc=iproc)
		
	save( runClusterParms, file=paramFilename )
	unlink( outFilename1 )
	nmax <- dim(runClusterParms$argsFRun$argsTwDEMCBatch$Zinit)[3]
	nProc=min(8,nmax)
	nProc=min(14,nmax)
	nProc=min(16,nmax)
	#	nP = min(8,dim(demcDrivers$Zinit[3])); copy2clip(as.character(GString("bsub -n ${nP} ~/bsubr_i.sh ${`getHamerDemcFilename('FGSD')`} iproc=41 nprocSingle=${nP} 'paramFile=\"parms/demcDrivers_1tFGSD.RData\"'")))
	#copy2clip(as.character(GString("~/bsubr_iLocal.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\\\"${paramFilenameAsom}\\\"' 'argsFRun=list(nGen=16)'" )))	#testing with 16 generations on pc026
	#copy2clip(as.character(GString("~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=16)'" )))	#testing with 16 generations on pc026
	copy2clip(as.character(GString("~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=4 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=16)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=16)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"'" )))

	load(outFilename1)
	#load(outFilename11)
	#load(file.path("..","parms","res_res_MDIHamer_FGSD1_1_1.RData"))
	res <- runClusterRes$res
	nRun1 <- (nrow(res$rLogLik)-1)*res$thin
	#to extend the run, copy the results file to param file (res will be used) and change nGen
	#file.copy(outFilename1,paramFilename,overwrite=TRUE);	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilename}\"' 'argsFRun=list(nGen=${`nRun1+1*1024`})'" )))
	
	.start=16
	#.start=145
	#windows(record=TRUE)
	par(mfrow=c(2,1))
	matplot( res$temp, type="l")
	abline(v=512/8*1:16, col="gray")
	matplot( twDEMCPopMeans(res$pAccept, ncol(res$temp), 1), type="l" )
	abline(h=0.1,col="gray")
	abline(v=512/8*1:16, col="gray")

	#.start=140
	matplot( res$rLogLik[-(1:.start),], type="l")
	rest <- thin(res, start=nRun1-512)
	matplot( rest$rLogLik, type="l", ylim=max(rest$rLogLik)*c(20,1.1))
	matplot( rest$parms["sA",,], type="l")
	matplot( twDEMCPopMeans(rest$pAccept, ncol(rest$temp), 1), type="l" )
	abline(h=0.1,col="gray")
	#windows(record=TRUE); rescoda <- as.mcmc.list(rest); plot(thinN(rescoda), smooth=FALSE)
	#
	#plot(thinN(rescoda), smooth=FALSE)
	#rest <- thin(res, start=(1024+512))
	#mtrace(ggplotDensity.twDEMC)
	p1 <- ggplotDensity.twDEMC( rest, poptDistr, pMin=0.1)
	#print(p1)

	#-------  second decrease temperature quite fast and observe retardation by acceptance rate
	# now do not use multistep Metropolis any more
	# tried with going down in 4*1024 steps with holding temp if acceptance rate drops to low 
	
	iproc=1;outFilename1 <- twResultFilename(paramFilename,iproc=iproc)
	#iproc=11;outFilename11 <- twResultFilename(paramFilenameB,iproc=iproc);load(outFilename11)
	load(outFilename1)
	#load(file.path("..","parms","res_res_MDIHamer_FGSD1_1_1.RData"))
	res <- runClusterRes$res
	res0 <- runClusterRes$argsFRun$argsTwDEMCBatch$Zinit
	runClusterRes$res <- runClusterRes$argsFRun$prevResRunCluster <- NULL	#else interferes with continued run, also save space
	
	#nRun0 <- (nrow(res0$rLogLik)-1)*res0$thin
	nRun1 <- (nrow(res$rLogLik)-1)*res$thin
	runClusterParms <- runClusterParms0
	runClusterParms$argsFRun$argsTwDEMCBatch <- within( runClusterRes$argsFRun$argsTwDEMCBatch,{
			Zinit<-res
			nBatch<-512	#	want to see output of progress and remove outlier chains 
			#nGen<-nRun1+32*1024	
			nGen<-nRun1+10*1024	
			nGenBurnin<-nRun1+8*1024	
			nPops<-ncol(res$temp)
			controlTwDEMC$minPCompAcceptTempDecr=0.18
			
	})
	runClusterParms$argsFRun$resFLogLikX<-character(0)
	runClusterParms$fRun <- twRunDEMC		# update functions to represent current dev state 
	runClusterParms$argsFRun$fLogLik <- ofb.hamer
	runClusterParms1 <- runClusterParms
	
	iproc=2; outFilename2 <- twResultFilename(paramFilename,iproc=iproc)
	#iproc=22; outFilename2 <- twResultFilename(paramFilename,iproc=iproc)
	save( runClusterParms, file=paramFilename )
	unlink( outFilename2 )

	#nProc=min(14,dim(Zinit)[3])
	nProc=min(16,dim(runClusterParms$argsFRun$argsTwDEMCBatch$Zinit$parms)[3])
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${`nRun1+16`},doRecordProposals=TRUE)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${nRun1}+1024*1)'" )))
	#copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' doRestart=TRUE" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"'" )))

	#--------	
	load(outFilename2)
	res <- runClusterRes$res
	res0 <- runClusterRes$argsFRun$argsTwDEMCBatch$Zinit
	runClusterRes$res <- runClusterRes$argsFRun$prevResRunCluster <- NULL	#else interferes with continued run, also save space
	nRun0 <- (nrow(res0$rLogLik)-1)*res0$thin
	nRun1 <- (nrow(res$rLogLik)-1)*res$thin
	
	.start=nRun0*2; .starti <- .start%/% res$thin
	.start=nRun1-3*512; .starti <- .start%/% res$thin
	#.start=res$nGenBurnin; .starti <- .start%/% res$thin
	#tmp<-apply(res$temp,1,max); .starti= min(which(tmp<80)); .start <- .starti*res$thin
	#tmp<-apply(res$temp,1,max); .starti= min(which(tmp<5)); .start <- .starti*res$thin
	#windows(record=TRUE)
	matplot( res$temp, type="l")
	matplot( res$temp[-(1:.starti),], type="l")
	matplot( twDEMCPopMeans(res$pAccept[-(1:.starti),], ncol(res$temp), 2), type="l" )
	abline(h=0.1,col="gray")
	matplot( res$rLogLik[-(1:.starti),], type="l")
	
	tmp<-apply(res$temp,1,max); .starti= min(which(tmp<80)); .start <- .starti*res$thin
	
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${`nRun1+5*1024`},nGenBurnin=${`res$nGenBurnin`})' doRestart=TRUE" )))
	
	rest0 <- rest <- thin(res, start=nRun1-4*512)
	rest0 <- rest <- thin(res, start=nRun1-3*512)
	rest0 <- rest <- thin(res, start=nRun1-256/mean(res$pAccept[nrow(res$pAccept),]))
	#rest0 <- thin(res, start=res$nGenBurnin)
	#mtrace(ggplotDensity.twDEMC)
	p0 <- ggplotDensity.twDEMC( rest0, poptDistr)
	#p0
	p0b <- ggplotDensity.twDEMC( rest0, poptDistr, doTransOrig=TRUE)
	#p0b

	tmp <- plotChainPopMoves(rest0)
	rest <- rest0
	matplot( rest$rLogLik[,], type="l")
	
	rest <- subChains( rest0, iPops=3)
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
	res <- resOrig <- runClusterRes$res
	#res <- subChains( resOrig, iPops=(1:8)[-c(3,7)])		#will cancel attribute batchCall
	res0 <- runClusterRes$argsFRun$argsTwDEMCBatch$Zinit
	runClusterRes$res <- runClusterRes$argsFRun$prevResRunCluster <- NULL	#else interferes with continued run, also save space
	nRun1 <- (nrow(res$rLogLik)-1)*res$thin
	nRun0 <- (nrow(res0$rLogLik)-1)*res0$thin

	iStep <- nRun1-nRun0	#nRun0 from outFilename1 above
	.burnin <- max( res$nGenBurnin,ceiling(nRun0+calcDEMCnGenBurnin(mean(res$temp[nRun0/res$thin,]), mean(res$temp[nRun1/res$thin,]), iStep)))
	.burnin <- nRun1 + 8*1024	
	
	runClusterParms <- runClusterRes
	runClusterParms$argsFRun$argsTwDEMCBatch <- within( runClusterParms$argsFRun$argsTwDEMCBatch,{
			Zinit<-res
			nBatch<-512	#	want to see output of progress and remove outlier chains 
			nGen<-.burnin+5*1024	
			nGenBurnin<-.burnin
			nPops <- ncol(res$temp)
		})
	runClusterParms$intResCompNames<-.internalComponents#intersect(.internalComponents, argsFLogLik$includeStreams)
	runClusterParms$fRun <- twRunDEMC		# update functions to represent current dev state 
	runClusterParms$argsFRun$fLogLik <- ofb.hamer
	iproc=3; outFilename3 <- twResultFilename(paramFilename,iproc=iproc)
	save( runClusterParms, file=paramFilename )
	unlink( outFilename3 )
	
	#nProc=min(15,dim(Zinit)[3])
	nProc=min(16,dim(res$parms)[3])
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${`nRun1+16`},doRecordProposals=TRUE)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${nRun1}+1024*1)'" )))
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"'" )))
	
	load(outFilename3)
	res <- runClusterRes$res
	res0 <- runClusterRes$argsFRun$argsTwDEMCBatch$Zinit
	runClusterRes$res <- runClusterRes$argsFRun$prevResRunCluster <- NULL	#else interferes with continued run, also save space
	nRun0 <- (nrow(res0$rLogLik)-1)*res0$thin
	nRun1 <- (nrow(res$rLogLik)-1)*res$thin
	
	.start=nRun0; .starti <- .start%/% res$thin
	#windows(record=TRUE)
	matplot( res$temp[-(1:.starti),], type="l")
	matplot( twDEMCPopMeans(res$pAccept[-(1:.starti),], ncol(res$temp), 10), type="l" )
	matplot( res$rLogLik[-(1:.starti),], type="l")
	copy2clip(as.character(GString("bsub -n ${nProc} ${bsubOptions} ~/bsubr_i.sh runCluster.R iproc=${iproc} nprocSingle=${nProc} 'paramFile=\"${paramFilenameAsom}\"' 'argsFRun=list(nGen=${`nRun1+5*1024`},nGenBurnin=${`res$nGenBurnin`})' doRestart=TRUE" )))
	
	rest <- rest0 <-  thin(res, start=nRun1-512)
	matplot( twDEMCPopMeans(rest0$pAccept[,], ncol(rest0$temp), 10), type="l" )
	matplot( rest0$rLogLik[,], type="l")
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
	#mtrace(ofb.hamer)
	tmp2 <- hamerPlot3.of( normpoptBest2, ofb.hamer, argsFLogLik=tmp$argsFLogLik)
	tmpDev <- dev.cur()
	windows(record=TRUE); tmpMai<-par()$mai;tmpMai=c(0.8,2,0.4,0.4); par(las=1, mai=tmpMai)
	barplot(-unlist(as.list(tmp2)), horiz = TRUE)
	barplot(attributes(tmp2)$logLikParms, horiz = TRUE)
	dev.set(tmpDev)
	
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
	pdev <- dev.cur(); windows(record=TRUE); par(mfrow=c(2,1))
	resOf <- hamerPlot3.of( normp, argsFLogLik=tmp$argsFLogLik) 
	unlist(as.list(resOf))
	sort(unlist(as.list(resOf)))
	sort(attr(resOf,"logLikParms"), decreasing=FALSE)

	#mtrace(hamerPlotResCols.resOf)
	#colnames(attributes(resOf)$out$amendm)
	pdev <- dev.cur(); windows(record=TRUE); par(mfrow=c(2,1))
	hamerPlotResCols.resOf( resOf, paste("decF",model$modMeta$colNames,sep="_"), "amendm" )
	hamerPlotResCols.resOf( resOf, paste("decG",model$modMeta$colNames,sep="_"), "amendm" )
	hamerPlotResCols.resOf( resOf, c("decF_s","decG_a"), "amendm" )
	hamerPlotResCols.resOf( resOf, c("propG"), "amendm" )
	hamerPlotResCols.resOf( resOf, paste("csums",model$modMeta$rowNames,sep="_"), "control" )
	hamerPlotResCols.resOf( resOf, paste("csums",c("F","A"),sep="_"), "control" )
	hamerPlotResCols.resOf( resOf, paste("csums",c("F","A"),sep="_"), "control" )
	hamerPlotResCols.resOf( resOf, paste(c("respT","decF_s"),sep="_"), "control" )
	dev.set(which = pdev)
	
	matplot( res$temp[-(1:.starti),], type="l")
	matplot( t(res$parms[c("sA","F0"),-(1:.starti),1]), type="l")
	matplot( t(res$parms[c("sA","F0"),-(1:1),1]), type="l")
	matplot( t(res$parms[c("sA","F0"),100:900,1]), type="l")
	matplot( t(res$parms[c("sA","F0"),200:600,5]), type="l")
	matplot( res$rLogLik[-(1:.starti),], type="l")
	matplot( res$rLogLik[-(1:1),], type="l")
	matplot( res$rLogLik[200:600], type="l")
	
	matplot( (res$parms[c("sA"),-(1:.starti),]), type="l")
	matplot( (res$parms[c("kS"),-(1:.starti),]), type="l")
	checkConvergenceGelman(rest)
	
	
	
}










