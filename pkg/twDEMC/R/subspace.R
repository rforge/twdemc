findSplit <- function(
	### determine the parameter and its value, where to split the parameter space
	aSample	##<< the sample (stacked chains)
	,nSplit = 4 ##<< how many points per parameter dimension to check for split
	,rVarCrit = 10^2	##<< minimla ratio of variances of a variable between right and left side of a splitting point
	,rAlphaSlopeCrit = base:::pi/3	##<< minimal ratio of angle between scaled slopes, defaults to a third of a half cirle
	,pSlopeMax=0.05					##<< confidence level for the slope, above which no difference in slopes is calculated 
	,iVars = 1:ncol(aSample)		##<< integer vector: index or parameter dimensions, i.e. variables, to check for split		
	,jVars = 1:ncol(aSample)		##<< integer vector: index or parameter dimensions, i.e. variables, to check for different scales of variance
){
	##details<< 
	## First it checks for different scales of variance in other variables. 
	if( is.character(iVars) ) iVars <- sapply(iVars, match, colnames(aSample))
	if( is.character(jVars) ) jVars <- sapply(jVars, match, colnames(aSample))
	foundSplit <- FALSE
	percentiles <- 1:(nSplit)/(nSplit+1)
	quantVar <- iVarLeft <- iVarRight <- matrix( NA_real_, nrow=length(iVars), ncol=nSplit )
	slope <- parVarLeft <- parVarRight <- array( NA_real_, dim=c(length(iVars),length(jVars),nSplit)
	, dimnames=list(iVar=colnames(aSample)[iVars], jVar=colnames(aSample)[jVars], iSub=NULL ))
	ssSubLeft <- ssSubRight <- vector(mode="list",length=length(iVars) )
	#i <- 1
	for( i in seq_along(iVars)){
		iVar <- iVars[i]
		p <- aSample[,iVar] 
		qp <- quantVar[i, ] <- quantile( p, percentiles )
		#iSub <- 1
		#iSub <- nSub
		ssSubLeft[[i]] <- lapply( 1:(nSplit), function(iSplit){
				ssSub <- aSample[ (p <= qp[iSplit]),]
			})
		ssSubRight[[i]] <- lapply( 1:(nSplit), function(iSplit){
				ssSub <- aSample[ (p >= qp[iSplit]),]
			})
		# ratio of variances of subsets left and right of splitting points
		#j <- 1
		for( j in seq_along(jVars)){
			jVar <- jVars[j]
			#iSplit=1
			if( jVar!=iVar){
				pVar <- rep( NA_real_, nSplit)
				for(iSplit in 1:nSplit){
					varLeft <- parVarLeft[i,j,iSplit] <- var(ssSubLeft[[i]][[iSplit]][,jVar])
					varRight <- parVarRight[i,j,iSplit] <- var(ssSubRight[[i]][[iSplit]][,jVar]) 
					#plot( ssSubLeft[[i]][[iSplit]][,iVar], ssSubLeft[[i]][[iSplit]][,jVar] )
					pVar[iSplit] <- varLeft /varRight
				}
				pVar[ pVar < 1 ] <- 1/pVar[pVar<1]
				if( max(pVar) > rVarCrit){
					iSplit <- which.max(pVar)
					res <- list( 
						split=structure( quantVar[i,iSplit], names=colnames(aSample)[iVar])
						, varName=colnames(aSample)[iVar]
						, perc=percentiles[iSplit] )
					return(res)
				}
			} # end if( jVar!=iVar)
		} # end for jVar variance
	} # for iVar
	
	##details<< 
	## Next it checks for different angles of the normalized slopes in relation of the parameters 
	#i<-1
	if( !foundSplit) for( i in seq_along(iVars)){
		# check difference in slope angles
		iVar <- iVars[i]
		# also need the variance of iVar for scaling the slope
		iPos <- match( iVar, jVars )
		iVarLeft[i, ] <- sapply( 1:(nSplit), function(iSplit){ var(ssSubLeft[[i]][[iSplit]][,iVar]) })
		iVarRight[i, ] <- sapply( 1:(nSplit), function(iSplit){ var(ssSubRight[[i]][[iSplit]][,iVar]) })
		#qp <- quantVar[i, ]
		#j <- 1
		for( j in seq_along(jVars)){
			jVar <- jVars[j]
			if( jVar!=iVar){
				#iSplit=4
				dAlphaSlope <- sapply( 1:nSplit, function(iSplit){
						#XXTODO names are not know to construct formula, lookup direct equation for
						#plot( ssSubLeft[[i]][[iSplit]][,iVar], ssSubLeft[[i]][[iSplit]][,jVar] )
						#plot( ssSubRight[[i]][[iSplit]][,iVar], ssSubRight[[i]][[iSplit]][,jVar] )
						lmLeft <- lm( as.formula(paste(colnames(aSample)[c(jVar,iVar)],collapse=" ~ ")), data=as.data.frame(ssSubLeft[[i]][[iSplit]][,c(jVar,iVar)]) )
						lmRight <- lm( as.formula(paste(colnames(aSample)[c(jVar,iVar)],collapse=" ~ ")), data=as.data.frame(ssSubRight[[i]][[iSplit]][,c(jVar,iVar)]) )
						#plot( ssSubLeft[[i]][[iSplit]][,iVar], ssSubLeft[[i]][[iSplit]][,jVar] ); abline(lmLeft)
						#plot( ssSubRight[[i]][[iSplit]][,iVar], ssSubRight[[i]][[iSplit]][,jVar] ); abline(lmRight)
						if( (anova(lmLeft)[["Pr(>F)"]][1] < pSlopeMax) || (anova(lmRight)[["Pr(>F)"]][1] < pSlopeMax) ){ 
							slopeLeft <- coef(lmLeft)[2] * sqrt(iVarLeft[i,iSplit]/parVarLeft[i,j,iSplit]) 
							slopeRight <- coef(lmRight)[2] * sqrt(iVarRight[i,iSplit]/parVarRight[i,j,iSplit])
							alphaRight <- asin( slopeLeft)
							alphaLeft <- asin( slopeRight )
							abs( alphaLeft - alphaRight )
						}else{ 0 }  # both slopes are non-significant
					})
				if( max(dAlphaSlope) > rAlphaSlopeCrit){
					iSplit <- which.max(dAlphaSlope)
					res <- list( 
						split=structure( quantVar[i,iSplit], names=colnames(aSample)[iVar])
						, varName=colnames(aSample)[iVar]
						, perc=percentiles[iSplit] )
					return(res)
				}
			} # end if( jVar!=iVar)
		} # end for jVar variance
	} # for iVar

	##value<< list with components  
	return( list(
		split=NA_real_	##<< named scalar: splitting value with the name of the parameter dimension that is to split
		,varName=NA_character_ ##<< name of the splitting variable
		,perc=NA_real_	##<< scalar: percentile of the splitting point
	))
	##end<< 
}
attr(findSplit,"ex") <- function(){
	data(den2dCorTwDEMC)
	ss1 <- ss <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
	#mtrace(findSplit)
	(res <- res1 <- findSplit(ss1))
	
	# successively find splits in subsets
	ss2 <- ss <- ss1[ ss1[,res$varName] > res$split,]
	res2 <- findSplit(ss2)
	ss3 <- ss <- ss2[ ss2[,res2$varName] > res2$split,]
	res3 <- findSplit(ss3)
	# NA indicates: no further split found
	
	# try different orders of the variables
	#mtrace(findSplit)
	(res <- findSplit(ss1, iVars=c(2,1) ))
	(res <- findSplit(ss1, iVars=c("b") ))
	plot( ss1[,1], ss1[,2], col=c("blue","red")[as.integer(ss1[,res1$varName] >= res1$split)+1])
}



getSubSpaces <- function( 
	### Recursively splitting the sample into subsamples and recording the upper and lower bounds
	aSample 			##<< numeric matrix with parameters in columns: sample or subsample that is devided
	,nSplit=4			##<< see \code{\link{findSplit}}, potentially modified due to argument \code{minPSub}
	,argsFSplit=list()	##<< further arguments to \code{\link{findSplit}}
	,pSub=1				##<< the fraction that a subSample comprises within a bigger sample (used in recursive calls) 
	,minPSub=0.05		##<< minimum fraction a subSample, below which the sample is not split further
){
	# search for a single splitting point
	nSplitMin <- min( nSplit, floor(pSub/minPSub)-1 )		# each subPopulation must comprise at least part minPSub, hence do not split in too many parts
	argsFSplitMod <- within(argsFSplit, nSplit<-nSplitMin )
	boMinPSub <- (nSplitMin < 1)
	if( !boMinPSub) resSplit <- do.call( findSplit, c(list(aSample=aSample), argsFSplitMod ))
	if( boMinPSub || is.na(resSplit$split) ){
		# when no splitting point was found, return the sample without parameter bounds
		##value<< list with an entry for each subspace. Entry is a list with 
		list( list( 
			sample = aSample		##<< numeric matrix: a subsample constrained to the subspace with col parameters
			, upperParBounds = c()	##<< list with each entry numeric scalar: upper parameter bounds
			, lowerParBounds = c()	##<< list with each entry numeric scalar: upper parameter bounds
			, pSub = pSub			##<< the proportion of the subSample to the overall Sample
		))	##end<<
	}else{
		# split the sample and recursively invoke \code{getSubSpaces}
		boLeft <- (aSample[,resSplit$varName] <= resSplit$split)
		sampleLeft <- aSample[boLeft, ]
		sampleRight <- aSample[!boLeft,]
		pSubLeft <- pSub*resSplit$perc
		pSubRight <- pSub*(1-resSplit$perc)
		spacesLeft <- getSubSpaces(sampleLeft, nSplit, argsFSplit, pSubLeft, minPSub)
		spacesRight <- getSubSpaces(sampleRight, nSplit, argsFSplit, pSubRight, minPSub)
		# add the splitting point as a new border to the results
		resLeft <- lapply( spacesLeft, function(spEntry){
				if( !(resSplit$varName %in% names(spEntry$upperParBounds)) )
					spEntry$upperParBounds <- c(resSplit$split, spEntry$upperParBounds)
				spEntry
			})
		resRight <- lapply( spacesRight, function(spEntry){
				if( !(resSplit$varName %in% names(spEntry$lowerParBounds)) )
					spEntry$lowerParBounds <- c(resSplit$split, spEntry$lowerParBounds)
				spEntry
			})
		c( resLeft, resRight )
	}
} 
attr(getSubSpaces,"ex") <- function(){
	data(den2dCorTwDEMC)
	aSample <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
	#mtrace(getSubSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.4)
	str(subSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.05)
	(tmp <- sapply( subSpaces, function(subSpace){nrow(subSpace$sample)})/nrow(aSample))
}



checkProblemsGelman1 <- function(
	### convergence between several populations
	aTwDEMC
	, criticalGelmanDiag=1.2
){
	
	resCoda <- as.mcmc.list(aTwDEMC)
	#str(resCoda)
	#plot( resCoda, smooth=FALSE )
	tmp2 <- gelman.diag(resCoda)
	list(
		hasProblem = (tmp2$mpsrf > criticalGelmanDiag)
		,resGelmanDiag = tmp2
	)
}

checkProblemsSpectral <- function(
	aTwDEMC
	, critSpecVarRatio=50
){
	aSamplePure <- stackChainsPop(aTwDEMC)[,-1,,drop=FALSE] 
	specVarRatio <- apply( aSamplePure, 3, function(aSamplePurePop){
			spec <- spectrum0.ar(aSamplePurePop)$spec
			varSample <- apply(aSamplePurePop, 2, var)
			spec/varSample
		})
	list(
		hasProblems = any(specVarRatio >= critSpecVarRatio)
		,specVarRatio = specVarRatio 
	)
}



divideTwDEMC <- function(
	### run twDEMC on subspaces
	aSample					##<< numeric array: rows:steps, col:parameters without logLik, 3: independent populations
	,nGen=512				##<< the number of generations for twDEMC
	, controlTwDEMC = list()	##<< list argument to \code{\link{twDEMC}} containing entry thin
	, ...					##<< further arguments to \code{\link{twDEMC}}
	, minPSub = 0.05		##<< see \code{\link{getSubSpaces}}
	, nChainsPop=4			##<< number of chains in subPopulations
	, m0 = calcM0twDEMC( ncol(aSample), nChains=nChainsPop )	##<< number of samples per chain to initialize subPopulations
	, nrow0 = 4000			##<< number of rows in initial overall sample
	, attachDetails=FALSE	##<< set TRUE to report upperParBounds, lowerParBounds, and pSubs per subPopulation
	, nChainPar=16			##<< number of chains to run in parallel, good choice is 2*nCpu
	, minNSamplesSub=32		##<< minimum number of records in a subspace sample, increase to avoid wrong estimation of weights
	, fCheckProblems=checkProblemsSpectral		
	, argsFCheckProblems=list()
	, dumpfileBasename="recover"
){
	##details<< 
	## the first column of aSample records the logDensity of the sample for consitency. 
	## It is not used and may be initialized to any value. 
	#samplePop <- aSample[,,1]
	if( length(dim(aSample))==2 ){
		warning("divideTwDEMC:third dimension missing of aSample missing - assuming 1 population.")
		aSample <- array(aSample, dim=c(dim(aSample),1), dimnames=c(dimnames(aSample), pops=NULL) )
	}
	nPops <- dim(aSample)[3]
	if( is.null(controlTwDEMC$thin) ) controlTwDEMC$thin <- 4
	thin <- controlTwDEMC$thin
	nGenThin <- (nGen %/% thin)*thin		# when dealing with thinned samples use the number of generations of last recorded thinned generation
	
	thinnedSample <- if( nrow(aSample) < nrow0 ) aSample else aSample[round(seq(1,nrow(aSample),length.out=nrow0)),,]
	nInit <- m0*nChainsPop		# number of samples to initialize population
	if( nInit > nrow(thinnedSample) ) stop(paste("divideTwDEMC: aSample or nGen has too few cases, need at least",nInit))
	minPSub <- max( minPSub, nInit/nrow(thinnedSample) )
	#mtrace(getSubSpaces)
	subSpaces <- apply( thinnedSample, 3, function(samplePop){
			#mtrace(getSubSpaces)
			getSubSpaces(samplePop[,-1], minPSub=minPSub)	# here omit the logDensity column
		})
	nSubPops <- sapply(subSpaces, length)	# number of subs per populations
	subSpacesFlat <- do.call( c, subSpaces )	# put all subPopulations of all populations on same level
	cumNSubPops <- cumsum(nSubPops)
	iSubPops <- lapply(1:nPops, function(iPop){ cumNSubPops[iPop]+1-(nSubPops[iPop]:1)})	# index of subs in flat version for each pop 			
	ZinitSubs <- abind( lapply( subSpacesFlat, function(ssEntry){
						#s0 <- ssEntry$sample[ round(seq(1,nrow(ssEntry$sample),length.out=nInit)), ]
						s0 <- ssEntry$sample[ sample.int(nrow(ssEntry$sample),nInit), ]	# here use random draws instead thinned interval
						Zinit <- array( t(s0), dim=c(ncol(s0),m0,nChainsPop), dimnames=list(parms=colnames(s0),steps=NULL,chains=NULL) )
					}), along=3)
	upperParBounds <- lapply( subSpacesFlat, function(ssEntry){ ssEntry$upperParBounds})
	lowerParBounds <- lapply( subSpacesFlat, function(ssEntry){ ssEntry$lowerParBounds})
	popQSubs <- lapply( subSpaces, sapply, "[[","pSub" ) # initial percentile of subspace of the overall space
	qSubs <- do.call(c,popQSubs)
	nGen0 <- pmax(m0*thin, ceiling(minNSamplesSub*thin/nChainsPop), nGen*qSubs)		# at minimum m0*thin generations to keep sufficient samples for extending the run
	nGen0Thin <- (ceiling(nGen0/thin)+2)*thin		# +2 to avoid extending runs with little shift 
	
	# in one call of twDEMC all chains have the same number of generations
	# thus, in order to save time, sort the populations by required number of generations and 
	# start several sequential calls to twDEMC differing by number of generations
	nSubs <- length(nGen0)
	nPopPar <- ceiling( nSubs / ceiling(nChainPar/nChainsPop))	# number of populations to run in parallel

	# in which twDEMCRun does subPopulation of index iis executed: iRunSubs
	iTwDEMCRun <- rep(1:(nSubs%/%nPopPar), each=nPopPar)
	nTwDEMCRun <- iTwDEMCRun[length(iTwDEMCRun)] 
	iTwDEMCRun <- c(iTwDEMCRun, rep(nTwDEMCRun, nSubs-length(iTwDEMCRun)) ) # append remaining populations to last twDEMCRun
	tmpOrd <- order(nGen0, decreasing =TRUE)
	iRunSubs <- sapply(1:nSubs, function(iSub){ iTwDEMCRun[ match(iSub,tmpOrd)] }) # iRun for sub of given index
	
	# mapping the sub-populations (subs) to runs and population within run
	# pop will be assigned below
	popSubPops <- do.call( c, lapply( 1:nPops, function(iPop){ rep(iPop,nSubPops[iPop])} ))		# each index gives the population that the subPopulatios is part of
	subInfo <- cbind( 
		pop=popSubPops			##<< population
		, run=iRunSubs			##<< index of run in resITwDEMC
		, rpop=NA_integer_		##<< population within run
		, nGen=NA_integer_		
		, nSamples=NA_integer_ 
	)
	resITwDEMC <- vector("list",nTwDEMCRun)
	#iRun <- 1
	for( iRun in 1:nTwDEMCRun ){
			iSubsRun <- which( subInfo[,"run"] == iRun )
			subInfo[,"rpop"][iSubsRun] <- tmp <- seq_along(iSubsRun)   
			#subInfo[,"rpop"][iSubsRun] <- tmp2 <- match(iSubsRun, which(iRunSubs==iRun)) 
			#if( !identical(tmp,tmp2) ) recover()
			nGenIRun <- max( nGen0Thin[ iSubsRun ] )
			#resTwDEMC <- 
			resITwDEMC[[iRun]] <- resTwDEMC <- twDEMC(
				ZinitSubs[,,as.vector(sapply( (iSubsRun-1)*nChainsPop, "+", (1:nChainsPop)))]
				, nPops=length(iSubsRun)
				, upperParBounds=upperParBounds[iSubsRun]
				, lowerParBounds=lowerParBounds[iSubsRun]
				, controlTwDEMC=controlTwDEMC
				, nGen=nGenIRun
				#, fLogDen=den2dCor
				, ...
			)
			subInfo[,"nGen"][iSubsRun] <- getNGen(resTwDEMC)
			subInfo[,"nSamples"][iSubsRun] <- getNSamples(resTwDEMC)
			#subInfo[iSubsRun,]
	}
	#pSubs1 <- subInfo[,"nGen"]/(nGen*pSubs0)		# sampling weights (>1 obtained more samples than required by quantile)
	
	##details<< 
	## A first estimate of the proportion of samples from different subspaces
	## are the initial percentiles qp. 
	## The proportion of the samples from different subspaces is estimated
	## by the proportions of integrated density of the subspaces.
	## These proportions are esimated by average density over all samples multiplied by an
	## estimate proportional to the volume: the initial quantiles.
	#iPop <- 1
	#iPop <- 2
	popLogMeanDensSubs <- lapply( 1:nPops, function(iPop){
			iSubs <- iSubPops[[iPop]]
			#iSub <- 2
			#iSub <- 14
			lw <- sapply(iSubs,function(iSub){
					#ssLogDen <- as.vector(ssSub[[iSub]][1,,])
					si <- subInfo[iSub, ]
					ssLogDen <- ssLogDenNew <- resITwDEMC[[ si["run"] ]]$rLogDen[,(si["rpop"]-1)*thin+(1:4) ]
					twLogSumExp(ssLogDen)-log(length(ssLogDen)) 
				})	# average the unnormalized densities
			# normalization is done below
			#lSumW <- twLogSumExp(lw) 
			#wSubs <- exp( lw-lSumW ) 			#w/sumW
			lw - max(lw)						#divide by a constant on exp scale to avoid numerical errors
		})	
	logMeanDensSubs <- do.call(c, popLogMeanDensSubs )
	#barplot(logMeanDensSubs)
	
	# estimate proportion of subspaces in the limiting distribution
	# it will be used to sample from the subspaces 
	pSubs <- do.call(c, lapply( 1:nPops, function(iPop){
			p2u <- exp(popLogMeanDensSubs[[iPop]])*popQSubs[[iPop]]		# two weights: quantile and density importance
			p2 <- p2u/sum(p2u)									# normalize to 1 within population 
		}))
	subPercChange <- pSubs/qSubs	# estimated proportion in limiting distribution to proportion of initial proportion 
	nSamplesSubsReq <- round(nGen/thin * pSubs)
	# because of rounding small differences in sum of sample numbers may occur 
	# for small deviations, adjust the number of samples in the largest sample
	for( iPop in 1:nPops ){
		dGen <- nGenThin/thin - sum(nSamplesSubsReq[iSubPops[[iPop]] ])  
		if( dGen != 0 ){
			if( abs(dGen)>length(iSubPops[[iPop]]) ) stop("divideTwDEMC: problem in calculating subspace sample numbers")
			iSubMax <- iSubPops[[iPop]][ which.max( nSamplesSubsReq[iSubPops[[iPop]] ]) ]
			nSamplesSubsReq[iSubMax] <- nSamplesSubsReq[iSubMax] + dGen
		} 
	}
	
	# extend chains of those populations that do not have enough samples yet
	nGenSubsTodo <- nSamplesSubsReq*thin - subInfo[,"nGen"]
	iSubsExt <- which(nGenSubsTodo > 0)
	prevIRunsExt <- subInfo[iSubsExt,"run"]
	uPrevIRunsExt <- unique(prevIRunsExt)
	#nGenSubsTodo[iSubsExt]
	#(iRun <- iRunsExt[1])
	
	nTwDEMCRunOld <- nTwDEMCRun
	resITwDEMC <- c( resITwDEMC, vector("list",length(uPrevIRunsExt)) )
	nTwDEMCRun <- length(resITwDEMC)
	
	for( iRunExt in  seq_along(unique(uPrevIRunsExt)) ){
		iRunPrev <- unique(uPrevIRunsExt)[iRunExt]
		iSubsRun <- iSubsExt[which( prevIRunsExt == iRunPrev )]
		nGenIRun <- max(nGenSubsTodo[iSubsRun])
		#twDEMC0 <- subChains( resITwDEMC[[iRun]], iPops= match( iSubsRun, which(iRunSubs==iRun)) ) 
		twDEMC0 <- subChains( resITwDEMC[[iRunPrev]], iPops= subInfo[ iSubsRun, "rpop"] )
		iRun <- nTwDEMCRunOld + iRunExt
		resTwDEMC <- resITwDEMC[[ iRun ]] <- twDEMC(twDEMC0
				, nGen=nGenIRun
				, nPops=length(iSubsRun)
				, upperParBounds=upperParBounds[iSubsRun]
				, lowerParBounds=lowerParBounds[iSubsRun]
				, controlTwDEMC=controlTwDEMC
				#, fLogDen=den2dCor 
				, ... 
			)
		## update subInfos
		subInfo[iSubsRun,c("run","rpop","nGen","nSamples")] <- 
			cbind(iRun, seq_along(iSubsRun), getNGen(resTwDEMC), getNSamples(resTwDEMC)	)
	}

	# do a subsampling
	# iPop=2
	subSamplePop <- function(iPop){
		iSubs  <- iSubPops[[iPop]] 
		#nSamplesSubsReq[iSubs]
		#(iSub <- iSubs[1])
		ssImpList <-  lapply( iSubs, function(iSub){
			if( nSamplesSubsReq[iSub] == 0) c() else{
				aTwDEMC <- resITwDEMC[[ subInfo[iSub,"run"] ]]
				aTwDEMCSub <- subChains( aTwDEMC, iPops=subInfo[iSub,"rpop"] )
				ss <- stackChains(aTwDEMCSub)
				#plot(as.mcmc.list(aTwDEMC), smooth=FALSE )
				#matplot( aTwDEMC$pAccept, type="l" )
				if( nSamplesSubsReq[iSub] == 1) 
					ss[ sample.int( nrow(ss),1),] 
				else {
					resFCheck <- do.call( fCheckProblems, c(list(aTwDEMCSub), argsFCheckProblems)) 
					if( resFCheck$hasProblems ){
						if( !is.null(dumpfileBasename) )
							if( dumpfileBasename == "recover") recover() else dump.frames(dumpfileBasename,TRUE)
						stop(paste("divideTwDEMC: checking function encountered problems. Dump in ",dumpfileBasename,sep=""))
					}
					sst <- ss[ round(seq(1,nrow(ss),length.out=nSamplesSubsReq[iSub])),]
				}
			}
		})
		ssImp <- do.call( rbind, ssImpList )
	}
	#mtrace(tmp.f)
	resl <- lapply( 1:nPops, subSamplePop)	# end lapply iPop: giving a single combined sample for each population
	#res <- ssImpPops <- abind(resl, rev.along=0)
	
	##value<< 
	## For each population, a list with entries
	res <- lapply(1:nPops,function(iPop){
			posPop <- which(popSubPops==iPop)
			list(
				sample=resl[[iPop]]					##<< numeric matrix: the sampled parameters (rows: cases, cols: 1: logDensity, others parameter dimensions
				,subPercChange=subPercChange[posPop]	##<< numeric vector: estimated proportion in limiting distribution to proportion of initial proportion 
				,upperParBounds=upperParBounds[posPop]	##<< list of named numeric vectors: the upper parameter bounds for the subspaces
				,lowerParBounds=lowerParBounds[posPop]	##<<  
				#,pSubs0=pSubs0[posPop]
				#,pSubs2=p2[posPop]
			)
		})
	##end<<
	
}
attr(divideTwDEMC,"ex") <- function(){
	data(den2dCorTwDEMC)
	aTwDEMC <- 	thin(den2dCorTwDEMC, start=300)
	aSample <- stackChainsPop(aTwDEMC)
	ss0 <- stackChains(aTwDEMC)
	#aSamplePure <- aSample[,-1,]
	#aSamplePurePop <- aSamplePure[,,1]
	#res <- divideTwDEMC(aSample, nGen=100, fLogDen=den2dCor, attachDetails=TRUE )
	#mtrace(divideTwDEMC)
	res <- divideTwDEMC(aSample, nGen=500, fLogDen=den2dCor )
	plot( b ~ a, as.data.frame(res[[1]]$sample), xlim=c(-0.5,2), ylim=c(-20,40) )
	#barplot(res[[1]]$subPercChange, names.arg=seq_along(res[[1]]$subPercChange) )
	#str(res[[1]]$lowerParBounds)
	ssImpPops1 <- ssImpPops <- abind( lapply( res, "[[", "sample"), rev.along=0 )
	plot(density( ssImpPops[,"a",1]));lines(density( ssImpPops[,"a",2]),col="green"); lines(density( ss0[,"a"]),col="blue")
	plot( b ~ a, as.data.frame(ss0), xlim=c(-0.5,2), ylim=c(-20,40) ); points(0.8, 0, col="red" )
	#plot( b ~ a, as.data.frame(ssImpPops[,,2]) ); points(xyMax[1], xyMax[2], col="red" )
	plot( b ~ a, as.data.frame(ssImpPops[,,2]), xlim=c(-0.5,2), ylim=c(-20,40) ); points(xyMax[1], xyMax[2], col="red" )
	plot( b ~ a, as.data.frame(ssImpPops[,,1]), xlim=c(-0.5,2), ylim=c(-20,40) ); points(xyMax[1], xyMax[2], col="red" )
	#plot(density( ssImpPops[,"b",1]));lines(density( ssImpPops[,"b",2]),col="green"); lines(density( ss0[,"b"]),col="blue")
	ssImpPops2 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops1[,,], nGen=500, fLogDen=den2dCor, attachDetails=TRUE ), "[[", "sample"), rev.along=0 )
	#mtrace(divideTwDEMC)
	#ssImpPops2 <- ssImpPops <- divideTwDEMC(ssImpPops1[,-1,], nGen=100, fLogDen=den2dCor, attachDetails=TRUE )
	ssImpPops3 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops2[,,], nGen=500, fLogDen=den2dCor ), "[[", "sample"), rev.along=0 )
	ssImpPops4 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops3[,,], nGen=500, fLogDen=den2dCor ), "[[", "sample"), rev.along=0 )
	ssImpPops5 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops4[,,], nGen=500, fLogDen=den2dCor, attachDetails=TRUE ), "[[", "sample"), rev.along=0 )
	ssImpPops6 <- ssImpPops <- abind( lapply( res <- divideTwDEMC(ssImpPops5[,,], nGen=500, fLogDen=den2dCor, attachDetails=TRUE ), "[[", "sample"), rev.along=0 )
	str(ssImpPops)
}
#mtrace(divideTwDEMC)
#twUtestF(divideTwDEMC, divertOutputFile=NULL)

setMethodS3("divideTwDEMCBatch","divideTwDEMCBatch", function(
	### append runs to result of \code{\link{divideTwDEMCBatch.default}}
	x						##<< result of former call to divideTwDEMCBatch
	, ...					##<< further arguments to \code{\link{divideTwDEMCBatch.default}}
){
	resArray <- divideTwDEMCBatch.default(x$sample,...)
	within(resArray,{
		wSubs <- rbind( x$wSubs, wSubs )
	})
})
		

setMethodS3("divideTwDEMCBatch","default", function( 
	### iteratively run batches using divideTwDEMCBatch
	x						##<< numeric array: rows:steps, col:parameters including logLik, 3: independent populations
	, ...					##<< further arguments to \code{\link{divideTwDEMC}}
	, nGen=2*512			##<< number of generations
	, nGenBatch=512			##<< number of generations within one batch
	, thinPastFac=0.2		##<< thinning the past between batches to speed up localization between 0 (no past) and 1 (keep entire past)
	, subPercChangeCrit=1.6	##<< if all subPercChange of all sub-populations are belwo this value in two last batches, may assume convergence and skip further batches 
){
	#divideTwDEMCBatch.default
	if( !is.null(thinPastFac) && thinPastFac == 1) thinPastFac <- NULL	# no need to thin
	if( !is.null(thinPastFac) && (length(thinPastFac) != 1 || !is.finite(thinPastFac) || thinPastFac<0 || thinPastFac>1) ) 
		stop("divideTwDEMCBatch: thinPastFac must a positive scalar between 0 and 1.")
	if( !is.null(thinPastFac) &&  ((nGen %% nGenBatch) != 0) )
		warning("divideTwDEMCBatch: nGenBatch is not a factor of nGen. Together with thinning the past, this may yield low sample number.")
	iBatch <- 0
	iN <- 0
	ss <- x
	nBatch <- ceiling(nGen/nGenBatch)
	nPops <- dim(x)[3]
	subPercChange <- matrix(vector("list",nBatch*nPops), nrow=nBatch )
	maxSubPercChange <- matrix( vector("numeric", nBatch*nPops ), nrow=nBatch, ncol=nPops )	# cols: populations
	for( iBatch in 1:nBatch ){
		cat(paste(iN," out of ",nGen," generations completed.     ",date(),"\n",sep=""))
		resBatch <- divideTwDEMC(  ss,nGen=min(nGenBatch, nGen-iN)	
		#,fLogDen=den2dCor 
		,...
		)
		subPercChange[iBatch,] <- subPercChangeCur <- lapply( resBatch, "[[", "subPercChange" )
		maxSubPercChange[iBatch,] <- maxSubPercChangeCur <- sapply( subPercChangeCur, max )
		ssCurrent <- abind(lapply(resBatch,"[[","sample"), rev.along=0)
		ssPast <-  if( !is.null(thinPastFac) ){
			iKeep <- round(seq(1,nrow(ss),length.out=nrow(ss)*thinPastFac))
			ss[iKeep,,]
		}else ss
		ss <- abind( ssPast, ssCurrent, along=1 )
		iN <- iN + nGenBatch
		# check if weights of sub-populations converged, then do not need further runs
		if( all(maxSubPercChange[iBatch-(0:1),] <= subPercChangeCrit) ) {
			cat("divideTwDEMCBatch: weights of sub-populations converged. Skip further generations.")
			cat("\n")
			break
		}
	}
	##value<< A list with components 
	res <- list(
		sample = ss					##<< numeric array same as argument \code{x} but thinned and new sample appended
		,subPercChange = subPercChange[1:iBatch,]	##<< matrix of numeric vectors giving the weights of sub-populations
		,maxSubPercChange = maxSubPercChange[1:iBatch,]	##<< numeric vector: calculated ratio of subspace weights, see argument \code{wSubFacMax}
		,resDivideTwDEMC = resBatch	##<< result of last call to \code{\link{divideTwDEMC}}
	)
	##end<<
	class( res ) <- c("divideTwDEMCBatch", class(res) )
	res
})
attr(divideTwDEMCBatch.default,"ex") <- function(){
	data(den2dCorTwDEMC)
	aTwDEMC <- 	thin(den2dCorTwDEMC, start=300)
	aSample <- stackChainsPop(aTwDEMC)
	#mtrace(divideTwDEMCBatch.default)
	ssRes1 <- ssRes <- divideTwDEMCBatch(aSample, nGen=512*6, fLogDen=den2dCor )
	#mtrace(divideTwDEMCBatch.divideTwDEMCBatch)
	#ssRes <- divideTwDEMCBatch(ssRes1, nGen=512*1, fLogDen=den2dCor )	#calling divideTwDEMCBatch.divideTwDEMCBatch
	#mtrace(getSubSpaces)
	#tmp <- getSubSpaces(ssRes1$sample[,-1,1])
	ssRes2 <- ssRes <- divideTwDEMCBatch(ssRes1, nGen=512*6, fLogDen=den2dCor )	#calling divideTwDEMCBatch.divideTwDEMCBatch
	
	#tmp <- divideTwDEMC( ssRes1$sample[,,2,drop=FALSE], nGen=512, fLogDen=den2dCor)
	#ssImpPops <- abind( tmp[[1]]$sample, tmp[[1]]$sample, rev.along=0) 
	
	#wSubsB <- ssRes$wSubs[,1]
	#wSubsB <- ssRes$wSubs[,2]
	wSubsB <- ssRes$wSubs
	sapply( lapply(wSubsB, quantile, probs = c(1, 0.25) ), function(entry){ entry[1] / entry[2] }) 
	
	ssImpPops <- ssRes$sample
	#ssImpPops <- ssRes1$sample
	plot(density( ssImpPops[,"a",1]));lines(density( ssImpPops[,"a",2]),col="green"); lines(density( ss0[,"a"]),col="blue")
	#ssImpPops <- ssRes1$sample
	plot( b ~ a, as.data.frame(ssImpPops[,,1]), xlim=c(-0.5,2), ylim=c(-20,40) ); points(0.8,0, col="red" )
	plot( b ~ a, as.data.frame(ssImpPops[,,2]), xlim=c(-0.5,2), ylim=c(-20,40) ); points(0.8,0, col="red" )
	plot( b ~ a, as.data.frame(ssImpPops[,,1])); points(0.8,0, col="red" )
	plot( b ~ a, as.data.frame(ssImpPops[,,2])); points(0.8,0, col="red" )
	str(ssImpPops)
}
#twUtestF(divideTwDEMCBatch)
