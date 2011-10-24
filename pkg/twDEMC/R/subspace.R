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
	ss1 <- ss <- pps[,-1]
	#mtrace(findSplit)
	(res <- findSplit(ss1))
		
	ss2 <- ss <- ss1[ ss1[,res$iVar] > res$splitValue,]
	res2 <- findSplit(ss2)
	ss3 <- ss <- ss2[ ss2[,res2$iVar] > res2$splitValue,]
	res3 <- findSplit(ss3)
	ss4 <- ss <- ss3[ ss3[,res3$iVar] > res3$splitValue,]
	res4 <- findSplit(ss4)
	
	#mtrace(findSplit)
	(res <- findSplit(ss1, iVars=c(2,1) ))
	(res <- findSplit(ss1, iVars=c(2) ))
	plot( ss[,1], ss[,2], col=c("blue","red")[as.integer(ss[,names(res)] >= res)+1])
}

getSubSpaces <- function( 
	### Splitting the sample into subsamples and recording the upper and lower bounds
	aSample 			##<< numeric matrix with parameters in columns: sample or subsample that is devided
	,nSplit=4			##<< see \code{\link{findSplit}}, potentially modified due to argument \code{minPSub}
	,argsFSplit=list()	##<< further arguments to \code{\link{findSplit}}
	,pSub=1				##<< the fraction that a subSample comprises within a bigger sample 
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
	aSample <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
	#mtrace(getSubSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.4)
	str(subSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.05)
	(tmp <- sapply( subSpaces, function(subSpace){nrow(subSpace$sample)})/nrow(aSample))
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
	, nChainPar=13			##<< number of chains to run in parallel, good choice is 2*nCpu
	, minNSamplesSub=32		##<< minimum number of records in a subspace sample, increase to avoid wrong estimation of weights
	, maxFacSubWeight		##<< maximum ratio of the weight of a subspace per average weight of subspaces XXTodo
){
	#samplePop <- aSample[,,1]
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
			getSubSpaces(samplePop, minPSub=minPSub)
		})
	nSubPops <- sapply(subSpaces, length)	# number of subs per populations
	subSpacesFlat <- do.call( c, subSpaces )	# put all subPopulations of all populations on same level
	cumNSubPops <- cumsum(nSubPops)
	iSubPops <- lapply(1:nPops, function(iPop){ cumNSubPops[iPop]+1-(nSubPops[iPop]:1)})	# index of subs in flat version for each pop 			
	ZinitSubs <- abind( lapply( subSpacesFlat, function(ssEntry){
						s0 <- ssEntry$sample[ round(seq(1,nrow(ssEntry$sample),length.out=nInit)), ]
						Zinit <- array( t(s0), dim=c(ncol(s0),m0,nChainsPop), dimnames=list(parms=colnames(s0),steps=NULL,chains=NULL) )
					}), along=3)
	upperParBounds <- lapply( subSpacesFlat, function(ssEntry){ ssEntry$upperParBounds})
	lowerParBounds <- lapply( subSpacesFlat, function(ssEntry){ ssEntry$lowerParBounds})
	pSubs0L <- lapply( subSpaces, sapply, "[[","pSub" ) # percentile of subspace of the overall space
	pSubs0 <- do.call(c,pSubs0L)
	nGen0 <- pmax(m0*thin, ceiling(minNSamplesSub*thin/nChainsPop), nGen*pSubs0)		# at minimum m0*thin generations to keep sufficient samples for extending the run
	nGen0Thin <- ceiling(nGen0/thin)*thin
	
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
	
	#iRun <- 1
	resITwDEMC <- lapply( 1:nTwDEMCRun, function(iRun){
			iSubsRun <- which( iRunSubs == iRun )
			nGenIRun <- max( nGen0Thin[ iSubsRun[1] ] )
			#resTwDEMC <- 
			resTwDEMC <- twDEMC(
				ZinitSubs[,,as.vector(sapply( (iSubsRun-1)*nChainsPop, "+", (1:nChainsPop)))]
				, nPops=length(iSubsRun)
				, upperParBounds=upperParBounds[iSubsRun]
				, lowerParBounds=lowerParBounds[iSubsRun]
				, controlTwDEMC=controlTwDEMC
				, nGen=nGenIRun
				, fLogDen=den2dCor 
			)
		})
	nGenIRuns <- sapply( resITwDEMC, function(r){getNGen(r)})
	nSamplesIRuns <- sapply( resITwDEMC, function(r){getNSamples(r)})
	nGenSubs <- sapply( 1:nSubs, function(iSub){ nGenIRuns[ iRunSubs[iSub] ]})
	nSamplesSubs <- sapply( 1:nSubs, function(iSub){ nSamplesIRuns[ iRunSubs[iSub] ]})
	pSubs1 <- nGenSubs/nGen0		# sampling weights (>1 obtained more samples than required by quantile)
	

	# get the samples of all subPopulations
	ssSub <- lapply( 1:nSubs, function(iSub){
			iRun <- iRunSubs[iSub]
			stackChains( subChains(resITwDEMC[[iRun]], iPops=match(iSub, which(iRunSubs==iRun)) ))
		}) 
	popSubPops <- do.call( c, lapply( 1:nPops, function(iPop){ rep(iPop,nSubPops[iPop])} ))		# each index gives the population that the subPopulatios is part of
	#iPop <- 1
	#iPop <- 2
	##details<< 
	## The samples of the subspaces are combined again after thinning the subsamples 
	## inversly proportional
	## to their contribution to the sum of unnormalized densities.
	wSubsL <- lapply( 1:nPops, function(iPop){
			iSubs <- iSubPops[[iPop]]
			lw <- sapply(iSubs,function(iSub){ twLogSumExp(ssSub[[iSub]][,1])-log(nSamplesSubs[iSub]) })		# average the unnormalized densities, providing log of the weights
			lSumW <- twLogSumExp(lw) 
			wSubs <- exp( lw-lSumW ) 			#w/sumW
		})
	wSubs <- do.call(c, wSubsL )
	
	
	# extend chains of those populations that do not have enough samples yet
	p2 <- do.call(c, lapply( 1:nPops, function(iPop){
			p2u <- wSubsL[[iPop]]*pSubs0L[[iPop]]		# two weights: quantile and density importance
			p2 <- p2u/sum(p2u)							# normalize to 1 within population 
		}))
	nSamplesSubsReq <- round(nGen/thin * p2)
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
	nGenSubsTodo <- nSamplesSubsReq*thin - nGenSubs
	iSubsExt <- which(nGenSubsTodo > 0)
	iRunsExt <- iRunSubs[iSubsExt]
	#nGenSubsTodo[iSubsExt]
	#(iRun <- iRunsExt[1])
	for( iRun in  unique(iRunsExt)){
		iSubsRun <- iSubsExt[which( iRunsExt == iRun )]
		nGenIRun <- max(nGenSubsTodo[iSubsRun])
		twDEMC0 <- subChains( resITwDEMC[[iRun]], iPops= match( iSubsRun, which(iRunSubs==iRun)) ) 
		resTwDEMC <- twDEMC(twDEMC0
				, nGen=nGenIRun
				, nPops=length(iSubsRun)
				, upperParBounds=upperParBounds[iSubsRun]
				, lowerParBounds=lowerParBounds[iSubsRun]
				, controlTwDEMC=controlTwDEMC
				, fLogDen=den2dCor 
			)
		## update the sample 
		#iiSub <- 1
		for( iiSub in seq_along(iSubsRun)){
			iSub <- iSubsRun[iiSub]			
			ssSub[[iSub]] <- stackChains( subChains(resTwDEMC, iPops=iiSub) )
		}
	}

	# do a subsampling
	# iPop=2
	tmp.f <- function(iPop){
		iSubs  <- iSubPops[[iPop]] 
		#nSamplesSubsReq[iSubs]
		#(iSub <- iSubs[8])
		ssImpList <-  lapply( iSubs, function(iSub){
				ss <- ssSub[[iSub]]
				if( nSamplesSubsReq[iSub] == 0) c() else
				if( nSamplesSubsReq[iSub] == 1) ss[ sample.int( nrow(ss),1),] else
					ss[ round(seq(1,nrow(ss),length.out=nSamplesSubsReq[iSub])),]
			})
		ssImp <- do.call( rbind, ssImpList )
	}
	#mtrace(tmp.f)
	resl <- lapply( 1:nPops, tmp.f)	# end lapply iPop: giving a single combined sample for each population
	ssImpPops <- abind(resl, rev.along=0)
	
	if( attachDetails ){
		attr(ssImpPops,"subSpaces") <- tmp <- lapply(1:nPops,function(iPop){
				posPop <- which(popSubPops==iPop)
				list(
					upperParBounds=upperParBounds[posPop]
					,lowerParBounds=lowerParBounds[posPop]
					,pSubs0=pSubs0[posPop]
					,wSubs=wSubs[posPop]
					,pSubs2=p2[posPop]
				)
			})
	}
	ssImpPops
}
attr(divideTwDEMC,"ex") <- function(){
	aTwDEMC <- 	thin(den2dCorTwDEMC, start=300)
	aSample <- stackChainsPop(aTwDEMC)[,-1,]
	#ssImpPops1 <- ssImpPops <- divideTwDEMC(aSample, nGen=100, fLogDen=den2dCor, attachDetails=TRUE )
	ssImpPops1 <- ssImpPops <- divideTwDEMC(aSample, nGen=500, fLogDen=den2dCor )
	plot( b ~ a, as.data.frame(ss0), xlim=c(-0.5,2), ylim=c(-20,40) ); points(xyMax[1], xyMax[2], col="red" )
	#plot( b ~ a, as.data.frame(ssImpPops[,,2]) ); points(xyMax[1], xyMax[2], col="red" )
	plot( b ~ a, as.data.frame(ssImpPops[,,2]), xlim=c(-0.5,2), ylim=c(-20,40) ); points(xyMax[1], xyMax[2], col="red" )
	plot( b ~ a, as.data.frame(ssImpPops[,,1]), xlim=c(-0.5,2), ylim=c(-20,40) ); points(xyMax[1], xyMax[2], col="red" )
	plot( b ~ a, as.data.frame(ssImpPops[,,2]), xlim=c(-2,3), ylim=c(-80,80) ); points(xyMax[1], xyMax[2], col="red" )
	plot( b ~ a, as.data.frame(ssImpPops[,,1]), xlim=c(-2,3), ylim=c(-80,80) ); points(xyMax[1], xyMax[2], col="red" )
	plot(density( ssImpPops[,"a",1]));lines(density( ssImpPops[,"a",2]),col="green"); lines(density( ss0[,"a"]),col="blue")
	#plot(density( ssImpPops[,"b",1]));lines(density( ssImpPops[,"b",2]),col="green"); lines(density( ss0[,"b"]),col="blue")
	ssImpPops2 <- ssImpPops <- divideTwDEMC(ssImpPops1[,-1,], nGen=500, fLogDen=den2dCor, attachDetails=TRUE )
	#mtrace(divideTwDEMC)
	#ssImpPops2 <- ssImpPops <- divideTwDEMC(ssImpPops1[,-1,], nGen=100, fLogDen=den2dCor, attachDetails=TRUE )
	ssImpPops3 <- ssImpPops <- divideTwDEMC(ssImpPops2[,-1,], nGen=500, fLogDen=den2dCor )
	ssImpPops4 <- ssImpPops <- divideTwDEMC(ssImpPops3[,-1,], nGen=500, fLogDen=den2dCor )
	ssImpPops5 <- ssImpPops <- divideTwDEMC(ssImpPops4[,-1,], nGen=500, fLogDen=den2dCor, attachDetails=TRUE )
	ssImpPops6 <- ssImpPops <- divideTwDEMC(ssImpPops5[,-1,], nGen=500, fLogDen=den2dCor, attachDetails=TRUE )
	str(ssImpPops)
}