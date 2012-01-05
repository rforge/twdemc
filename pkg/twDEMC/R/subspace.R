findSplit <- function(
	### determine the parameter and its value, where to split the parameter space
	aSample	##<< the sample (stacked chains)
	,nSplit = 4 ##<< how many points per parameter dimension to check for split
	,rVarCrit = 3^2	##<< minimal ratio of variances of a variable between right and left side of a splitting point
	,rAlphaSlopeCrit = base:::pi/4	##<< minimal ratio of angle between scaled slopes, defaults to a third of a half cirle
	,pSlopeMax=0.05					##<< confidence level for the slope, above which no difference in slopes is calculated 
	,iVars = 1:ncol(aSample)		##<< integer vector: index or parameter dimensions, i.e. variables, to check for split		
	,jVarsVar = 1:ncol(aSample)		##<< integer vector: index or parameter dimensions, i.e. variables, to check for different scales of variance
	,jVarsSlope = 1:ncol(aSample)	##<< integer vector: index or parameter dimensions, i.e. variables, to check for differences in angle of normalized slopes
	,isBreakEarly = TRUE			##<< if TRUE and argument \code{checkSlopesFirst} is given, then check slope angles for variables given in checkSlopesFirst first and if break is found, do not evaluate all the other slope angles
	,checkSlopesFirst = data.frame(ivar=1,j1AlphaSlope=1,jAlphaSlope=2)[FALSE,]	##<< data.frame with entries ivar, j1AlphaSlope and j2AlphaSlope as returned in entry resD in result of ressplit
	,maxNSample=128					##<< if given a value, then aSample is thinned before to given number of records (for efficiently calculating variances)
	,calcVarParallel=TRUE			##<< having less than 30 parameters it is often faster not to parallelize variance calculation 
	,debugSequential=FALSE			##<< by default calcualation is distributed to sup-processes for parralel execution, set to TRUE for non-distributed execution
){
	##details<< 
	## First it checks for different scales of variance in other variables.
	# i and j are indexes within iVars and jVarsVar which are indexes in colnames(aSample)
	cNames <- colnames(aSample)
	if( is.character(iVars) ) iVars <- as.numeric(sapply(iVars, match, cNames))
	if( is.character(jVarsVar) ) jVarsVar <- as.numeric(sapply(jVarsVar, match, cNames))
	if( is.character(jVarsSlope) ) jVarsSlope <- as.numeric(sapply(jVarsSlope, match, cNames))
	jVarsVar <- union(jVarsSlope, jVarsVar)	# make sure to check variances for all slope variables, jVarsSlope must be first so that indexing j works in parVarRight and parVarLeft
	if( is.finite(maxNSample) && (maxNSample < nrow(aSample)) ) aSample <- aSample[ round(seq(1,nrow(aSample),length.out=maxNSample)),]
	foundSplit <- FALSE
	percentiles <- 1:(nSplit)/(nSplit+1)
	# value of the variable (var x split)
	quantVar <- iVarLeft <- iVarRight <- matrix( NA_real_, nrow=length(iVars), ncol=nSplit )
	# variance of the left subset and variance of the right subset (splitVar, varVar, split)
	parVarLeft <- parVarRight <- array( NA_real_, dim=c(length(iVars),length(jVarsVar),nSplit)
		, dimnames=list(iVar=cNames[iVars], jVar=cNames[jVarsVar], iSplit=NULL ))
	# slope of the left subset and variance of the right subset (splitVar, aVar1, aVar2, split)
	# results of pVar and slope for best split per variable combination (varj x varj x (val,pVar,perc))
	pVars <- array( NA_real_, dim=c(length(iVars),length(jVarsVar),3), dimnames=list(iVar=cNames[iVars],jVar=cNames[jVarsVar],c("varVal","perc","val")) )
	dAlphaSlopes <- array( NA_real_, dim=c(length(iVars),length(jVarsSlope),length(jVarsSlope),3), dimnames=list(iVar=cNames[iVars],jVar1=cNames[jVarsSlope],jVar2=cNames[jVarsSlope],c("varVal","perc","val")) )
	# left and right samples var -> split -> sample
	ssSubLeft <- ssSubRight <- vector(mode="list",length=length(iVars) )
	# result details
	resD <- data.frame(iVar=iVars				##<< index of the splitting variable
		, jPVar=NA_integer_						##<< index of the variable with maxium difference of variance when splitting on variable iVar 
		, pVar=NA_real_							##<< max( abs(proportion of the variances))
		, j1AlphaSlope=NA_integer_				##<< index of the first variable with maxium difference in slopes
		, j2AlphaSlope=NA_integer_				##<< index of the second variable with maxium difference in slopes
		, dAlphaSlope=NA_real_					##<< maximum diffrence in slope angle
	)
	
	# calculate all quantiles and subsamples
	#i <- 1
	for( i in seq_along(iVars)){
		iVar <- iVars[i]
		p <- aSample[,iVar] 
		qp <- quantVar[i, ] <- quantile( p, percentiles )
		#iSub <- 1
		#iSub <- nSub
		ssSubLeft[[i]] <- ssiSubLeft <- lapply( 1:(nSplit), function(iSplit){
				ssSub <- aSample[ (p <= qp[iSplit]),]
			})
		ssSubRight[[i]] <- ssiSubRight <- lapply( 1:(nSplit), function(iSplit){
				ssSub <- aSample[ (p >= qp[iSplit]),]
			})
	}
	
	# invoke calculating variances in parallel manner
	# ?sfFArgsApplyLB
	F_ARGS <- function(i){ list(
			iVar=iVars[i]
			,qp=quantVar[i,]
			,ssiSubLeft=ssSubLeft[[i]]
			,ssiSubRight=ssSubRight[[i]]
		)}
	resVarIL <- sfFArgsApplyLB( length(iVars), F_ARGS, .fPVarI		,jVarsVar=jVarsVar, nSplit=nSplit, cNames=cNames, percentiles=percentiles		,debugSequential=(!calcVarParallel || debugSequential)	)
	#resVarI <- .fPVarI(iVar=iVar,p=p,qp=qp,ssiSubLeft=ssiSubLeft,ssiSubRight=ssiSubRight, jVarsVar=jVarsVar, nSplit=nSplit, cNames=cNames, percentiles=percentiles)
	
	# reformat the results
	# i=2
	for( i in seq_along(iVars)){
		resVarI <- resVarIL[[i]]
		parVarLeft[i,,] <- resVarI$parVarLeft
		parVarRight[i,,] <- resVarI$parVarRight
		pVars[i,,] <- resVarI$pVars
		#pVars[i,(colnames(pVars) == cNames[i]),"val"] <- NA  # variance ratio for splitting variable does not make sense and confuses selection of highest pVar	
		
		jmax <- which.max( pVars[i,,"val"])		
		resD[i,"jPVar"] <- jVarsVar[jmax]
		resD[i,"pVar"] <- pVars[i,jmax,"val"]
	} # for iVar
	
	#maxPVar <- suppressWarnings( apply(pVars[,,"val"],1, max , na.rm=TRUE ) )
	#sortI <- order(maxPVar, decreasing = TRUE)		
	#sortJ <- order(pVars[sortI[1],,"val"], decreasing = TRUE)
	oi <- order(resD$pVar, decreasing=TRUE)
	imax <- oi[1]
	if( resD$pVar[imax] > rVarCrit ){
		varName=cNames[ iVars[imax] ]
		pVarsMax <- pVars[imax, match( resD$jPVar[imax],jVarsVar),]
		res <- list( 
			split=structure( pVarsMax["varVal"], names=varName)
			, varName=varName
			, perc=as.numeric(pVarsMax["perc"]) 
			#, iVars=resD$iVar[oi] 
			#, jVarsVar=resD$jPVar[oi]
			, resD=resD
		)
		return(res)
	}	
	
	##details<< 
	## Next it checks for different angles of the normalized slopes in relation of the parameters
	F_ARGS <- function(i){ list(
			qp=quantVar[i,]
			,ssiSubLeft=ssSubLeft[[i]]
			,ssiSubRight=ssSubRight[[i]]
			,pariVarLeft=adrop(parVarLeft[i,,,drop=FALSE],1)	# only drop the i dimension, avoid dropping only one jVar or only one iSplit
			,pariVarRight=adrop(parVarRight[i,,,drop=FALSE],1)
		)}

	if( isBreakEarly && is.data.frame(checkSlopesFirst) ){
		checkSlopesFirst <- checkSlopesFirst[ 
			(checkSlopesFirst$iVar %in% iVars) &  
			is.finite(checkSlopesFirst$j1AlphaSlope) & 
			is.finite(checkSlopesFirst$j2AlphaSlope)
			,]	# do not check the ones that are not in iVars or have non-finite jAlphaSlope
		ind <- match( checkSlopesFirst$iVar, iVars )
		if( (nRowResD <- nrow(checkSlopesFirst)) >0 ){
			# first check the correlations given in checkSlopesFirst
			# mtrace(.fCheckAlphaSlopeSplit)
			#iResD=3
			F_ARGS3 <- function(iResD){ c(
					list( j1=match( checkSlopesFirst$j1AlphaSlope[iResD], jVarsSlope)
						,j2=match( checkSlopesFirst$j2AlphaSlope[iResD], jVarsSlope) )
					,F_ARGS( ind[iResD] )
				)}
			resCheckL <- sfFArgsApplyLB( nRowResD, F_ARGS3, .fCheckAlphaSlopeSplit
					, jVarsSlope=jVarsSlope, nSplit=nSplit, cNames=cNames, pSlopeMax=pSlopeMax, percentiles=percentiles 
					, debugSequential=debugSequential
				)
			resCheck <- do.call(rbind, resCheckL)
			if( !is.null(resCheck) ){
				resD[ match(checkSlopesFirst$iVar, resD$iVar), c("j1AlphaSlope","j2AlphaSlope","dAlphaSlope")] <- 
					cbind(checkSlopesFirst[,c("j1AlphaSlope","j2AlphaSlope")], resCheck[,"val"] )
				# search the maximum and if found critical slope, return early
				oiAlpha <- order(resD$dAlphaSlope, decreasing=TRUE)
				imax <- oiAlpha[1]
				if( resD$dAlphaSlope[imax] > rAlphaSlopeCrit ){
					varName=cNames[ iVars[imax] ]
					resCheckMax <- resCheck[imax, ] 
					res <- list( 
						split=structure( resCheckMax["varVal"], names=varName)
						, varName=varName
						, perc=as.numeric(resCheckMax["perc"]) 
						#, iVars=resD$iVar[oi] 
						#, jVarsVar=resD$jPVar[oi]
						, resD=resD
					)
					return(res)
				}
			}# resCheck not null
		}# nrow(checkSlopesFirst)>0
	} # if breakearly
	#i<-1
	if( length(jVarsSlope) >= 2 ){
		#mtrace(.fDAlphaSlopeI)
		resDAlphaSlopesIL <- sfFArgsApplyLB( length(iVars), F_ARGS, .fDAlphaSlopeI
			, jVarsSlope=jVarsSlope, nSplit=nSplit, cNames=cNames, pSlopeMax=pSlopeMax, percentiles=percentiles
			, fCheckAlphaSlopeSplit=.fCheckAlphaSlopeSplit	# must be transported to foreign process
			,debugSequential=debugSequential
		)
		#dAlphaSlopes[i,,,] <- .fDAlphaSlopeI(qp=qp, ssiSubLeft=ssiSubLeft, ssiSubRight=ssiSubRight, pariVarLeft=pariVarLeft, pariVarRight=pariVarRight, jVarsSlope=jVarsSlope, nSplit=nSplit, cNames=cNames, pSlopeMax=pSlopeMax, percentiles=percentiles)
		
		for( i in seq_along(iVars)){
			iVar <- iVars[i]
			dAlphaSlopes[i,,,] <- resDAlphaSlopesIL[[i]]
			x <- dAlphaSlopes[i,,,"val"]
			j12max <- arrayInd( which.max(x), dim(x) )
			resD[i,c("j1AlphaSlope","j2AlphaSlope")] <- jVarsSlope[ j12max ]  
			resD[i,"dAlphaSlope"] <- x[j12max[1],j12max[2] ]
			dAlphaSlopesMax <- 	dAlphaSlopes[i, match(resD$j1AlphaSlope[i],jVarsSlope), match(resD$j2AlphaSlope[i],jVarsSlope),] 
		} # for iVar
		oiAlpha <- order(resD$dAlphaSlope, decreasing=TRUE)
		imax <- oiAlpha[1]
		if( resD$dAlphaSlope[imax] > rAlphaSlopeCrit ){
			varName=cNames[ iVars[imax] ]
			dAlphaSlopeMax <- dAlphaSlopes[imax, match(resD$j1AlphaSlope[imax],jVarsSlope), match(resD$j2AlphaSlope[imax],jVarsSlope),] 
			res <- list( 
				split=structure( dAlphaSlopeMax["varVal"], names=varName)
				, varName=varName
				, perc=as.numeric(dAlphaSlopeMax["perc"]) 
				#, iVars=resD$iVar[oi] 
				#, jVarsVar=resD$jPVar[oi]
				, resD=resD
			)
			return(res)
		}
	}
	
	##value<< list with components  
	return( list(
			split=NA_real_	##<< named scalar: splitting value with the name of the parameter dimension that is to split
			,varName=NA_character_ ##<< name of the splitting variable
			,perc=NA_real_	##<< scalar: percentile of the splitting point
			#, iVars=resD$iVar[oi] ##<< integer vector: argument \code{iVars} ordered by maximum proportion of variances or difference in slopes 
			#, jVarsVar=resD$jPVar[oi]##<< integer vector: argument \code{jVarsVar} ordered by proportion of varainces or difference in slopes for iVars[1]
			, resD=resD		##<< dataframe of result details giving for each split variable the indices of variable with maximum proportion of variances and and the variable with maximum angle in correlation  
		))
	##end<< 
}
attr(findSplit,"ex") <- function(){
	#twUtestF(findSplit) # there are unit tests for this function
	data(den2dCorEx)
	ss1 <- ss <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
	#mtrace(findSplit)
	(res <- res0 <- findSplit(ss1))	# returns before checking slope angles
	(res <- res1 <- findSplit(ss1, rVarCrit=Inf))	# find slope angles
	
	# successively find splits in subsets
	ss2 <- ss <- ss1[ ss1[,res$varName] > res$split,]
	res2 <- findSplit(ss2)
	ss3 <- ss <- ss2[ ss2[,res2$varName] > res2$split,]
	res3 <- findSplit(ss3)
	
	# try different orders of the variables
	(res <- findSplit(ss1, iVars=c("b") )) # NA indicates: no further split found
	
	#visualize the split
	plot( ss1[,"a"], ss1[,"b"], col=c("blue","red")[as.integer(ss1[,res0$varName] >= res0$split)+1])
}


.fPVarI <- function(
	### remote part of calculating proportion of variance between two subsamples
	iVar,qp,ssiSubLeft,ssiSubRight, jVarsVar, nSplit, cNames, percentiles
){
	# variance of the left subset and variance of the right subset (splitVar, varVar, split)
	parVarLeft <- parVarRight <- array( NA_real_, dim=c(length(jVarsVar),nSplit)
		, dimnames=list(jVar=cNames[jVarsVar], iSplit=NULL ))
	# results of pVar  for best split per variable combination (varj x (val,pVar,perc))
	pVars <- array( NA_real_, dim=c(length(jVarsVar),3), dimnames=list(jVar=cNames[jVarsVar],c("varVal","perc","val")) )
	# ratio of variances of subsets left and right of splitting points
	#j <- 1
	for( j in seq_along(jVarsVar)){
		jVar <- jVarsVar[j]
		#iSplit=1
		pVar <- rep( NA_real_, nSplit)
		for(iSplit in 1:nSplit){
			varLeft <- parVarLeft[j,iSplit] <- var(ssiSubLeft[[iSplit]][,jVar])
			varRight <- parVarRight[j,iSplit] <- var(ssiSubRight[[iSplit]][,jVar]) 
			#plot( ssSubLeft[[i]][[iSplit]][,iVar], ssSubLeft[[i]][[iSplit]][,jVar] )
			#plot( ssSubRight[[i]][[iSplit]][,iVar], ssSubRight[[i]][[iSplit]][,jVar] )
			pVar[iSplit] <- varLeft /varRight
			#pVars[i,j,] <- c( val=percentiles[iSplit], pVar=pVar[iSplit])
		}
		pVar[ pVar < 1 ] <- 1/pVar[pVar<1]
		iSplit <- which.max(pVar)	# get the best split
		pVars[j,c("varVal","val","perc")] <- c( varVal=qp[[iSplit]], val=pVar[[iSplit]], perc=percentiles[[iSplit]] )	# store value and variance ratio
		if( jVar==iVar) pVars[j,"val"] <- NA  # variance ratio of the splitting variable does ot make sense and confuses selection of maximum variance
		# but still calculated variances in parVarLeft and parVarRight are used 
	} # end for jVar variance
	##value<< 
	list(
		parVarLeft = parVarLeft			##<< numeric matrix( j x iSplit ): variance of the left side of the split of variable jvar 
		,parVarRight = parVarRight		##<< numeric matrix( j x iSplit ): variance of the right side of the split of variable jvar
		,pVars=pVars					##<< numeric matrix( j x c("varVal","val","perc") ): split value, variance ratio and subspace percentile
		)
	##end<<
} # f1Var

.tmp.f <- function(){
	# useful commands inside 
	tmpj<-jVar2
	jVar2<-jVar1
	jVar1<-tmpj
	tmpj<-j1
	j1<-j2
	j2<-tmpj
	plot( ssiSubLeft[[iSplit]][,jVar1], ssiSubLeft[[iSplit]][,jVar2] )
	abline(lmLeft)
	points( ssiSubRight[[iSplit]][,jVar1], ssiSubRight[[iSplit]][,jVar2], col="blue" )
	abline(lmRight, col="blue")
	plot( ssiSubRight[[iSplit]][,jVar1], ssiSubRight[[iSplit]][,jVar2] )
	abline(lmRight)
}

.fCheckAlphaSlopeSplit <- function(
	### calculating angle of slopes between two subsamples for one variable combination
	j1,j2,qp, ssiSubLeft, ssiSubRight, pariVarLeft, pariVarRight, jVarsSlope, nSplit, cNames, pSlopeMax, percentiles
){
	#iSplit=4
	jVar1 <- jVarsSlope[j1]
	jVar2 <- jVarsSlope[j2]
	tmpf <- function(iSplit){
		# for debugging commands see .tmp.f above
		lmLeft <- lm( as.formula(paste(cNames[c(jVar2,jVar1)],collapse=" ~ ")), data=as.data.frame(ssiSubLeft[[iSplit]][,c(jVar2,jVar1)]) )
		lmRight <- lm( as.formula(paste(cNames[c(jVar2,jVar1)],collapse=" ~ ")), data=as.data.frame(ssiSubRight[[iSplit]][,c(jVar2,jVar1)]) )
		if( (anova(lmLeft)[["Pr(>F)"]][1] < pSlopeMax) || (anova(lmRight)[["Pr(>F)"]][1] < pSlopeMax) ){ 
			slopeLeft <- coef(lmLeft)[2] * sqrt(pariVarLeft[j1,iSplit]/pariVarLeft[j2,iSplit]) 
			slopeRight <- coef(lmRight)[2] * sqrt(pariVarRight[j1,iSplit]/pariVarRight[j2,iSplit])
			alphaRight <- asin( slopeLeft)
			alphaLeft <- asin( slopeRight )
			abs( alphaLeft - alphaRight )
		}else{ 0 }  # both slopes are non-significant
	}
	#mtrace(tmpf)
	dAlphaSlope <- sapply( 1:nSplit, tmpf )
	iSplit <- which.max(dAlphaSlope)
	#plot( ssSubLeft[[i]][[iSplit]][,iVar], ssSubLeft[[i]][[iSplit]][,jVar] ); abline(lmLeft)
	#plot( ssSubRight[[i]][[iSplit]][,iVar], ssSubRight[[i]][[iSplit]][,jVar] ); abline(lmRight)
	c( varVal=qp[[iSplit]], val=dAlphaSlope[[iSplit]], perc=percentiles[[iSplit]] )	# store split value and differenc ein slope angles					
}

.fDAlphaSlopeI <- function(
	### remote part of calculating angle of slopes between two subsamples for all variable combinations 
	qp, ssiSubLeft, ssiSubRight, pariVarLeft, pariVarRight, jVarsSlope, nSplit, cNames, pSlopeMax, percentiles, fCheckAlphaSlopeSplit=.fCheckAlphaSlopeSplit
){
	diAlphaSlopes <- array( NA_real_, dim=c(length(jVarsSlope),length(jVarsSlope),3), dimnames=list(jVar1=cNames[jVarsSlope],jVar2=cNames[jVarsSlope],c("varVal","perc","val")) )
	for( j1 in seq_along(jVarsSlope)){
		for( j2 in j1:length(jVarsSlope) ){	# do not repeat regressions so start from j1
			if( j1!=j2){
				diAlphaSlopes[j1,j2,c("varVal","val","perc")] <- fCheckAlphaSlopeSplit(j1,j2
					,qp=qp, ssiSubLeft=ssiSubLeft, ssiSubRight=ssiSubRight, pariVarLeft=pariVarLeft, pariVarRight=pariVarRight, jVarsSlope=jVarsSlope, nSplit=nSplit, cNames=cNames, pSlopeMax=pSlopeMax, percentiles=percentiles)				
			} # end if( jVar1!=jVar2)
		} #end for jVar2
	} # end for jVar1
	### numeric array (j1 x j2 x iSplit ) of difference in angle of normalized slopes
	diAlphaSlopes	
} 




.parBoundsEnvelope <- function(
	### get the parameter bounds that encompasses all given parBounds 
	popsParBounds	##<< list of populations with entries upperParBounds and upperParBounds both named numeric vectors 
){
	#stop(".parBoundsEnvelope: not implemented yet.")
	#parName <- colnames(pop1$parms)[1]
	ub <- lb <- list()
	#pop <- popsParBounds[[2]]
	pnames <- unique( do.call(c, c( 
				lapply(popsParBounds, function(pop){ names(pop$upperParBounds) }) 
				,lapply(popsParBounds, function(pop){ names(pop$lowerParBounds) }) 
			)))
	#parName <- pnames[1]
	for( parName in pnames ){
		ubPop <- lapply( popsParBounds, function(pop){ pop$upperParBounds[parName] })
		lbPop <- lapply( popsParBounds, function(pop){ pop$lowerParBounds[parName] })
		if( all(sapply(ubPop, function(x){ !is.null(x) && is.finite(x)})) ) 
			ub[parName] <- max(unlist(ubPop)) 
		if( all(sapply(lbPop, function(x){ !is.null(x) && is.finite(x)})) ) 
			lb[parName] <- min(unlist(lbPop)) 
	}
	##value<< list with entries
	list( ##describe<<
		upperParBounds = unlist(ub)		##<< named numeric vector of upper parameter bounds
		,lowerParBounds = unlist(lb)	##<< named numeric vector of lower parameter bounds
	) ##end<<
}
attr(.parBoundsEnvelope,"ex") <- function(){
	data(den2dCorEx)
	mc0 <- den2dCorEx$mcSubspaces0
	#mtrace(.parBoundsEnvelope)
	.parBoundsEnvelope( mc0$pops[1:4] )
	.parBoundsEnvelope( mc0$pops[3:4] )
}







getSubSpaces <- function( 
	### Recursively splitting the sample into subsamples and recording the upper and lower bounds
	aSample 			##<< numeric matrix with parameters in columns: sample or subsample that is devided (do not include rLogLik in first column)
	,nSplit=4			##<< see \code{\link{findSplit}}, potentially modified due to argument \code{minPSub}
	,isBreakEarly = TRUE	##<< if TRUE and argument \code{checkSlopesFirst} is given, then check slope angles for variables given in checkSlopesFirst first and if break is found, do not evaluate all the other slope angles
	,isBreakEarlySubs = TRUE ##<< same as argument \code{isBreakEarly} for recursive call from within 
	,checkSlopesFirst = data.frame(ivar=1,j1AlphaSlope=1,jAlphaSlope=2)[FALSE,]	##<< data.frame with entries ivar, j1AlphaSlope and j2AlphaSlope as returned in entry resD in result of ressplit
	,argsFSplit=list()	##<< further arguments to \code{\link{findSplit}}
	,pSub=1				##<< the fraction that a subSample comprises within a bigger sample (used in recursive calls) 
	,minPSub=0.05		##<< minimum fraction a subSample, below which the sample is not split further
	, upperParBounds = numeric(0)
	### named numeric vectors: giving upper parameter bounds: lowerBound < par <= upperBound
	### for exploring subspaces of the limiting distribution, see details
	, lowerParBounds = numeric(0)  ##<< similar to upperParBounds: sample > bound
	,splits=numeric(0)	##<< named numeric numeric vector: history of splitting points
	,maxNSample=128					##<< if given a value, then aSample is thinned before to given number of records (for efficiently calculating variances)
	, verbose=FALSE		##<< if TRUE then prints on every terminating subspace
){
	if( 0 != length(maxNSample) && is.finite(maxNSample) && (nrow(aSample) > maxNSample) ){
		aSample <- aSample[ sample.int( nrow(aSample), maxNSample),]		
	}
	##details<< 
	# uses interval: lower < val <= upper
	# search for a single splitting point
	nSplitMin <- min( nSplit, floor(pSub/minPSub)-1 )		# each subPopulation must comprise at least part minPSub, hence do not split in too many parts
	if( is.null(argsFSplit) ) argsFSplit <- list()
	argsFSplitMod <- within(argsFSplit, nSplit<-nSplitMin )
	boMinPSub <- (nSplitMin < 1)
	#mtrace(findSplit)
	if( !boMinPSub) 
		resSplit <- do.call( findSplit, c(list(aSample=aSample, isBreakEarly=isBreakEarly, checkSlopesFirst=checkSlopesFirst), argsFSplitMod ))
	if( boMinPSub || is.na(resSplit$split) ){
		if( verbose ) print(paste("getSubSpaces: pSub=",signif(pSub,2),"splits=",paste(names(splits),signif(splits,2),sep=":",collapse=", ")))
		# when no splitting point was found, return the sample without parameter bounds
		##value<< a list with entries
		list( 
			spaces=list(list(	##<< a list with an entry for each subspace. Each Entry is a list with entries \itemize{
				##describe<< 
				sample = aSample		##<< numeric matrix: a subsample constrained to the subspace lb < val <= ub with col parameters
				, upperParBounds = upperParBounds	##<< list with each entry numeric scalar: upper parameter bounds
				, lowerParBounds = lowerParBounds	##<< list with each entry numeric scalar: upper parameter bounds
				, pSub = pSub			##<< the proportion of the subSample to the overall Sample
				, splits=splits		##<< named numeric vector of splitting points
				##end<< 
			))
			#, iVars=NA	##<< integer vector: the indices of variables to check for splits. If \code{isBreakEarly=FALSE} its has been updated. NA if no split was found. 
			#, jVarsVar=NA	##<< integer vector: the indices of variables to check scales of variances given a split in iVar. If \code{isBreakEarly=FALSE} its has been updated. NA if no split was found.
			, resD=data.frame(ivar=1,jPVar=1,pVar=NA_real_, j1AlphaSlope=1,jAlphaSlope=2, dAlphaSlope=NA_real_)[FALSE,] ##<< result details of \code{\link{findSplit}}
		)	##end<< 
	}else{
		# split the sample and recursively invoke \code{getSubSpaces}
		#argsFSplit$iVars <- resSplit$iVars
		#argsFSplit$jVarsVar <- resSplit$jVarsVar
		boLeft <- (aSample[,resSplit$varName] <= resSplit$split)
		sampleLeft <- aSample[boLeft, ]
		sampleRight <- aSample[!boLeft,]
		pSubLeft <- pSub*resSplit$perc
		pSubRight <- pSub*(1-resSplit$perc)
		upperParBoundsLeft <- upperParBounds;  vectorElements(upperParBoundsLeft) <- resSplit$split
		lowerParBoundsRight <- lowerParBounds; vectorElements(lowerParBoundsRight) <- resSplit$split
		splitsNew <- c(splits, resSplit$split)
		# adjust the parBounds
		
		#spacesLeft <- getSubSpaces(sampleLeft, nSplit=nSplit, isBreakEarly=TRUE, argsFSplit=argsFSplit, pSub=pSubLeft, minPSub=minPSub)$spaces
		#spacesRight <- getSubSpaces(sampleRight, nSplit=nSplit, isBreakEarly=TRUE, argsFSplit=argsFSplit, pSub=pSubRight, minPSub=minPSub)$spaces
		#even when provided iVars for first level, need to check on subspaces for different variables agin
		resLeft <- spacesLeft <- getSubSpaces(sampleLeft, nSplit=nSplit, isBreakEarly=isBreakEarlySubs, isBreakEarlySubs=isBreakEarlySubs, checkSlopesFirst=resSplit$resD, argsFSplit=argsFSplit, pSub=pSubLeft, minPSub=minPSub, splits=splitsNew, upperParBounds=upperParBoundsLeft, lowerParBounds=lowerParBounds)$spaces
		resRight <- spacesRight <- getSubSpaces(sampleRight, nSplit=nSplit, isBreakEarly=isBreakEarlySubs, isBreakEarlySubs=isBreakEarlySubs, checkSlopesFirst=resSplit$resD, argsFSplit=argsFSplit, pSub=pSubRight, minPSub=minPSub, splits=splitsNew, upperParBounds=upperParBounds, lowerParBounds=lowerParBoundsRight)$spaces
		# add the splitting point as a new border to the results
		# XXadd return iVars and jVarsVar
		.tmp.f <- function(){
			resLeft <- lapply( spacesLeft, function(spEntry){
					posVar <- match(resSplit$varName, names(spEntry$upperParBounds) )
					spEntry$upperParBounds <- if( is.finite(posVar) )
						c(spEntry$upperParBounds[posVar], spEntry$upperParBounds[-posVar]) # make this variable first position in parBounds
					else
						spEntry$upperParBounds <- c(resSplit$split, spEntry$upperParBounds)
					spEntry
				})
			resRight <- lapply( spacesRight, function(spEntry){
					posVar <- match(resSplit$varName, names(spEntry$lowerParBounds) )
					spEntry$lowerParBounds <- if( is.finite(posVar) )
							c(spEntry$lowerParBounds[posVar], spEntry$lowerParBounds[-posVar])
						else
							spEntry$lowerParBounds <- c(resSplit$split, spEntry$lowerParBounds)
					spEntry
				})
		}
		res <- list(
			spaces = c( resLeft, resRight)	
			#, iVars=resSplit$iVars	 
			#, jVarsVar=resSplit$jVarsVar
			,resD=resSplit$resD
		)
	}
} 
attr(getSubSpaces,"ex") <- function(){
	data(den2dCorEx)
	aSample <- stackChains(thin(concatPops(den2dCorEx$mcBulk), start=300))[,-1]
	#mtrace(getSubSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.4 )
	str(subSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.05, argsFSplit=list(debugSequential=TRUE))
	subSpaces <- getSubSpaces(aSample, minPSub=0.05)
	(tmp <- sapply( subSpaces$spaces, function(subSpace){nrow(subSpace$sample)})/nrow(aSample)) # percentiles
	(tmp <- sapply( subSpaces$spaces, "[[", "upperParBounds" )) # bounds
	#visualize the splits
	ss1 <- do.call(rbind, lapply( seq_along(subSpaces$spaces), function(i){cbind(spaceInd=i, subSpaces$spaces[[i]]$sample)}))
	plot(b~a, as.data.frame(ss1), col=rainbow(length(subSpaces$spaces))[spaceInd] )
	plot(b~a, as.data.frame(ss1), col=rainbow(length(subSpaces$spaces))[spaceInd], ylim=c(-5000,+5000) )
	#twUtestF(getSubSpaces)	# there are unit tests for this function
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

.tmp.f <- function(){
	# TODO update implementation for Pops and Blocks
setMethodS3("divideTwDEMCBatch","default", function( 
	### iteratively run batches using divideTwDEMCBatch
	x						##<< numeric array: rows:steps, col:parameters including logLik, 3: independent populations
	, ...					##<< further arguments to \code{\link{divideTwDEMC}}
	, nGen=2*512			##<< number of generations
	, nGenBatch=512			##<< number of generations within one batch
	, thinPastFac=0.2		##<< thinning the past between batches to speed up localization between 0 (no past) and 1 (keep entire past)
	, subPercChangeCrit=1.6	##<< if all subPercChange of all sub-populations are belwo this value in two last batches, may assume convergence and skip further batches 
	, argsFSplitPop=vector("list",dim(x)[3])	##<< for each population: list of arguments  passed \code{\link{getSubSpaces}} and further to \code{\link{findSplit}}, e.g. for passing order of variables to check in \code{iVars} and \code{jVarsVar}
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
	# if iVars and jVarsVar for all populations have been supplied, then do not need to check all splits on first batch
		#argsFSplit<- argsFSplitPop[[1]]
	firstIsCheckSplitVars <- any( sapply(argsFSplitPop, function(argsFSplit){ is.null(argsFSplit$iVars) || is.null(argsFSplit$jVarsVar) }) )
	for( iBatch in 1:nBatch ){
		# need to check on every batch, because variable ordering might be different in nested spaces
		#isCheckAllSplitVars <- (iBatch %% iCheckAllSplitVars == 0) || (firstIsCheckSplitVars && iBatch==1) 
		isCheckAllSplitVars <- (firstIsCheckSplitVars && iBatch==1) 
		cat(paste(iN," out of ",nGen," generations completed.     ",date(),"\n",sep=""))
		resBatch <- divideTwDEMC(  ss,nGen=min(nGenBatch, nGen-iN), isBreakEarly=!isCheckAllSplitVars, argsFSplitPop=argsFSplitPop	
		#,fLogDen=den2dCor 
		,...
		)
		if(isCheckAllSplitVars) for( iPop in 1:nPops ){
			argsFSplitPop[[iPop]]$iVars <- resBatch[[iPop]]$iVars	# update variable order
			argsFSplitPop[[iPop]]$jVarsVar <- resBatch[[iPop]]$jVarsVar	# update variable order
		}
		subPercChange[iBatch,] <- subPercChangeCur <- lapply( resBatch, "[[", "subPercChange" )
		maxSubPercChange[iBatch,] <- maxSubPercChangeCur <- sapply( subPercChangeCur, max )
		ssCurrent <- abind(lapply(resBatch,"[[","sample"), rev.along=0)
		ssPast <-  if( !is.null(thinPastFac) ){
			iKeep <- round(seq(1,nrow(ss),length.out=nrow(ss)*thinPastFac))
			ss[iKeep,,,drop=FALSE]
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
	data(den2dCorEx)
	aTwDEMC <- 	thin(den2dCorTwDEMC, start=300)
	aSample <- stackChainsPop(aTwDEMC)
	#mtrace(divideTwDEMCBatch.default)
	ssRes1 <- ssRes <- divideTwDEMCBatch(aSample, nGen=512*3, fLogDen=den2dCor )
	#ssRes <- divideTwDEMCBatch(ssRes1, nGen=512*1, fLogDen=den2dCor )	#calling divideTwDEMCBatch.divideTwDEMCBatch
	#mtrace(getSubSpaces)
	#tmp <- getSubSpaces(ssRes1$sample[,-1,1])
	#mtrace(divideTwDEMCBatch.divideTwDEMCBatch)
	ssRes2 <- ssRes <- divideTwDEMCBatch(ssRes1, nGen=512*6, fLogDen=den2dCor )	#calling divideTwDEMCBatch.divideTwDEMCBatch
	
	#tmp <- divideTwDEMC( ssRes1$sample[,,2,drop=FALSE], nGen=512, fLogDen=den2dCor)
	#ssImpPops <- abind( tmp[[1]]$sample, tmp[[1]]$sample, rev.along=0) 
	
	#wSubsB <- ssRes$wSubs[,1]
	#wSubsB <- ssRes$wSubs[,2]
	ssRes$maxSubPercChange
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

twRunDivideTwDEMCBatch <- function(
	### wrapper to divideTwDEMC compliant to runCluster.R
	argsDivideTwDEMCBatch	 
	### Arguments passed to DivideTwDEMCBatch -> DivideTwDEMC -> twDEMCBlockInt -> fLogDen.
	### It is updated by \dots.
	### After update it must contain entries Zinit and fLogDen
	### It is further searched for entries nPops, and argsFLogDen, and logDenCompX. The latter are initialized to defaults  \code{1,list(),character(0)} respectively if not found.   
	,...				 	##<< further arguments passed to DivideTwDEMCBatch -> DivideTwDEMC -> twDEMCBlockInt -> fLogDen
	,prevResRunCluster=NULL	##<< results of call to twRunDEMC, argument required to be called from runCluster.R
	,restartFilename=NULL	##<< name of the file to store restart information, argument required to be called from runCluster.R 
){
	#update argsDEMC to args given 
	argsDEMC <- argsDivideTwDEMCBatch
	
	# update the restartFilename argument
	.dots <- list(...)
	argsDEMC[ names(.dots) ] <- .dots
	#for( argName in names(.dots) ) argsDEMC[argName] <- .dots[argName]
	if( !is.null(restartFilename)) argsDEMC$restartFilename <- restartFilename
	
	# do the actual call
	res <- do.call( divideTwDEMCBatch, argsDEMC )
}
}


divideTwDEMCPop <- function( 
	### splits a specific populations of x into subpopulations for given subSpaces 
	pop			##<< numeric array (nStep, nParm, nChain): an object of class aTwDEMCPops returned by \code{\link{twDEMCBlockInt}}
	,subSpaces	##<< a list (nPop): each entry listing subspaces for the given population with entries lowerParBounds, upperParBounds, and splits
	,isCheckBorderConsistency = TRUE	##<< if FALSE the check for border consistency is omitted.
){
	nSub <- length(subSpaces)
	if( nSub < 2){ 
		warning("divideTwDEMCPop: called with less than 2 subspaces.")
		return( list(pop) )
	}
	nChainPop <- dim(pop$parms)[3]
	# list of parms for each chain in pop
	aSampleChain <- lapply(1:nChainPop, function(iChain){ adrop(pop$parms[,,iChain ,drop=FALSE],3) })
	#iSub <- nSub
	#iSub <- 6
	newPops <- lapply(1:nSub, function(iSub){
		subSpace <- subSpaces[[iSub]]
		# check that the subspace borders are consistent with the population
		if( isCheckBorderConsistency && (length(pop$splits) != 0)){
			if( !isTRUE(all.equal(subSpace$splits[ 1:length(pop$splits)] , pop$splits)) )
				stop("divideTwDEMCPop: encountered subspaces with inconsitent splits to pop")
			for( pName in names(pop$lowerParBounds) ){
				if( !(pName %in% names(subSpace$lowerParBounds)) || 
					!(subSpace$lowerParBounds[pName] >= pop$lowerParBounds[pName]) )  
						stop("divideTwDEMCPop: encountered subspaces with lower lower parameter bound than pop")
			}
			for( pName in names(pop$upperParBounds) ){
				if( !(pName %in% names(subSpace$upperParBounds)) || 
					!(subSpace$upperParBounds[pName] <= pop$upperParBounds[pName]) )  
					stop("divideTwDEMCPop: encountered subspaces with higher upper parameter bound than pop")
			}
			#.getParBoundsPops(list(pop, subSpace))
		}
		
		#iChain <- nChainPop
		.subSamplesChain <- lapply( 1:nChainPop, function(iChain){
			aSample <- aSampleChain[[iChain]]
			#i <- length(subSpace$lowerParBounds) 
			boLower <- lapply( seq_along(subSpace$lowerParBounds), function(iBound){
					aSample[,names(subSpace$lowerParBounds)[iBound] ] > subSpace$lowerParBounds[iBound]
				})
			boUpper <- lapply( seq_along(subSpace$upperParBounds), function(iBound){
					aSample[,names(subSpace$upperParBounds)[iBound] ] <= subSpace$upperParBounds[iBound]
				})
			boMat <- abind( c(boLower,boUpper), along=2 )
			boKeep <- apply(boMat,1,all)
			.subsetTwDEMCPop( pop, iKeep=boKeep, iChain=iChain)
		})	# iChain
		# collect the subsamples of all chains into new chains of equal length
		ss <- combineTwDEMCPops(.subSamplesChain)
		ssu <- .unstackPopsTwDEMCPops(ss$pop, nChainPop)	# here may loose a few samples due to not a multiple of nChains
		newPop <- .concatChainsTwDEMCPops(ssu)	# combine chains into one array
		# do the new parameter bounds
		newPop$lowerParBounds <- subSpace$lowerParBounds	
		newPop$upperParBounds <- subSpace$upperParBounds
		newPop$splits <- subSpace$splits
	#if( iSub==9) recover()
		.tmp.f <- function(){
			# add initial lower and upper bounds to unbounded parameters
			for( pName in pop$lowerParBounds ){
				if( !(pName %in% names(newPop$lowerParBounds)) ) 
					newPop$lowerParBounds <- c( pop$lowerParBounds[pName], newPop$lowerParBounds )
			}
			for( pName in pop$upperParBounds ){
				if( !(pName %in% names(newPop$upperParBounds)) ) 
					newPop$upperParBounds <- c( pop$upperParBounds[pName], newPop$upperParBounds )
			}
		}
		newPop$spaceInd <- pop$spaceInd		# inherit the reference to space replicate
		newPop
	}) # iSub
	### populations with lBound < val <= uBound
	newPops
}
attr(divideTwDEMCPop,"ex") <- function(){
	data(den2dCorEx)
	x <- den2dCorEx$mcBulk
	newPops <- divideTwDEMCPop(x$pops[[1]], den2dCorEx$subspaces0[[1]]$spaces )
	length(newPops)
	str(newPops[[1]])
	#lapply(newPops,"[[","splits")
}


divideTwDEMCPops <- function( 
	### splits all populations of x into subpopulations for given subSpaces 
	x			##<< an object of class aTwDEMCPops returned by \code{\link{twDEMCBlockInt}}
	,subSpacesPop	##<< a list (nPop): each entry listing subspaces for the given population
){
	nPop <- getNPops(x)
	if( length(subSpacesPop) != nPop ) 
		stop("divideTwDEMCPops: argument subSpacesPop must list subspaces for all populations in x")
	res <- x
	newPops <- lapply( 1:nPop, function(iPop){
			#mtrace(divideTwDEMCPop)
			resPopsI <- divideTwDEMCPop( x$pops[[iPop]], subSpacesPop[[iPop]]$spaces )		
		})
	res$pops <- do.call( c, newPops )
	### argument x with pops splitted into subspaces lBound < val <= uBound
	res
}
attr(divideTwDEMCPops,"ex") <- function(){
	data(den2dCorEx)
	#x <- den2dCorEx$mcSubspaces0
	x <- den2dCorEx$mcBulk
	#sapply( x$pops, "[[", "spaceInd")
	#subSpaces <- do.call( c, lapply(den2dCorEx$subspaces0,"[[","spaces"))
	subSpacesPop <- den2dCorEx$subspaces0
	#mtrace(divideTwDEMCPops)
	xNew <- divideTwDEMCPops(x, subSpacesPop )
	sapply(den2dCorEx$subspaces0, function(subSpacesPop){ length(subSpacesPop$spaces) } )
	getNPops(xNew)
	getSpacesPop(xNew)
	str(xNew$pops[[9]])
	#lapply(xNew$pops,"[[","splits")
}




