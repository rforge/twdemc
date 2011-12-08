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
	data(den2dCorTwDEMC)
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
	,splitHist=NA_real_[FALSE]	##<< named numeric numeric vector: history of splitting points
){
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
		print(paste("getSubSpaces: pSub=",signif(pSub,2),"splits=",paste(names(splitHist),signif(splitHist,2),sep=":",collapse=", ")))
		# when no splitting point was found, return the sample without parameter bounds
		##value<< a list with entries
		list( 
			spaces=list(list(	##<< a list with an entry for each subspace. Each Entry is a list with entries \itemize{
				##describe<< 
				sample = aSample		##<< numeric matrix: a subsample constrained to the subspace with col parameters
				, upperParBounds = c()	##<< list with each entry numeric scalar: upper parameter bounds
				, lowerParBounds = c()	##<< list with each entry numeric scalar: upper parameter bounds
				, pSub = pSub			##<< the proportion of the subSample to the overall Sample
				, splits=splitHist		##<< named numeric vector of splitting points
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
		splitHistNew <- c(splitHist, resSplit$split)
		#spacesLeft <- getSubSpaces(sampleLeft, nSplit=nSplit, isBreakEarly=TRUE, argsFSplit=argsFSplit, pSub=pSubLeft, minPSub=minPSub)$spaces
		#spacesRight <- getSubSpaces(sampleRight, nSplit=nSplit, isBreakEarly=TRUE, argsFSplit=argsFSplit, pSub=pSubRight, minPSub=minPSub)$spaces
		#even when provided iVars for first level, need to check on subspaces for different variables agin
		spacesLeft <- getSubSpaces(sampleLeft, nSplit=nSplit, isBreakEarly=isBreakEarlySubs, isBreakEarlySubs=isBreakEarlySubs, checkSlopesFirst=resSplit$resD, argsFSplit=argsFSplit, pSub=pSubLeft, minPSub=minPSub, splitHist=splitHistNew)$spaces
		spacesRight <- getSubSpaces(sampleRight, nSplit=nSplit, isBreakEarly=isBreakEarlySubs, isBreakEarlySubs=isBreakEarlySubs, checkSlopesFirst=resSplit$resD, argsFSplit=argsFSplit, pSub=pSubRight, minPSub=minPSub, splitHist=splitHistNew)$spaces
		# add the splitting point as a new border to the results
		# XXadd return iVars and jVarsVar
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
		res <- list(
			spaces = c( resLeft, resRight)	
			#, iVars=resSplit$iVars	 
			#, jVarsVar=resSplit$jVarsVar
			,resD=resSplit$resD
		)
	}
} 
attr(getSubSpaces,"ex") <- function(){
	data(den2dCorTwDEMC)
	aSample <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
	#mtrace(getSubSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.4 )
	str(subSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.05, argsFSplit=list(debugSequential=TRUE))
	subSpaces <- getSubSpaces(aSample, minPSub=0.05)
	(tmp <- sapply( subSpaces$spaces, function(subSpace){nrow(subSpace$sample)})/nrow(aSample)) # percentiles
	(tmp <- sapply( subSpaces$spaces, "[[", "upperParBounds" )) # bounds
	#visualize the splits
	ss1 <- do.call(rbind, lapply( seq_along(subSpaces$spaces), function(i){cbind(iSpace=i, subSpaces$spaces[[i]]$sample)}))
	plot(b~a, as.data.frame(ss1), col=rainbow(length(subSpaces$spaces))[iSpace] )
	plot(b~a, as.data.frame(ss1), col=rainbow(length(subSpaces$spaces))[iSpace], ylim=c(-5000,+5000) )
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


#XXTODO: adapt divideTwDEMC to changed getSubSpaces that returns now updated iVars and jVarsVar 
divideTwDEMC <- function(
	### run twDEMCBlock on subspaces
	aSample					##<< numeric array: rows:steps, col:parameters including rLogLik in first column, 3: independent populations
	,nGen=512				##<< the number of generations for twDEMCBlock
	, controlTwDEMC = list()	##<< list argument to \code{\link{twDEMCBlock}} containing entry thin
	, ...					##<< further arguments to \code{\link{twDEMCBlock}}
	, minPSub = 0.1			##<< passed to \code{\link{getSubSpaces}}
	, subPercChangeCrit=1.6	##<< if all subPercChange of all sub-populations are belwo this value in two last batches, may assume convergence and skip further batches 
	, isBreakEarly = FALSE	##<< passed to \code{\link{getSubSpaces}}
	, argsFSplitPop=vector("list",dim(aSample)[3])	##<< for each population: list of arguments  passed \code{\link{getSubSpaces}} and further to \code{\link{findSplit}}, e.g. for passing order of variables to check in \code{iVars} and \code{jVarsVar}
	, nChainPop=4			##<< number of chains in subPopulations
	, m0 = calcM0twDEMC( ncol(aSample), nChainPop=nChainPop )	##<< number of samples per chain to initialize subPopulations
	, nrow0 = 4000			##<< number of rows in initial overall sample
	, attachDetails=FALSE	##<< set TRUE to report upperParBounds, lowerParBounds, and pSubs per subPopulation
	, nChainPar=16			##<< number of chains to run in parallel, good choice is 2*nCpu
	, minNSamplesSub=32		##<< minimum number of records in a subspace sample, increase to avoid wrong estimation of weights
	, fCheckProblems=checkProblemsSpectral	##<< function applied to twDEMCBlock results of each subspace	
	, argsFCheckProblems=list() ##<< further arguments to argsFCheckProblems
	, dumpfileBasename="recover"
	, subSpacesPop=list()	##<< may provide previous results of \code{\link{getSubSpaces}} per population to save computing time
){
	##details<< 
	## the first columns of aSample records the logDensity of the sample for consitency. 
	## It is not used and may be initialized to any value. 
	#samplePop <- aSample[,,1]
	if( length(dim(aSample))==2 ){
		warning("divideTwDEMC:third dimension missing of aSample missing - assuming 1 population.")
		aSample <- array(aSample, dim=c(dim(aSample),1), dimnames=c(dimnames(aSample), pops=NULL) )
	}
	nPop <- dim(aSample)[3]
	if( is.null(controlTwDEMC$thin) ) controlTwDEMC$thin <- 4
	thin <- controlTwDEMC$thin
	nGenThin <- (nGen %/% thin)*thin		# when dealing with thinned samples use the number of generations of last recorded thinned generation
	
	thinnedSample <- if( nrow(aSample) < nrow0 ) aSample else aSample[round(seq(1,nrow(aSample),length.out=nrow0)),,]
	nInit <- m0*nChainPop		# number of samples to initialize population
	if( nInit > nrow(thinnedSample) ) stop(paste("divideTwDEMC: aSample or nGen has too few cases, need at least",nInit))
	minPSub <- max( minPSub, nInit/nrow(thinnedSample) )
	#mtrace(getSubSpaces)
	#iPop=2
	nBlock <-  attr(aSample, "nBlock")
	subSpacesL <- if( 0==length(subSpacesPop) ) lapply( 1:dim(thinnedSample)[3], function(iPop){
			samplePop <- thinnedSample[,,iPop]
			#mtrace(getSubSpaces)
			#mtrace(findSplit)
			getSubSpaces(samplePop[,-(1:nBlock)], minPSub=minPSub, isBreakEarly=isBreakEarly, argsFSplit=argsFSplitPop[[iPop]])	# here omit the logDensity column
		}) else subSpacesPop
	subSpaces <- lapply( subSpacesL, "[[" ,"spaces" )
	nSubPops <- sapply(subSpaces, length)	# number of subs per populations
	subSpacesFlat <- do.call( c, subSpaces )	# put all subPopulations of all populations on same level
	nSub <- length(subSpacesFlat)			# number of total subs
	iSubs <- 1:nSub						# set of sub indices
	#popSub <- sapply( 1:nPop, function(iPop){ rep,each=nSubPops)	
	iSubsPop <- {				# index of subs in flat version for each pop
		cumNSubPops <- cumsum(nSubPops)
		lapply(1:nPop, function(iPop){ cumNSubPops[iPop]+1-(nSubPops[iPop]:1)})	
	}

	# calculating initial quantiles and number of genrations
	popQSubs <- lapply( subSpaces, sapply, "[[","pSub" ) # initial percentile of subspace of the overall space
	qSubs <- do.call(c,popQSubs)
	nGen0 <- pmax(m0*thin, ceiling(minNSamplesSub*thin/nChainPop), nGen*qSubs)		# at minimum m0*thin generations to keep sufficient samples for extending the run
	nGen0Thin <- (ceiling(nGen0/thin)+2)*thin		# +2 to avoid extending runs with little shift 
	
	# setting up initial populations (each sub of overall populations becomes one population)
	subs0 <- vector("list",nSub)	
	#iSub <- nSub
	for( iSub in iSubs){
		subSpacesI <- subSpacesFlat[[iSub]]
		s0 <- subSpacesI$sample[sample.int(nrow(subSpacesI$sample),nInit), ]
		.tmp <- matrix( 1:(m0*nChainPop), ncol=nChainPop )
		#iChainPop <- 4
		Zinit <- abind( lapply(1:nChainPop, function(iChainPop){s0[.tmp[,iChainPop],]}), rev.along=0)
		names(dimnames(Zinit)) <- c("steps","parms","chains")
		subs0[[iSub]]$parms <- Zinit
		subs0[[iSub]]$upperParBounds <- subSpacesI$upperParBounds
		subs0[[iSub]]$lowerParBounds <- subSpacesI$lowerParBounds
	}
	
	resTwDEMC <- resTwDEMC0 <- twDEMCBlockInt(pops=subs0, nGen=nGen, controlTwDEMC=controlTwDEMC, ...)
	##details<< 
	## A first estimate of the proportion of samples from different subspaces
	## are the initial percentiles qp. 
	## The proportion of the samples from different subspaces is estimated
	## by the proportions of integrated density of the subspaces.
	## These proportions are esimated by average density over all samples multiplied by an
	## estimate proportional to the volume: the initial quantiles.
	#iPop <- 1
	#iPop <- 2
	popLogMeanDensSubs <- lapply( 1:nPop, function(iPop){
			iSubs <- iSubsPop[[iPop]]
			#iSub <- 2
			#iSub <- nSub
			lw <- sapply(iSubs,function(iSub){  # average the unnormalized densities
					ss <- stackChains(resTwDEMC$pops[[iSub]]$logDen)
					ssLogDen <- rowSums(ss)
					twLogSumExp(ssLogDen, shiftUpperBound=TRUE)-log(length(ssLogDen)) 
				})	
			lw - max(lw, na.rm=TRUE)			#divide by a constant on exp scale to avoid numerical errors
		})	
	logMeanDensSubs <- do.call(c, popLogMeanDensSubs )
	#barplot(logMeanDensSubs)
	
	# estimate proportion of subspaces in the limiting distribution
	# it will be used to sample from the subspaces 
	pSubs <- do.call(c, lapply( 1:nPop, function(iPop){
			p2u <- exp(popLogMeanDensSubs[[iPop]])*popQSubs[[iPop]]		# two weights: quantile and density importance
			p2 <- p2u/sum(p2u)									# normalize to 1 within population 
		}))
	subPercChange <- pSubs/qSubs	# ratio of estimated proportion in limiting distribution to proportion of initial proportion 
	nSamplesSubsReq <- round(nGen/thin * pSubs)
	# because of rounding small differences in sum of sample numbers may occur 
	# for small deviations, adjust the number of samples in the largest sample
	# add argument nGen to pops
	for( iPop in 1:nPop ){
		dGen <- nGenThin/thin - sum(nSamplesSubsReq[iSubsPop[[iPop]] ])  
		if( dGen != 0 ){
			if( abs(dGen)>length(iSubsPop[[iPop]]) ) stop("divideTwDEMC: problem in calculating subspace sample numbers")
			iSubMax <- iSubsPop[[iPop]][ which.max( nSamplesSubsReq[iSubsPop[[iPop]] ]) ]
			nSamplesSubsReq[iSubMax] <- nSamplesSubsReq[iSubMax] + dGen
		} 
	}
	
	# extend the previous runs
	nSamplesAdd <- pmax(0, nSamplesSubsReq-getNSamples(resTwDEMC0))
	resTwDEMC <- resTwDEMC1 <- if( max(nSamplesAdd) == 0) resTwDEMC0 else
		#mtrace(twDEMCBlock.twDEMCPops)
		#mtrace(twDEMCBlockInt)
		resTwDEMC <- twDEMCBlock( resTwDEMC0, nGen=nSamplesAdd*thin, controlTwDEMC=controlTwDEMC, ...  )
	#all( getNSamples(resTwDEMC) >= nSamplesSubsReq)

	# do a subsampling
	# iPop=2
	subSamplePop <- function(iPop){
		iSubs  <- iSubsPop[[iPop]] 
		#nSamplesSubsReq[iSubs]
		#(iSub <- iSubs[length(iSubs)])
		ssImpList <-  lapply( iSubs, function(iSub){
			if( nSamplesSubsReq[iSub] == 0) c() else{
				aTwDEMCSub <- concatPops(subPops( resTwDEMC, iPops=iSub ))
				ss <- stackChains(aTwDEMCSub)
				ssp <- ss[,-(1:getNBlocks(aTwDEMCSub))]
				# see .tmp.f below for commands to trace error
				if( nSamplesSubsReq[iSub] == 1) 
					ss[ sample.int( nrow(ss),1),] 
				else {
					resFCheck <- do.call( fCheckProblems, c(list(aTwDEMCSub), argsFCheckProblems))
		if( iSub==21) recover()
					if( resFCheck$hasProblems ){
						if( !is.null(dumpfileBasename) )
							if( dumpfileBasename == "recover"){
								cat("encountered convergence problems in a subSpace. calling recover() \n ")
								recover()
							}  else dump.frames(dumpfileBasename,TRUE)
						stop(paste("divideTwDEMC: checking function encountered problems. Dump in ",dumpfileBasename,sep=""))
					}
					sst <- ss[ round(seq(1,nrow(ss),length.out=nSamplesSubsReq[iSub])),]
				}
			}
		})
		ssImp <- do.call( rbind, ssImpList )
	}
	#mtrace(tmp.f)
	resl <- lapply( 1:nPop, subSamplePop)	# end lapply iPop: giving a single combined sample for each population
	#res <- ssImpPops <- abind(resl, rev.along=0)
	
	##value<< 
	## For each population, a list with entries
	res <- lapply(1:nPop,function(iPop){
			iSubs <- iSubsPop[[iPop]]
			list(
				sample=resl[[iPop]]					##<< numeric matrix: the sampled parameters (rows: cases, cols: 1: logDensity, others parameter dimensions
				,subPercChange=subPercChange[iSubs]	##<< numeric vector: estimated proportion in limiting distribution to proportion of initial proportion 
				#,upperParBounds=subSpacesFlat[iSubs] upperParBounds[iSubs]	##<< list of named numeric vectors: the upper parameter bounds for the subspaces
				#,lowerParBounds=lowerParBounds[iSubs]	##<< dito upperParBounds
				,iVars=subSpacesL[[iPop]]$iVars		##<< reordered index or parameter dimensions to check for split, see \code{\link{findSplit}}
				,jVarsVar=subSpacesL[[iPop]]$jVarsVar		##<< reordered index or parameter dimensions to check for split, see \code{\link{findSplit}}
				#,pSubs0=pSubs0[posPop]
				#,pSubs2=p2[posPop]
			)
		})
	##end<<
}
attr(divideTwDEMC,"ex") <- function(){
	data(den2dCorTwDEMC)
	aTwDEMC <- 	thin(den2dCorTwDEMC, start=300)
	plot( b ~ a, as.data.frame(stackChains(subChains(aTwDEMC,iPops=1)$parms)), xlim=c(-0.5,2), ylim=c(-20,40) )
	aSample <- stackChainsPop(aTwDEMC)
	ss0 <- stackChains(aTwDEMC)
	#aSamplePure <- aSample[,-1,]
	#aSamplePurePop <- aSamplePure[,,1]
	#mtrace(divideTwDEMC)
	#res <- divideTwDEMC(aSample, nGen=128, dInfos=list(list(fLogDen=den2dCor)),  isBreakEarly=FALSE, debugSequential=TRUE )
	res <- divideTwDEMC(aSample, nGen=256, dInfos=list(list(fLogDen=den2dCor)),  isBreakEarly=FALSE, debugSequential=TRUE,  controlTwDEMC=list(DRGamma=0.1) )
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

.debug.divideTwDEMC <- function(){
	tmp <- concatPops(subPops(resTwDEMC0,iPops=iSub))
	ss1 <- stackChains(tmp)
	plot(as.mcmc.list(tmp), smooth=FALSE )
	matplot( tmp$pAccept[,1,], type="l" )
	
	windows(record=TRUE); plot(as.mcmc.list(aTwDEMCSub), smooth=FALSE )
	matplot( aTwDEMCSub$pAccept[,1,], type="l" )
	matplot( aTwDEMCSub$logDen[,1,], type="l" )
	plot(ss[,"a"], ss[,"b"], col=rev(heat.colors(100))[ twRescale(ss[,1],c(20,100)) ] )
	
	resFCheck
	
	tmp <- findSplit( ss[,-1], isBreakEarly=FALSE, rVarCrit = 2^2 )
	mtrace(getSubSpaces)
	tmp2 <- getSubSpaces( ss[,-1], isBreakEarly=FALSE, argsFSplit=list(rVarCrit = 2^2, rAlphaSlopeCrit=base:::pi/2),minPSub=minPSub, pSub=pSubs[iSub] )
	lapply(tmp2$spaces, "[[", "upperParBounds" )
	colnames(ss)[-1][tmp$iVars]

	iSubsCheckSubspaces <- which(  (pSubs > 2*minPSub) & (subPercChange > subPercChangeCrit) )
	#iSub <- iSubsCheck[length(iSubsCheck)]
	if( iSub %in% iSubsCheckSubspaces ){
		subSubSpaces <- getSubSpaces( ssp, isBreakEarly=FALSE, argsFSplit=argsFSplitPop[[ popSub]],minPSub=minPSub, pSub=pSubs[iSub] )
	}
	
	
	pSub > 2*minPSub
}


setMethodS3("divideTwDEMCBatch","divideTwDEMCBatch", function(
	### append runs to result of \code{\link{divideTwDEMCBatch.default}}
	x						##<< result of former call to divideTwDEMCBatch
	, ...					##<< further arguments to \code{\link{divideTwDEMCBatch.default}}
){
	.dots <- list(...)
	# extract iVars and jVarsVar per population from x$resDivideTwDEMC
	if( is.null(.dots$argsFSplitPops) ){
		 nPops = dim(x$sample)[3]
		 .dots$argsFSplitPop <- vector("list",nPops)
		 for( iPop in 1:nPops ){
			.dots$argsFSplitPop[[iPop]]$iVars <- x$resDivideTwDEMC[[iPop]]$iVars	# update variable order
			.dots$argsFSplitPop[[iPop]]$jVarsVar <- x$resDivideTwDEMC[[iPop]]$jVarsVar	# update variable order
		 }
	}
	resArray <- do.call( divideTwDEMCBatch.default, c(list(x$sample), .dots) )
	within(resArray,{
		subPercChange <- rbind( x$subPercChange, subPercChange )
		maxSubPercChange <- c(x$maxSubPercChange, maxSubPercChange)
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
	data(den2dCorTwDEMC)
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
	### Arguments passed to DivideTwDEMCBatch -> DivideTwDEMC -> twDEMCInt -> fLogDen.
	### It is updated by \dots.
	### After update it must contain entries Zinit and fLogDen
	### It is further searched for entries nPops, and argsFLogDen, and logDenCompX. The latter are initialized to defaults  \code{1,list(),character(0)} respectively if not found.   
	,...				 	##<< further arguments passed to DivideTwDEMCBatch -> DivideTwDEMC -> twDEMCInt -> fLogDen
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

divideTwDEMCPops <- function( 
	### splits all populations of x into subpopulations for given subSpaces 
	x			##<< an object of class aTwDEMCPops returned by \code{\link{twDEMCBlockInt}}
	,subSpacesPop	##<< a list (nPop): each entry listing subspaces for the given population
){
	nPop <- getNPops(x)
	if( length(subSpacesPop) != nPop ) 
		stop("divideTwDEMCPops: argument subSpacesPop must list subspaces for all populations in x")
	
	res <- x
	res$pops <- newPops <- lapply( 1:nPop, function(iPop){
		res <- divideTwDEMCPop( x$pops[[iPop]], subSpacesPop[[iPop]] )		
	})
	### argument x with pops splitted into subspaces
	res
}

divideTwDEMCPop <- function( 
	### splits a specific populations of x into subpopulations for given subSpaces 
	pop			##<< numeric array (nStep, nParm, nChain): an object of class aTwDEMCPops returned by \code{\link{twDEMCBlockInt}}
	,subSpaces	##<< a list (nPop): each entry listing subspaces for the given population
){
	nSub <- length(subSpaces$spaces)
	#iSub <- nSub
	#iSub <- 6
	newPops <- lapply(1:nSub, function(iSub){
		subSpace <- subSpaces$spaces[[iSub]]
		nChainPop <- dim(pop$parms)[3]
		#iChain <- nChainPop
		.subSamplesChain <- lapply( 1:nChainPop, function(iChain){
			aSample <- adrop(pop$parms[,,iChain ,drop=FALSE],3)
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
		ss <- .stackChainsTwDEMCPops(.subSamplesChain)
		ssu <- .unstackChainsTwDEMCPops(ss, nChainPop)	# here may loose a few samples due to not a multiple of nChains
		newPop <- .concatChainsTwDEMCPops(ssu)
		newPop$lowerParBounds <- subSpace$lowerParBounds	
		newPop$upperParBounds <- subSpace$upperParBounds	
		newPop
	}) # iSub
}
attr(divideTwDEMCPop,"ex") <- function(){
	data(den2dCorTwDEMC)
	x <- den2dCorTwDEMCPops
	newPops <- divideTwDEMCPop(x$pops[[1]], den2dCorSubSpaces[[1]] )
}


