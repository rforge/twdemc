.calcRhat2 <- function(
	### calculate squared RHat Gelman Diagnostics for a single dimension for several chains
	parmsi			##<< numeric matrix steps x chains 
	,n=nrow(parmsi)	##<< number of steps: maybe passed for efficiency in repeated calculations 
	,m=ncol(parmsi)	##<< number of chains: maybe passed for efficiency in repeated calculations
){
	muyj <- sapply(1:m, function(j){ mean(parmsi[,j]) }) # mean of each chain
	muy <- mean(muyj)	# mean across chains
	B <- n/(m-1)*sum((muyj-muy)^2) # between chain variance
	s2j <- sapply(1:m, function(j){sum((parmsi[,j]-muyj[j])^2) })/(n-1) #squared sums of each chain
	W <- sum(s2j)/m	# within chain variance
	VarPlus=(n-1)/n*W + 1/n*B
	#c( B=B, W=W, VarPlus=VarPlus, Rhat=round(sqrt(VarPlus/W),3), effSize=min(n*m,n*m*VarPlus/B) )
	VarPlus/W
}

checkConvergenceGelman <- function(
	### Gelman RHat criterion applied to twDEMC result assuming all chains one population.
	res					##<< see return value of \code{\link{twDEMCInt}} ($parms (d x nStep x nChain) )
	,burninFrac=0.5 	##<< fraction of the chain to be discarded (default 0.5)
	,rHatMin = 1.1		##<< rHat criterion, upper bound that is regarded as convergence 
){
	# checkConvergenceGelman
	
	##details<< 
	## see Gelman04 (twutz:Gelman04_3#Inference_and_assessing_convergence)
	
	##seealso<<  
	## \code{\link{twDEMCInt}}
	
	##details<< 
	## There are several methods to get diagnostics for a twDEMC run. \itemize{
	## \item{ the Gelman criterion: this method  } 
	## \item{ the theorectical minimum logDen-Value for significant model difference : \code{\link{getRLogDenQuantile}}  } 
	##}
	
	l <- dim(res$parms)[2]
	parms <- res$parms[ ,max(1,ceiling(l*burninFrac)):l, ]	# later part of all the chains
	n <- dim(parms)[2]	# number of steps
	m <- dim(parms)[3]   # number of chains
	rl <- sqrt(apply( parms, 1, .calcRhat2, n=n, m=m  ))
	names(rl) <- rownames(parms)
	res <- all(rl <= rHatMin )
	attr(res,"rHat") <- rl
	res
	### all rl <= criterion for each chain
}
#mtrace( checkConvergenceGelman )


checkConvergenceGelmanPops <- function(
	### Gelman RHat criterion applied to each population and between populations
	aTwDEMC				##<< see return value of \code{\link{twDEMCInt}} ($parms (d x nStep x nChain) )
	,burninFrac=aTwDEMC$nGenBurnin/getNGen(aTwDEMC) 	##<< fraction of the chain to be discarded (default 0.5)
	,rHatMin = 1.1		##<< rHat criterion, upper bound that is regarded as convergence 
){
	# checkConvergenceGelmanPops
	thinned <- thin(aTwDEMC,start=ceiling(getNGen(aTwDEMC)*burninFrac)  )
	nPops <- getNPops(thinned)
	#popChains <- matrix(1:nChains,ncol=nPops)
	resPops <- lapply(1:nPops,function(iPop){subChains(thinned,iPops=iPop)})
	popCrit <- lapply( resPops, checkConvergenceGelman, burninFrac=0, rHatMin=rHatMin )
	#if( !all(unlist(popCrit)) ) return(FALSE)
	
	#stop(".checkConvergenceGelmanPops not implemented yet.")
	#tmp <- do.call(combinePops,resPops) #does not stack within population
	tmp <- stackChainsPop( thin(thinned,start=ceiling(getNGen(thinned)*burninFrac)), varInRows=TRUE)
	#tmp <- stackChains(thin(res,start=ceiling(getNGen(res)*burninFrac)))

	l <- dim(thinned$parms)[2]
	res2 <- thinned$parms[ ,ceiling(l*burninFrac):l, ]	# second half of all the chains
	n <- dim(res2)[2]	# number of steps
	m <- dim(res2)[3]   # number of chains
	rl <- sapply( 1:dim(res2)[1], function(vn){	#over all variables (estimands)
			muyj <- sapply(1:m, function(j){ mean(res2[vn,,j]) }) # mean of each chain
			muy <- mean(muyj)	# mean across chains
			B <- n/(m-1)*sum((muyj-muy)^2) # between chain variance
			s2j <- sapply(1:m, function(j){sum((res2[vn,,j]-muyj[j])^2) })/(n-1) #squared sums of each chain
			W <- sum(s2j)/m	# within chain variance
			VarPlus=(n-1)/n*W + 1/n*B
			#c( B=B, W=W, VarPlus=VarPlus, Rhat=round(sqrt(VarPlus/W),3), effSize=min(n*m,n*m*VarPlus/B) )
			VarPlus/W
		})
	names(rl) <- rownames(res2)
	thinned <- all(rl <= (rHatMin)^2 )
	attr(thinned,"rHat2") <- rl
	thinned
	### all rl <= criterion for each chain
}
#mtrace( checkConvergenceGelmanPops )
attr(checkConvergenceGelmanPops,"ex") <- function(){
	data(twdemcEx1)
	aTwDEMC <- twdemcEx1
	checkConvergenceGelmanPops(twdemcEx1)
}

getRLogDenQuantile <- function(
	### Quantile of logDensity below which models are significantly different from the best model, i.e. parameterization
	stackedSample				##<< numeric matrix: first column logDensity, other columns free parameters, see \code{\link{stackChains.twDEMC}}
	,maxLogDen=max(stackedSample[,1])	##<< maximum logDen Density
	,df=ncol(stackedSample)-1	##<< degress of freedom: number of fitted parameters
	,perc=0.95					##<< percentile of significance
){
	##seealso<<  
	## \code{\link{checkConvergenceGelman}}
	## \code{\link{twDEMCInt}}
	
	##details<<
	## See Hilborn97 for explanation of Density ratio test for nested models.
	x2 <- qchisq(perc, df=df )
	maxLogDen -x2/2
	### numeric scalar: minimum LogDensity below which models are significantly different 
}

checkConvergenceTrend <- function(
	### checks whether the first and last fifth mean of populations differ significantly 
	resB			##<< the twDEMC to examine
	, iChains = rep(1:ncol(resB$temp), each=ncol(resB$rLogDen)%/%ncol(resB$temp))
){
	iGen <- structure( cbind( floor(c(0,1/5)*nrow(resB$rLogDen))+1, floor(c(4/5,1)*nrow(resB$rLogDen)) )
		,dimnames= list(pop=NULL,part=c("start","end")) )
	nPops <- getNPops(resB)
	# stack rLogDen for each population
	#rLogDenStart <- popApplyTwDEMC( resB$rLogDen[iGen[,1],], nPops=1, as.vector )
	rLogDenStart <- popApplyTwDEMC( resB$rLogDen[iGen[,1],], nPops=nPops, as.vector )
	rLogDenEnd <- popApplyTwDEMC( resB$rLogDen[iGen[,2],], nPops=nPops, as.vector )
	res <- sapply( 1:nPops, function(i){ t.test(rLogDenStart[,i],rLogDenEnd[,i],alternative="greater")$p.value })
	res
}
attr(checkConvergenceTrend,"ex") <- function(){
	data(twdemcEx1)
	# p.value for differenc in means for each population
	(res <- checkConvergenceTrend(twdemcEx1))
	# none of them has converged on a 5% level
	res < 0.05
	#discard burnin 
	(res <- checkConvergenceTrend(thin(twdemcEx1,start=40)))
	# none of them has converged on a 5% level
	res < 0.05
}

