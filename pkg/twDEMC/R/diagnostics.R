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
	## \item{ the theorectical minimum logLik-Value for significant model difference : \code{\link{getRLogLikQuantile}}  } 
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


.checkConvergenceGelmanPops <- function(
	### Gelman RHat criterion applied to each population and between populations
	res					##<< see return value of \code{\link{twDEMCInt}} ($parms (d x nStep x nChain) )
	,burninFrac=0.5 	##<< fraction of the chain to be discarded (default 0.5)
	,rHatMin = 1.1		##<< rHat criterion, upper bound that is regarded as convergence 
){
	stop(".checkConvergenceGelmanPops not implemented yet.")
	# checkConvergenceGelmanPops
	res <- thin(res,start=ceiling(calcNGen(res)*burninFrac)  )
	nPops <- getNPops(res)
	#popChains <- matrix(1:nChains,ncol=nPops)
	resPops <- lapply(1:nPops,function(iPop){subChains(res,iPops=iPop)})
	popCrit <- lapply( resPops, checkConvergenceGelman, burninFrac=0, rHatMin=rHatMin )
	if( !all(unlist(popCrit)) ) return(FALSE)
	
	#tmp <- do.call(combinePops,resPops) #does not stack within population
	tmp <- stackChains(thin(res,start=ceiling(calcNGen(res)*burninFrac)))

	l <- dim(res$parms)[2]
	res2 <- res$parms[ ,ceiling(l*burninFrac):l, ]	# second half of all the chains
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
	res <- all(rl <= (rHatMin)^2 )
	attr(res,"rHat2") <- rl
	res
	### all rl <= criterion for each chain
}
#mtrace( checkConvergenceGelman )

getRLogLikQuantile <- function(
	### Quantile of logLikelihood below which models are significantly different from the best model, i.e. parameterization
	stackedSample				##<< numeric matrix: first column log-Likelihood, other columns free parameters, see \code{\link{stackChains.twDEMC}}
	,maxLogLik=max(stackedSample[,1])	##<< maximum logLik Likelihood
	,df=ncol(stackedSample)-1	##<< degress of freedom: number of fitted parameters
	,perc=0.95					##<< percentile of significance
){
	##seealso<<  
	## \code{\link{checkConvergenceGelman}}
	## \code{\link{twDEMCInt}}
	
	##details<<
	## See Hilborn97 for explanation of Likelihood ratio test for nested models.
	x2 <- qchisq(perc, df=df )
	maxLogLik -x2/2
	### numeric scalar: minimum Log-Likelihood below which models are significantly different 
}

