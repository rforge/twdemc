checkConvergenceGelman <- function(
	### Gelman RHat criterion applied to twDEMC result assuming all chains one population.
	res					##<< see return value of \code{\link{twDEMCInt}} ($parms (d x nStep x nChain) )
	,burninFrac=0.5 	##<< fraction of the chain to be discarded (default 0.5)
	,rHatMin = 1.1		##<< rHat criterion, upper bound that is regarded as convergence 
){
	# checkConvergenceGelman
	##details<< 
	## see Gelman04 (twutz:Gelman04_3#Inference_and_assessing_convergence)
	l <- dim(res$parms)[2]
	res2 <- res$parms[ ,max(1,ceiling(l*burninFrac)):l, ]	# later part of all the chains
	n <- dim(res2)[2]	# number of steps
	m <- dim(res2)[3]   # number of chains
	rl <- sqrt(sapply( 1:dim(res2)[1], function(vn){	#over all variables (estimands)
			muyj <- sapply(1:m, function(j){ mean(res2[vn,,j]) }) # mean of each chain
			muy <- mean(muyj)	# mean across chains
			B <- n/(m-1)*sum((muyj-muy)^2) # between chain variance
			s2j <- sapply(1:m, function(j){sum((res2[vn,,j]-muyj[j])^2) })/(n-1) #squared sums of each chain
			W <- sum(s2j)/m	# within chain variance
			VarPlus=(n-1)/n*W + 1/n*B
			#c( B=B, W=W, VarPlus=VarPlus, Rhat=round(sqrt(VarPlus/W),3), effSize=min(n*m,n*m*VarPlus/B) )
			VarPlus/W
		}))
	names(rl) <- rownames(res2)
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
