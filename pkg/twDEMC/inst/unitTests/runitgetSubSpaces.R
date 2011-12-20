.setUp <- function(){
	data(den2dCorEx)
	#attach(twLinreg1)
	#attach(twdemcEx1)
}

.tearDown <- function(){
	#detach( twdemcEx1 )
	#detach( twLinreg1 )
}

test.getSubSpaces04 <- function(){
	aSample <- stackChains(thin(concatPops(den2dCorEx$mcBulk), start=300))[,-1]
	#mtrace(getSubSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.4 )
	str(subSpaces)
	checkEquals( length(subSpaces$spaces),2,"expected two subspaces with minPSub=0.4")
}

test.getSubSpaces005 <- function(){
	aSample <- stackChains(thin(concatPops(den2dCorEx$mcBulk), start=300))[,-1]
	#mtrace(getSubSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.05, verbose=TRUE )
	str(subSpaces)
	checkTrue( length(subSpaces$spaces) > 2," with minPSub=0.05 there should be more than 2 subspaces")
	checkTrue( all( sapply(subSpaces$spaces, "[[", "pSub") >= 0.05), " percentiles of all subspaces should be >= 0.05")
}

test.getSubSpacesOfSubspace <- function(){
	aSample0 <- stackChains(thin(concatPops(den2dCorEx$mcBulk), start=300))[,-1]
	subSpaces0 <- getSubSpaces(aSample0, minPSub=0.4 )
	pop0 <- subSpaces0$spaces[[1]]
	aSample <- pop0$sample
	#mtrace(getSubSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.2, splits=pop0$splits, upperParBounds=pop0$upperParBounds, lowerParBounds=pop0$lowerParBounds )
	str(subSpaces)
	#pop <- subSpaces$spaces[[1]]
	for( pop in subSpaces$spaces ){
		checkEquals( pop0$splits , pop$splits[ 1:length(pop0$splits)])
		#pName <- names(pop0$upperParBounds)[1]
		for( pName in names(pop0$upperParBounds) ){
			checkTrue( pop0$upperParBounds[pName] >= pop$upperParBounds[pName] )
		}
		#pName <- names(pop0$lowerParBounds)[1]
		for( pName in names(pop0$lowerParBounds) ){
			checkTrue( pop0$lowerParBounds[pName] <= pop$lowerParBounds[pName] )
		}
	}
	(splitsPop <- lapply(subSpaces$spaces,"[[","splits"))
}

test.splitMerge <- function(){
	mc0 <- subPops(den2dCorEx$mcBulk,iSpace=1)	# one population
	n0 <- getNSamples(mc0)
	
	#------- split the single population along a
	split1 <- c(a=-2)
	spaces1 <- list(
		low1 = list(
			splits = split1
			,upperParBounds=split1
		)
		,upp1 = list(
			splits = split1
			,lowerParBounds=split1
		)
	)
	mc1 <- mc0
	mc1$pops <- divideTwDEMCPop(mc0$pops[[1]], spaces1)
	#.getParBoundsPops(mc1$pops)
	checkEquals( 2, length(mc1$pops) )
	checkEquals( spaces1$low1, mc1$pops[[1]][ names(spaces1$low1) ]  )
	checkEquals( 0,  length(mc1$pops[[1]]$lowerParBounds))
	checkEquals( spaces1$upp1, mc1$pops[[2]][ names(spaces1$upp1) ]  )
	checkEquals( 0,  length(mc1$pops[[2]]$upperParBounds))
	# merge again
	.nS <- getNSamples(mc1)
	#c( n0, sum(.nS) )
	pSubs1 <- .nS/sum(.nS)
	resMerge <- .mergePopTwDEMC( mc1$pops, 1, pSubs1 )
	checkEquals( 1, length(resMerge$pops) )
	checkEquals( 0,  length(resMerge$pops[[1]]$upperParBounds))
	checkEquals( 0,  length(resMerge$pops[[1]]$lowerParBounds))
	
	#----- split the upper population along a again 
	# merge pop 1 with pop2 only
	split2 <- c(a=-1)
	spaces2 <- list(
		low2 = list(
			splits = split1
			,upperParBounds=split1
		)
		,upp2Low = list(
			splits = c(split1,split2)
			,lowerParBounds=c(split1)
			,upperParBounds=c(split2)
		)
		,upp2Upp = list(
			splits = c(split1,split2)
			,lowerParBounds=c(split2)
		)
	)
	mc2 <- mc0
	mc2$pops <- divideTwDEMCPop(mc0$pops[[1]], spaces2)
	#.getParBoundsPops(mc2$pops)
	checkEquals( 3, length(mc2$pops) )
	checkEquals( spaces2$low2, mc2$pops[[1]][ names(spaces2$low2) ]  )
	checkEquals( 0,  length(mc2$pops[[1]]$lowerParBounds))
	checkEquals( spaces2$upp2Low, mc2$pops[[2]][ names(spaces2$upp2Low) ]  )
	checkEquals( spaces2$upp2Upp, mc2$pops[[3]][ names(spaces2$upp2Upp) ]  )
	checkEquals( 0,  length(mc2$pops[[3]]$upperParBounds))
	#-- do the same by starting from upp1
	mc2b <- mc1
	newPops <- divideTwDEMCPop(mc1$pops[[2]], spaces2[c("upp2Low","upp2Upp")])
	mc2b$pops <- c(mc1$pops[-2], newPops )
	#.getParBoundsPops(newPops)
	checkEquals( .getParBoundsPops(mc2$pops),  .getParBoundsPops(mc2b$pops))
	# merge first to second pop
	.nS <- getNSamples(mc2)
	#c(n0, sum(.nS) )		# may have lost a few samples in subSpacing
	pSubs2 <- .nS/sum(.nS)
	resMerge <- .mergePopTwDEMC( mc2$pops, 1, pSubs2 )
	#.getParBoundsPops(resMerge$pops)
	checkEquals( 2, length(resMerge$pops) )
	checkEquals( 0,  length(resMerge$pops[[1]]$lowerParBounds))
	checkEquals( split2,  (resMerge$pops[[1]]$upperParBounds))
	checkEquals( split2,  resMerge$pops[[1]]$splits)
	checkEquals( split2,  (resMerge$pops[[2]]$lowerParBounds))
	checkEquals( 0,  length(resMerge$pops[[2]]$upperParBounds))
	checkEquals( split2,  resMerge$pops[[2]]$splits)
	
	#----- split the upp2low population along b again
	# merge splitted pop1 with 2 and 3, but not with 4
	split3 <- c(b=0)
	spaces3 <- list(
		low2 = list(
			splits = split1
			,upperParBounds=split1
		)
		,upp2Low3Low = list(
			splits = c(split1,split2,split3)
			,lowerParBounds=c(split1)
			,upperParBounds=c(split2,split3)
		)
		,upp2Low3Upp = list(
			splits = c(split1,split2,split3)
			,lowerParBounds=c(split1,split3)
			,upperParBounds=c(split2)
		)
		,upp2Upp = list(
			splits = c(split1,split2)
			,lowerParBounds=c(split2)
		)
	)
	mc3 <- mc0
	mc3$pops <- divideTwDEMCPop(mc0$pops[[1]], spaces3)
	#.getParBoundsPops(mc3$pops)
	checkEquals( 4, length(mc3$pops) )
	# do the same by starting from upp2Low
	mc3b <- mc2
	newPops <- divideTwDEMCPop(mc2$pops[[2]], spaces3[c("upp2Low3Low","upp2Low3Upp")])
	#.getParBoundsPops(newPops)
	mc3b$pops <- c(mc2$pops[1], newPops, mc2$pops[3] )
	checkEquals( .getParBoundsPops(mc3$pops),  .getParBoundsPops(mc3b$pops))
	#-- merge first to second and third pop
	.nS <- getNSamples(mc3)
	#c(n0, sum(.nS) )		# may have lost a few samples in subSpacing
	pSubs3 <- .nS/sum(.nS)
	resMerge <- .mergePopTwDEMC( mc3$pops, 1, pSubs3 )
	#.getParBoundsPops(resMerge$pops)
	checkEquals( 3, length(resMerge$pops) )
	checkEquals( 0,  length(resMerge$pops[[1]]$lowerParBounds))
	checkEquals( c(split2,split3),  (resMerge$pops[[1]]$upperParBounds))
	checkEquals( c(split2,split3),  resMerge$pops[[1]]$splits)
	checkEquals( split3,  (resMerge$pops[[2]]$lowerParBounds))
	checkEquals( split2,  (resMerge$pops[[2]]$upperParBounds))
	checkEquals( c(split2,split3),  resMerge$pops[[2]]$splits)
	checkEquals( split2,  (resMerge$pops[[3]]$lowerParBounds))
	checkEquals( 0,  length(resMerge$pops[[3]]$upperParBounds))
	checkEquals( c(split2),  resMerge$pops[[3]]$splits)
	#-- merge second to third pop (do not touch split of 1 pop)
	resMerge <- .mergePopTwDEMC( mc3$pops, 2, pSubs3 )
	#.getParBoundsPops(resMerge$pops)
	checkEquals( 3, length(resMerge$pops) )
	checkEquals( 0,  length(resMerge$pops[[1]]$lowerParBounds))
	checkEquals( c(split1),  (resMerge$pops[[1]]$upperParBounds))
	checkEquals( c(split1),  resMerge$pops[[1]]$splits)
	checkEquals( split1,  (resMerge$pops[[2]]$lowerParBounds))
	checkEquals( split2,  (resMerge$pops[[2]]$upperParBounds))
	checkEquals( c(split1,split2),  resMerge$pops[[2]]$splits)
	checkEquals( split2,  (resMerge$pops[[3]]$lowerParBounds))
	checkEquals( 0,  length(resMerge$pops[[3]]$upperParBounds))
	checkEquals( c(split1,split2),  resMerge$pops[[3]]$splits)




}

