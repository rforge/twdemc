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
	subSpaces <- getSubSpaces(aSample, minPSub=0.05 )
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
	subSpaces <- getSubSpaces(aSample, minPSub=0.05, splits=pop0$splits, upperParBounds=pop0$upperParBounds, lowerParBounds=pop0$lowerParBounds )
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
	splitsPop <- lapply(subSpaces,)
}

