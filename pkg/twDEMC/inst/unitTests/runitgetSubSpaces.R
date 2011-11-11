.setUp <- function(){
	data(den2dCorDivideTwDEMC)
	data(den2dCorTwDEMC)
	#attach(twLinreg1)
	#attach(twdemcEx1)
}

.tearDown <- function(){
	#detach( twdemcEx1 )
	#detach( twLinreg1 )
}

test.getSubSpaces04 <- function(){
	aSample <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
	#mtrace(getSubSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.4 )
	str(subSpaces)
	checkEquals( length(subSpaces$spaces),2,"expected two subspaces with minPSub=0.4")
}

test.getSubSpaces005 <- function(){
	aSample <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
	#mtrace(getSubSpaces)
	subSpaces <- getSubSpaces(aSample, minPSub=0.05 )
	str(subSpaces)
	checkTrue( length(subSpaces$spaces) > 2," with minPSub=0.05 there should be more than 2 subspaces")
	checkTrue( all( sapply(subSpaces$spaces, "[[", "pSub") >= 0.05), " percentiles of all subspaces should be >= 0.05")
}
