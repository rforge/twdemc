.setUp <- function(){
	data(den2dCorTwDEMC)
	#attach(twLinreg1)
	#attach(twdemcEx1)
}

.tearDown <- function(){
	#detach( twdemcEx1 )
	#detach( twLinreg1 )
}

test.findSplitBreakEarly <- function(){
	ss1 <- ss <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
	#mtrace(findSplit)
	# a is the first split
	(res <- res1 <- findSplit(ss1))
	checkEquals(res$varName, "a", msg="expected to find split in variable a")
	checkInterval( res$split, -0.9, -0.6, msg="split in a not in expected interval" )	
	checkEquals( res$perc, 0.8, msg="split in a not in expected percentile" )	
	checkEquals( res$jVars, c(1,2), msg="ordering of variables should not change." )
	
	# here b is the first found split
	(res <- res2 <- findSplit(ss1, iVars=c("b","a"), rVarCrit=2 ))
	checkEquals(res$varName, "b", msg="expected to find split in variable a")
	checkEquals( res$perc, 0.8, msg="split in a not in expected percentile" )	
	checkEquals( res$iVars, c(2,1), msg="ordering of variables should not change." )	
	checkEquals( res$jVars, c(1,2), msg="ordering of variables should not change." )	
	
	# correlation: b is the first found split
	(res <- res1 <- findSplit(ss1, rVarCrit=Inf))
	checkEquals(res$varName, "b", msg="expected to find split in variable a")
	checkEquals( res$perc, 0.6, msg="split in a not in expected percentile" )	
	checkEquals( res$jVars, c(1,2), msg="ordering of variables should not change." )
	
}

test.findSplitNonBreakEarly <- function(){
	ss1 <- ss <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
	#mtrace(findSplit)
	(res <- res1 <- findSplit(ss1, isBreakEarly=FALSE))
	checkEquals(res$varName, "a", msg="expected to find split in variable a")
	checkInterval( res$split, -0.9, -0.6, msg="split in a not in expected interval" )	
	checkEquals( res$perc, 0.8, msg="split in a not in expected percentile" )	
	checkEquals( res$jVars, c(2,1), msg="ordering of should have changed." )
	
	# here b is found first, but not broken early - so better break in a is found
	(res <- res2 <- findSplit(ss1, iVars=c("b","a"), rVarCrit=2, isBreakEarly=FALSE ))
	checkEquals(res$varName, "a", msg="expected to find split in variable a")
	checkEquals( res$perc, 0.8, msg="split in a not in expected percentile" )	
	checkEquals( res$iVars, c(1,2), msg="ordering of should have changed." )	
	checkEquals( res$jVars, c(2,1), msg="ordering of should have changed." )
	
	# b is the first changed slope break point found, here in addition the order of variables should change
	(res <- res1 <- findSplit(ss1, rVarCrit=Inf, isBreakEarly=FALSE))
	checkEquals(res$varName, "b", msg="expected to find split in variable a")
	checkEquals( res$perc, 0.6, msg="split in a not in expected percentile" )	
	checkEquals( res$iVars, c(2,1), msg="ordering of variables should change." )
	
}

