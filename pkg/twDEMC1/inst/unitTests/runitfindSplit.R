.setUp <- function(){
	data(den2dCorTwDEMC)
	#attach(twLinreg1)
	#attach(twdemcEx1)
}

.tearDown <- function(){
	#detach( twdemcEx1 )
	#detach( twLinreg1 )
}

test.findSplit <- function(){
	ss1 <- ss <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
	#mtrace(findSplit)
	# a is the first split
	(res <- res1 <- findSplit(ss1))
	checkEquals(res$varName, "a", msg="expected to find split in variable a")
	checkInterval( res$split, -0.9, -0.6, msg="split in a not in expected interval" )	
	checkEquals( res$perc, 0.8, msg="split in a not in expected percentile" )	
	checkEquals( res$resD$iVar, c(1,2), msg="expected same order of resD$iVar as argument iVars all variables back." )
	checkEquals( res$resD$jPVar, c(2,1), msg="ordering resD$jPVar should have changed." )
	
	# here a is the first found split but order or resD is flipped
	(res <- res2 <- findSplit(ss1, iVars=c("b","a"), rVarCrit=2 ))
	checkEquals(res$varName, "a", msg="expected to find split in variable a")
	checkInterval( res$split, -0.9, -0.6, msg="split in a not in expected interval" )	
	checkEquals( res$perc, 0.8, msg="split in a not in expected percentile" )	
	checkEquals( res$resD$iVar, c(2,1), msg="expected same order of resD$iVar as argument iVars all variables back." )
	checkEquals( res$resD$jPVar, c(1,2), msg="ordering resD$jPVar should have changed." )
	
	# correlation: b is the first found split
	(res <- res1 <- findSplit(ss1, rVarCrit=Inf))
	checkEquals(res$varName, "b", msg="expected to find split in variable a")
	checkEquals( res$perc, 0.2, msg="split in a not in expected percentile" )	
	checkTrue( res$resD$dAlphaSlope[2] > 1.2, msg="expected high angle between slopes" )	
}

test.findSplitBreakEarly <- function(){
	ss1 <- ss <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
	#mtrace(findSplit)
	res <- res1 <- findSplit(ss1, rVarCrit=Inf)	# calculate slopeAngles
	res <- res2 <- findSplit(ss1, rVarCrit=Inf, isBreakEarly=TRUE, checkSlopesFirst=res1$resD )
	checkEquals(res1,res2,"isBreakEarly should with two variables give the same")
}

