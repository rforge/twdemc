.setUp <- function(){
	data(den2dCorEx)
	#attach(twLinreg1)
	#attach(twdemcEx1)
}

.tearDown <- function(){
	#detach( twdemcEx1 )
	#detach( twLinreg1 )
}

test.findSplit <- function(){
	ss1 <- ss <- stackChains(thin(den2dCorEx$mcBulk, start=300))[,-1]
	#mtrace(findSplit)
	# a is the first split
	(res <- res1 <- findSplit(ss1))
	checkEquals(res$varName, "a", msg="res1: expected to find split in variable a")
	checkInterval( res$split, -1.2, -0.4, msg=paste("res1: split in a not in expected interval, res$split=",res$split,sep="") )	
	checkEquals( res$perc, 0.8, msg="res1: split in a not in expected percentile" )	
	checkEquals( res$resD$iVar, c(1,2), msg="res1: expected same order of resD$iVar as argument iVars all variables back." )
	checkEquals( res$resD$jPVar, c(2,1), msg="res1: ordering resD$jPVar should have changed." )
	
	# here a is the first found split but order or resD is flipped
	(res <- res2 <- findSplit(ss1, iVars=c("b","a"), rVarCrit=2 ))
	checkEquals(res$varName, "a", msg="res2: expected to find split in variable a")
	checkInterval( res$split, -1.2, -0.5, msg=paste("res2:split in a not in expected interval res$split=",res$split,sep="") )	
	checkEquals( res$perc, 0.8, msg="res2:split in a not in expected percentile" )	
	checkEquals( res$resD$iVar, c(2,1), msg="res2:expected same order of resD$iVar as argument iVars all variables back." )
	checkEquals( res$resD$jPVar, c(1,2), msg="res2:ordering resD$jPVar should have changed." )
	
	# correlation: b is the first found split
	(res <- res3 <- findSplit(ss1, rVarCrit=Inf))
	checkEquals(res$varName, "b", msg="res3: expected to find split in variable b")
	#checkEquals( res$perc, 0.6, msg="split in b in a not in expected percentile" )	
	checkTrue( res$resD$dAlphaSlope[2] > 1.0, msg="res3: expected high angle between slopes" )	
}

test.findSplitBreakEarly <- function(){
	ss1 <- ss <- stackChains(thin(den2dCorEx$mcBulk, start=300))[,-1]
	#mtrace(findSplit)
	res <- res1 <- findSplit(ss1, rVarCrit=Inf)	# calculate slopeAngles
	res <- res2 <- findSplit(ss1, rVarCrit=Inf, isBreakEarly=TRUE, checkSlopesFirst=res1$resD )
	checkEquals(res1,res2,"isBreakEarly should with two variables give the same")
}

