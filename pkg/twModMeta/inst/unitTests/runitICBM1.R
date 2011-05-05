.setUp <-function () {
	data(delta14Catm)
}

.tearDown <- function () {}


test.calcLagged14CSeries <- function(){
	input <- data.frame(yr=1991:1996,leaf=1:6)
	res <- calcLagged14CSeries(inputYr=input$yr, inputValue=input$leaf, delta14Catm=delta14Catm, lag=0)
	checkEqualsNumeric( delta2iR14C(delta14Catm$delta14C[ delta14Catm$yr==1991])*1, res[1])
	checkEqualsNumeric( delta2iR14C(delta14Catm$delta14C[ delta14Catm$yr==1996])*6, res[6])
	
	#delta14Catm[ delta14Catm$yr %in% 1991:1999, ]
	#mtrace(.calcLagged14CSeries)
	res <- calcLagged14CSeries(inputYr=input$yr, inputValue=input$leaf, delta14Catm=delta14Catm, lag=2)
	checkEqualsNumeric( delta2iR14C(delta14Catm$delta14C[ delta14Catm$yr==1991-2])*1, res[1])
	checkEqualsNumeric( delta2iR14C(delta14Catm$delta14C[ delta14Catm$yr==1996-2])*6, res[6])
	
	res <- calcLagged14CSeries(inputYr=input$yr, inputValue=input$leaf, delta14Catm=delta14Catm, lag=2, iR14CStandard=1)
	checkEqualsNumeric( delta2iR14C(delta14Catm$delta14C[ delta14Catm$yr==1991-2])*1/c14Constants$iR14CStandard, res[1])
	checkEqualsNumeric( delta2iR14C(delta14Catm$delta14C[ delta14Catm$yr==1996-2])*6/c14Constants$iR14CStandard, res[6])
}

test.decay <- function(){
	# unit is tC/ha
	mm <- modMetaICBMDemo()
	parms0 <- list(
		tLagLeaf=1	# other is recent-C, which is not accounted in the Howland study
		,tLagRoot=5
		,kY=0.8
		,kO=0.006
		,h=0.2
		,cY=0.2		
	)
	input0 <- list(
		leaf = cbind( yr=c(1900,2000), obs=0 )
		,root = cbind( yr=c(1900,2000), obs=0 )
	)
	
	yr0 <- 1950
	iRNew <- delta2iR14C(delta14Catm$delta14C[delta14Catm$yr==yr0])
	tvrOld <- 1000 #1/parms0$kO	#1000 yr old carbon
	iROld <- decayIR14C( yr=yr0, iR0=delta2iR14C(delta14Catm$delta14C[1]), yr0=1950-tvrOld )	# near 1 (standard of old wood)
	
	#mtrace(initStateModMeta)
	x0 <- initStateICBMDemo( xc12=101.0*c(parms0$cY,(1-parms0$cY)),iR=matrix(c(iRNew,iROld),ncol=1,dimnames=list(NULL,"c14")) )
	
	#mtrace(solveICBMDemo)
	#mtrace(derivICBMDemo)
	res <- solveICBMDemo(
		x0=x0,	times=yr0:2007	
		,parms=parms0
		,input=input0
		,useRImpl=TRUE
	)
	#colnames(res)
	matplot(res[,"time"], res[,c("Y_c12","Y_c14","O_c12","O_c14","R_c12","R_c14")], type="l" )
	matplot(res[,"time"], res[,c("respY_c12","respY_c14","respO_c12","respO_c14")], type="l" )
	matplot(res[,"time"], res[,c("F14C_Y","F14C_O","respF14CT")], type="l" )
	#matplot(res[,"time"], res[,c("respF14C_Y","respF14C_O")], type="l" ) # same as pool
	
	# test if the DLL-compiled version gives the same results
	res2 <- solveICBMDemo(
		x0=x0,	times=yr0:2007	
		,parms=parms0
		,input=input0
		,useRImpl=FALSE
	)
	checkTrue( max(res2-res, na.rm=TRUE) < 1e-8 )		
	
}

test.constInput <- function(){
	# unit is tC/ha
	mm <- modMetaICBMDemo()
	parms0 <- list(
		tLagLeaf=1	# other is recent-C, which is not accounted in the Howland study
		,tLagRoot=5
		,kY=0.5
		,kO=0.006
		,h=0.2
		,cY=0.2		
	)
	yr0 <- 1950
	yrEnd <- 2007
	times <- yr0:yrEnd
	input1 <- list(
		leaf = cbind( yr=c(1900,2000), obs=2.3 )
		,root = cbind( yr=c(1900,2000), obs=2.5 )
	)
	
	iRNew <- delta2iR14C(delta14Catm$delta14C[delta14Catm$yr==yr0])
	tvrOld <- 1000 #1/parms0$kO	#1000 yr old carbon
	iROld <- decayIR14C( yr=yr0, iR0=delta2iR14C(delta14Catm$delta14C[1]), yr0=1950-tvrOld )	# near 1 (standard of old wood)
	
	#mtrace(initStateModMeta)
	x0 <- initStateICBMDemo( xc12=101.0*c(parms0$cY,(1-parms0$cY)),iR=matrix(c(iRNew,iROld),ncol=1,dimnames=list(NULL,"c14")) )
	
	#mtrace(derivICBMDemo)
	res <- solveICBMDemo(
		x0=x0,	times=times
		,parms=parms0
		,input=input1
	)
	#colnames(res)
	matplot(res[,"time"], res[,c("inputLeaf_c12","inputLeaf_c14","inputRoot_c12","inputRoot_c14")], type="l" )
	matplot(res[,"time"], res[,c("Y_c12","Y_c14","O_c12","O_c14")], type="l" )
	matplot(res[,"time"], res[,c("respY_c12","respY_c14","respO_c12","respO_c14")], type="l" )
	matplot(res[,"time"], res[,c("F14C_Y","F14C_O","respF14CT")], type="l" )
	lines( delta2iR14C(delta14C)/mm$iR14CStandard ~ yr, data=delta14Catm, col="blue")	
}





