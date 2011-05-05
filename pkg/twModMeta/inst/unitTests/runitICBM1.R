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

test.calcSteadyHcY_ICBM1 <- function(){
	p0 <- list(
		kY=1
		,kO=1/30
		,iY=480
		,Ctot=1010
		,h=NA
		,cY=NA
	)
	p1 <- with(p0, {calcSteadyHcY_ICBM1(Ctot=Ctot,iY=iY,parms=p0)})
	#p0[c("h","cY")]
	#attach(p1)
	with(p1,{
		checkEqualsNumeric( iY, kY*cY*Ctot )	
		checkEqualsNumeric( h*iY, kO*(1-cY)*Ctot )	
	})
	#detach()
	p2 <- with(p1, {calcSteadyK_ICBM1(Ctot=Ctot,iY=iY,parms=p1)})
	checkEqualsNumeric( unlist(p1[c("kY","kO")]), unlist(p2[c("kY","kO")]) )
}

test.decay <- function(){
	# unit is tC/ha
	mm <- modMetaICBM1()
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
	
	#mtrace(initStateSoilMod)
	x0 <- initStateICBM1( xc12=101.0*c(parms0$cY,(1-parms0$cY)),iR=matrix(c(iRNew,iROld),ncol=1,dimnames=list(NULL,"c14")) )
	
	#mtrace(solveICBM1)
	#mtrace(derivICBM1)
	res <- solveICBM1(
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
	
	res <- solveICBM1(
		x0=x0,	times=yr0:2007	
		,parms=parms0
		,input=input0
		,useRImpl=FALSE
	)
	
}

test.constInput <- function(){
	# unit is tC/ha
	mm <- modMetaICBM1()
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
	
	#mtrace(initStateSoilMod)
	x0 <- initStateICBM1( xc12=101.0*c(parms0$cY,(1-parms0$cY)),iR=matrix(c(iRNew,iROld),ncol=1,dimnames=list(NULL,"c14")) )
	
	#mtrace(derivICBM1)
	res <- solveICBM1(
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



test.steadyStateHowland <- function(){
	data(Howland14C)
	
	# unit is tC/ha
	mm <- modMetaICBM1()
	parms0 <- list(
		tLagLeaf=1	# other is recent-C, which is not accounted in the Howland study
		,tLagRoot=5
		,kY=2.4
		,kO=0.1
		,h=0.2
		,cY=0.2		
	)
	obs <- Howland14C$obsNutrientSite
	Ctot <- obs$somStock[1,2]
	input <- lapply(Howland14C$litter, "[", ,1:2, drop=FALSE)
	sumInput <- sum(sapply(input,"[",2))
	parms0 <- calcSteadyHcY_ICBM1(Ctot=Ctot,iY=sumInput,parms=parms0)
	yr0 <- 1940
	yrEnd <- 2007
	times <- yr0:yrEnd
	
	iRNew <- delta2iR14C(delta14Catm$delta14C[delta14Catm$yr==yr0])
	tvrOld <- 1/parms0$kO	#1000 yr old carbon
	iROld <- decayIR14C( yr=yr0, iR0=delta2iR14C(delta14Catm$delta14C[1]), yr0=1950-tvrOld )	# near 1 (standard of old wood)
	
	#mtrace(initStateSoilMod)
	x0 <- initStateICBM1( xc12=as.vector(Ctot*c(parms0$cY,(1-parms0$cY))),iR=matrix(c(iRNew,iROld),ncol=1,dimnames=list(NULL,"c14")) )
	
	#mtrace(derivICBM1)
	#mtrace(solveICBM1)
	system.time(res <- resR <- resSolve <- solveICBM1(
		x0=x0,	times=times
		,parms=parms0
		,input=input
		,useRImpl=TRUE
	))
	#colnames(res)
	matplot(res[,"time"], res[,c("inputLeaf_c12","inputLeaf_c14","inputRoot_c12","inputRoot_c14")], type="l" )
	matplot(res[,"time"], res[,c("Y_c12","Y_c14","O_c12","O_c14","cStock")], type="l" )
	matplot(res[,"time"], res[,c("respY_c12","respY_c14","respO_c12","respO_c14")], type="l" )
	matplot(res[,"time"], res[,c("respF14CT", "F14CT", "F14C_Y","F14C_O")], type="l" )
	lines( delta2iR14C(delta14C)/mm$iR14CStandard ~ yr, data=delta14Catm, col="gray")

	#must be steady state
	checkEqualsNumeric( rep(0, nrow(res)-1), diff(res[,"cStock"]) )
	
	#mtrace(solveICBM1)
	system.time(res <- resc <- solveICBM1(
		x0=x0,	times=yr0:2007	
		,parms=parms0
		,input=input
		,useRImpl=FALSE
	))
	checkEqualsNumeric( resR, resc, tol=1e-5 )
	
	profile.f <- function(){
		Rprof()
		for( i in 1:10 ) solveICBM1(
				x0=x0,	times=times
				,parms=parms0
				,input=input
				,useRImpl=TRUE
		)
		Rprof(NULL)
		head(summaryRprof()$by.self, n=12)
		
		system.time({ Rprof()
		for( i in 1:500 ) solveICBM1(
				x0=x0,	times=times
				,parms=parms0
				,input=input
				,useRImpl=FALSE
			)
		Rprof(NULL)})
		head(summaryRprof()$by.self, n=12)
		#speedup of about 50x
	}
	
}


