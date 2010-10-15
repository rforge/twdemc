.setUp <-function () {
	data(delta14Catm)
	data(Howland14C)
	data(HowlandParameterPriors)
}

.tearDown <- function () {}


test.ofHowland <- function(){
	parms0 <- HowlandParameterPriors$parms0
	poptnames <- c("h","cY")
	poptDistr <- twConstrainPoptDistr(poptnames, HowlandParameterPriors$parDistr)
	model <- list(
		modMeta=modMeta.ICBM1()
		,fInitState=initState.howland.ICBM1SteadyState
		,fSolve=solve.ICBM1
	)
	#mtrace(deriv.ICBM1)
	#mtrace(of.howlandSteady)
	resOf <- of.howlandSteady(poptDistr$mu,model=model,poptDistr=poptDistr)
	res <- attr(resOf,"out")
	#colnames(res)
	matplot(res[,"time"], res[,c("inputLeaf_c12","inputLeaf_c14","inputRoot_c12","inputRoot_c14")], type="l" )
	matplot(res[,"time"], res[,c("Y_c12","Y_c14","O_c12","O_c14","cStock")], type="l" )
	matplot(res[,"time"], res[,c("respY_c12","respY_c14","respO_c12","respO_c14")], type="l" )
	matplot(res[,"time"], res[,c("F14C_Y","F14C_O","F14CT","respF14CT")], type="l" )
	lines( delta2iR14C(delta14C)/mm$iR14CStandard ~ yr, data=delta14Catm, col="blue")
	
	profile.f <- function(){
		Rprof()
		for( i in 1:10 ) of.howlandSteady(poptDistr$mu,model=model,poptDistr=poptDistr, useRImpl=TRUE)
		Rprof(NULL)
		head(summaryRprof()$by.self, n=12)
		
		system.time({Rprof()
		for( i in 1:500 ) of.howlandSteady(poptDistr$mu,model=model,poptDistr=poptDistr)
		Rprof(NULL)})
		head(summaryRprof()$by.self, n=12)
	}
	
	resCl <- sfClusterCall(of.howlandSteady,poptDistr$mu,model=model,poptDistr=poptDistr)
	
}






