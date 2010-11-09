.setUp <-function () {
	data(delta14Catm)
	data(Howland14C)
	data(HowlandParameterPriors)
}

.tearDown <- function () {}


test.ofHowland <- function(){
	parms0 <- HowlandParameterPriors$parms0
	poptnames <- c("kY","kO")
	poptDistr <- twConstrainPoptDistr(poptnames, HowlandParameterPriors$parDistr)
	model <- list(
		modMeta=modMetaICBM1()
		,fInitState=initState.howland.ICBM1SteadyState
		,fSolve=solveICBM1
	)
	#mtrace(derivICBM1)
	#mtrace(of.howlandSteadyRootConstr)
	resOf <- of.howlandSteadyRootConstr(poptDistr$mu,model=model,poptDistr=poptDistr)
	res <- attr(resOf,"out")
	#colnames(res)
	matplot(res[,"time"], res[,c("inputLeaf_c12","inputLeaf_c14","inputRoot_c12","inputRoot_c14")], type="l" )
	matplot(res[,"time"], res[,c("Y_c12","Y_c14","O_c12","O_c14","cStock")], type="l" )
	matplot(res[,"time"], res[,c("respY_c12","respY_c14","respO_c12","respO_c14")], type="l" )
	matplot(res[,"time"], res[,c("F14C_Y","F14C_O","F14CT","respF14CT")], type="l" )
	lines( delta2iR14C(delta14C)/model$modMeta$iR14CStandard ~ yr, data=delta14Catm, col="gray")
	
	testClusterSetup <- function(){
		# see wheter all necessary information is loaded in remote processes
		argsFLogLik <- list(
			model=model
			,poptDistr=poptDistr		
			,obs=Howland14C$obsNutrientSite
			,input=Howland14C$litter
			, parms=HowlandParameterPriors$parms0           ##<< default parameters to the model
			, fCalcBiasedInput=meanInput    ##<< function(input,padj,...){obs} possibility to account for bias and to optimize bias parameters 
		#, fTransOrigPopt=transOrigPopt.default  ##<< function that translates parameters from normal to original scale
		)
		resCl <- do.call( sfClusterCall, c(list(of.howlandSteadyRootConstr,poptDistr$mu),argsFLogLik) )
		checkEqualsNumeric( resOf, resCl[[1]] )
	}
	
	profile.f <- function(){
		Rprof()
		for( i in 1:10 ) of.howlandSteadyRootConstr(poptDistr$mu,model=model,poptDistr=poptDistr, useRImpl=TRUE)
		Rprof(NULL)
		head(summaryRprof()$by.self, n=12)
		
		system.time({Rprof()
		lapply(1:500, function(x){of.howlandSteadyRootConstr(poptDistr$mu,model=model,poptDistr=poptDistr)})
		Rprof(NULL)})
		head(summaryRprof()$by.self, n=12)
	
		#using sfRemoteWrapper and exporting before 
		argsFLogLik <- list(
			remoteFun=of.howlandSteadyRootConstr
			,model=model
			,poptDistr=poptDistr		
			,obs=Howland14C$obsNutrientSite
			,input=Howland14C$litter
			, parms=HowlandParameterPriors$parms0           ##<< default parameters to the model
			, fCalcBiasedInput=meanInput    ##<< function(input,padj,...){obs} possibility to account for bias and to optimize bias parameters 
		)
		sfExport("argsFLogLik")
		resCl <- sfClusterCall( sfRemoteWrapper, normpopt=poptDistr$mu, remoteFunArgs=substitute(argsFLogLik) )
		system.time({Rprof()
				sfLapply(1:500, sfRemoteWrapper, normpopt=poptDistr$mu, remoteFunArgs=substitute(argsFLogLik) )
				Rprof(NULL)})
		head(summaryRprof()$by.self, n=12)
	}
	
	
}






