# same as steady 3 but including bias in difference between respiration and litter fall (root input)

.tmp.init <- function(){
	setupClusterHowlandDev(pkgDir = ".")
	parms0 <- HowlandParameterPriors$parms0
	model <- list(
		modMeta=modMetaICBM1()
		,fInitState=initState.howland.ICBM1SteadyState
		,fSolve=solveICBM1
	)
	
	argsFLogLik <- argsFLogLikRemoteFun <- list(
		model=model
		#remoteFun=of.howlandSteadyRootConstr		# will not work with twDEMC
		#,poptDistr=poptDistr			# must set when parameters are known
		,obs=Howland14C$obsNutrientSite
		,input=Howland14C$litter
		, parms=parms0           		##<< default parameters to the model
		, fCalcBiasedInput=meanInput    ##<< function(input,padj,...){obs} possibility to account for bias and to optimize bias parameters
		, fCalcIROLayer=calcIROLayer	##<< function to calculate iRofO-Layer
		, fCalcSteadyPars=calcSteadyHcY_ICBM1
	)
	#sfExport("argsFLogLik")			# done after assigning poptDistr	
	
	windows(width=4.4,height=3.4,pointsize=10, record=TRUE)
	par( las=1 )					#also y axis labels horizontal
	par(mar=c(2.0,3.3,0,0)+0.3 )  #margins
	par(tck=0.02 )				#axe-tick length inside plots             
	par(mgp=c(1.1,0.2,0) )  #positioning of axis title, axis labels, axis
	
}






mdi.klagLitter <- function(){
	#deprecated from steady3-fit does not find the global maximum:
	#normpopt <- c(structure(c(1.26184949077747, 0.0270660705723808, 1.08291416737787,6.52418312193126), .Names = c("kY", "kO", "tLagLeaf", "tLagRoot")),biasDiffRespLitterfall=0)
	# from twDEMC later on:
	normpopt <- structure(c(0.228439444501037, -3.58235519122273, 0.0356533561112984, 1.87008926040014, 38.9385098325832), .Names = c("kY", "kO", "tLagLeaf","tLagRoot", "biasDiffRespLitterfall"))
	
	#using sfRemoteWrapper and exporting before 
	poptnames <- names(normpopt)
	poptDistr <- argsFLogLik$poptDistr <- twConstrainPoptDistr(poptnames, HowlandParameterPriors$parDistr)
	sfExport("argsFLogLik")
	
	#mtrace(of.howlandSteadyRootConstr)
	#resOf <- sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik)	
	fOpt <- function(normpopt){
		sum(sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik ))
	}
	#fOpt(normpopt)
	#resOpt <- optim(normpopt, fOpt, method="Nelder-Mead", hessian = TRUE, control=list(maxit=1000, fnscale=-1))
	resOpt <- optim(normpopt, fOpt, method="BFGS", control=list(fnscale=-1), hessian = TRUE)
	transOrigPopt(poptDistr$mu, HowlandParameterPriors$parDistr$trans[poptnames])
	transOrigPopt(normpopt, HowlandParameterPriors$parDistr$trans[poptnames])
	(tmp <- transOrigPopt(resOpt$par, HowlandParameterPriors$parDistr$trans[poptnames]))
	#copy2clip(deparse(tmp))
	#structure(c(1.31004509263524, 0.0276374254632283, 1.0829503483462,6.65177437410762, 38.9279591345289), .Names = c("kY", "kO", "tLagLeaf","tLagRoot", "biasDiffRespLitterfall"))


	
	argsFLogLik2 <- argsFLogLik
	#tmp <- argsFLogLik2$remoteFun; mtrace(tmp); argsFLogLik2$remoteFun<-tmp
	resOf <- sfRemoteWrapper( normpopt=resOpt$par, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik2)
	#resOf <- sfRemoteWrapper( normpopt=c(cY=logit(cYOpt),h=logit(hOpt)), remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik2)
	sort(resOf)
	sort(attr(resOf,"logLikParms"))
	
	res <- attr(resOf,"out")
	#colnames(res)
	matplot(res[,"time"], res[,c("inputLeaf_c12","inputLeaf_c14","inputRoot_c12","inputRoot_c14")], type="l" )
	matplot(res[,"time"], res[,c("Y_c12","Y_c14","O_c12","O_c14","cStock")], type="l", ylim=c(0,1100) )
	matplot(res[,"time"], res[,c("respY_c12","respY_c14","respO_c12","respO_c14")], type="l" )
	plotHowlandFM( res, attr(resOf,"obs"))
	
	#------ explore posterior with MCMC using prior
	tmp.fcovarPrior <- function(){
		covMat <- poptDistr$sigma    
		.nPops=3
		Zinit <- initZtwDEMCNormal( resOpt$par, covMat, nChains=4*.nPops, nPops=.nPops)
		#resMC <- twDEMCBatch( Zinit, nGen=4*5, debugSequential=TRUE, fLogLik=of.howlandSteadyRootConstr, argsFLogLik=argsFLogLik, nPops=.nPops )
		.remoteDumpfileBasename=file.path("tmp","dumpRemote")
		#resMC <- twDEMCBatch( Zinit, nGen=500, fLogLik=of.howlandSteadyRootConstr, argsFLogLik=argsFLogLik, nPops=.nPops, remoteDumpfileBasename=.remoteDumpfileBasename )
		resMC <- twDEMCBatch( Zinit, nGen=500, fLogLik=of.howlandSteadyRootConstr, argsFLogLik=argsFLogLik, nPops=.nPops )
		matplot(resMC$pAccept, type="l")
		plot(as.mcmc.list(resMC))
		resMC <- twDEMCBatch( resMC, nGen=1000, doRecordProposals=TRUE )
		#resMC <- twDEMCBatch( resMC, nGen=3000, doRecordProposals=TRUE )	# for integrating over nuisance tLagRoot and tLagLeaf
		resMCB <- thin(resMC, start=400)
		#save(resMCB,file=file.path("tmp","resMCB_steady4.RData"))
		plotThinned(as.mcmc.list(resMCB))
		matplot( resMCB$rLogLik, type="l" )
		plotChainPopMoves(resMCB)
	}
	load(file=file.path("tmp","resMCB_steady4.RData"))	#resMCB
	
	resMCBO <- transOrigPopt(resMCB,poptDistr=poptDistr$trans)
	plot(as.mcmc.list(resMCBO))
	
	sampleN <- sample <-  stackChains(resMCB)
	tmp <- sampleN[ which.max(sampleN[,1]),-1]
	#copy2clip(deparse(tmp))
	#structure(c(0.228439444501037, -3.58235519122273, 0.0356533561112984,1.87008926040014, 38.9385098325832), .Names = c("kY", "kO", "tLagLeaf","tLagRoot", "biasDiffRespLitterfall"))


	
	minLogLik <- quantile(sampleN[,1], probs=c(0.05) )	# empirical 95%
	sampleN0 <- sample0 <- sampleN[ sampleN[,1] >= minLogLik, ]
	minLogLik2 <- getRLogLikQuantile(sampleN) 	# theoretical criterion
	sampleN02 <- sample02 <- sampleN[ sampleN[,1] >= minLogLik2, ]
	sample[,-1] <- transOrigPopt(sampleN[,-1],  poptDistr=poptDistr$trans)
	sample0[,-1] <- transOrigPopt(sampleN0[,-1],  poptDistr=poptDistr$trans)
	sample02[,-1] <- transOrigPopt(sampleN02[,-1],  poptDistr=poptDistr$trans)
	cor(sampleN02[,-1])
	hist(sampleN02[,1])
	
	#inspect loglik-suface
	Ys <- stackChains(resMC$Y)
	Ys0 <- Ys0O <- Ys[Ys[,1]>=minLogLik,]
	Ys0O[,poptnames] <- transOrigPopt(Ys0[,poptnames], poptDistr$trans[poptnames])
	
	library(lattice)
	# round numbers to see something in levelplot else points get too small
	#sampleSig <- apply(sample[ (sample[,"rLogLik"] >= max(sample[,"rLogLik"]-1.9) & (sample[,"h"]<0.2)), ],2,function(var){
	smp <- sampleN02
	smp <- sample02
	smp <- sampleN0
	smp <- sample0
	#sampleSig <- apply(smp[ (smp[,"rLogLik"] >= max(smp[,"rLogLik"]-1.9)), ],2,function(var){
	smp <- cbind(smp, tvrY=1/smp[,"kY"], tvrO=1/smp[,"kO"])
	plotConditional2D(smp,"kY","kO")	#in plotTwDEMC.R
	plotConditional2D(smp,"tvrY","tvrO")
	plotConditional2D(smp,"tLagLeaf","tLagRoot")
	plotConditional2D(smp,"biasDiffRespLitterfall","kO")
	
	smp <- Ys0O	
	smp <- cbind(smp, tvrY=1/smp[,"kY"], tvrO=1/smp[,"kO"])
	plotConditional2D(smp,"tvrY","tvrO","somOFM")	
	plotConditional2D(smp,"tLagLeaf","tLagRoot","somOFM")	# almost no constrain
	plotConditional2D(smp,"tvrY","tvrO","respFM")			# almost no constrain
	plotConditional2D(smp,"tLagLeaf","tLagRoot","respFM")	# constraines tLagRoot above 6
	plotConditional2D(smp,"tvrY","tvrO","parms")			# almost no constrain
	plotConditional2D(smp,"tLagLeaf","tLagRoot","parms")	# constraines tLagRoot above 6

	plotConditional2D(smp,"biasDiffRespLitterfall","tLagRoot","parms")
	plotConditional2D(smp,"biasDiffRespLitterfall","tLagRoot","respFM")
	plotConditional2D(smp,"biasDiffRespLitterfall","tLagRoot","somOFM")
	
	p1 <- ggplotDensity.twDEMC(resMC, poptDistr=poptDistr)
	#print(p1)
	p2 <- ggplotDensity.twDEMC(resMC, poptDistr=poptDistr, doTransOrig=TRUE)
	#print(p2)
	twPairs(sampleN0)
	
	
	
}


