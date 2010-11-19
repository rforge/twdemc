# same as steady 4 but stochastic year to year variability of litter fall

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
		#, fCalcBiasedInput=meanInput    ##<< function(input,padj,...){obs} possibility to account for bias and to optimize bias parameters
		#, fCalcBiasedInput=meanInputFluctuating    ##<< function(input,padj,...){obs} possibility to account for bias and to optimize bias parameters
		, fCalcBiasedInput=NULL
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

.generateFluctuatingInput <- function(
){
	set.seed(0815)
	obsLitter = Howland14C$litter
	padj=HowlandParameterPriors$parms0
	filename=file.path("data","inputFluctuating.RData")
	
	#padj <- HowlandParameterPriors$parms0 
	inputFluctuating <- inputF <- input <- meanInputFluctuating(obsLitter, padj=padj)	# in of_steadyState.R
	# root input and associated sd will be adjusted in of.howlandSteadyRootConstr by a factor 
	#inputF <- inputadj
	plot(obs~times,data=inputF$leaf)
	plot(obs~times,data=inputF$root)
	plot(inputF$root[,"obs"]+inputF$leaf[,"obs"]~inputF$leaf[,1])
	mean(inputF$root[,"obs"]+inputF$leaf[,"obs"])
	plot(inputF$root[,"obs"]~inputF$leaf[,"obs"])
	
	save(inputFluctuating,file=filename)

	# generate entire series
	inputFluctuatingEnsemble <- lapply(1:30, function(i){
			inputFluctuating <- inputF <- input <- meanInputFluctuating(obsLitter, padj=padj)	# in of_steadyState.R
		})
	filename1=file.path("data",paste("inputFluctuatingEnsemble",".RData",sep=""))
	save(inputFluctuatingEnsemble,file=filename1)
	
	#second: try to decrease sd of litterfall
	relErr <- 0.05
	obsLitter$leaf[,"sdObs"] <- mean(obsLitter$leaf[,"obs"])*relErr 
	inputFluctuating <- inputF <- input <- meanInputFluctuating(obsLitter, padj=padj)	# in of_steadyState.R
	plot(obs~times,data=inputF$leaf)
	
	filename2=file.path("data",paste("inputFluctuating_",relErr*100,".RData",sep=""))
	save(inputFluctuating,file=filename2)

	# generate entire series
	inputFluctuatingEnsemble <- lapply(1:30, function(i){
		inputFluctuating <- inputF <- input <- meanInputFluctuating(obsLitter, padj=padj)	# in of_steadyState.R
	})
	filename3=file.path("data",paste("inputFluctuatingEnsemble",relErr*100,".RData",sep=""))
	save(inputFluctuatingEnsemble,file=filename3)
}


mdi.kLagLitterVar <- function(){
	#.generateFluctuatingInput()
	#from steady4-fit:
	normpopt <- structure(c(1.31004509263524, 0.0276374254632283, 1.0829503483462,6.65177437410762, 38.9279591345289), .Names = c("kY", "kO", "tLagLeaf","tLagRoot", "biasDiffRespLitterfall"))
	poptnames <- names(normpopt)
	poptDistr <- argsFLogLik$poptDistr <- twConstrainPoptDistr(poptnames, HowlandParameterPriors$parDistr)
	popt <- transOrigPopt(normpopt,poptDistr$trans)

	#.generateFluctuatingInput
	load(file.path("data","inputFluctuating.RData"))	#inputFluctuating, generated above
	#load(file.path("data","inputFluctuating_5.RData"))	#inputFluctuating, generated above
	argsFLogLik$input <- inputFluctuating
	sfExport("argsFLogLik")

	# gradient method will not work with stochastic nature of litter inputs, which vary between runs
	# instead of generating stochastic litter input, use several previously generated scenarios 
	
	#mtrace(of.howlandSteadyRootConstr)
	#resOf <- sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik)	
	fOpt <- function(normpopt){
		sum(sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik ))
	}
	fOpt <- function(normpopt, argsFLogLik){
		sum(do.call( of.howlandSteadyRootConstr, c(list(normpopt=normpopt),argsFLogLik)))
	}
	#fOpt(normpopt)
	resOpt <- optim(normpopt, fOpt, method="Nelder-Mead", hessian = TRUE, control=list(maxit=1000, fnscale=-1), argsFLogLik=argsFLogLik)
	#resOpt <- optim(normpopt, fOpt, method="BFGS", control=list(fnscale=-1), hessian = TRUE)
	transOrigPopt(poptDistr$mu, HowlandParameterPriors$parDistr$trans[poptnames])
	transOrigPopt(normpopt, HowlandParameterPriors$parDistr$trans[poptnames])
	(tmp <- transOrigPopt(resOpt$par, HowlandParameterPriors$parDistr$trans[poptnames]))
	#copy2clip(deparse(tmp))
	#structure(c(1.07109167581197, 0.0266326047938156, 1.06776201540149,6.64386347044198, 47.1583648738795), .Names = c("kY", "kO", "tLagLeaf","tLagRoot", "biasDiffRespLitterfall"))

	
	argsFLogLik2 <- argsFLogLik
	#tmp <- argsFLogLik2$remoteFun; mtrace(tmp); argsFLogLik2$remoteFun<-tmp
	#mtrace(of.howlandSteadyRootConstr)
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
	
	#------ determine by ensemble
	#load(file.path("data","inputFluctuatingEnsemble5.RData"))	#inputFluctuatingEnsemble, generated above
	load(file.path("data","inputFluctuatingEnsemble.RData"))	#inputFluctuatingEnsemble, generated above
	i=1
	sfExport("of.howlandSteadyRootConstr")
	parEnsNorm <- t(sfSapply( seq_along(inputFluctuatingEnsemble), function(i, argsFLogLik, inputFluctuatingEnsemble, normpopt, fOpt){
			argsFLogLik$input <- inputFluctuatingEnsemble[[i]]
			#sfExport("argsFLogLik")
			resOpt <- optim(normpopt, fOpt, method="Nelder-Mead", hessian = TRUE, control=list(maxit=1000, fnscale=-1), argsFLogLik=argsFLogLik)
			resOpt$par
		},argsFLogLik=argsFLogLik, inputFluctuatingEnsemble=inputFluctuatingEnsemble, normpopt=resOpt$par, fOpt=fOpt))
	parEns <- transOrigPopt(parEnsNorm, HowlandParameterPriors$parDistr$trans[poptnames])
	hist(parEns[,"tLagRoot"])
	plot( tLagRoot ~ biasDiffRespLitterfall, data=parEns)
	qplot( 1/kO, tLagRoot, data=as.data.frame(parEns))
	
	
	#------ explore posterior with MCMC using prior
	# make it a stochastic model: at each step generate new input
	argsFLogLikFluct <- within( argsFLogLik,{
			input=Howland14C$litter
			fCalcBiasedInput=meanInputFluctuating
		})
	sfExport("argsFLogLikFluct")
	
	#mtrace(of.howlandSteadyRootConstr)
	resOf <- sfRemoteWrapper( normpopt=resOpt$par, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLikFluct)
	.remoteDumpfileBasename=file.path("tmp","dumpRemote")
	resOfCl <- sfClusterCall( sfRemoteWrapper, normpopt=resOpt$par, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLikFluct, remoteDumpfileBasename=.remoteDumpfileBasename)
	
	tmp.fcovarPrior <- function(){
		covMat <- poptDistr$sigma    
		.nPops=3
		Zinit <- initZtwDEMCNormal( resOpt$par, covMat, nChains=4*.nPops, nPops=.nPops)
		#mtrace(of.howlandSteadyRootConstr)
		#resMC <- twDEMCBatch( Zinit, nGen=4*5, debugSequential=TRUE, fLogLik=of.howlandSteadyRootConstr, argsFLogLik=argsFLogLikFluct, nPops=.nPops, intResCompNames="parms" )
		.remoteDumpfileBasename=file.path("tmp","dumpRemote")
		resMC <- twDEMCBatch( Zinit, nGen=500, fLogLik=of.howlandSteadyRootConstr, argsFLogLik=argsFLogLikFluct, nPops=.nPops, intResCompNames="parms", remoteDumpfileBasename=.remoteDumpfileBasename )
		matplot(resMC$pAccept, type="l")
		plot(as.mcmc.list(resMC))
		#resMC <- twDEMCBatch( resMC, nGen=1000, doRecordProposals=TRUE )
		resMC <- twDEMCBatch( resMC, nGen=6000, doRecordProposals=TRUE )	# for integrating over nuisance tLagRoot and tLagLeaf
		resMCB <- thin(resMC, start=500)
		try(dir.create("tmp"))
		save(resMCB,file=file.path("tmp","resMCB_steady5.RData"))
		plotThinned(as.mcmc.list(resMCB))
		matplot( resMCB$rLogLik, type="l" )
		plotChainPopMoves(resMCB)
	}
	tmp.fDebug <- function(){
		load(paste(.remoteDumpfileBasename,".rda",sep=""))
		debugger(get(.remoteDumpfileBasename))
		#mtrace(remoteFun); eval(body)
	}
	load(file.path("tmp","resMCB_steady5.RData"))	#resMCB
	
	resMCBO <- transOrigPopt(resMCB,poptDistr=poptDistr$trans)
	plot(as.mcmc.list(resMCBO))
	
	sampleN <- sample <-  stackChains(resMCB)
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
	
	p1 <- ggplotDensity.twDEMC(resMCB, poptDistr=poptDistr)
	#print(p1)
	p2 <- ggplotDensity.twDEMC(resMCB, poptDistr=poptDistr, doTransOrig=TRUE)
	#print(p2)
	twPairs(sampleN0)		# no big correlations
	
	
	
}

