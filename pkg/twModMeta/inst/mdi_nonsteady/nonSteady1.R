# repeat steady2 (estimating kY and kO directly and calculating h and cY from steady state)
# with allowing initial input of old carbon to be lower than steady state
# the root input is adjusted to compensate for the missing respiration due to less input

.tmp.init <- function(){
	setupClusterHowlandDev(pkgDir = ".")
	# replace prior of kY by narrow version kYNarrow
	parms0 <- HowlandParameterPriors$parms0; parms0$kY <- parms0$kYNarrow; 
	parDistrKYNarrow <- lapply( HowlandParameterPriors$parDistr, function(comp){ comp["kY"] <- as.vector(comp["kYNarrow"]); comp }); # as.vector needed, else name is copied too
	# assume some accumulation
	times0 <- 1950:2007
	parms0$dO = as.numeric(Howland14C$obsNutrientSite$somStock[1,"obs"]*1/5 /(Howland14C$obsTowerSite$somStock[1,"times"]-times0[1]))	# increase rate so that over period 1/5 is accumulated 
	model <- list(
		modMeta=modMetaICBM1()
		,fInitState=initState.howland.ICBM1SteadyState
		,fSolve=solveICBM1
	)
	
	argsFLogLik <- argsFLogLikRemoteFun <- list(
		model=model
		#remoteFun=of.howlandSteadyRootConstr		# will not work with twDEMC
		#,poptDistr=poptDistr			# must set when parameters are known
		, times=times0
		,obs=Howland14C$obsNutrientSite
		,input=Howland14C$litter
		, parms=parms0           		##<< default parameters to the model
		, fCalcBiasedInput=meanInput    ##<< function(input,padj,...){obs} possibility to account for bias and to optimize bias parameters
		, fCalcIROLayer=calcIROLayer	##<< function to calculate iRofO-Layer
		, fCalcSteadyPars=calcRelaxSteadyHcYC0_ICBM1	##<< adjust h,C, and initial Ctot
	)
	#sfExport("argsFLogLik")			# done after assigning poptDistr	
	
	windows(width=4.4,height=3.4,pointsize=10, record=TRUE)
	par( las=1 )					#also y axis labels horizontal
	par(mar=c(2.0,3.3,0,0)+0.3 )  #margins
	par(tck=0.02 )				#axe-tick length inside plots             
	par(mgp=c(1.1,0.2,0) )  #positioning of axis title, axis labels, axis
}

mdi.kYkO <- function(){
	#using sfRemoteWrapper and exporting before 
	poptnames <- c("kY","kO")
	poptDistr <- argsFLogLik$poptDistr <- twConstrainPoptDistr(poptnames, parDistrKYNarrow)
	sfExport("argsFLogLik")
	
	pnorm <- transNormPopt(unlist(parms0), parDistrKYNarrow$trans[names(parms0)])
	normpopt <- pnorm[poptnames]
	#mtrace(of.howlandSteadyRootConstr)
	#resOf <- sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik)
	#mtrace(calcRelaxSteadyHcYC0_ICBM1); argsFLogLik$fCalcSteadyPars=calcRelaxSteadyHcYC0_ICBM1
	#mtrace(calcRelaxSteadyHcY_ICBM1)
	fOpt <- function(normpopt){
		sum(sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik ))
	}
	fOpt(normpopt)
	#resOpt <- optim(normpopt, fOpt, method="Nelder-Mead", hessian = TRUE, control=list(maxit=1000, fnscale=-1))
	resOpt <- optim(normpopt, fOpt, method="BFGS", control=list(fnscale=-1), hessian = TRUE)
	1/transOrigPopt(poptDistr$mu, parDistrKYNarrow$trans[poptnames])
	1/transOrigPopt(normpopt, parDistrKYNarrow$trans[poptnames])
	1/transOrigPopt(resOpt$par, parDistrKYNarrow$trans[poptnames])
	
	
	argsFLogLik2 <- argsFLogLik
	#argsFLogLik2$useRImpl=TRUE
	#tmp <- argsFLogLik2$remoteFun; mtrace(tmp); argsFLogLik2$remoteFun<-tmp
	#tmp <- argsFLogLik2$model$fSolve; mtrace(tmp); argsFLogLik2$model$fSolve<-tmp
	#mtrace(derivICBM1)
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
}

mdi.kYkOC0 <- function(){	
	#using sfRemoteWrapper and exporting before 
	poptnames <- c("kY","kO","Ctot0")
	argsFLogLik$fCalcSteadyPars=calcRelaxSteadyHcY_ICBM1
	
	poptDistr <- argsFLogLik$poptDistr <- twConstrainPoptDistr(poptnames, parDistrKYNarrow)
	sfExport("argsFLogLik")
	
	pnorm <- transNormPopt(unlist(parms0), parDistrKYNarrow$trans[names(parms0)])
	normpopt <- pnorm[poptnames]
	#mtrace(of.howlandSteadyRootConstr)
	#resOf <- sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik)
	#mtrace(calcRelaxSteadyHcYC0_ICBM1); argsFLogLik$fCalcSteadyPars=calcRelaxSteadyHcYC0_ICBM1
	#mtrace(calcRelaxSteadyHcY_ICBM1); argsFLogLik$fCalcSteadyPars=calcRelaxSteadyHcY_ICBM1
	#mtrace(calcRelaxSteadyHcY_ICBM1)
	#argsFLogLik$useRImpl=TRUE
	#mtrace(of.howlandSteadyRootConstr)
	fOpt <- function(normpopt){
		sum(sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik ))
	}
	fOpt(normpopt)
	#resOpt <- optim(normpopt, fOpt, method="Nelder-Mead", hessian = TRUE, control=list(maxit=1000, fnscale=-1))
	resOpt <- optim(normpopt, fOpt, method="BFGS", control=list(fnscale=-1), hessian = TRUE)
	1/transOrigPopt(poptDistr$mu, parDistrKYNarrow$trans[poptnames])
	1/transOrigPopt(normpopt, parDistrKYNarrow$trans[poptnames])
	1/transOrigPopt(resOpt$par, parDistrKYNarrow$trans[poptnames])
	
	
	argsFLogLik2 <- argsFLogLik
	#argsFLogLik2$useRImpl=TRUE
	#tmp <- argsFLogLik2$remoteFun; mtrace(tmp); argsFLogLik2$remoteFun<-tmp
	#tmp <- argsFLogLik2$model$fSolve; mtrace(tmp); argsFLogLik2$model$fSolve<-tmp
	#mtrace(derivICBM1)
	resOf <- sfRemoteWrapper( normpopt=resOpt$par, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik2)
	#resOf <- sfRemoteWrapper( normpopt=c(cY=logit(cYOpt),h=logit(hOpt)), remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik2)
	sort(resOf)
	sort(attr(resOf,"logLikParms"))
	
	res <- attr(resOf,"out")
	#colnames(res)
	matplot(res[,"time"], res[,c("inputLeaf_c12","inputLeaf_c14","inputRoot_c12","inputRoot_c14")], type="l" )
	matplot(res[,"time"], res[,c("Y_c12","Y_c14","O_c12","O_c14","cStock")], type="l", ylim=c(0,1100) )
	matplot(res[,"time"], res[,c("Y_c12","O_c12","cStock")], type="l", ylim=c(0,1200) )
	matplot(res[,"time"], res[,c("decY_c12","respY_c12")], type="l" )
	matplot(res[,"time"], res[,c("respY_c12","respY_c14","respO_c12","respO_c14")], type="l" )
	plotHowlandFM( res, attr(resOf,"obs"))
}
mdi.kYkOC0dO <- function(){	
	#using sfRemoteWrapper and exporting before 
	poptnames <- c("kY","kO","Ctot0","dO")
	argsFLogLik$fCalcSteadyPars=calcRelaxSteadyHcY_ICBM1
	
	poptDistr <- argsFLogLik$poptDistr <- twConstrainPoptDistr(poptnames, parDistrKYNarrow)
	sfExport("argsFLogLik")
	
	pnorm <- transNormPopt(unlist(parms0), parDistrKYNarrow$trans[names(parms0)])
	normpopt <- pnorm[poptnames]
	#mtrace(of.howlandSteadyRootConstr)
	#resOf <- sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik)
	#mtrace(calcRelaxSteadyHcYC0_ICBM1); argsFLogLik$fCalcSteadyPars=calcRelaxSteadyHcYC0_ICBM1
	#mtrace(calcRelaxSteadyHcY_ICBM1); argsFLogLik$fCalcSteadyPars=calcRelaxSteadyHcY_ICBM1
	#mtrace(calcRelaxSteadyHcY_ICBM1)
	#argsFLogLik$useRImpl=TRUE
	#mtrace(of.howlandSteadyRootConstr)
	fOpt <- function(normpopt){
		sum(sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik ))
	}
	fOpt(normpopt)
	#resOpt <- optim(normpopt, fOpt, method="Nelder-Mead", hessian = TRUE, control=list(maxit=1000, fnscale=-1))
	resOpt <- optim(normpopt, fOpt, method="BFGS", control=list(fnscale=-1), hessian = TRUE)
	1/transOrigPopt(poptDistr$mu, parDistrKYNarrow$trans[poptnames])
	1/transOrigPopt(normpopt, parDistrKYNarrow$trans[poptnames])
	1/transOrigPopt(resOpt$par, parDistrKYNarrow$trans[poptnames])
	
	
	argsFLogLik2 <- argsFLogLik
	#argsFLogLik2$useRImpl=TRUE
	#tmp <- argsFLogLik2$remoteFun; mtrace(tmp); argsFLogLik2$remoteFun<-tmp
	#tmp <- argsFLogLik2$model$fSolve; mtrace(tmp); argsFLogLik2$model$fSolve<-tmp
	#mtrace(derivICBM1)
	resOf <- sfRemoteWrapper( normpopt=resOpt$par, remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik2)
	#resOf <- sfRemoteWrapper( normpopt=c(cY=logit(cYOpt),h=logit(hOpt)), remoteFun=of.howlandSteadyRootConstr, remoteFunArgs=argsFLogLik2)
	sort(resOf)
	sort(attr(resOf,"logLikParms"))
	
	res <- attr(resOf,"out")
	#colnames(res)
	matplot(res[,"time"], res[,c("inputLeaf_c12","inputLeaf_c14","inputRoot_c12","inputRoot_c14")], type="l" )
	matplot(res[,"time"], res[,c("Y_c12","Y_c14","O_c12","O_c14","cStock")], type="l", ylim=c(0,1100) )
	matplot(res[,"time"], res[,c("Y_c12","O_c12","cStock")], type="l", ylim=c(0,1200) )
	matplot(res[,"time"], res[,c("decY_c12","respY_c12")], type="l" )
	matplot(res[,"time"], res[,c("respY_c12","respY_c14","respO_c12","respO_c14")], type="l" )
	plotHowlandFM( res, attr(resOf,"obs"))
	
	#------ explore posterior with MCMC using prior
	tmp.fcovarPrior <- function(){
		covMat <- poptDistr$sigma    
		.nPops=3
		Zinit <- initZtwDEMCNormal( resOpt$par, covMat, nChains=4*.nPops, nPops=.nPops)
		resMC <- twDEMCBatch( Zinit, nGen=500, fLogLik=of.howlandSteadyRootConstr, argsFLogLik=argsFLogLik, nPops=.nPops )
		matplot(resMC$pAccept, type="l")
		plot(as.mcmc.list(resMC))
		#resMC <- twDEMCBatch( resMC, nGen=1000, doRecordProposals=TRUE )
		#resMC <- twDEMCBatch( resMC, nGen=2500, doRecordProposals=TRUE )
		plot(as.mcmc.list(resMC))
		matplot( resMC$rLogLik[-(1:10),], type="l" )
		#resMC <- twDEMCBatch( resMC, nGen=2000 )
		#plot(as.mcmc.list(resMC))
		resMCB <- thin(resMC, start=200)
		plot(as.mcmc.list(resMCB))
		matplot( resMCB$rLogLik, type="l" )
		plotChainPopMoves(resMCB)
		save(resMCB,file=file.path("tmp","resMCB_nonSteady1.RData"))
	}
	load(file.path("tmp","resMCB_nonSteady1.RData"))	#resMCB
	
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
	
	twPairs(sampleN0)		# no big correlations
	p1 <- ggplotDensity.twDEMC(resMCB, poptDistr=poptDistr)
	#print(p1)
	p2 <- ggplotDensity.twDEMC(resMCB, poptDistr=poptDistr, doTransOrig=TRUE)
	#print(p2)
	
	
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
	sampleSig <- apply(smp,2,function(var){
			grain <- diff(range(var))/60
			round(var/grain)*grain
		})
	levelplot(rLogLik~kY*kO, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	levelplot(rLogLik~tvrY*tvrO, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	
	smp <- Ys0O	
	smp <- cbind(smp, tvrY=1/smp[,"kY"], tvrO=1/smp[,"kO"])	
	sampleSig <- apply(smp,2,function(var){
			grain <- diff(range(var))/60
			round(var/grain)*grain
		})
	colnames(sampleSig)
	levelplot(somOFM~tvrY*tvrO, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	levelplot(respFM~tvrY*tvrO, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	levelplot(parms~tvrY*tvrO, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	
	library(Rcmdr)
	#sampleSig <- sample[ (sample[,"rLogLik"] >= max(sample[,"rLogLik"]-1.9)& (sample[,"h"]<0.2)), ]
	sampleSig <- sample[ sample[,"rLogLik"] >= max(sample[,"rLogLik"]-1.9), ]
	ds <- as.data.frame(sampleSig)
	scatter3d(ds$h, ds$rLogLik, ds$cY
		, surface=FALSE
		,bg="white", axis.scales=TRUE, grid=TRUE, ellipsoid=FALSE, xlab="h" 
		,ylab="rLogLik", zlab="cY"
		, point.col=rev(heat.colors(100))[round(rescale(ds$rLogLik,to=c(1,100)))]
	)
	
	ds <- as.data.frame(sampleN0)
	ds <- as.data.frame(sample0)
	ds <- as.data.frame(sampleN02)
	ds <- as.data.frame(sample02)

	Ys
	ds <- as.data.frame(Ys0)
	#ds <- as.data.frame(Ys0O)
	colnames(ds)
	tmp.var <- "rLogLik"
	tmp.var <- "parms"
	#tmp.var <- "respCum"  # range e-24, practically zero
	tmp.var <- "respFM"
	tmp.var <- "somOFM"
	scatter3d(ds$h, ds[[tmp.var]], ds$cY
		, surface=FALSE
		,bg="white", axis.scales=TRUE, grid=TRUE, ellipsoid=FALSE, xlab="h" 
		,ylab=tmp.var, zlab="cY"
		, point.col=rev(heat.colors(100))[round(rescale(ds[[tmp.var]],to=c(1,100)))]
	)
	
	plot(sample0[,1]~sample0[,"h"])
	plot(sample0[,1]~sample0[,"cY"])
	
	#------- display k values
	smp <- sample02
	smp <- sample0
	sample0k <- cbind( rLogLik=smp[,1], calcSteadyK_ICBM1(
		Ctot = argsFLogLik$obs$somStock[1,2] 
		,iY = sum(sapply(argsFLogLik$input,"[",1,2))
		,cY = smp[,"cY"]
		,h = smp[,"h"]
	))
	str(sample0k)
	sample0tvr <- cbind( rLogLik=smp[,1], 1/calcSteadyK_ICBM1(
			Ctot = argsFLogLik$obs$somStock[1,2] 
			,iY = sum(sapply(argsFLogLik$input,"[",1,2))
			,cY = smp[,"cY"]
			,h = smp[,"h"]
		))
	colnames(sample0tvr)[-1] <- c("tvrY","tvrO")
	str(sample0tvr)
	
	
	sampleSig <- apply(sample0k,2,function(var){
			grain <- diff(range(var))/60
			round(var/grain)*grain
		})
	levelplot(rLogLik~kY*kO, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))

	sampleSig <- apply(sample0tvr,2,function(var){
			grain <- diff(range(var))/60
			round(var/grain)*grain
		})
	levelplot(rLogLik~tvrY*tvrO, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
}


