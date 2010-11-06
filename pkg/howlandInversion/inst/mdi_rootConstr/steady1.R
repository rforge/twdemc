.tmp.f <- function(){
	setupClusterHowlandDev(pkgDir = ".")
	parms0 <- HowlandParameterPriors$parms0
	model <- list(
		modMeta=modMetaICBM1()
		,fInitState=initState.howland.ICBM1SteadyState
		,fSolve=solveICBM1
	)
	
	argsFLogLik <- argsFLogLikRemoteFun <- list(
		model=model
		#remoteFun=of.howlandSteady		# will not work with twDEMC
		#,poptDistr=poptDistr			# must set when parameters are known
		,obs=Howland14C$obsNutrientSite
		,input=Howland14C$litter
		, parms=parms0           		##<< default parameters to the model
		, fCalcBiasedInput=meanInput    ##<< function(input,padj,...){obs} possibility to account for bias and to optimize bias parameters
		, fCalcIROLayer=calcIROLayer	##<< function to calculate iRofO-Layer
		, fCalcSteadyPars=calcSteadyK_ICBM1
	)
	#sfExport("argsFLogLik")			# done after assigning poptDistr	
	
	windows(width=4.4,height=3.4,pointsize=10, record=TRUE)
	par( las=1 )					#also y axis labels horizontal
	par(mar=c(2.0,3.3,0,0)+0.3 )  #margins
	par(tck=0.02 )				#axe-tick length inside plots             
	par(mgp=c(1.1,0.2,0) )  #positioning of axis title, axis labels, axis
	
}




mdi.cY <- function(){
	#using sfRemoteWrapper and exporting before 
	#poptnames <- c("h","cY")
	poptnames <- c("cY")
	poptDistr <- argsFLogLik$poptDistr <- twConstrainPoptDistr(poptnames, HowlandParameterPriors$parDistr)
	sfExport("argsFLogLik")
	
	argsFLogLik2 <- argsFLogLik
	#tmp <- argsFLogLik2$remoteFun; mtrace(tmp); argsFLogLik2$remoteFun<-tmp
	#resOf <- sfRemoteWrapper( normpopt=c(cY=logit(cYopt),h=logit(hOpt)), remoteFunArgs=argsFLogLik2)
	resOf <- sfRemoteWrapper( normpopt=poptDistr$mu["cY"], remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik2)
	
	#sfRemoteWrapper(normpopt=poptDistr$mu, remoteFun=of.howlandSteady, remoteFunArgs=substitute(argsFLogLik))
	resCl <- sfClusterCall( sfRemoteWrapper, normpopt=poptDistr$mu, remoteFun=of.howlandSteady, remoteFunArgs=substitute(argsFLogLik) )
	head(resCl[[1]])
	
	
	pnorm <- transNormPopt(unlist(parms0), HowlandParameterPriors$parDistr$trans[names(parms0)])
	cYs <- seq(0.5,1,length.out=101)[-c(1,101)]
	#cYs <- seq(0.99,1,length.out=101)[-c(1,101)]
	cYi <- cYs[1]
	fOpt <- function(cYi, remoteFun=of.howlandSteady){
		normpopt <- c(cY=logit(cYi))
		sum(sfRemoteWrapper( normpopt=normpopt, remoteFun=remoteFun, remoteFunArgs=substitute(argsFLogLik) ))
	}
	tmp <- sfSapply(cYs, fOpt, remoteFun=of.howlandSteady)
	plot( tmp ~ cYs)
	cYOpt <- optimize(fOpt, interval=c(0.6,0.8), maximum = TRUE)$maximum
	resOf <- sfRemoteWrapper( normpopt=c(cY=cYOpt), remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik)
	sort(resOf)
	res <- attr(resOf,"out")
	#colnames(res)
	matplot(res[,"time"], res[,c("inputLeaf_c12","inputLeaf_c14","inputRoot_c12","inputRoot_c14")], type="l" )
	matplot(res[,"time"], res[,c("Y_c12","Y_c14","O_c12","O_c14","cStock")], type="l" )
	matplot(res[,"time"], res[,c("respY_c12","respY_c14","respO_c12","respO_c14")], type="l" )
	matplot(res[,"time"], res[,c("F14C_Y","F14C_O","F14CT","respF14CT")], type="l" )
	lines( delta2iR14C(delta14C)/model$modMeta$iR14CStandard ~ yr, data=delta14Catm, col="gray")
	points( argsFLogLik$obs$somOFM[,"times"],argsFLogLik$obs$somOFM[,"obs"], col="red")
	points( argsFLogLik$obs$respFM[,"times"],argsFLogLik$obs$respFM[,"obs"], col="blue")
}
	
mdi.cY_h_1dInner <- function(){
	#using sfRemoteWrapper and exporting before 
	poptnames <- c("h","cY")
	poptDistr <- argsFLogLik$poptDistr <- twConstrainPoptDistr(poptnames, HowlandParameterPriors$parDistr)
	sfExport("argsFLogLik")
	
	pnorm <- transNormPopt(unlist(parms0), HowlandParameterPriors$parDistr$trans[names(parms0)])
	hOrig=hs[1]; cYi=0.2
	fOptInner <- function(cYi,hOrig){
		normpopt <- c(cY=logit(cYi),h=logit(hOrig))
		sum(sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteady, remoteFunArgs=substitute(argsFLogLik) ))
	}
	sfExport("fOptInner")
	sfExport("of.howlandSteady")
	fOpt <- function(h){
		optimize(fOptInner, interval=c(0.05,0.9), hOrig=h, maximum = TRUE)$objective
	}
	tmp.f <- function(){
		hs <- seq(0.1,0.5,length.out=11)[-c(1,11)]
		tmp <- sapply(hs, fOpt)
 	}
	hs <- seq(0,0.2,length.out=31)[-c(1,31)]
	#hs <- seq(0,1,length.out=31)[-c(1,31)]
	tmp <- sfSapply(hs, fOpt)
	plot(tmp~hs)
	#hOpt.res <- optimize(fOpt, c(0.02,0.08), maximum=TRUE )
	hOpt.res <- optimize(fOpt, c(0.02,0.14), maximum=TRUE )
	hOpt.res$objective
	hOpt <- hOpt.res$maximum
	cYOpt.res <- optimize(fOptInner, interval=c(0.05,0.9), hOrig=hOpt, maximum = TRUE)
	cYOpt.res$objective
	cYOpt <- cYOpt.res$maximum
	c(h=hOpt,cY=cYOpt)
	
	argsFLogLik2 <- argsFLogLik
	#tmp <- argsFLogLik2$remoteFun; mtrace(tmp); argsFLogLik2$remoteFun<-tmp
	resOf <- sfRemoteWrapper( normpopt=c(cY=logit(cYOpt),h=logit(hOpt)), remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik2)
	sort(resOf)
	sort(attr(resOf,"logLikParms"))
	res <- attr(resOf,"out")
	#colnames(res)
	matplot(res[,"time"], res[,c("inputLeaf_c12","inputLeaf_c14","inputRoot_c12","inputRoot_c14")], type="l" )
	matplot(res[,"time"], res[,c("Y_c12","Y_c14","O_c12","O_c14","cStock")], type="l" )
	matplot(res[,"time"], res[,c("respY_c12","respY_c14","respO_c12","respO_c14")], type="l" )
	
	matplot(res[,"time"], res[,c("F14C_Y","F14C_O","F14CT","respF14CT")], type="l", xlab="Nutrient Site, Time (yr)", ylab="" )
	mtext( "14C (Fraction modern)", side=2, line=2.5, las=0)
	#matplot(res[,"time"], res[,c("F14C_Y","F14C_O","respF14CT")], type="l" )
	lines( res[,"time"], attr(resOf,"iROLayer"), col="orange")
	#lines( res[,"time"], res[,"F14C_Y"], col="orange")
	#lines( res[,"time"], res[,"F14T"], col="orange")
	lines( delta2iR14C(delta14C)/model$modMeta$iR14CStandard ~ yr, data=delta14Catm, col="gray")
	points( argsFLogLik$obs$somOFM[,"times"],argsFLogLik$obs$somOFM[,"obs"], col="orange")
	points( argsFLogLik$obs$respFM[,"times"],argsFLogLik$obs$respFM[,"obs"], col="blue")
	legend( "topright", inset=0.05, xjust=1, yjust=1, bg="white"
		,legend=c( "RespSOM","O-Layer","Young","Old","SOM","Atmosphere")
		,lty=(c(4,1,1,2,3,1))
		,pch=(c(1,1,NA,NA,NA,NA))
		,col=(c("blue","orange","black","red","green","gray"))
	)
	
}

mdi.cY_h <- function(){
	#using sfRemoteWrapper and exporting before 
	poptnames <- c("h","cY")
	poptDistr <- argsFLogLik$poptDistr <- twConstrainPoptDistr(poptnames, HowlandParameterPriors$parDistr)
	sfExport("argsFLogLik")
	
	pnorm <- transNormPopt(unlist(parms0), HowlandParameterPriors$parDistr$trans[names(parms0)])
	normpopt <- pnorm[poptnames]
	fOpt <- function(normpopt){
		sum(sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteady, remoteFunArgs=substitute(argsFLogLik) ))
	}
	#fOpt(normpopt)
	resOpt <- optim(normpopt, fOpt, method="Nelder-Mead", hessian = TRUE, control=list(maxit=1000, fnscale=-1))
	#resOpt <- optim(normpopt, fOpt, method="BFGS", control=list(fnscale=-1), hessian = TRUE)
	transOrigPopt(resOpt$par, HowlandParameterPriors$parDistr$trans[names(resOpt$par)])
	
	argsFLogLik2 <- argsFLogLik
	#tmp <- argsFLogLik2$remoteFun; mtrace(tmp); argsFLogLik2$remoteFun<-tmp
	resOf <- sfRemoteWrapper( normpopt=resOpt$par, remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik2)
	#resOf <- sfRemoteWrapper( normpopt=c(cY=logit(cYOpt),h=logit(hOpt)), remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik2)
	sort(resOf)
	sort(attr(resOf,"logLikParms"))
	
	res <- attr(resOf,"out")
	#colnames(res)
	matplot(res[,"time"], res[,c("inputLeaf_c12","inputLeaf_c14","inputRoot_c12","inputRoot_c14")], type="l" )
	matplot(res[,"time"], res[,c("Y_c12","Y_c14","O_c12","O_c14","cStock")], type="l" )
	matplot(res[,"time"], res[,c("respY_c12","respY_c14","respO_c12","respO_c14")], type="l" )
	plotHowlandFM( res, attr(resOf,"obs"))
	

	#------ explore posterior with MCMC using Hessian
	#tmp.fcovarHessian <- function(){
		covMat <- solve( -resOpt$hessian )		#solve(resOpt$hessian)*resOpt
		.nPops=3
		Zinit <- initZtwDEMCNormal( resOpt$par, covMat, nChains=4*.nPops, nPops=.nPops)
		#resMC <- twDEMCBatch( Zinit, nGen=4*5, fLogLik=of.howlandSteady, argsFLogLik=argsFLogLik, nPops=.nPops, debugSequential=TRUE )
		resMC <- twDEMCBatch( Zinit, nGen=500, fLogLik=of.howlandSteady, argsFLogLik=argsFLogLik, nPops=.nPops )
		plot(as.mcmc.list(resMC))
		resMC <- twDEMCBatch( resMC, nGen=1000 )
		plot(as.mcmc.list(resMC))
		resMCB <- thin(resMC,start=300)
		matplot( resMCB$rLogLik, type="l" )
		plot(as.mcmc.list(resMCB))
		sample <- stackChains(resMCB)
		cor(sample[,-1])
		hist(sample[,1])
	#}
	
	#------ explore posterior with MCMC using prior
	tmp.fcovarPrior <- function(){
		covMat <- poptDistr$sigma    
		.nPops=3
		Zinit <- initZtwDEMCNormal( resOpt$par, covMat, nChains=4*.nPops, nPops=.nPops)
		resMC <- twDEMCBatch( Zinit, nGen=500, fLogLik=of.howlandSteady, argsFLogLik=argsFLogLik, nPops=.nPops )
		matplot(resMC$pAccept, type="l")
		plot(as.mcmc.list(resMC))
		resMC <- twDEMCBatch( resMC, nGen=1000, doRecordProposals=TRUE )
		#resMC <- twDEMCBatch( resMC, nGen=2500, doRecordProposals=TRUE )
		plot(as.mcmc.list(resMC))
		matplot( resMC$rLogLik[-(1:10),], type="l" )
		#resMC <- twDEMCBatch( resMC, nGen=2000 )
		#plot(as.mcmc.list(resMC))
		resMCB <- thin(resMC, start=300)
		plot(as.mcmc.list(resMCB))
		matplot( resMCB$rLogLik, type="l" )
		plotChainPopMoves(resMCB)
	}
	
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
	
	#-- start chains only in optimum
	tmp.f <- function(){
		Zinit2 <- initZtwDEMCSub( sampleN02[,-1], nChains=4*.nPops, nPops=.nPops) 
		resMC2 <- twDEMCBatch( Zinit2, nGen=500, fLogLik=of.howlandSteady, argsFLogLik=argsFLogLik, nPops=.nPops )
		matplot( resMC2$rLogLik, type="l" )
		resMC2 <- twDEMCBatch( resMC2, nGen=1000 )
		matplot( resMC2$rLogLik, type="l" )
		#moves out again
	}
	#inspect loglik-suface
	tmp.f <- function(){
		Ys <- stackChains(resMC$Y)
		Ys0 <- Ys0O <- Ys[Ys[,1]>=minLogLik,]
		Ys0O[,poptnames] <- transOrigPopt(Ys0[,poptnames], poptDistr$trans[poptnames])
		#mtrace(of.howlandSteady)
		tmp <- do.call( of.howlandSteady, c( list(normpopt=c(h=20,cY=15)), argsFLogLik))
		sort(tmp)
	}
	
	library(lattice)
	# round numbers to see something in levelplot else points get too small
	#sampleSig <- apply(sample[ (sample[,"rLogLik"] >= max(sample[,"rLogLik"]-1.9) & (sample[,"h"]<0.2)), ],2,function(var){
	smp <- sampleN02
	smp <- sample02
	smp <- sampleN0
	smp <- sample0
	#sampleSig <- apply(smp[ (smp[,"rLogLik"] >= max(smp[,"rLogLik"]-1.9)), ],2,function(var){
	sampleSig <- apply(smp,2,function(var){
			grain <- diff(range(var))/60
			round(var/grain)*grain
		})
	levelplot(rLogLik~h*cY, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	
	smp <- Ys0O	
	sampleSig <- apply(smp,2,function(var){
			grain <- diff(range(var))/60
			round(var/grain)*grain
		})
	colnames(sampleSig)
	levelplot(somOFM~h*cY, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	levelplot(respFM~h*cY, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	levelplot(parms~h*cY, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	
	
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


