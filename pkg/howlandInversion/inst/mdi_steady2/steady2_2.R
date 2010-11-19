# estimating kY and kO directly and calculating h and cY from steady state

.tmp.init <- function(){
	setupClusterHowlandDev(pkgDir = ".")
	parms0 <- HowlandParameterPriors$parms0
	obs <- Howland14C
	obsSite <- Howland14C$obsNutrientSite
	# iL + iR + biasIR = resp 
	parms0$biasLitterRoot <- mean(obsSite$respCum[,"obs"]) - mean(obs$litter$leaf[,"obs"]) - mean(obs$litter$root[,"obs"]) 
	model <- list(
		modMeta=modMetaICBM1()
		,fInitState=initState.howland.ICBM1SteadyState
		,fSolve=solveICBM1
	)
	
	argsFLogLik <- argsFLogLikRemoteFun <- list(
		model=model
		#remoteFun=of.howlandSteady		# will not work with twDEMC
		#,poptDistr=poptDistr			# must set when parameters are known
		,obs=obsSite
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

mdi.kYkO <- function(){
	# test for repeating former version, see bias scenarios below
	#using sfRemoteWrapper and exporting before 
	poptnames <- c("kY","kO")
	poptDistr <- argsFLogLik$poptDistr <- twConstrainPoptDistr(poptnames, HowlandParameterPriors$parDistr)
	sfExport("argsFLogLik")
	
	pnorm <- transNormPopt(unlist(parms0), HowlandParameterPriors$parDistr$trans[names(parms0)])
	normpopt <- pnorm[poptnames]
	#mtrace(of.howlandSteady)
	resOf <- sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik)	
	fOpt <- function(normpopt){
		sum(sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik ))
	}
	#fOpt(normpopt)
	#resOpt <- optim(normpopt, fOpt, method="Nelder-Mead", hessian = TRUE, control=list(maxit=1000, fnscale=-1))
	resOpt <- optim(normpopt, fOpt, method="BFGS", control=list(fnscale=-1), hessian = TRUE)
	1/transOrigPopt(poptDistr$mu, HowlandParameterPriors$parDistr$trans[poptnames])
	1/transOrigPopt(normpopt, HowlandParameterPriors$parDistr$trans[poptnames])
	1/transOrigPopt(resOpt$par, HowlandParameterPriors$parDistr$trans[poptnames])
	
	
	argsFLogLik2 <- argsFLogLik
	#tmp <- argsFLogLik2$remoteFun; mtrace(tmp); argsFLogLik2$remoteFun<-tmp
	#mtrace(of.howlandSteady)
	resOf <- sfRemoteWrapper( normpopt=resOpt$par, remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik2)
	#resOf <- sfRemoteWrapper( normpopt=c(cY=logit(cYOpt),h=logit(hOpt)), remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik2)
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
		resMC <- twDEMCBatch( Zinit, nGen=500, fLogLik=of.howlandSteady, argsFLogLik=argsFLogLik, nPops=.nPops )
		matplot(resMC$pAccept, type="l")
		plot(as.mcmc.list(resMC))
		#resMC <- twDEMCBatch( resMC, nGen=1000, doRecordProposals=TRUE )
		#resMC <- twDEMCBatch( resMC, nGen=2500, doRecordProposals=TRUE )
		plot(as.mcmc.list(resMC))
		matplot( resMC$rLogLik[-(1:10),], type="l" )
		#resMC <- twDEMCBatch( resMC, nGen=2000 )
		#plot(as.mcmc.list(resMC))
		resMCB <- thin(resMC, start=200)
		try(dir.create("tmp"))
		save(resMCB,file=file.path("tmp","resMCB_steady2_2noBias.RData"))
		
		plot(as.mcmc.list(resMCB))
		matplot( resMCB$rLogLik, type="l" )
		plotChainPopMoves(resMCB)
	}
	load(file.path("tmp","resMCB_steady2_2_noBias.RData"))	#resMCB
	
	
	sampleN <- sample <-  stackChains(resMCB)
	minLogLik <- quantile(sampleN[,1], probs=c(0.05) )	# empirical 95%
	sampleN0 <- sample0 <- sampleN[ sampleN[,1] >= minLogLik, ]
	minLogLik2 <- getRLogLikQuantile(sampleN) 	# theoretical criterion
	sampleN02 <- sample02 <- sampleN[ sampleN[,1] >= minLogLik2, ]
	sample[,-1] <- transOrigPopt(sampleN[,-1],  poptDistr=poptDistr$trans)
	sample0[,-1] <- transOrigPopt(sampleN0[,-1],  poptDistr=poptDistr$trans)
	sample02[,-1] <- transOrigPopt(sampleN02[,-1],  poptDistr=poptDistr$trans)
	
	#inspect loglik-suface
	Ys <- stackChains(resMC$Y)
	Ys0 <- Ys0O <- Ys[Ys[,1]>=minLogLik,]
	Ys0O[,poptnames] <- transOrigPopt(Ys0[,poptnames], poptDistr$trans[poptnames])

	twPairs(sampleN0)		# no big correlations
	p1 <- ggplotDensity.twDEMC(resMCB, poptDistr=poptDistr)
	#print(p1)
	p2 <- ggplotDensity.twDEMC(resMCB, poptDistr=poptDistr, doTransOrig=TRUE)
	#print(p2)
}

mdi.kYkObiasLitter <- function(){
	poptnames <- c("kY","kO","biasLitterRoot")
	poptnames <- c("kY","kO","biasLitterRoot","biasLitterLeaf")
	poptDistr <- argsFLogLik$poptDistr <- twConstrainPoptDistr(poptnames, HowlandParameterPriors$parDistr)
	sfExport("argsFLogLik")
	
	pnorm <- transNormPopt(unlist(parms0), HowlandParameterPriors$parDistr$trans[names(parms0)])
	normpopt <- pnorm[poptnames]
	#mtrace(of.howlandSteady)
	resOf <- sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik)	
	fOpt <- function(normpopt){
		sum(sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik ))
	}
	#fOpt(normpopt)
	#resOpt <- optim(normpopt, fOpt, method="Nelder-Mead", hessian = TRUE, control=list(maxit=1000, fnscale=-1))
	normpopt["biasLitterRoot"] <- -100		# test if correct bias is retrieved
	resOpt <- optim(normpopt, fOpt, method="BFGS", control=list(fnscale=-1), hessian = TRUE)	# slightly depends on initial conditions
	transOrigPopt(poptDistr$mu, HowlandParameterPriors$parDistr$trans[poptnames])
	transOrigPopt(normpopt, HowlandParameterPriors$parDistr$trans[poptnames])
	transOrigPopt(resOpt$par, HowlandParameterPriors$parDistr$trans[poptnames])
	
	argsFLogLik2 <- argsFLogLik
	#tmp <- argsFLogLik2$remoteFun; mtrace(tmp); argsFLogLik2$remoteFun<-tmp
	#mtrace(of.howlandSteady)
	resOf <- sfRemoteWrapper( normpopt=resOpt$par, remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik2)
	#resOf <- sfRemoteWrapper( normpopt=c(cY=logit(cYOpt),h=logit(hOpt)), remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik2)
	sort(resOf)
	sort(attr(resOf,"logLikParms"))
	
	res <- attr(resOf,"out")
	#colnames(res)
	matplot(res[,"time"], res[,c("inputLeaf_c12","inputLeaf_c14","inputRoot_c12","inputRoot_c14")], type="l" )
	matplot(res[,"time"], res[,c("Y_c12","Y_c14","O_c12","O_c14","cStock")], type="l", ylim=c(0,1100) )
	matplot(res[,"time"], res[,c("respY_c12","respY_c14","respO_c12","respO_c14")], type="l" )
	plotHowlandFM( res, attr(resOf,"obs"))
	plot( rowSums(res[,c("respY_c12","respO_c12")]) ~ res[,"time"], type="l" )
	points( obs ~ times, data=as.data.frame(obsSite$respCum), col="maroon" )
	
	#------ explore posterior with MCMC using prior
	tmp.fcovarPrior <- function(){
		set.seed(0815)
		covMat <- poptDistr$sigma    
		.nPops=3
		Zinit <- initZtwDEMCNormal( resOpt$par, covMat, nChains=4*.nPops, nPops=.nPops)
		resMC <- twDEMCBatch( Zinit, nGen=500, fLogLik=of.howlandSteady, argsFLogLik=argsFLogLik, nPops=.nPops )
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
		save(resMCB,file=file.path("tmp","resMCB_steady2_2.RData"))
	}
	load(file.path("tmp","resMCB_steady2_2.RData"))	#resMCB
	
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
	
	#replace rate by turnvoer time
	sample0tvr <- sample0
	sample0tvr[,c("kY","kO")] <- 1 / exp(sampleN0[,c("kY","kO")])
	colnames(sample0tvr)[ colnames(sample0tvr) %in% c("kY","kO")] <- c("tvrY","tvrO") 
	
	#inspect loglik-suface
	Ys <- stackChains(resMCB$Y)
	Ys0 <- Ys0O <- Ys[Ys[,1]>=minLogLik,]
	Ys0O[,poptnames] <- transOrigPopt(Ys0[,poptnames], poptDistr$trans[poptnames])
	Ys0O[,poptnames] <- transOrigPopt(Ys0[,poptnames], poptDistr$trans[poptnames])
	YsTvr <- cbind(Ys0O, tvrY=1/Ys0O[,"kY"], tvrO=1/Ys0O[,"kO"])
	
	twPairs(sampleN0)		# no big correlations
	p1 <- ggplotDensity.twDEMC(resMCB, poptDistr=poptDistr)
	#print(p1)
	#mtrace(ggplotDensity.twDEMC)
	p2 <- ggplotDensity.twDEMC(resMCB, poptDistr=poptDistr, doTransOrig=TRUE, doDispLogLik=FALSE)
	#print(p2)

	#-------- fit the marginal density
	library(ks)
	x <- sample0tvr[,c("tvrY","tvrO")]
	H <- Hpi(x, binned=TRUE)
	fhat2 <- kde(x, H=H, compute.cont=TRUE, approx.cont=TRUE)
	plot(fhat2, display="filled.contour", color.palette=function(n){rev(heat.colors(n))})

	tmp.f.ks <- function(){
		# using package ks: most fancy (3D plot) take long and density is too smooth
		# library(ks)
		x2 <- sample0tvr[,c("tvrY","tvrO")]
		# takes long
		H2 <- Hpi(x2[1:100,], pilot="samse")
		# takes long
		#fhat3e <- kde(x2, H=H3, compute.cont=TRUE, approx.cont=TRUE)
		fhat2 <- kde(x2, H=H2)
		plot(fhat2, "filled.contour", color.palette=function(n){rev(heat.colors(n))})
		
		# using package ks: most fancy and robost, takes long
		x3 <- sample0tvr[,c("tvrY","tvrO","biasLitterRoot")]
		# takes long
		H3 <- Hpi(x3[1:100,], pilot="samse")
		# takes long
		#fhat3e <- kde(x2, H=H3, compute.cont=TRUE, approx.cont=TRUE)
		fhat3 <- kde(x3, H=H3)
		plot(fhat3)
		
		x3.gr <- cut2(sampleN02[,"biasLitterLeaf"],g=3)
		bo <- sample.int(length(x3.gr),100)
		H4 <- Hkda(x3[1:100,], x3.gr[1:100], bw="plugin")
		fhat4 <- kda.kde(x3, x3.gr, Hs=H4, compute.cont=TRUE)
		plot(fhat4)   
		
	}
	tmp.f.np <- function(){
		# seem to work faster with a fixed default bandwidth
		# estimating adaptive bandwidth on a subsample seems ok within minutes
		library(np)
		ds <- as.data.frame(sample0tvr)
		# cross validation not feasible with large dataset
		bw0 <- npudensbw(formula=~tvrY+tvrO, data=ds, bwtype="fixed", bwmethod="normal-reference", tol=.01, ftol=.01)
		#bw1 <- npudensbw(formula=~tvrY, data=ds, bwtype="generalized_nn", tol=.01, ftol=.01)
		#bw2 <- npudensbw(formula=~tvrY, data=ds, bwtype="adaptive_nn", tol=.1, ftol=.1)
	
		ds100 <- ds[sample.int(nrow(ds),300),]	# tradeoff: larger sample sizes become really slow
		bw1 <- npudensbw(formula=~tvrY+tvrO, data=ds100, bwtype="generalized_nn", tol=.01, ftol=.01)
		bw4 <- npudensbw(formula=~tvrY+tvrO+biasLitterRoot+biasLitterLeaf, data=ds100, bwtype="generalized_nn", tol=.01, ftol=.01)
		
		dsPred <- model.frame(~tvrY+tvrO,ds)
		#est <- npudens(bws=bw1, tdat=dsPred, edat=dsPred)		
		tmpf <- function(x,y, bws){
			npudens(bws,edat=data.frame(tvrY=x,tvrO=y))$dens
		}
		#mtrace(twPlot2DFun)
		tmp <- twPlot2DFun.contour( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws=bw0), xdiv=80)
		tmp2 <- twPlot2DFun.contour( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws=bw1), xdiv=80)
		tmp3 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmp2, xdiv=nrow(tmp2))
		# note the three humps corresponding to three chains
	}
	
	#--------- estimate a response surface of rLoglik
	#mtrace(plotMarginal2D)
	tmpMarg <- plotMarginal2D( smp, "tvrY", "tvrO" )
	tmp <- plotConditional2D( smp, "tvrY", "tvrO" )

	tmp.npreg <- function(){
		# uses MAD mean absolute deviation: close to integrating marginal distribution (tau=0.5)
		require(np)
		form=rLogLik~tvrY+tvrO
		#ds <- as.data.frame(YsTvr) # need to use sample for integrating as the number of points represents p(X)
		ds <- as.data.frame(sample0tvr)
		ds100 <- ds[sample.int(nrow(ds),200),]	# tradeoff: larger sample sizes become really slow
		#bw0 <- npregbw(formula=form, ds100, bwmethod="normal-reference" )
		#bw3 <- bw30 <- npregbw(formula=form, ds100, bwtype="generalized_nn" ,tol=.1, ftol=.1 )
		bw3 <- bw30 <- npregbw(formula=form, ds100, bwtype="adaptive_nn" ,tol=.1, ftol=.1 )
		# refine (but does not change)
		#bw3 <- npregbw(bws=bw30, formula=form, ds100, bwtype="adaptive_nn", tol=0.01, ftol=0.01 )
		#bw3$bw <- bw30$bw*2
		#plot(bw32)
		dsPred <- model.frame(form,ds)
		#est <- npreg(bws=bw3, tdat=dsPred, edat=dsPred)		
		lMin <- min(ds$rLogLik)
		tmpf <- function(x,y, bws, minLogLik=lMin){
			tmp <- npreg(bws=bws,txdat=dsPred[,-1], tydat=dsPred[,1], exdat=data.frame(tvrY=x,tvrO=y))
			ifelse( tmp$mean>minLogLik, tmp$mean, as.numeric(NA))
		}
		#tmpi<-sample.int(nrow(ds),30);bws=bw3;x=ds$tvrY[tmpi];y=ds$tvrO[tmpi]
		#tmpf(2,20,bw3)
		#mtrace(twPlot2DFun)
		#tmp <- twPlot2DFun.contour( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws=bw0), xdiv=80)
		#tmp2 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=bw3$bw*8;tmp}), xdiv=60, col=rev(heat.colors(20)) )
		tmp2 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=c(200,200);tmp}), xdiv=60, col=rev(heat.colors(20)) )
		#tmp2 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws=bw3), xdiv=60, col=rev(heat.colors(20)) )
		#tmp2 <- twPlot2DFun.contour( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws=bw3), xdiv=80, key.title=title(sub="Log-Like-\nlihood\n"), color.palette=function(n){rev(heat.colors(n))})
		tmp3 <- twPlot2DFun.contour( ds$tvrY, ds$tvrO, tmp2, xdiv=nrow(tmp2), xlab="tvrY",ylab="tvrO", key.title=title(sub="Log-Like-\nlihood\n\n"), color.palette=function(n){rev(heat.colors(n))})
		tmp3 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmp2, xdiv=nrow(tmp2), col=rev(heat.colors(20)) )
		# note the three humps corresponding to three chains		
	}
	tmp.loess <- function(){
		#loess fit: relies on gaussian errors: ok to approximate all dimensions to find optimum but do not use for marginals
		form = rLogLik ~ tvrY+tvrO+biasLitterLeaf+biasLitterRoot
		ds <- as.data.frame(YsTvr)
		modmat <- model.matrix(form, ds)[,-1]
		#rs1 <- loess(form, as.data.frame(sample0tvr))	#span=0.75 3/4 of the data used in each local fit
		rs2 <- loess(form, ds, span=(5^ncol(modmat))/nrow(ds) )
		plot( fitted(rs2) ~ ds$rLogLik )
		# conditional LogLik: at median litter-Bias
		mLeaf <- median(dsNew$biasLitterLeaf)
		mRoot <- median(dsNew$biasLitterRoot)
		lMin <- min(ds$rLogLik)
		tmpf <- function(tvrY,tvrO, biasLitterLeaf=mLeaf, biasLitterRoot=mRoot, mod=rs2, minLogLik=lMin){
			dfNew <- data.frame(tvrY=tvrY,tvrO=tvrO,biasLitterRoot=biasLitterLeaf, biasLitterLeaf=biasLitterRoot)
			tmp <- predict(mod, dfNew )
			ifelse( tmp>minLogLik, tmp, NA)
		}
		#mtrace(tmpf)
		tmp2 <- twPlot2DFun.contour( ds$tvrY, ds$tvrO, tmpf, xdiv=80, key.title=title(sub="Log-Like-\nlihood\n"), color.palette=function(n){rev(heat.colors(n))})
		tmp3 <- twPlot2DFun.contour( ds$tvrY, ds$tvrO, tmp2, xdiv=nrow(tmp2), key.title=title(sub="Log-Like-\nlihood\n"), color.palette=function(n){rev(heat.colors(n))})
		tmp3 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmp2, xdiv=nrow(tmp2), col=rev(heat.colors(20)) )
		# note the three humps corresponding to three chains
	
		# find the optimum solution
	}

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
	parOrig <- par()
	levelplot(rLogLik~kY*kO, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	levelplot(rLogLik~tvrY*tvrO, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	par(parOrig)

	#mtrace(plotMarginal2D)
	tmp <- plotMarginal2D( smp, "tvrY", "tvrO" )
	tmp <- plotConditional2D( smp, "tvrY", "tvrO" )
	
	smp2 <- Ys0O	
	smp2 <- cbind(smp, tvrY=1/smp[,"kY"], tvrO=1/smp[,"kO"])
	sampleSig2 <- apply(smp2,2,function(var){
			grain <- diff(range(var))/60
			round(var/grain)*grain
		})
	colnames(sampleSig)
	levelplot(rLogLik~tvrY*tvrO, data=as.data.frame(sampleSig), col.regions=rev(heat.colors(100)))
	levelplot(rLogLik~tvrY*tvrO, data=as.data.frame(sampleSig2), col.regions=rev(heat.colors(100)))
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



