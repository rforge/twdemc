# including varying litter inputs 

.tmp.init <- function(){
	# for plotting 
	source("../twDEMC/R/plotTwDEMC.R")
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
		#, fCalcBiasedInput=meanInput    ##<< function(input,padj,...){obs} possibility to account for bias and to optimize bias parameters
		, fCalcBiasedInput=meanInputFluctuating	
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

mdi.kYkObiasLitter <- function(){
	#poptnames <- c("kY","kO","biasLitterRoot","biasLitterLeaf")
	poptnames <- c("kY","kO","biasLitterRoot","biasLitterLeaf","tLagLeaf","tLagRoot")
	poptDistr <- argsFLogLik$poptDistr <- twConstrainPoptDistr(poptnames, HowlandParameterPriors$parDistr)
	sfExport("argsFLogLik")
	
	pnorm <- transNormPopt(unlist(parms0), HowlandParameterPriors$parDistr$trans[names(parms0)])
	normpopt <- pnorm[poptnames]
	#mtrace(of.howlandSteady)
	resOf <- sfRemoteWrapper( normpopt=normpopt, remoteFun=of.howlandSteady, remoteFunArgs=argsFLogLik)
	
	tmp.f <- function(){
		# optimization makes no sense with stochastic model
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
		#plot( rowSums(res[,c("respY_c12","respO_c12")]) ~ res[,"time"], type="l" )
		#points( obs ~ times, data=as.data.frame(obsSite$respCum), col="maroon" )
	}
	
	#------ explore posterior with MCMC using prior
	tmp.fcovarPrior <- function(){
		set.seed(0815)
		covMat <- poptDistr$sigma    
		.nPops=3
		Zinit <- initZtwDEMCNormal( resOpt$par, covMat, nChains=4*.nPops, nPops=.nPops)
		#options(error=recover)
		resMC <- twDEMCBatch( Zinit, nGen=500, fLogLik=of.howlandSteady, argsFLogLik=argsFLogLik, nPops=.nPops
			#,debugSequential=TRUE
		)
		matplot(resMC$pAccept, type="l")
		plot(as.mcmc.list(resMC))
		resMC <- twDEMCBatch( resMC, nGen=5000, doRecordProposals=TRUE )
		plot(as.mcmc.list(resMC))
		matplot( resMC$rLogLik[-(1:10),], type="l" )
		#resMC <- twDEMCBatch( resMC, nGen=2000 )
		#plot(as.mcmc.list(resMC))
		resMCB <- thin(resMC, start=200, newThin = 16)
		plot(as.mcmc.list(resMCB))
		matplot( resMCB$rLogLik, type="l" )
		plotChainPopMoves(resMCB)
		save(resMCB,file=file.path("tmp","resMCB_steady2_4.RData"))
	}
	load(file.path("tmp","resMCB_steady2_4.RData"))	#resMCB
	
	sampleN <- sample <-  stackChains(resMCB)
	minLogLik <- quantile(sampleN[,1], probs=c(0.05) )	# empirical 95%
	sampleN0 <- sample0 <- sampleN[ sampleN[,1] >= minLogLik, ]
	minLogLik2 <- getRLogLikQuantile(sampleN) 	# theoretical criterion
	c(minLogLik,minLogLik2)
	sampleN02 <- sample02 <- sampleN[ sampleN[,1] >= minLogLik2, ]
	sample[,-1] <- transOrigPopt(sampleN[,-1],  poptDistr=poptDistr$trans)
	sample0[,-1] <- transOrigPopt(sampleN0[,-1],  poptDistr=poptDistr$trans)
	sample02[,-1] <- transOrigPopt(sampleN02[,-1],  poptDistr=poptDistr$trans)
	
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
	
	windows(width=6.4,height=4.4,pointsize=10, record=TRUE)
	twPairs(sampleN0)		# no big correlations
	p1 <- ggplotDensity.twDEMC(resMCB, poptDistr=poptDistr)
	#print(p1)
	#mtrace(ggplotDensity.twDEMC)
	p2 <- ggplotDensity.twDEMC(resMCB, poptDistr=poptDistr, doTransOrig=TRUE, doDispLogLik=FALSE)
	#print(p2)

	
	#--------- estimate a response surface of rLoglik
	#mtrace(plotMarginal2D)
	tmpMarg <- plotMarginal2D( sample0tvr, "tvrY", "tvrO" )
	#tmp <- plotConditional2D( smp, "tvrY", "tvrO" )

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
		#tmp2 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=bw3$bw*8;tmp}), xdiv=60, col=rev(heat.colors(20)) )
		tmp1 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=c(200,200);tmp}), xdiv=60, col=rev(heat.colors(20)) )
		tmp2 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmpf, argsFUN=list(bws={tmp<-bw3;tmp$bw=c(300,300);tmp}), xdiv=60, col=rev(heat.colors(20)) )
		tmp3 <- twPlot2DFun.contour( ds$tvrY, ds$tvrO, tmp2, xdiv=nrow(tmp2), xlab="tvrY",ylab="tvrO", key.title=title(sub="Log-Like-\nlihood\n\n"), color.palette=function(n){rev(heat.colors(n))})
		tmp3 <- twPlot2DFun( ds$tvrY, ds$tvrO, tmp2, xdiv=nrow(tmp2), col=rev(heat.colors(20)) )
		# note the three humps corresponding to three chains		
	}
}



