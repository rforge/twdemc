#source(file.path(system.file(package="twKinrespTutorial"),"kinresp0.R"))

.tmp.f <- function(){
	#for development of package only:
	library(nlme)
	data(dsKinrespTut)
	source(file.path("R","kinrespModel1.R"))	# provided with package
	source(file.path("inst","kinresp0.R"))
}
ds24 <- subset(dsKinrespTut, experiment==17 & replicate==2 )
plot(resp~time,data=ds24,col=ds24$replicate)
	
#fit for initializing the run
fitLS6 <- gnls( resp ~ kinrespModel(x0,r0,mumax,time), ds24
	,params=x0+r0+mumax~1
	,start=c(x0=140, r0=0.1, mumax=0.24)
	,weights=varPower(fixed=0.5)
)
coef(fitLS6)
plot((resid(fitLS6,type="pearson"))~ds24$time)	# seems ok for weighted residuals
confint(fitLS6)

#------------ using package FME
# more general
.tmp.f <- function(){
	library(FME)
	#?modFit
	fResidualsFME  <- function(pars,ds,...) {
		# here we need to return vector of residuals 
		predResp <- kinrespModelVec( pars, kinrespModel, time=ds$time, ...)
		#(ds$resp-predResp)^2
		(ds$resp-predResp)
	}
	fCostResid <- function(pars,ds,...){
		resid <- fResidualsFME(pars,ds,...)
		sum(resid^2)
	}
	
	P <- modFit(f=fResidualsFME,p=coef(fitLS), ds=ds24)
	covarProp <- summary(P)$cov.scaled * 2.4^2/3
	
	#?modMCMC
	resMC1 <- modMCMC(f=fResidualsFME
		,p=P$par
		,jump=covarProp
		,niter=30000
		,lower=c(0,0,0)
		,var0=P$ms
		, updatecov=100
		#,ntrydr=2	# adaptive jumping distribution
		,wvar0=1
		,burninlength = 2000
		,ds=ds24
	)
	summary(resMC1)
	plot(resMC1)	# there is still drift in the fit
	pairs(resMC1)
	
	library("misc3d")
	tmp.s <- sample(nrow(resMC1$par),5000)
	sampleDf1 <- as.data.frame(resMC1$par[tmp.s,])
	sampleDf1$logLik <- apply(resMC1$par[tmp.s,], 1,fCostResid, ds=ds24 )
	sampleDf1$logLik[!is.finite(sampleDf1$logLik)] <- NA
	hist(sampleDf1$logLik)
	plot3d(sampleDf1$x0, sampleDf1$r0, sampleDf1$mumax
		,col=(heat.colors(100))[round(rescale(sampleDf1$logLik,to=c(1,100)))]
		,xlab="x0", ylab="r0", zlab="mumax"
	)
	
}


#-------------- using package twDEMC
# optimized for usage on a computer cluster, restarting, debugging
# integration with simulated annealing to avoid local minima
# no treatment of nuisance parameters (yet)
library(twDEMC)
#sfInit(parallel=TRUE, cpus=4); sfSource(file.path("R","kinrespModel1.R"))
Zinit <- initZtwDEMCNormal( coef(fitLS6), fitLS6$varBeta, nChain=2*4, nPops=2 )

# need to specify uncertainties of respiration
# take from gnls fit with weights=varPower(fixed=0.5)
ds24$sdResp <- attr(fitLS6$residuals, "std")

fCostMC2 <- function(
		### log-Likelihood of observations
		par		##<< parameters: numeric vector (x0,r0,mumax)
		, ds	##<< dataset with columns resp and time, and sdResp
		, ...
){
	predResp <- kinrespModelVec( par, kinrespModel, time=ds$time, ...)
	obs <- cbind( -1/2*sum(((ds$resp-predResp)/ds$sdResp)^2) )
	### log-Likelihood
}
print(fCostMC2( coef(fitLS), ds=ds24 ))	#test
	
#mtrace(twDEMCInt)
resMC2 <- twDEMCBatch( Zinit, nGen=1024, nPops=2
	, fLogLik=fCostMC2, argsFLogLik=list(ds=ds24)
	, nGenBurnin=500
)
plot(as.mcmc.list(resMC2))	# still trend in chain: still burnin, not yet converged to limiting distribution
resMC2 <- twDEMCBatch( resMC2, nGen=2*1024 )	#extend the run
plot(as.mcmc.list(resMC2))
plot(as.mcmc.list(thin(resMC2,start=700)))

res2 <- stackChains(thin(resMC2,start=700))
# 3d plot results
#library("Rcmdr")
#tmp.s <- sample(nrow(res2),min(nrow(res2),400))
#tmp.s <- sample(nrow(res2),min(nrow(res2),800))	# sample 800 rows
tmp.s <- sample(nrow(res2),min(nrow(res2),1200))
sampleDf <- as.data.frame(res2[tmp.s,])
plot( denL <- density(sampleDf$rLogLik))
fl <- approxfun(denL$x,denL$y)
lScaled <- rescale(sampleDf$rLogLik*fl(sampleDf$rLogLik),to=c(1,30))
#tmpGrid <- seq( min(sampleDf$rLogLik),max(sampleDf$rLogLik),length.out=81)
#lines( fl(tmpGrid)~tmpGrid, type="l" )
plot3d(sampleDf$x0, sampleDf$r0, sampleDf$mumax
	,col=rev(heat.colors(100))[round(rescale(sampleDf$rLogLik,to=c(1,100)))]
	,xlab="x0", ylab="r0", zlab="mumax"
)
plot( res2[,1]~res2[,"mumax"])	#marginal distribution
tmp <- optim( coef(fitLS6), fCostMC2, ds=ds24, control=list(maxit=2000,fnscale=-1))
tmp$par			# classical optimization

# marginal distribution of x0
#qplot( res2[,"x0"], geom="density", xlab="x0" )
#qplot( res2[,"mumax"], geom="density", xlab="mumax" )

#----------- +Gaussian prior for x0 
.tmp.f <- function(){
	# assume we get a value for x0 by a different method, e.g. fumigation extration
	# how do the results for mumax change?
	
	parms0 <- c(x0=570)
	parms0.sd <- parms0*0.08	# relative error of 8%
	invCovarParms0 <- diag(1/(parms0.sd)^2,nrow = length(parms0))
	dimnames(invCovarParms0) <- list( names(parms0), names(parms0) )
	
	
	fCostMC3 <- function(
		### Log-Likelihood of sum of sqared differences between y and pred
		parms		##<< parameters: numeric vector (x0,r0,mumax)
		, parms0
		, invCovarParms0=NULL
		, ds	##<< dataset with columns resp and time, and sdResp
	, ...
	){
		predResp <- kinrespModelVec( parms, kinrespModel, time=ds$time, ...)
		resObs1 <- -1/2*sum(((ds$resp-predResp)/ds$sdResp)^2)
		
		iPrior <- match(names(parms0), names(parms) )
		tmp.diffParms <- parms[iPrior] - parms0
		resPrior <- -1/2*as.numeric(t(tmp.diffParms) %*% invCovarParms0[iPrior,iPrior] %*% tmp.diffParms) 
		
		c(obs=resObs1, parms=resPrior)
		### Log-Likelihood of observations and prior
	}
	
	#mtrace(fCostMC3)
	fCostMC3( coef(fitLS), parms0=parms0, invCovarParms0=invCovarParms0, ds=ds24 ) #test
	
	#mtrace(twDEMCInt)
	resMC3 <- twDEMCBatch( Zinit, nGen=2*1024, nPops=2
		, fLogLik=fCostMC3, argsFLogLik=list(ds=ds24, parms0=parms0, invCovarParms0=invCovarParms0)
		, nGenBurnin=500
	)
	plotThinned(as.mcmc.list(resMC3))
	
	tmp3 <- stackChains(thin(resMC3,start=nGenBurnin))
	
	tmp.prior <- data.frame(x=seq(400,700,length.out=80))
	tmp.prior$y <- c(0,diff(pnorm( tmp.prior$x, mean=parms0["x0"], sd=parms0.sd["x0"]))/diff(tmp.prior$x))
	
	p1 <- qplot( res2[,"x0"], geom="density", xlab="x0" )
	p1 <- p1 + geom_density(aes(x=tmp3[,"x0"]), color="blue")
	p1 <- p1 + geom_line(aes(x=x,y=y), data=tmp.prior, color="magenta")
	print(p1)
	
	p2 <- qplot( res2[,"mumax"], geom="density", xlab="mumax" )
	p2 = p2 + geom_density(aes(x=tmp3[,"mumax"]), color="blue")
	print(p2)
	
	p3 <- qplot( res2[,"r0"], geom="density", xlab="r0" )
	p3 = p3 + geom_density(aes(x=tmp3[,"r0"]), color="blue")
	print(p3)
}

#----------- +Log-normal prior for mumax 
# Assume we have prior knowledge about mumax (from a different experiment, or from theory)
# How do the results for x0 change?

# define the prior
parms0 <- c(mumax=0.17)		
parms0.sd <- parms0*0.5	# relative error of 50%: vague prior
parms0.sigma <- log(1+(parms0.sd/parms0)^2)   
parms0.mu <- log(parms0)-1/2*parms0.sigma^2

# vizalize the prior
xr <- qlnorm( c(0.01,0.99), meanlog=parms0.mu, sdlog=parms0.sigma)
xgrid <- seq(xr[1],xr[2],length.out=30)
dx <- dlnorm(xgrid, meanlog=parms0.mu, sdlog=parms0.sigma)
plot( dx ~ xgrid, xlab="mumax", ylab="prior probability density")

fCostMC4 <- function(
	### Log-Likelihood of sum of sqared differences between y and pred and lognormal prior for mumax
	parms		##<< parameters: numeric vector (x0,r0,mumax)
	, mu,sigma	##<< priors for mumax-lognormal distribution
	, ds		##<< dataset with columns resp, time, and sdResp
	, ...
){
	predResp <- kinrespModelVec( parms, kinrespModel, time=ds$time, ...)
	resObs1 <- -1/2*sum(((ds$resp-predResp)/ds$sdResp)^2)
	
	mumax <- parms["mumax"]
	resPrior <- if( mumax <= 0 ) -Inf else {
		logMumax <- log(mumax)
		# note the formula for log-likelihood of lognormal distribution 		
		-(1/2*((logMumax-mu)^2)/sigma^2 + logMumax)	# take care to multiply by Jacobian (+logMumax)
	}
	c(obs=resObs1, parms=as.numeric(resPrior))
	### Log-Likelihood of observations and prior
}

#mtrace(fCostMC4)
fCostMC4( coef(fitLS6), mu=parms0.mu, sigma=parms0.sigma, ds=ds24 ) #test

#mtrace(twDEMCInt)
resMC4 <- twDEMCBatch( Zinit, nGen=2*1024, nPops=2
	, fLogLik=fCostMC4, argsFLogLik=list(mu=parms0.mu, sigma=parms0.sigma, ds=ds24)
	, nGenBurnin=500
)
plotThinned(as.mcmc.list(resMC4))
plotThinned(as.mcmc.list(thin(resMC4,start=resMC4$nGenBurnin)))

tmp4 <- stackChains(thin(resMC4,start=resMC4$nGenBurnin))

tmp.prior <- data.frame(x=xgrid,y=dx)
p1 <- ggplot( as.data.frame(res2),aes(x=mumax), xlab="mumax") 
p1 <- p1 + geom_line(aes(x=x,y=y, color="prior"), data=tmp.prior)
p1 <- p1 + geom_density(aes(x=tmp4[,"mumax"], color="posterior with prior"))
p1 <- p1 + geom_density( aes(color="posterior without prior") )
p1l <- p1 + scale_colour_manual("", c("black","blue","magenta") ) 
print(p1l)

p2 <- ggplot( as.data.frame(res2), aes(x=x0)) 
p2 <- p2 + geom_density(aes(color="implicit" ))
p2 <- p2 + geom_density(aes(x=tmp4[,"x0"], color="lognormal mumax"))
p2 <- p2 + scale_colour_manual(name="prior", c("blue","magenta") ) 
print(p2)
#significantly different

p3 <- ggplot( as.data.frame(res2), aes(x=r0)) 
p3 <- p3 + geom_density(aes(color="implicit" ))
p3 <- p3 + geom_density(aes(x=tmp4[,"r0"], color="lognormal mumax"))
p3 <- p3 + scale_colour_manual("prior", c("blue","magenta") ) 
print(p3)

#----------- +Logit-normal prior for r0 
# Assume we have prior knowledge about mumax (from a different experiment, or from theory)
# How do the results for x0 change?

parms0 <- twCoefLogitnormMLE( mle=0.007, quant=0.015)

xr <- qlogitnorm( c(0.01,0.99), mu=parms0[1], sigma=parms0[2])
xgrid <- seq(xr[1],xr[2],length.out=50)
dx <- dlogitnorm(xgrid, mu=parms0[1], sigma=parms0[2] )
plot( dx ~ xgrid)

fCostMC4 <- function(
	### Log-Likelihood of sum of sqared differences between y and pred and lognormal prior for mumax
	parms		##<< parameters: numeric vector (x0,r0,mumax)
	, parms0	##<< priors for r0-logitnormal distribution
	, ds		##<< dataset with columns resp, time, and sdResp
	, ...
){
	predResp <- kinrespModelVec( parms, kinrespModel, time=ds$time, ...)
	resObs1 <- -1/2*sum(((ds$resp-predResp)/ds$sdResp)^2)
	
	r0 <- parms["r0"]
	resPrior <- if( r0 <= 0 | r0 >= 1) -Inf else {
			mu <- parms0[1]
			sigma <- parms0[2]
			## XX implement log-Likelihood for logitnormal distribution
			0
		}
	c(obs=resObs1, parms=as.numeric(resPrior))
	### Log-Likelihood of observations and prior
}

#mtrace(fCostMC4)
fCostMC4( coef(fitLS6), parms0=parms0, ds=ds24 ) #test

#mtrace(twDEMCInt)
resMC4 <- twDEMCBatch( Zinit, nGen=2*1024, nPops=2
	, fLogLik=fCostMC4, argsFLogLik=list(parms0=parms0, ds=ds24)
	, nGenBurnin=500
)
plotThinned(as.mcmc.list(resMC4))

tmp4 <- stackChains(thin(resMC4,start=resMC4$nGenBurnin))

tmp.prior <- data.frame(x=xgrid,y=dx)
p3 <- qplot( res2[,"r0"], geom="density", xlab="r0" )
p3 = p3 + geom_density(aes(x=tmp4[,"r0"]), color="blue")
p3 = p3 + geom_line(aes(x=x,y=y),data=tmp.prior, color="maroon")
print(p3)

p1 <- qplot( res2[,"mumax"], geom="density", xlab="mumax" )
p1 <- p1 + geom_density(aes(x=tmp4[,"mumax"]), color="blue")
print(p1)

p2 <- qplot( res2[,"x0"], geom="density", xlab="x0" )
p2 = p2 + geom_density(aes(x=tmp4[,"x0"]), color="blue")
print(p2)


p1 <- ggplot( as.data.frame(res2),aes(x=r0)) 
p1 <- p1 + geom_line(aes(x=x,y=y, color="prior"), data=tmp.prior)
p1 <- p1 + geom_density(aes(x=tmp4[,"r0"], color="posterior with prior"))
p1 <- p1 + geom_density( aes(color="posterior without prior") )
p1l <- p1 + scale_colour_manual("", c("black","blue","magenta") ) 
print(p1l)

p2 <- ggplot( as.data.frame(res2), aes(x=x0)) 
p2 <- p2 + geom_density(aes(color="implicit" ))
p2 <- p2 + geom_density(aes(x=tmp4[,"x0"], color="logitnormal r0"))
p2 <- p2 + scale_colour_manual(name="prior", c("blue","magenta") ) 
print(p2)
#significantly different

p3 <- ggplot( as.data.frame(res2), aes(x=mumax)) 
p3 <- p3 + geom_density(aes(color="implicit" ))
p3 <- p3 + geom_density(aes(x=tmp4[,"mumax"], color="logitnormal r0"))
p3 <- p3 + scale_colour_manual("prior", c("blue","magenta") ) 
print(p3)



#----------- Multiple data streams
# Assume we get have other observational data by an substrate monod-curve
