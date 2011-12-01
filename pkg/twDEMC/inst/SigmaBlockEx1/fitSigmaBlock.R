data(twLinreg1)

denSigmaEx1Obs <- function(
	### unnormalized log density for observations for given parameters
	theta				 ##<< named numeric vector a,b,logSigma2
	,twLinreg=twLinreg1  ##<< list with components xval and obs 
){
	pred <- twLinreg$fModel(theta[1:2],twLinreg$xval)
	resid <- pred - twLinreg$obs
	varObs <- exp(theta[3])
	-1/2 * sum(resid^2)/varObs
}
attr(denSigmaEx1Obs,"ex") <- function(){
	data(twLinreg1)
	#mtrace(denSigmaEx1Obs)
	logSigma2 <- 2*log( mean(twLinreg1$sdObs) )
	denSigmaEx1Obs( c( twLinreg1$theta0, logSigma2=logSigma2 ) )
	
	#likelihood profile of a and b conditional on given logSigma2
	lmDummy <- with(twLinreg1, lm( obs ~ xval, weights=1/sdObs^2))		# results without priors
	(.expTheta <- coef(lmDummy))
	(.expCovTheta <- vcov(lmDummy))		# a is very weak constrained, negative covariance between a an b
	(.expSdTheta <- structure(sqrt(diag(.expCovTheta)), names=c("a","b")) )
	nGrid <- 50
	aGrid <- .expTheta[1] + .expSdTheta[1]*seq(-2, 2, length.out=nGrid)	
	bGrid <- .expTheta[2] + .expSdTheta[2]*seq(-2, 2, length.out=nGrid)
	abGrid <- expand.grid(a=aGrid,b=bGrid)
	Zinit <- cbind(abGrid, logSigma2=logSigma2)
	logLikab <- apply(Zinit, 1, denSigmaEx1Obs  )
	logLikabm <- matrix(logLikab, nrow=nGrid)
	logLikabm[ logLikabm < max(logLikabm)-3 ] <- NA	# do not color regions with very low likelihood
	image( aGrid, bGrid, logLikabm, col=rev(heat.colors(80)) )
}


denSigmaEx1Sigma <- function(
	### unnormalized log density for observations uncertainty logSigma2 for given parameters
	theta				 ##<< named numeric vector a,b,logSigma2
	,twLinreg=twLinreg1 ##<< list with components fModel, xval and obs
){
	#see http://www.wutzler.net/reh/index.php?title=twutz:Gelman04_1#normal_distribution_with_known_mean_but_unknown_variance
	pred <- twLinreg$fModel(theta[1:2],twLinreg$xval)
	resid <- pred - twLinreg$obs
	n <- length(pred)
	logSigma2 <- theta[3]
	sigma2 <- exp(logSigma2)
	-1/2*( n*logSigma2  + sum(resid^2)/sigma2 )
}
attr(denSigmaEx1Sigma,"ex") <- function(){
	data(twLinreg1)
	#mtrace(denSigmaEx1Sigma)
	logSigma2 <- 2*log( mean(twLinreg1$sdObs) )
	denSigmaEx1Sigma( c( twLinreg1$theta0, logSigma2=logSigma2 ) )
	
	#likelihood profile of varSigma conditional on the true parameters
	Zinit <- cbind( matrix(twLinreg1$thetaTrue,nrow=41,ncol=2,byrow=TRUE, dimnames=list(case=NULL,par=names(twLinreg1$thetaTrue)))
	, logSigm2=logSigma2*seq(0.8,1.2,length.out=41))
	head(Zinit)
	logLikSigma <- apply( Zinit, 1, denSigmaEx1Sigma )
	bo <- logLikSigma > max(logLikSigma)-3
	plot( logLikSigma[bo] ~ exp(Zinit[bo,3]/2) )
	abline(h=max(logLikSigma)-2)
}

denSigmaEx1SigmaInvX <- function(
	### unnormalized log density for observations uncertainty logSigma2 for given parameters using invchisq distribution
	theta				 ##<< named numeric vector a,b,logSigma2
	,twLinreg=twLinreg1 ##<< list with components fModel, xval and obs
){
	#see http://www.wutzler.net/reh/index.php?title=twutz:Gelman04_1#normal_distribution_with_known_mean_but_unknown_variance
	# this time calculating the density by a scaled inverse Chi-Square distribution
	pred <- twLinreg$fModel(theta[1:2],twLinreg$xval)
	resid <- pred - twLinreg$obs
	n <- length(pred)
	v <- sum( resid^2 ) / n
	logSigma2 <- theta[3]
	sigma2 <- exp(logSigma2)	
	dinvchisq( sigma2, n, v, log=TRUE )
}
attr(denSigmaEx1SigmaInvX,"ex") <- function(){
	data(twLinreg1)
	#mtrace(denSigmaEx1SigmaInvX)
	logSigma2 <- 2*log( mean(twLinreg1$sdObs) )
	(tmp2 <- denSigmaEx1SigmaInvX( c( twLinreg1$theta0, logSigma2=logSigma2 ) ))
	
	#likelihood profile of varSigma conditional on the true parameters
	Zinit <- cbind( matrix(twLinreg1$thetaTrue,nrow=41,ncol=2,byrow=TRUE, dimnames=list(case=NULL,par=names(twLinreg1$thetaTrue)))
		, logSigma2=logSigma2*seq(0.8,1.2,length.out=41))
	head(Zinit)
	logLikSigma <- apply( Zinit, 1, denSigmaEx1SigmaInvX )
	bo <- logLikSigma > max(logLikSigma)-3
	plot( logLikSigma[bo] ~ exp(Zinit[bo,3]/2) )
	abline(h=max(logLikSigma)-2)
}



logSigma2 <- 2*log( mean(twLinreg1$sdObs) )		# expected sigma2
logSigma2 <- logSigma2+2	# bad starting conditions: assume bigger variance
varLogSigma2=2#var( 2*log(twLinreg1$sdObs) )		# variance for initialization

lmDummy <- with(twLinreg1, lm( obs ~ xval, weights=1/sdObs^2))		# results without priors
(.expTheta <- coef(lmDummy))
(.expCovTheta <- vcov(lmDummy))		# a is very weak constrained, negative covariance between a an b
(.expSdTheta <- structure(sqrt(diag(.expCovTheta)), names=c("a","b")) )


.nPop=2
.nChainPop=4
ZinitPops <- with(twLinreg1, initZtwDEMCNormal( c(theta0,logSigma2=logSigma2), diag(c(.expSdTheta^2,varLogSigma2)), nChainPop=.nChainPop, nPop=.nPop))
#dim(ZinitPops)
#head(ZinitPops[,,1])
theta0 <- ZinitPops[nrow(ZinitPops),,1]

do.call( denSigmaEx1Obs, c(list(theta=theta0)) )
do.call( denSigmaEx1Sigma, c(list(theta=theta0)) )

.nGen=128
#.nGen=16
#mtrace(twDEMCBlockInt)
resa <-  concatPops( resBlock <- twDEMCBlock( ZinitPops, nGen=.nGen, 
		dInfos=list(
			dObs=list(fLogDen=denSigmaEx1Obs, argsFLogDen=argsFLogDen1)
			,dSigma=list(fLogDen=denSigmaEx1Sigma, argsFLogDen=argsFLogDen2)
		)
		,blocks=list(
			bObs=list(dInfoPos="dObs", compPos=c("a","b"))
			,bSigma=list(dInfoPos="dSigma", compPos=c("logSigma2"))
		)
		,nPop=.nPop
		,controlTwDEMC=list(thin=4)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	))
#str(resBlock$pops[[1]])
#str(resa)
plot( as.mcmc.list(resa), smooth=FALSE )
res <- thin(resa, start=64)
plot( rescoda <- as.mcmc.list(res), smooth=FALSE )

matplot(res$pAccept[,"bObs",], type="l")
matplot(res$pAccept[,"bSigma",], type="l")
matplot(res$logDen[,"dObs",], type="l")
matplot(res$logDen[,"dSigma",], type="l")
matplot(resa$parms[,"logSigma2",], type="l")
matplot(res$parms[,"logSigma2",], type="l")
matplot(exp(0.5*res$parms[,"logSigma2",]), type="l")

#str(summary(rescoda))
suppressWarnings({	
		summary(rescoda)$statistics[,"Mean"]
		summary(rescoda)$statistics[,"SD"]
		(.popmean <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"Mean"]}))
		(.popsd <- lapply(list(p1=1:4,p2=5:8),function(i){summary(rescoda[i])$statistics[,"SD"]}))
	})

# 1/2 orders of magnitude around prescribed sd for theta
.pop=1
for( .pop in seq(along.with=.popsd) ){
	checkMagnitude( .expSdTheta, .popsd[[.pop]][1:2] )
}

# check that thetaTrue is in 95% interval 
.pthetaTrue <- sapply(1:2, function(.pop){
		pnorm(.expTheta, mean=.popmean[[.pop]][1:2], sd=.popsd[[.pop]][1:2])
	})
checkInterval( .pthetaTrue ) 



#---------- now replace metropolis sampling of sigma by direct parameter sampling



