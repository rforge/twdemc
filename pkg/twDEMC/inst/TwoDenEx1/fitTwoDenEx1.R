.denSparce <- function(
	### density of sparce observations
	theta
	,twTwoDemEx=twTwoDenEx1
){
	pred <- twTwoDemEx$fModel(theta, xSparce=twTwoDemEx$xSparce, xRich=twTwoDemEx$xRich)
	misfit <- twTwoDemEx$obs$y1 - pred$y1
	-1/2 * sum((misfit/twTwoDemEx$sdObs$y1)^2)
}

.denRich <- function(
	### density of sparce observations
	theta=thetaEff
	,twTwoDemEx=twTwoDenEx1
){
	pred <- twTwoDemEx$fModel(theta, xSparce=twTwoDemEx$xSparce, xRich=twTwoDemEx$xRich)
	misfit <- twTwoDemEx$obs$y2 - pred$y2
	-1/2 * sum((misfit/twTwoDemEx$sdObs$y2)^2)
}

.denBoth <- function(
	### density of sparce observations
	theta
	,twTwoDemEx=twTwoDenEx1
	,weights=c(1,1)			# weights for the two data streams
){
	pred <- twTwoDemEx$fModel(theta, xSparce=twTwoDemEx$xSparce, xRich=twTwoDemEx$xRich)
	misfit <- twTwoDemEx$obs$y1 - pred$y1
	d1 <- -1/2 * sum((misfit/twTwoDemEx$sdObs$y1)^2)
	misfit <- twTwoDemEx$obs$y2 - pred$y2
	d2 <- -1/2 * sum((misfit/twTwoDemEx$sdObs$y2)^2)
	db <- c(y1=d1,y2=d2)*weights
	db
}
.denBoth(twTwoDenEx1$thetaEff)



#---------  fit only the short term observations
.nPop=2
.nChainPop=4
ZinitPops <- with(twTwoDenEx1, initZtwDEMCNormal( thetaEff, diag((thetaEff*0.3)^2), nChainPop=.nChainPop, nPop=.nPop))
plot(density(ZinitPops[,"a",]))
plot(density(ZinitPops[,"b",]))
#dim(ZinitPops)
#head(ZinitPops[,,1])
theta0 <- ZinitPops[nrow(ZinitPops),,1]

.nGen=512
#.nGen=16
#mtrace(twDEMCBlockInt)
resa1 <- resa <- concatPops( resBlock <- twDEMCBlock( ZinitPops, nGen=.nGen, 
		dInfos=list(dRich=list(fLogDen=.denRich))
		,nPop=.nPop
		,controlTwDEMC=list(thin=4)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	))
plot( as.mcmc.list(resa), smooth=FALSE)
res1 <- res <- thin(resa, start=200)
plot( as.mcmc.list(res), smooth=FALSE)
apply(stackChains(res)[,-1],2,quantile, probs=c(0.025,0.5,0.975) )
twTwoDenEx1$thetaEff
# wide cf intervals a biased downwards, b biased upwards

#---------  fit both data streams in one density
resa <- resa2 <-  concatPops( resBlock <- twDEMCBlock( ZinitPops, nGen=.nGen, 
		dInfos=list(dBoth=list(fLogDen=.denBoth))
		,nPop=.nPop
		,controlTwDEMC=list(thin=4)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	))
plot( as.mcmc.list(resa), smooth=FALSE)
#mtrace(subset.twDEMC)
res2 <- res <- thin(resa, start=100)
plot( as.mcmc.list(res), smooth=FALSE)
apply(stackChains(res)[,-1],2,quantile, probs=c(0.025,0.5,0.975) )
twTwoDenEx1$thetaEff
# much more constrained but a biased downwards, b biased upwards

#----------- fit b to shortterm and a to longterm observations
.nGen=1024
resa <- resa3 <- concatPops( resBlock <- twDEMCBlock( 
		res2$parms
		, nGen=.nGen 
		,dInfos=list(
			dSparce=list(fLogDen=.denSparce)
			,dRich=list(fLogDen=.denRich)
		)
		,blocks = list(
			a=list(dInfoPos="dSparce", compPos="a")
			,b=list(dInfoPos="dRich", compPos="b")
		)
		,nPop=.nPop
		,controlTwDEMC=list(thin=8)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	))
plot( as.mcmc.list(resa), smooth=FALSE)
#mtrace(subset.twDEMC)
res3 <- res <- thin(resa, start=100)
plot( as.mcmc.list(res), smooth=FALSE)
apply(stackChains(res)[,-(1:2)],2,quantile, probs=c(0.025,0.5,0.975) )
twTwoDenEx1$thetaEff
# a biased downwards, b biased upwards

plot(density(stackChains(res2)[,"a"]))
lines( density(stackChains(res3)[,"a"]), col="red")
plot(density(stackChains(res2)[,"b"]))
lines( density(stackChains(res3)[,"b"]), col="red")


