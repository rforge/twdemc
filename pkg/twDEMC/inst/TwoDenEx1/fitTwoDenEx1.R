.denSparce <- function(
	### density of sparce observations
	theta
	,twTwoDenEx=twTwoDenEx1
	,theta0=twTwoDenEx$thetaTrue
	,...
){
	theta0[names(theta)] <- theta
	pred <- twTwoDenEx$fModel(theta0, xSparce=twTwoDenEx$xSparce, xRich=twTwoDenEx$xRich, ...)
	misfit <- twTwoDenEx$obs$y1 - pred$y1
	-1/2 * sum((misfit/twTwoDenEx$sdObs$y1)^2)
}

.denRich <- function(
	### density of sparce observations
	theta=theta0
	,twTwoDenEx=twTwoDenEx1
){
	pred <- twTwoDenEx$fModel(theta, xSparce=twTwoDenEx$xSparce, xRich=twTwoDenEx$xRich)
	misfit <- twTwoDenEx$obs$y2 - pred$y2
	-1/2 * sum((misfit/twTwoDenEx$sdObs$y2)^2)
}

.denBoth <- function(
	### density of sparce observations
	theta
	,twTwoDenEx=twTwoDenEx1
	,weights=c(1,1)			# weights for the two data streams
){
	pred <- twTwoDenEx$fModel(theta, xSparce=twTwoDenEx$xSparce, xRich=twTwoDenEx$xRich)
	misfit <- twTwoDenEx$obs$y1 - pred$y1
	d1 <- -1/2 * sum((misfit/twTwoDenEx$sdObs$y1)^2)
	misfit <- twTwoDenEx$obs$y2 - pred$y2
	d2 <- -1/2 * sum((misfit/twTwoDenEx$sdObs$y2)^2)
	db <- c(y1=d1,y2=d2)*weights
	db
}
.denBoth(twTwoDenEx1$thetaTrue)


#--------- fit only a to the long term observations
.nPop=2
.nChainPop=4
ZinitPopsA <- with(twTwoDenEx1, initZtwDEMCNormal( thetaMean["a"], diag((thetaMean["a"]*0.3)^2,nrow=1), nChainPop=.nChainPop, nPop=.nPop))
#plot(density(ZinitPopsA[,"a",]))
#dim(ZinitPopsA)
#head(ZinitPopsA[,,1])
theta0 <- adrop(ZinitPopsA[nrow(ZinitPopsA),,1 ,drop=FALSE],3)

.nGen=512
#.nGen=16
#mtrace(twDEMCBlockInt)

# with correct b
resa1 <- resa <- concatPops( resBlock <- twDEMCBlock( ZinitPopsA, nGen=.nGen, 
		dInfos=list(dSparce=list(fLogDen=.denSparce, argsFLogDen=list(theta0=c(a=1,b=2))))
		,nPop=.nPop
		,controlTwDEMC=list(thin=4)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	))
plot( as.mcmc.list(resa), smooth=FALSE)

# with biased b
resa1 <- resa <- concatPops( resBlock <- twDEMCBlock( ZinitPopsA, nGen=.nGen, 
		dInfos=list(dSparce=list(fLogDen=.denSparce, argsFLogDen=list(theta0=c(a=1,b=1.6))))
		,nPop=.nPop
		,controlTwDEMC=list(thin=4)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	))
plot( as.mcmc.list(resa), smooth=FALSE)

#---------  fit only the short term observations to both parameters
thetaMean <- twTwoDenEx1$thetaTrue
.nPop=2
.nChainPop=4
ZinitPops <- with(twTwoDenEx1, initZtwDEMCNormal( thetaMean, diag((thetaMean*0.3)^2), nChainPop=.nChainPop, nPop=.nPop))
#plot(density(ZinitPops[,"a",]))
#plot(density(ZinitPops[,"b",]))
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
ss <- stackChains(res)
twTwoDenEx1$thetaTrue
(thetaBest <- thetaBest1 <- ss[ which.max(ss[,1]), -1])
(.qq <- apply(ss[,-1],2,quantile, probs=c(0.025,0.5,0.975) ))
# both biased downwards

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
ss <- stackChains(res)
twTwoDenEx1$thetaTrue
(thetaBest <- thetaBest2 <- ss[ which.max(ss[,1]), -1])
(.qq <- apply(ss[,-1],2,quantile, probs=c(0.025,0.5,0.975) ))
# still biased downwards

# clearly not effective to match short data stream
pred <- pred2 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
pred1 <- with( twTwoDenEx1, fModel(thetaBest1, xSparce=xSparce, xRich=xRich) )
plot( pred$y1 ~ twTwoDenEx1$obs$y1 ); abline(0,1)
points( pred1$y1 ~ twTwoDenEx1$obs$y1, col="blue")

#----------- fit b to shortterm and a to longterm observations
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
ss <- stackChains(res)
twTwoDenEx1$thetaTrue
(thetaBest <- thetaBest3 <- ss[ which.max(rowSums(ss[,1:2])), -(1:2)])
(.qq <- apply(ss[,-(1:2)],2,quantile, probs=c(0.025,0.5,0.975) ))
# a a bit biased upwards (effect of b)

plot(density(stackChains(res2)[,"a"]), xlim=c(0.8,1.2))
lines( density(stackChains(res3)[,"a"]), col="red")

pred <- pred3 <- with( twTwoDenEx1, fModel(thetaBest3, xSparce=xSparce, xRich=xRich) )
plot( pred2$y1 ~ twTwoDenEx1$obs$y1, ylim=range(c(pred1$y1,pred2$y1,pred3$y1)) ); abline(0,1)
points( pred1$y1 ~ twTwoDenEx1$obs$y1, col="green")
points( pred3$y1 ~ twTwoDenEx1$obs$y1, col="red")


plot( ss[,"dSparce"] ~ ss[,"a"] )
plot( ss[,"dRich"] ~ ss[,"a"] )
plot( rowSums(ss[,1:2]) ~ ss[,"a"] )
plot( ss[,"dRich"] ~ ss[,"b"] )
plot( ss[,"b"] ~ ss[,"a"], col=rev(heat.colors(100))[round(twRescale(ss[,"dRich"],c(1,100)))])
plot( ss[,"b"] ~ ss[,"a"], col=rev(heat.colors(100))[round(twRescale(ss[,"dSparce"],c(1,100)))])

lDen <- twRescale(ss[,"dSparce"])+twRescale(ss[,"dRich"])
bo <- (lDen > quantile(lDen,0.1)) 
denCol <- rgb( 
	twRescale(ss[bo,"dRich"])
	,0
	,twRescale(ss[bo,"dSparce"])
)
plot( ss[bo,"b"] ~ ss[bo,"a"], col=denCol)
# red (rich) increases with higher b along the valley, blue (sparce) constaines a

plot3d( ss[bo,"a"], ss[bo,"b"], lDen[bo], col=denCol )
