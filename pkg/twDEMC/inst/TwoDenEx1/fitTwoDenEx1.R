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

pred <- pred1 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
plot( pred$y2 ~ twTwoDenEx1$obs$y2 ); abline(0,1) # seems very good
plot( pred$y1 ~ twTwoDenEx1$obs$y1 ); abline(0,1) # far off 

#---------  fit both data streams in one unweighted density
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


pred <- pred2 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
plot( pred$y2 ~ twTwoDenEx1$obs$y2 ); abline(0,1) # not as good but ok
plot( pred$y1 ~ twTwoDenEx1$obs$y1 ); abline(0,1) # still no relation 

#----------- 2b try to infer by weights
res <- res2
lDen <- stackChains(res$resLogDen)
medianDen <- apply(lDen,2,function(lDenI){ median(lDenI - max(lDenI)) })
weights <- wMedian <- medianDen/min(medianDen)

resa <- resa2b <-  concatPops( resBlock <- twDEMCBlock( ZinitPops, nGen=.nGen, 
		dInfos=list(dBoth=list(fLogDen=.denBoth, argsFLogDen=list(weights=weights)))
		,nPop=.nPop
		,controlTwDEMC=list(thin=4)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	))
plot( as.mcmc.list(resa), smooth=FALSE)
#mtrace(subset.twDEMC)
res2b <- res <- thin(resa, start=100)
plot( as.mcmc.list(res), smooth=FALSE)
ss <- stackChains(res)
twTwoDenEx1$thetaTrue
(thetaBest <- thetaBest2 <- ss[ which.max(ss[,1]), -1])
(.qq <- apply(ss[,-1],2,quantile, probs=c(0.025,0.5,0.975) ))
# still biased downwards

pred <- pred2b <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
plot( pred$y2, twTwoDenEx1$obs$y2 ); abline(0,1) # mismatch evident
plot( pred$y1, twTwoDenEx1$obs$y1 ); abline(0,1) # at least some relation 

#XX check same magnitude of resLogDen
res <- res2b
lDen <- stackChains(res$resLogDen)
lDenDiff <- apply(lDen,2,function(lDenI){ lDenI - max(lDenI) })
q10lDenDiff <- apply(lDenDiff,2,quantile,0.1)
lDenDiffM <- melt(lDenDiff)
lDenDiffM2 <- lDenDiffM[
	((lDenDiffM$resComp=="y1") & (lDenDiffM$value >= q10lDenDiff[1])) |
	((lDenDiffM$resComp=="y2") & (lDenDiffM$value >= q10lDenDiff[2])) 
		,]
boxplot( value ~ resComp, lDenDiffM2 )
if( require(ggplot2) ){
	p1 <- ggplot(lDenDiffM2, aes(y=value, x=resComp)) + geom_boxplot() +
		opts(axis.title.x=theme_blank())
}

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

pred <- pred3 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
plot( pred$y2 ~ twTwoDenEx1$obs$y2 ); abline(0,1) # here the mismatch becomes clear
plot( pred$y1 ~ twTwoDenEx1$obs$y1 ); abline(0,1) # not super but relation existing 


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
plot( ss[bo,"b"] ~ ss[bo,"a"], col=denCol, xlab="a", ylab="b")
# red (rich) increases with higher b along the valley, blue (sparce) constaines a

plot3d( ss[bo,"a"], ss[bo,"b"], lDen[bo], col=denCol )


#----------- ggplot of violin plots of the distributions 
require(ggplot2)
.tmp.f <- function(){
	# creating violin plots with ggplot2
	#http://markmail.org/search/?q=list%3Ar-project+violin+plots+ggplot#query:list%3Ar-project%20violin%20plots%20ggplot+page:1+mid:ka6kymn6agvvdakg+state:results
	p <- rep(c(rep("condition_a", 4), rep("condition_b", 4)), 2) 
	q <- c(rep("grp_1", 8), rep("grp_2", 8)) 
	r <- sample(1:5, 16, rep = T) 
	d <- data.frame(p, q, r) 
	ggplot(d, aes(x = r)) + 
		geom_ribbon( aes(ymax = ..density.., ymin = -..density..), stat = "density") + 
		facet_grid(q ~ p) + 
		coord_flip()
	
	ggplot(d, aes(x = r, fill = p))	+ 
		geom_ribbon(alpha = 0.2, aes(ymax =	..density.., ymin = -..density..), stat = "density") + 
		facet_wrap( ~ q) +
		coord_flip()
}

.nSample <- 128
dfDen <- rbind(
	cbind( data.frame( scenario="R", {tmp <- stackChains(res1)[,-1]; tmp[ seq(1,nrow(tmp),length.out=.nSample),] }))
	,cbind( data.frame( scenario="RS", {tmp <- stackChains(res2)[,-1]; tmp[ seq(1,nrow(tmp),length.out=.nSample),] }))
	,cbind( data.frame( scenario="RSw", {tmp <- stackChains(res2b)[,-1]; tmp[ seq(1,nrow(tmp),length.out=.nSample),] }))
	,cbind( data.frame( scenario="D", {tmp <- stackChains(res3)[,-c(1:2)]; tmp[ seq(1,nrow(tmp),length.out=.nSample),] }))
)
dfDenM <- melt(dfDen)
str(dfDenM)

StatDensityNonZero <- StatDensity$proto( calculate<- function(.,data, scales, dMin=0.001, ...){
		res <- StatDensity$calculate(data,scales,...)
		# append a zero density at the edges
		if( res$density[1] > 0) res <- rbind( data.frame(x=res$x[1]-diff(res$x[1:2]),density=0,scaled=0,count=0),res)
		if( res$density[nrow(res)] > 0) res <- rbind( res, data.frame(x=res$x[nrow(res)]+diff(res$x[nrow(res)-c(1,0)]),density=0,scaled=0,count=0))
		bo <- (res$density < dMin) & (c(res$density[-1],0) < dMin) & (c(0,res$density[-length(res$density)]) < dMin)
		res$density[ bo ] <- NA
		res$count[ bo ] <- NA
		res$scaled[ bo ] <- NA
		res
	}, required_aes=c("x"))
#with(StatDensityNonZero, mtrace(calculate)) 
#with(StatDensityNonZero, mtrace(calculate,F)) 

stat_densityNonZero <- function (
	### constructs a new StatPrior statistics based on aesthetics x and parName
	mapping = NULL, data = NULL, geom = "line", position = "stack", 
	adjust = 1,	...
){ 
	StatDensityNonZero$new(mapping = mapping, data = data, geom = geom, 
		position = position, adjust = adjust, ...)
}

#pgm <- geom_ribbon( alpha=0.8, aes(ymax = ..density.., ymin = -..density..), stat = "density")
#pgm <- geom_ribbon(  aes(ymax = ..density.., ymin = -..density..), stat = "density")
#pgm <- stat_densityNonZero(  aes(ymax = -..density..),  size=1, geom="line")
pgm <- geom_line(aes(y = +..scaled..), stat="densityNonZero", size=1)
pgm2 <- geom_line(aes(y = -..scaled..), stat="densityNonZero", size=1) 
optsm <- opts(axis.title.x = theme_blank(), axis.text.y = theme_blank(), axis.ticks.y = theme_blank() ) 
pa <- ggplot(dfDen, aes(x = a, colour=scenario, linetype=scenario)) + pgm + pgm2 + optsm + scale_y_continuous('a')
pb <- ggplot(dfDen, aes(x = b, colour=scenario, linetype=scenario)) + pgm + pgm2 + optsm + scale_y_continuous('b')
pb
windows(width=7, height=3)
grid.newpage()
pushViewport( viewport(layout=grid.layout(2,1)))
print(pa , vp = viewport(layout.pos.row=1,layout.pos.col=1))	
print(pb + opts(legend.position = "none") , vp = viewport(layout.pos.row=2,layout.pos.col=1))	

#----------- ggplot predictive posterior
#str(pred1)
dfPred <- rbind(
	 data.frame( scenario="R", variable="y1", value = pred1$y1 )
	,data.frame( scenario="R", variable="y2", value = pred1$y2 )
	,data.frame( scenario="RS", variable="y1", value = pred2$y1 )
	,data.frame( scenario="RS", variable="y2", value = pred2$y2 )
	,data.frame( scenario="RSw", variable="y1", value = pred2b$y1 )
	,data.frame( scenario="RSw", variable="y2", value = pred2b$y2 )
	,data.frame( scenario="D", variable="y1", value = pred3$y1 )
	,data.frame( scenario="D", variable="y2", value = pred3$y2 )
)
dfPred$observations <- NA_real_
dfPred$observations[dfPred$variable=="y1"] <- twTwoDenEx1$obs$y1
dfPred$observations[dfPred$variable=="y2"] <- twTwoDenEx1$obs$y2
str(dfPred)

p1 <- ggplot(dfPred, aes(x=value, y=observations, colour=scenario) ) + geom_point() + 
	facet_wrap( ~ variable, scales="free") +
	opts(axis.title.x=theme_blank() ) +
	geom_abline(colour="black") +
	c()
p1


