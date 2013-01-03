data(twTwoDenEx1)
require(reshape)
sfInit(parallel=TRUE, cpus=2)
set.seed(0815)
# moved densities to denSparceRichBoth.R



.updateTwoDenA <- function( 
	### update paramter A conditional on observations and b
	theta				##<< numeric vector: current state in parameter space with components a and b 
	,argsFUpdateBlock=list()	##<< additional information from .updateBlocksTwDEMC
	,twTwoDenEx=twTwoDenEx1
){
	# subtract the effect of parameter b
	obsOffset <- twTwoDenEx$obs$y1 - theta[2]*mean(twTwoDenEx$xRich)/10
	# fit a linear model without offset to obtain mean and sd of the slope, i.e. parameter theta1
	lm1 <- lm( obsOffset ~ twTwoDenEx$xSparce -1 )
	meanA <- coef(lm1)[1]
	sdA <- sqrt(vcov(lm1)[1,1])
	# do a random draw
	newA <- rnorm(1,mean=meanA,sd=sdA)
	##value<< list with components
	list(	##describe<<
		accepted=1			##<< boolean scalar: if step was accepted
		, xC=c(a=newA)	##<< numeric vector: components of position in parameter space that are being updated
	)	##end<<
	#})
}
attr(.updateTwoDenA,"ex") <- function(){
	data(twTwoDenEx1)
	.updateTwoDenA( twTwoDenEx1$thetaTrue )
}
.tmp.f <- function(){
	plot( obsOffset ~ twTwoDenEx$xSparce )
	abline(lm1)
	xGrid <- meanA + 2*sdA*seq(-1,1,length.out=31)
	plot( dnorm(xGrid, mean=meanA, sd=sdA) ~ xGrid )
}


.updateTwoDenB <- function( 
	### update paramter A conditional on observations and b
	theta				##<< numeric vector: current state in parameter space with components a and b 
	,argsFUpdateBlock=list()	##<< additional information from .updateBlocksTwDEMC
	,twTwoDenEx=twTwoDenEx1
){
	thresholdCovar <- if( 0 != length(argsFUpdateBlock$argsFLogDen) && 0 != length(argsFUpdateBlock$argsFLogDen$thresholdCovar)){
		argsFUpdateBlock$argsFLogDen$thresholdCovar
	} else 0
	# subtract the effect of parameter b
	obsOffset <- twTwoDenEx$obs$y2 - theta[1]*pmax(0,twTwoDenEx$xSparce[1])
	# fit a linear model without offset to obtain mean and sd of the slope, i.e. parameter theta1
	lm1 <- lm( obsOffset ~ I(twTwoDenEx$xRich-thresholdCovar) -1 )
	meanB <- coef(lm1)[1]
	sdB <- sqrt(vcov(lm1)[1,1])
	# do a random draw
	newB <- rnorm(1,mean=meanB,sd=sdB)
	##value<< list with components
	list(	##describe<<
		accepted=1			##<< boolean scalar: if step was accepted
		, xC=c(b=newB)	##<< numeric vector: components of position in parameter space that are being updated
	)	##end<<
	#})
}
attr(.updateTwoDenB,"ex") <- function(){
	data(twTwoDenEx1)
	#trace(.updateTwoDenB, recover )
	.updateTwoDenB( twTwoDenEx1$thetaTrue, argsFUpdateBlock=list(argsFLogDen=list(thresholdCovar=0.3)) )
}
.tmp.f <- function(){
	plot( obsOffset ~ twTwoDenEx$xRich )
	abline(lm1)
	xGrid <- meanB + 2*sdB*seq(-1,1,length.out=31)
	plot( dnorm(xGrid, mean=meanB, sd=sdB) ~ xGrid )
}


thresholdCovar = 0.3	# the true value
thresholdCovar = 0		# the effective model that glosses over this threshold
thetaMean <- twTwoDenEx1$thetaTrue

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
		dInfos=list(
                dSparce=list(fLogDen=denSparce, argsFLogDen=list(theta0=c(a=1,b=2),thresholdCovar=thresholdCovar)))
		,nPop=.nPop
		,controlTwDEMC=list(thin=4)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	))
plot( as.mcmc.list(resa), smooth=FALSE)

# with biased b
resa1 <- resa <- concatPops( resBlock <- twDEMCBlock( ZinitPopsA, nGen=.nGen, 
		dInfos=list(dSparce=list(fLogDen=denSparce, argsFLogDen=list(theta0=c(a=1,b=1.6),thresholdCovar=thresholdCovar)))
		,nPop=.nPop
		,controlTwDEMC=list(thin=4)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	))
plot( as.mcmc.list(resa), smooth=FALSE)

#---------  fit only the short term observations to both parameters
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
		dInfos=list(dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar)))
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
		dInfos=list(dBoth=list(fLogDen=denBoth, argsFLogDen=list(thresholdCovar=thresholdCovar)))
		,nPop=.nPop
		,controlTwDEMC=list(thin=4)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	))
plot( as.mcmc.list(resa), smooth=FALSE)
#mtrace(subset.twDEMC)
resRS <-res2 <- res <- thin(resa, start=100)
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
res <- resRS
lDen <- stackChains(res$resLogDen)
#plot(density(lDen[lDen>-1000]))
lDenM <- melt(lDen)
boxplot( value ~ resComp, lDenM[lDenM$resComp != "parms",] )

madDen <- apply(lDen,2,function(lDenI){ median(abs(lDenI - median(lDenI))) })
madDen <- apply(lDen,2, median)     # actually wrong: the spread (difference to best) is important
weights <- wMedian <- ifelse( madDen, min(madDen[madDen != 0])/madDen, 0)
#weights <- c(1,1,0.1)

resa <- resa2b <-  concatPops( resBlock <- twDEMCBlock( ZinitPops, nGen=.nGen, 
		dInfos=list(dBoth=list(fLogDen=denBoth, argsFLogDen=list(weights=weights, thresholdCovar=thresholdCovar)))
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
lDen2 <- stackChains(res$resLogDen)
lDenM2 <- melt(lDen2)
boxplot( value ~ resComp, lDenM2[lDenM2$resComp != "parms",] )

lDenDiff <- apply(lDen2,2,function(lDenI){ abs(lDenI - median(lDenI)) })
apply(lDenDiff,2, median)
q10lDenDiff <- apply(lDenDiff,2,quantile,0.1)
lDenDiffM <- melt(lDenDiff)
lDenDiffM2 <- lDenDiffM[
	((lDenDiffM$resComp=="y1") & (lDenDiffM$value >= q10lDenDiff[1])) |
	((lDenDiffM$resComp=="y2") & (lDenDiffM$value >= q10lDenDiff[2])) 
		,]
boxplot( value ~ resComp, lDenDiffM )
boxplot( value ~ resComp, lDenDiffM2 )
if( require(ggplot2) ){
	p1 <- ggplot(lDenDiffM2, aes(y=value, x=resComp)) + geom_boxplot() +
		theme(axis.title.x=element_blank()) + coord_flip()
}

#----------- fit b to shortterm and a to longterm observations using direct sampling
# takes quite long
#untrace(.updateTwoDenB)
#trace(.updateTwoDenB,recover)
resa <- resa3a <- concatPops( resBlock <- twDEMCBlock( 
		res2$parms
		, nGen=.nGen 
		,dInfos=list(
			dBoth=list(fLogDen=denBoth, argsFLogDen=list(thresholdCovar=thresholdCovar))
		)
		,blocks = list(
			a=list(compPos="a", fUpdateBlock=.updateTwoDenA, requiresUpdatedDen=FALSE)
			,b=list(compPos="b", fUpdateBlock=.updateTwoDenB, requiresUpdatedDen=FALSE)
		)
		,nPop=.nPop
		,controlTwDEMC=list(thin=4)		
		,debugSequential=TRUE
		,doRecordProposals=TRUE
	))
plot( as.mcmc.list(resa), smooth=FALSE)
#mtrace(subset.twDEMC)
res3a <- res <- thin(resa, start=100)
plot( as.mcmc.list(res), smooth=FALSE)
ss <- stackChains(res)
twTwoDenEx1$thetaTrue
(thetaBest <- thetaBest3 <- ss[ which.max(ss[,1]), -1])
(.qq <- apply(ss[,-(1)],2,quantile, probs=c(0.025,0.5,0.975) ))
# a a bit biased upwards (effect of b)

pred <- pred3a <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich) )
plot( pred$y2 ~ twTwoDenEx1$obs$y2 ); abline(0,1) # here the mismatch becomes clear
plot( pred$y1 ~ twTwoDenEx1$obs$y1 ); abline(0,1) # not super but relation existing 

#----------- fit b to shortterm and a to longterm observations using Metropolis
resa <- resa3 <- concatPops( resBlock <- twDEMCBlock( 
		res3a$parms
		, nGen=.nGen*2 
		,dInfos=list(
			dSparce=list(fLogDen=denSparce, argsFLogDen=list(thresholdCovar=thresholdCovar))
			,dRich=list(fLogDen=denRich, argsFLogDen=list(thresholdCovar=thresholdCovar))
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

plot(density(stackChains(res3a)[,"a"]), xlim=c(0.8,1.2))
lines( density(stackChains(res3)[,"a"]), col="red")

pred <- pred3 <- with( twTwoDenEx1, fModel(thetaBest, xSparce=xSparce, xRich=xRich,thresholdCovar=thresholdCovar) )
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

if( require(rgl) )
    plot3d( ss[bo,"a"], ss[bo,"b"], lDen[bo], col=denCol )


#----------- ggplot of violin plots of the distributions 
require(ggplot2)
require(grid)

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
	,cbind( data.frame( scenario="DG", {tmp <- stackChains(res3a)[,-1]; tmp[ seq(1,nrow(tmp),length.out=.nSample),] }))
	,cbind( data.frame( scenario="DM", {tmp <- stackChains(res3)[,-c(1:2)]; tmp[ seq(1,nrow(tmp),length.out=.nSample),] }))
)
dfDenM <- melt(dfDen)
str(dfDenM)

StatDensityNonZero <- proto( ggplot2:::StatDensity ,{
    calculate <- function(.,data, scales, dMin=0.001, ...){
		res <- ggplot2:::StatDensity$calculate(data,scales,...)
		# append a zero density at the edges
		if( res$density[1] > 0) res <- rbind( data.frame(x=res$x[1]-diff(res$x[1:2]),density=0,scaled=0,count=0),res)
		if( res$density[nrow(res)] > 0) res <- rbind( res, data.frame(x=res$x[nrow(res)]+diff(res$x[nrow(res)-c(1,0)]),density=0,scaled=0,count=0))
		bo <- (res$density < dMin) & (c(res$density[-1],0) < dMin) & (c(0,res$density[-length(res$density)]) < dMin)
		res$density[ bo ] <- NA
		res$count[ bo ] <- NA
		res$scaled[ bo ] <- NA
		res
	}
    required_aes=c("x")
})
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
optsm <- theme(axis.title.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank() ) 
pa <- ggplot(dfDen, aes(x = a, colour=scenario, linetype=scenario)) + pgm + pgm2 + optsm + scale_y_continuous('a') +
   geom_vline( aes(xintercept= a, col="true", linetype="true"), data=as.data.frame(t(twTwoDenEx1$thetaTrue)))    
pb <- ggplot(dfDen, aes(x = b, colour=scenario, linetype=scenario)) + pgm + pgm2 + optsm + scale_y_continuous('b') + 
        geom_vline( aes(xintercept= b, col="true", linetype="true"), data=as.data.frame(t(twTwoDenEx1$thetaTrue)))    
#pb
windows(width=7, height=3)
grid.newpage()
pushViewport( viewport(layout=grid.layout(2,1)))
print(pa , vp = viewport(layout.pos.row=1,layout.pos.col=1))	
print(pb + theme(legend.position = "none") , vp = viewport(layout.pos.row=2,layout.pos.col=1))	

# XX TODO print true values

#----------- ggplot predictive posterior
	scenarios <- c("R","RS","RSw","DG","DM")
	nScen <- length(scenarios)
	# infer quantiles of predictions
	.nSample=128
	y1M <- array( NA_real_, dim=c(.nSample, length(twTwoDenEx1$obs$y1),nScen), dimnames=list(sample=NULL,iObs=NULL,scenario=scenarios)  )
	y2M <- array( NA_real_, dim=c(.nSample, length(twTwoDenEx1$obs$y2),nScen), dimnames=list(sample=NULL,iObs=NULL,scenario=scenarios) )
	resScen <- list( res1, res2, res2b, res3a, res3 ); names(resScen) <- scenarios
	#scen <- "R"
	for( scen in scenarios ){
		ss <- stackChains(resScen[[scen]])[,-(1:getNBlocks(resScen[[scen]]))]
		ssThin <- ss[round(seq(1,nrow(ss),length.out=.nSample)),]
		#i <- .nSample
		for( i in 1:.nSample){
			pred <-  twTwoDenEx1$fModel(ssThin[i,], xSparce=twTwoDenEx1$xSparce, xRich=twTwoDenEx1$xRich, thresholdCovar=thresholdCovar) 
			y1M[i,,scen] <- pred$y1
			y2M[i,,scen] <- pred$y2
		}
	}
	predQuantilesY1 <- lapply( scenarios, function(scen){
		t(apply( y1M[,,scen],2, quantile, probs=c(0.025,0.5,0.975) ))
	}); names(predQuantilesY1) <- scenarios
	predQuantilesY2 <- lapply( scenarios, function(scen){
			t(apply( y2M[,,scen],2, quantile, probs=c(0.025,0.5,0.975) ))
		}); names(predQuantilesY2) <- scenarios
	str(predQuantilesY1)
	rm( y1M, y2M)	# free space, we only need the summary
	
	#str(pred1)
	dfPred <- rbind(
		 cbind( data.frame( scenario="R", variable="y1", value = pred1$y1 ), predQuantilesY1[["R"]] )
		,cbind( data.frame( scenario="R", variable="y2", value = pred1$y2 ), predQuantilesY2[["R"]] )
		,cbind( data.frame( scenario="RS", variable="y1", value = pred2$y1 ), predQuantilesY1[["RS"]] )
		,cbind( data.frame( scenario="RS", variable="y2", value = pred2$y2 ), predQuantilesY2[["RS"]] )
		,cbind( data.frame( scenario="RSw", variable="y1", value = pred2b$y1 ), predQuantilesY1[["RSw"]] )
		,cbind( data.frame( scenario="RSw", variable="y2", value = pred2b$y2 ), predQuantilesY2[["RSw"]] )
		,cbind( data.frame( scenario="DG", variable="y1", value = pred3a$y1 ), predQuantilesY1[["DG"]] )
		,cbind( data.frame( scenario="DG", variable="y2", value = pred3a$y2 ), predQuantilesY2[["DG"]] )
		,cbind( data.frame( scenario="DM", variable="y1", value = pred3$y1 ), predQuantilesY1[["DM"]] )
		,cbind( data.frame( scenario="DM", variable="y2", value = pred3$y2 ), predQuantilesY2[["DM"]] )
	)
	dfPred$observations <- NA_real_
	dfPred$observations[dfPred$variable=="y1"] <- twTwoDenEx1$obs$y1
	dfPred$observations[dfPred$variable=="y2"] <- twTwoDenEx1$obs$y2
	colnames(dfPred)[ match(c("2.5%","50%","97.5%"),colnames(dfPred)) ] <- c("lower","median","upper")
	str(dfPred)

.tmp.f <- function(){
	dfPred <- rbind(
		data.frame( scenario="R", variable="y1", value = pred1$y1 )
		,data.frame( scenario="R", variable="y2", value = pred1$y2 )
		,data.frame( scenario="RS", variable="y1", value = pred2$y1 )
		,data.frame( scenario="RS", variable="y2", value = pred2$y2 )
		,data.frame( scenario="RSw", variable="y1", value = pred2b$y1 )
		,data.frame( scenario="RSw", variable="y2", value = pred2b$y2 )
		,data.frame( scenario="DG", variable="y1", value = pred3a$y1 )
		,data.frame( scenario="DG", variable="y2", value = pred3a$y2 )
		,data.frame( scenario="DM", variable="y1", value = pred3$y1 )
		,data.frame( scenario="DM", variable="y2", value = pred3$y2 )
	)
	dfPred$observations <- NA_real_
	dfPred$observations[dfPred$variable=="y1"] <- twTwoDenEx1$obs$y1
	dfPred$observations[dfPred$variable=="y2"] <- twTwoDenEx1$obs$y2
	str(dfPred)

	p1 <- ggplot(dfPred, aes(x=value, y=observations, colour=scenario) ) +
		#geom_errorbarh(aes(xmax = upper, xmin = lower)) +
		geom_point() +
		facet_wrap( ~ variable, scales="free") +
		opts(axis.title.x=element_blank() ) +
		geom_abline(colour="black") +
		c()
	p1
}

windows(width=7, height=3)
p2 <- ggplot(dfPred, aes(x=median, y=observations, colour=scenario) ) +
	#geom_errorbarh(aes(xmax = upper, xmin = lower)) +
	geom_point() +
	facet_wrap( ~ variable, scales="free") +
	opts(axis.title.x=element_blank() ) +
	geom_abline(colour="black") +
	c()
p2




