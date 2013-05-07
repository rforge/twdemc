modDenBoth <- function(
	theta
	,twTwoDenEx1=twTwoDenEx1 ##<< list with components rich and Sparse 
){
	theta[1] * x[,1] +
}

denBoth <- function(
	theta
	,twTwoDenEx1=twTwoDenEx1 ##<< list with components rich and Sparse 
){
	obsRich <- with( twTwoDenEx1$rich, logDenGaussian(theta, fModel=fModel, obs=obs, invCovar=varObs, xval=xval ))
	obsSparse <- with( twTwoDenEx1$Sparse, logDenGaussian(theta, fModel=fModel, obs=obs, invCovar=varObs, xval=xval ))
	obsRich["obs"] + obsSparse["obs"]
	### the misfit: scale *( t(tmp.diffObs) %*% invCovar %*% tmp.diffObs + t(tmp.diffParms) %*% invCovarTheta %*% tmp.diffParms )
}

data( twTwoDenEx1 )
denBoth( twTwoDenEx1$rich$theta0, twTwoDenEx1 )

attach( twTwoDenEx1 )




.nPops=2
.nChainPop=4
ZinitPops <- initZtwDEMCNormal( rich$theta0, diag(rich$varTheta), nChainPop=.nChainPop, nPops=.nPops)
#dim(ZinitPops)
#head(ZinitPops[,,1])
pops <- pops0 <- list(
	pop1 <- list(
		Zinit = ZinitPops[1:3,,1:4,drop=FALSE]	# the first population with less initial conditions
		,nGen=8
	),
	pop2 <- list(
		Zinit = ZinitPops[,,5:8,drop=FALSE]	# the first population with less initial conditions
		,nGen=100
		,T0=10
	)
)
#tmp <- .checkPop(pops[[1]])
# both fLogDen compare the same model against the same observations but use different priors
dInfoDefault <- list(
	fLogDen=logDenGaussian
)
dInfos0 <- dInfox <- list(
	rich = within(dInfoDefault,{
			argsFLogDen = list(
				fModel=rich$fModel,		### the model function, which predicts the output based on theta 
				obs=rich$obs,			### vector of data to compare with
				invCovar=rich$varObs,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
				xval=rich$xval
			)
		})
	,Sparse = within(dInfoDefault,{
			argsFLogDen = list(
				fModel=Sparse$fModel,		### the model function, which predicts the output based on theta 
				obs=Sparse$obs,			### vector of data to compare with
				invCovar=Sparse$varObs,		### do not constrain by data, the inverse of the Covariance of obs (its uncertainty)
				xval=Sparse$xval
			)
		})
	,both = list( fLogDen=denBoth, argsFLogDen = list(twTwoDenEx1=twTwoDenEx1) )
)
#mtrace(logDenGaussian)
do.call( logDenGaussian, c(list(theta=rich$theta0), dInfos0$rich$argsFLogDen))
do.call( logDenGaussian, c(list(theta=Sparse$theta0), dInfos0$Sparse$argsFLogDen))
do.call( denBoth, c(list(theta=rich$theta0), dInfos0$both$argsFLogDen))


#fit only to rich data
blocks <- list(
	rich <- list(compPos=c("a","b"), dInfoPos="rich")
)
.nGen=100
.thin=4
#mtrace(.updateBlockTwDEMC)
#mtrace(.updateBlocksTwDEMC)
#mtrace(.updateIntervalTwDEMCPar)
#mtrace(twDEMCBlockInt)
resRichOnly <- res <- twDEMCBlock( pops=pops, dInfos=dInfos0, blocks=blocks, nGen=.nGen, controlTwDEMC=list(thin=.thin) )
str(res$pops[[2]])















detach( twTwoDenEx1 )