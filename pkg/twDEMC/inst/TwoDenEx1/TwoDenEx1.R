TwoDen1 <- function(
	### Example of log-Unnormalized Density function with changing correlation and scales in variances 
	theta	##<< named numeric vector with components "a" and "b"
	,muA=0.8	##<< parameter for the log-normal distribution of component a
	,sigmaA=1	
){
	a <- theta[1]
	lda <- ((a-muA)/sigmaA)^2
	sdb <- exp(-2.8*a)*80
	mub <- (0.8-a)^2*20
	ldb <- ((b-mub)/sdb)^2
	#-1/2 * ( lda + max(ldb,0.9*lda))  # may not be a density 	
	c( den = as.numeric(-1/2 * ( lda + ldb )) )			   # also not a density because sigmab depends on a and does not factor out
}
attr(den2dCor,"ex") <- function(){
	#gridlogx <- seq(log(0.1),log(+4),length.out=91)
	#gridx <- exp(gridlogx)
	gridx <- a <- seq(-0.5,2,length.out=91)
	#plot( lda ~ a)
	gridy <- seq(-20,+40,length.out=91)
	gridX <- expand.grid(gridx, gridy)
	den2dCor(c(0.8,0.8))
	luden <- apply( gridX, 1, den2dCor ) 
	mLuden <- matrix(luden,nrow=length(gridx))
	#plot( rowSums(mLuden) ~ gridx )
	imax <- which( matrix(luden,nrow=length(gridx))==max(luden), arr.ind=TRUE)
	#c( gridx[ imax[1] ], gridy[ imax[2] ] )

	image( gridx, gridy,  mLuden, col = rev(heat.colors(100)), xlab="a", ylab="b" )
	xyMax <- c(x=gridx[ imax[1] ], y=gridy[ imax[2] ])
	
	image( gridx, gridy,  matrix(exp(luden),nrow=length(gridx)), col = rev(heat.colors(100)), xlab="a", ylab="b" )
	points( gridx[ imax[1] ], gridy[ imax[2] ]  )
	
	q20 <- quantile(luden,0.2) 
	plot(density(luden[luden>q20]))
	head(sort(luden,dec=TRUE))

	### todo: normalizing constant: 
	
	##------------------ do an MCMC run
	(.expTheta <- c(a=0,b=0) )
	(.expCovTheta <- diag(c(a=2,b=2)) )		
	.nPops=2
	Zinit <- initZtwDEMCNormal( .expTheta, .expCovTheta, nChains=4*.nPops, nPops=.nPops)
	#mtrace(twDEMCBlockInt)
	
	den2dCorTwDEMC <- twDEMCBatch(Zinit, nGen=500, fLogDen=den2dCor, nPops=.nPops )
	den2dCorTwDEMC <- twDEMCBatch(den2dCorTwDEMC, nGen=1000)
	
	plot( thinN(as.mcmc.list(den2dCorTwDEMC)))
	matplot( den2dCorTwDEMC$pAccept, type="l" )
	pps <- pps0 <- stackChains(thin(den2dCorTwDEMC,start=300))
	ss <- ss0 <- pps[,-1]
	#plot( ss[,1], ss[,2] )
	plot( ss[,1], ss[,2], ylim=c(-40,80) )
	plot( density(ss[,1]) )
	plot( ecdf( ss[,1] ) )
	plot( ecdf( ss[,2] ) )
}
