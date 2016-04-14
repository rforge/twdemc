den2dCor <- function(
	### Example of log-Unnormalized Density function with changing correlation and scales in variances 
	theta	##<< named numeric vector with components "a" and "b"
    ,intermediates=list()
    ,logDensityComponents=numeric(0)
	,muA=-1.5,sigmaA=1.2	##<< parameter for the log-normal distribution of component a
	,cb=1/4	##<< multiplier for density in b
	#,cb=1e-900	##<< multiplier for density in b
	,cDen=0
){
    ##description<<
    ## a ~ N( theta[1], sigmaA )
    ## b ~ N( fMu(a), fSigma(a) )
    ##
    ## With fMu(a) quadratic
    ## With fSigma(a) exponentially decreasing with increasing a
    #
	a <- theta[1]
	b <- theta[2]
	ssqA <- ((a-muA)/sigmaA)^2
    #mub <- (0.8-a)^2*20
    mub <- (a+0.4)^6
    #sdb <- exp(-2.8*a)*80
    sdb <- exp(-1.5*(a))*200
    ssqB <- ((b-mub)/sdb)^2 
	#-1/2 * ( lda + max(ldb,0.9*lda))  # may not be a density
	c( den = as.numeric(-1/2*( ssqA + ssqB ))-log(sdb*sqrt(2*pi))-log(sigmaA*sqrt(2*pi)) )			   
}
attr(den2dCor,"ex") <- function(){
    gridx <- a <- seq(-0.5,2,length.out=91)
    gridx <- a <- seq(-1,2,length.out=91)
    gridx <- a <- seq(-0.5,2,,length.out=91)
    gridx <- a <- seq(-1,3,,length.out=91)
    gridx <- a <- seq(-2,2,,length.out=91)
    gridy <- seq(-20,+40,length.out=91)
    gridy <- seq(-180,+180,length.out=91)
    gridX <- expand.grid(gridx, gridy)
	den2dCor(c(0.8,0.8))
    muA <- -1.5
    sigmaA <- 1.2
	luden <- apply( gridX, 1, den2dCor, muA=muA, sigmaA=sigmaA ) 
	mLuden <- matrix(luden,nrow=length(gridx))
	#plot( rowSums(mLuden) ~ gridx )
	imax <- which( matrix(luden,nrow=length(gridx))==max(luden), arr.ind=TRUE)
	#c( gridx[ imax[1] ], gridy[ imax[2] ] )

	#image( gridx, gridy,  mLuden, col = rev(heat.colors(100)), xlab="a", ylab="b" )
	#xyMax <- c(x=gridx[ imax[1] ], y=gridy[ imax[2] ])
	
	image( gridx, gridy,  matrix(exp(luden),nrow=length(gridx)), col = rev(heat.colors(100)), xlab="a", ylab="b" )
	points( gridx[ imax[1] ], gridy[ imax[2] ]  )
	
    pA <- rowSums(exp(mLuden)) 
    plot( pA/max(pA)~ gridx)
    pX <- dnorm(gridx, muA, sigmaA)
    lines( pX/max(pX) ~ gridx)
    abline(v=muA)
    
	q20 <- quantile(luden,0.2) 
	plot(density(luden[luden>q20]))
	head(sort(luden,dec=TRUE))
}


denBanana <- function(
        theta
        ,intermediates=list()
        ,logDensityComponents=numeric(0)
        ,mleA=0.5       ##<< center of parameter a
        ,coefA= twCoefLogitnormMLE(mleA,0)[1,]	   ##<< alternative way of specifying logitnormal parameters for a,
            ##<< numeric vector (2): first mu, second sigma
        ,sdB=0.1
){
    muB <- ((theta[1]-mleA)*2)^2
    logdenA <- suppressWarnings(log(dlogitnorm(theta[1], coefA[1], coefA[2])))
    logdenB <- dnorm(theta[2], muB, sdB, log = TRUE)
    ans <- c(logDenBanana=as.vector(logdenA + logdenB))
    if( is.na(ans) ){     ans[] <- -Inf   }
    ans
}
attr(denBanana,"ex") <- function(){
    #denBanana( c(0.8,0) )
    mleA <- 0.5
    coefA <- twCoefLogitnormMLE(mleA,0)[1,]	   
    gridx <- a <- seq(0,1,length.out=91)[-c(1,91)]
    gridy <- seq(-0.2,1.1,length.out=51)
    gridX <- expand.grid(gridx, gridy)
    #
    luden <- apply( gridX, 1, denBanana, coefA=coefA)
    mLuden <- matrix(luden,nrow=length(gridx))
    imax <- which( matrix(luden,nrow=length(gridx))==max(luden), arr.ind=TRUE)[1,]
    #c( gridx[ imax[1] ], gridy[ imax[2] ] )
    
    image( gridx, gridy,  matrix(exp(luden),nrow=length(gridx)), col = rev(heat.colors(100)), xlab="a", ylab="b" )
    points( gridx[ imax[1] ], gridy[ imax[2] ]  )
    
    # the denstiy of the two variables is multiplicative
    # hence the marginal of parameter a, is the specified density
    pA <- rowSums(exp(mLuden)) 
    plot( pA/max(pA)~ gridx)
    pX <- dlogitnorm(gridx, coefA[1], coefA[2])
    lines( pX/max(pX) ~ gridx)
}
