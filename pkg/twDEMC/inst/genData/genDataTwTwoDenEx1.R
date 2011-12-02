.mod2Pred <- function(
	### example model giving two predictions that can be compared to different observations
	theta	##<< model parameters a and b
	, xInAnnual	  ##<< numeric vector of annual GPP
	, xRespCov	  ##<< covariate for respiration
	, meanRespCov=quantile(xRespCov,0.2)	##<< annual aggregate of the covariate
){
	assAnn <- theta[1]*xInAnnual
	respAnn <- assAnn * theta[2]*meanRespCov
	allocAnn <- assAnn - respAnn
	
	assDailyMean <- assAnn/365
	respDaily <- assDailyMean * theta[2]*xRespCov
   ##value<< list with model predictions
   list( ##describe<<
	      allocAnn = allocAnn
   		, respDaily = respDaily
	) ##end<<
}


set.seed(0815)
twTwoDenEx1 <- within( list(), {
		fModel <- .mod2Pred
		thetaA <- 0.8
		xInAnnual <- c(1000,runif(9,min=800,max=1200))  # only 10 observations 10 years, but high precision
		xRespCov <- rlnorm(1000, meanlog=0-(0.5^2)/2, sdlog=0.5)  	 # 1000 observations
		# parameter b is an effective parameter 
		thetaTrue <- cbind( a=thetaA, b=rlnorm(length(xRespCov), meanlog=log(0.5), sdlog=0.4) )
		thetaEff <- apply(thetaTrue,2,mean)
		obsTrue <- list(
			allocAnn = .mod2Pred( thetaEff, xInAnnual=xInAnnual,xRespCov=xRespCov)$allocAnn
			,respDaily = sapply( seq_along(xRespCov), function(i){ .mod2Pred(thetaTrue[i,],xInAnnual=xInAnnual[1],xRespCov=xRespCov[i])$respDaily})
		)
		obsEff <- .mod2Pred(thetaEff,xInAnnual=xInAnnual[1],xRespCov=xRespCov)
		sdObsTrue <- sdObs <- list( 
			 allocAnn=mean(obsTrue$allocAnn)*0.01
			,respDaily=mean(obsTrue$respDaily)*0.01 )
		obs <- list(
			allocAnn = obsTrue$allocAnn + rnorm(length(obsTrue$allocAnn),sd=sdObsTrue$allocAnn) 
			,respDaily = obsTrue$respDaily + rnorm(length(obsTrue$respDaily),sd=sdObsTrue$respDaily) )
		sdObs$respDaily <- sd(resid(lm(obs$respDaily ~ xRespCov)))	# error includes the effects of varying true b 
	})

.tmp.f.plot <- function(){
	attach(twTwoDenEx1)
	with(twTwoDenEx1,{
		c( xInAnnual[1], mean(xRespCov) )
		plot( density( thetaTrue[,"b"]))
		abline(v=thetaEff["b"])
		plot(obsTrue$allocAnn ~ xInAnnual)
		abline(lm0 <- lm(obsTrue$allocAnn ~ xInAnnual))
		coef(lm0)
		
		plot(obsEff$respDaily ~ xRespCov)
		lmEff <- lm(obsEff$respDaily ~ xRespCov) 
		
		plot(obsTrue$respDaily ~ xRespCov)
		abline(lm1 <- lm(obsTrue$respDaily ~ xRespCov))
		abline(lmEff, col="maroon")
		coef(lm1)
		confint(lm1)
		(lm2 <- lm(I(obsTrue$respDaily-10) ~ xRespCov -1 ))
		thetaEff
		
		
		
		plot(obs$allocAnn ~ xInAnnual)
		abline(lm0o <- lm(obs$allocAnn ~ xInAnnual))
		coef(lm0o)
		confint(lm0o)
		
		plot(obs$respDaily ~ xRespCov)
		abline(lm1o <- lm(obs$respDaily ~ xRespCov))
		coef(lm1o)
		confint(lm1o)
		(lm2 <- lm(I(obs$respDaily-10) ~ xRespCov -1 ))
		thetaEff
	}) # with
	detach(twTwoDenEx1)
}

save(twTwoDenEx1,file="data/twTwoDenEx1.RData")

