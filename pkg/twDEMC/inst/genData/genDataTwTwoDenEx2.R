.mod2Pred <- function(
	### example model giving two predictions that can be compared to different observations
	theta	##<< model parameters a and b
	, xSparce	##<< numeric vector of sparce input/output relationship 
	, xRich		##<< numeric vector of rich input/output relationship
){
   ##details<< model output y1 represents a longterm observations
   ## It is based on longterm average of xRich instead of detailed values
   ## \cr Model output y1 represents a short measurement campaing. 
   ## During this campaing xSparce does not vary but detailed measurements of xRich are utilized 

   ##value<< list with model predictions
   list( ##describe<<
	   	# with the following line, bias of theta[2] fully leaks into estimate of a
		#y1 = as.numeric(theta[1]*xSparce + theta[2]*mean(xRich))	##<< theta[1]*x1 + theta[2]*8
		y1 = as.numeric(theta[1]*xSparce + theta[2]*mean(xRich)/10)	##<< theta[1]*x1 + theta[2]*mean(xRich)/10
		,y2 = as.numeric(theta[1]*xSparce[1] + theta[2]*xRich)		##<< theta[1]*1 + theta[1]*1 + theta[2]*xRich 
	) ##end<<
}


set.seed(0815)
twTwoDenEx1 <- within( list(), {
		fModel <- .mod2Pred
		thetaA <- 100
		xSparce <- c(1,runif(9,min=0.5,max=1.5))  # only 10 observations
		xRich <- runif(1000,min=70,max=100)  	 # 1000 observations	
		# parameter b is an effective parameter 
		thetaTrue <- cbind( a=thetaA, b=rlnorm(length(xRich), meanlog=log(2), sdlog=0.4) )
		thetaEff <- apply(thetaTrue,2,mean)
		obsTrue <- list(
			y1 = .mod2Pred( thetaEff, xSparce=xSparce,xRich=xRich)$y1
			,y2 = sapply( seq_along(xRich), function(i){ .mod2Pred(thetaTrue[i,],xSparce=xSparce[1],xRich=xRich[i])$y2})
		)
		sdObsTrue <- sdObs <- list( 
			 y1=mean(obsTrue$y1)*0.01
			,y2=mean(obsTrue$y2)*0.01 )
		obs <- list(
			y1 = obsTrue$y1 + rnorm(length(obsTrue$y1),sd=sdObsTrue$y1) 
			,y2 = obsTrue$y2 + rnorm(length(obsTrue$y2),sd=sdObsTrue$y2) )
		sdObs$y2 <- sd(resid(lm(obs$y2 ~ xRich)))	# error includes the effects of varying true b 
	})

.tmp.f.plot <- function(){
	attach(twTwoDenEx1)
	with(twTwoDenEx1,{
		c( xSparce[1], mean(xRich) )
		plot( density( thetaTrue[,"b"]))
		abline(v=thetaEff["b"])
		plot(obsTrue$y1 ~ xSparce)
		abline(lm0 <- lm(obsTrue$y1 ~ xSparce))
		coef(lm0)
		plot(obsTrue$y2 ~ xRich)
		abline(lm1 <- lm(obsTrue$y2 ~ xRich))
		abline(thetaEff, col="maroon")
		coef(lm1)
		confint(lm1)
		(lm2 <- lm(I(obsTrue$y2-10) ~ xRich -1 ))
		thetaEff
		
		
		plot(obs$y1 ~ xSparce)
		abline(lm0o <- lm(obs$y1 ~ xSparce))
		coef(lm0o)
		confint(lm0o)
		
		plot(obs$y2 ~ xRich)
		abline(lm1o <- lm(obs$y2 ~ xRich))
		coef(lm1o)
		confint(lm1o)
		(lm2 <- lm(I(obs$y2-10) ~ xRich -1 ))
		thetaEff
	}) # with
	detach(twTwoDenEx1)
}

save(twTwoDenEx1,file="data/twTwoDenEx1.RData")

