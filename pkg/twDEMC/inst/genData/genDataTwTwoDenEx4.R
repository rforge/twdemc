.mod2Pred <- function(
	### example model giving two predictions that can be compared to different observations
	theta	##<< model parameters a and b
	, xSparce	##<< numeric vector of sparce input/output relationship 
	, xRich		##<< numeric vector of rich input/output relationship
	, thresholdCovar=0	##<< model structural deficiency
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
		,y2 = as.numeric(theta[1]*xSparce[1] + theta[2]*pmax(0,xRich-thresholdCovar) ) 		 
	) ##end<<
}


set.seed(0815)
twTwoDenEx1 <- within( list(), {
		fModel <- .mod2Pred
		thetaA <- 1
		xSparce <- c(1,runif(9,min=0.5,max=1.5))  # only 10 observations
		xRich <- runif(1000,min=.7,max=1)  	 # 1000 observations	
		# parameter b is an effective parameter 
		thetaTrue <- c( a=thetaA, b=2 )
		obsTrue <- .mod2Pred( thetaTrue, xSparce=xSparce,xRich=xRich, thresholdCovar=0.3)
		obsBiased <- .mod2Pred( thetaTrue, xSparce=xSparce,xRich=xRich)		
		sdObsTrue <- sdObs <- list( 
			 y1=mean(obsTrue$y1)*0.06
			,y2=mean(obsTrue$y2)*0.02 )
		obs <- list(
			y1 = obsTrue$y1 + rnorm(length(obsTrue$y1),sd=sdObsTrue$y1) 
			,y2 = obsTrue$y2 + rnorm(length(obsTrue$y2),sd=sdObsTrue$y2) )
	})

.tmp.f.plot <- function(){
	attach(twTwoDenEx1)
	with(twTwoDenEx1,{
		c( xSparce[1], mean(xRich) )
		plot(obsTrue$y1 ~ xSparce)
		abline(lm0 <- lm(obsTrue$y1 ~ xSparce))
		coef(lm0)

		plot(obsTrue$y2 ~ xRich)
		abline(lm1 <- lm(obsTrue$y2 ~ xRich))
		coef(lm1)
		confint(lm1)
		
		plot(obs$y1 ~ xSparce)
		abline(lm0o <- lm(obs$y1 ~ xSparce))
		coef(lm0o)
		confint(lm0o)
		
		plot( obsBiased$y2 ~ xRich, col="green")
		(lm1b <- lm(obsBiased$y2 ~ xRich))
		confint(lm1b)
	
		plot(obs$y2 ~ xRich)
		abline(lm1o <- lm(obs$y2 ~ xRich))
		coef(lm1o)
		confint(lm1o)
		
	}) # with
	detach(twTwoDenEx1)
}

save(twTwoDenEx1,file="data/twTwoDenEx1.RData")

