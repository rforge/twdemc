# moved modTwTwoDenEx1 to R/denSparseRichBoth.R so that docu is generated

set.seed(0815)
twTwoDenEx1 <- within( list(), {
		fModel <- modTwTwoDenEx1
		thetaA <- 1
		xSparse <- c(1,runif(9,min=0.5,max=1.5))  # only 10 observations
		xRich <- runif(1000,min=.7,max=1)  	 # 1000 observations	
		# parameter b is an effective parameter 
		thetaTrue <- c( a=thetaA, b=2 )
		obsTrue <- modTwTwoDenEx1( thetaTrue, xSparse=xSparse,xRich=xRich
			, thresholdCovar=0.3)
		obsBiased <- modTwTwoDenEx1( thetaTrue, xSparse=xSparse,xRich=xRich)		
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
		c( xSparse[1], mean(xRich) )
		plot(obsTrue$y1 ~ xSparse)
		abline(lm0 <- lm(obsTrue$y1 ~ xSparse))
		coef(lm0)

		plot(obsTrue$y2 ~ xRich)
		abline(lm1 <- lm(obsTrue$y2 ~ xRich))
		coef(lm1)
		confint(lm1)
		
		plot(obs$y1 ~ xSparse)
		abline(lm0o <- lm(obs$y1 ~ xSparse))
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

