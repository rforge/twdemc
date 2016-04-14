.fixtureLinReg1 <- function(){
    set.seed(0815)
    fxLinReg <- within( list(),{
            fModel <- linearModelFunctionExample
            xval <- runif(30,min=5,max=10)
            thetaTrue = c(a=10,b=5)	# the parameter vector
            sdObs = 3*xval^0.8
            obs = fModel(thetaTrue, xval) + rnorm(length(xval), sd=sdObs)		### vector of data to compare with
            invCovar = diag(1/sdObs^2,nrow = length(sdObs))		### the inverse of the Covariance of obs (its uncertainty)
            sdTheta= thetaTrue*0.05	# 5% relativeerror
            invCovarTheta = diag(1/(sdTheta)^2,nrow = length(sdTheta))
            theta0 = thetaTrue + rnorm(length(thetaTrue),sd=sdTheta)
        })
    fxLinReg$blocks <- list(met1=blockSpec(,,new("MetropolisBlockUpdater",
                fLogDensity=function(theta, intermediates, logDensityComponents, fxLinReg){
                    pred <- linearModelFunctionExample(theta, fxLinReg$xval)
                    obs <- fxLinReg$obs
                    -1/2*sum( ((pred-obs)/fxLinReg$sdObs)^2 )
                },
                argsFLogDensity=list(fxLinReg=fxLinReg),
                logDensityComponentNames = c("obs1")))
    )
    fxLinReg
}
attr(.fixtureLinReg1,"ex") <- function(){
    fxLinReg <- .fixtureLinReg1()
    with( fxLinReg,{
            plot( xval, obs )
            abline( thetaTrue )
            abline( theta0, col="gray" )
        })
    
}
