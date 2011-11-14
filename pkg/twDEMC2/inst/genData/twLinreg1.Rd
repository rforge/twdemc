     \name{twLinreg1}
     \docType{data}
     \title{Artificial data and parameters for DEMC-Example}
     \alias{twLinreg1}
     \description{
Artificial data and parameters for DEMC-Example.     
     }
     \usage{twLinreg1}
     \format{A list components \describe{
		\item{fModel}{\code{\link{dummyTwDEMCModel}} }
		\item{xval}{random numbers \code{runif(30,min=5,max=10)} }
		\item{obs}{artifical data: \code{fModel(thetaTrue, xval) + rnorm(length(xval), sd=sdObs)} }
		\item{sdObs}{standard deviation for each observations}
		\item{invCovar}{inverse of the covariance matrix: \code{diag(1/sdObs^2,nrow = length(sdObs))} }
		\item{thetaTrue}{Model parameters that were used to generate the artifical data}
		\item{sdTheta}{Standard devation of prior for model parameters, here 5 percent relative error: \code{thetaTrue*0.05} }
		\item{invCovarTheta}{inverse of the covariance matrix: \code{diag(1/(sdTheta)^2,nrow = length(sdTheta))}}
		\item{theta0}{Starting values, something off the true parameters: \code{thetaTrue + rnorm(length(thetaTrue),sd=sdTheta)} }
		}
     }
     \examples{
data(twLinreg1)     
attach(twLinreg1)
plot( xval, obs ) 
abline( thetaTrue )     
abline( theta0, col="gray" )     
     }
     \source{Thomas Wutzler}
     \keyword{datasets}
