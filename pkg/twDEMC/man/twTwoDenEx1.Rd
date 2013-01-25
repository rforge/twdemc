     \name{twTwoDenEx1}
     \docType{data}
     \title{Example model for artificial inverse experiment with two densities }
     \alias{twTwoDenEx1}
     \description{
        TODO     

     See test unit \code{twDEMC} testcase {ofMultiIntermediate}.    
     }
     \usage{twTwoDenEx1}
\format{ A list with components \describe{
    \item{ fModel }{ data generating model outputs \code{\link{modTwTwoDenEx1}}  }
    \item{ thetaTrue }{ true parameter vector} 
    \item{ sdObsTrue }{ Covariance matrix used to generate Multinormal samples } 
    \item{ obsTrue }{ observations generated with \code{thresholdCovar=0.3} } 
    \item{ obsBiased }{ biased observations generated with \code{thresholdCovar=0} } 
    \item{ sdObs }{ standard deviation for each observation } 
    \item{ xSparce }{ covariate to generate observations } 
    \item{ xRich }{ covariate to generate observations } 
}}
     \examples{
data(twTwoDenEx1)     
names(twTwoDenEx1)
     }
     \source{Thomas Wutzler}
     \keyword{datasets}
