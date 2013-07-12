     \name{den2dCorEx}
     \docType{data}
     \title{Example model for artificial inverse experiment with changing correlation and scales in variances }
     \alias{den2dCorEx}
     \description{
        TODO     

     See test unit \code{divideTwDEMCStep}    
     }
     \usage{den2dCorEx}
\format{ A list with components \describe{
    \item{ \code{\link{den2dCor}} }{ the model  }
    \item{ thetaPrior }{ true parameter vector} 
    \item{ covarTheta }{ Covariance matrix used to generate Multinormal samples } 
    \item{ mcBulk }{ MCMC result without subspace sampling } 
    \item{ subspaces0 }{ intial partitioning of the parameter space including the initial samples } 
    \item{ mcSubspaces0 }{ first step MCMC result with subspace sampling } 
}}
     \examples{
data(den2dCorEx)     
names(den2dCorEx)
names(den2dCorEx$mcSubspaces0)
plot(as.mcmc.list(concatPops(den2dCorEx$mcSubspaces0)), smooth=FALSE)
sample <- stackChains(den2dCorEx$mcSubspaces0)
sample0 <- sample[ sample[,1] >= quantile(sample[,1],0.05), ]
plot(sample0[,-1], col=rev(heat.colors(100))[round(twMisc::twRescale(sample0[,1],to=c(1,100)))] )
     }
     \source{Thomas Wutzler}
     \keyword{datasets}
