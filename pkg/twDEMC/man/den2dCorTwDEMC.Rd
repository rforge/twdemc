     \name{den2dCorTwDEMC}
     \docType{data}
     \title{Example DEMC result}
     \alias{den2dCorTwDEMC}
     \description{
Result of running DEMC-Example without subspaces for density with different scales of variance.     
     }
     \usage{den2dCorTwDEMC}
     \format{A list with components described by \code{\link{twDEMCInt}}  }
     \examples{
data(den2dCorTwDEMC)     
str(den2dCorTwDEMC)
rownames(den2dCorTwDEMC)
plot(as.mcmc.list(den2dCorTwDEMC), smooth=FALSE)
sample <- stackChains(den2dCorTwDEMC)
sample0 <- sample[ sample[,1] >= quantile(sample[,1],0.05), ]
plot(sample0[,-1], col=rev(heat.colors(100))[round(rescale(sample0[,1],to=c(1,100)))] )
     }
     \source{Thomas Wutzler}
     \keyword{datasets}
