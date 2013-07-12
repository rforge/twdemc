     \name{twdemcEx1}
     \docType{data}
     \title{Example DEMC result}
     \alias{twdemcEx1}
     \description{
Result of running DEMC-Example of \code{\link{twLinreg1}}.     
     }
     \usage{twdemcEx1}
     \format{A list with components described by \code{\link{twDEMCBlockInt}}  }
     \examples{
data(twdemcEx1)     
str(twdemcEx1)
rownames(twdemcEx1)
plot(as.mcmc(twdemcEx1))
sample <- stackChains(twdemcEx1)
sample0 <- sample[ sample[,1] >= quantile(sample[,1],0.05), ]
plot(sample0[,-1], col=rev(heat.colors(100))[round(twMisc::twRescale(sample0[,1],to=c(1,100)))] )
     }
     \source{Thomas Wutzler}
     \keyword{datasets}
