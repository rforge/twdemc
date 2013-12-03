\name{plot2DKDMarginals}
\alias{plot2DKDMarginals}
\title{plot2DKDMarginals}
\description{Evaluates the Kernel-density regression on a grid and displays a figure.}
\usage{plot2DKDMarginals(form, ds, dims = 40, bandwidth = NULL, nSample = 200, 
    bwtype = "adaptive_nn", tol = 0.1, ftol = 0.1, bws = NULL, 
    ...)}
\arguments{
  \item{form}{formula of regression e.g. rLogDen~dO+tvrO}
  \item{ds}{dataframe or matrix with column names corresponding to variables in the formula}
  \item{dims}{vector of spaces across predictor space where}
  \item{bandwidth}{see adaptive bandwidth of npreg.
integer vector: bandwidth for each dimension
about Records to include for estimating bandwidth.
integer vector: bandwidth for each dimension
If only one value is supplied, then it is repeated for each predictor
integer vector: bandwidth for each dimension}
  \item{nSample}{The larger the sample, the longer it takes
the sample size to estimate the bandwidth and the marginal}
  \item{bwtype}{parameters for \code{\link{npregbw}}}
  \item{tol}{
}
  \item{ftol}{
}
  \item{bws}{may directly supply a bandwidht object}
  \item{\dots}{further arguments to \code{\link{plot.twApply2DMesh}}}
}
\details{Uses package snowfall to split calculation to subprocesses.
Requires package np loaded in cluster (sfLibrary(np))}
\value{Value of \code{\link{twApply2DMesh}}. 
- may be used with \code{\link{plot.twApply2DMesh}} from package twMisc - the side effect of this method.
Additional entry \code{bandwidth}, the result of bandwidth calculation \code{\link{npregbw}}.}

\author{Thomas Wutzler}




\examples{
.tmp.f <- function(){
    #wrapped inside function because it takes long
    data(twdemcEx1)
    sample <- stackChains(twdemcEx1)
    sample0 <- sample[ sample[,1] >= quantile(sample[,1],0.05), ]
    # just plot the sample
    plot( sample0[,-1], col=rev(heat.colors(20))[round(twRescale(sample0[,1],c(1,20)))] )
    
    #simple plot
    #mtrace(plotMarginal2D)
    plotMarginal2D( sample0, "a", "b", minN=1 )	# actually not a marginal
    
    # kernel density regression estimate
    #mtrace(plot2DKDMarginals)
    rl <- sample0[,1]
    set.seed(0815)
    sample0b <- cbind(sample0, rLogDen2 = rl+rnorm(length(rl),sd=diff(range(rl))/4 ))
    plot( rLogDen2 ~ dInfo1, sample0b)
    plot( sample0b[,c("a","b")], col=rev(heat.colors(20))[round(twRescale(sample0b[,"rLogDen2"],c(1,20)))] )
    #refit disturbed sample
    tmp2 <- tmp3 <- plot2DKDMarginals( rLogDen2 ~ a+b, sample0b, dim=20  )
    tmp2$bandwidth$x
    # the bandwidth seems to be underestimated for the second dimension, repeat with higher bandwidth
    #tmp3 <- plot2DKDMarginals( rLogDen2 ~ a+b, sample0b, dim=20, bandwidth=30  )
    #refine plot
    plot(tmp3, xlab="Intercept a", ylab="Slope b", zlab="LogDensity", contour=TRUE )
}
}
