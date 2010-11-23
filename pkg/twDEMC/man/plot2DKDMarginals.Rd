\name{plot2DKDMarginals}
\alias{plot2DKDMarginals}
\title{plot2DKDMarginals}
\description{Evaluates the Kernel-density regression on a grid and displays a figure.}
\usage{
plot2DKDMarginals(form, ds, dim = 40, bandwidth = NULL, nSampleBandwidth = 200, bwtype = "adaptive_nn", tol = 0.1, ftol = 0.1)
}
\arguments{
  \item{form}{formula of regression e.g. rLogLik~dO+tvrO}
  \item{ds}{dataframe or matrix with column names corresponding to variables in the formula}
  \item{dim}{vector of spaces across predictor space where}
  \item{bandwidth}{see adaptive bandwidth of npreg.
about Records to include for estimating bandwidth.
If only one value is supplied, then it is repeated for each predictor

integer vector: bandwidth for each dimension}
  \item{nSampleBandwidth}{The larger the sample, the longer it takes

the sample size to estimate the bandwidth}
  \item{bwtype}{parameters for \code{\link{npregbw}}}
  \item{tol}{
}
  \item{ftol}{
}
}
\details{Uses package snowfall to split calculation to subprocesses.
Requires package np loaded in cluster (sfLibrary(np))}
\value{Matrix of predictions for the grid.}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(twdemcEx1)
sample <- stackChains(twdemcEx1)
sample0 <- sample[ sample[,1] >= quantile(sample[,1],0.05), ]
#simple plot
plotMarginal2D( sample0, "a", "b", minN=1 )	# actually not a marginal

# kernel density regression estimate
#mtrace(plot2DKDMarginals)
tmp2 <- plot2DKDMarginals( rLogLik ~ a+b, sample0, dim=20  )

#refine plot
xy <- apply(do.call(cbind,dimnames(tmp2)),2,as.numeric  )
twPlot2DFun.contour( xy[,1], xy[,2], tmp2, xdiv=nrow(tmp2)
	, xlab="Intercept a"
	, ylab="Slope b"
	, key.title=title(sub="Log-Like-\nlihood\n\n")
	, color.palette=function(n){rev(heat.colors(n))
})
}
