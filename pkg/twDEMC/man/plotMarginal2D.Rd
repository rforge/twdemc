\name{plotMarginal2D}
\alias{plotMarginal2D}
\title{plotMarginal2D}
\description{Plot the 2D marginal Likelihood of a sample}
\usage{
plotMarginal2D(smp, xCol = 2, yCol = 3, intCol = 1, grains = 18, FUN = mean, argsFUN = list(na.rm = TRUE), col = rev(heat.colors(20)), minN = 7, ...)
}
\arguments{
  \item{smp}{numeric matrix: first column log-Likelihood, other columns free parameters, see \code{\link{stackChains.twDEMC}}}
  \item{xCol}{index (column number or column name) of column for x ordinate}
  \item{yCol}{index (column number or column name) of column for y ordinate}
  \item{intCol}{index (column number or column name) of column of LogLikelihood values}
  \item{grains}{vector of length 1 or 2 giving the number of groups for x and y classes respectively}
  \item{FUN}{the function applied over log-Likelihoods in the classes}
  \item{argsFUN}{additional arguments to FUN}
  \item{col}{vector of colors}
  \item{minN}{minimum number of items in classpixel to be plotted}
  \item{\dots}{additional arguments to twPlot2D}
}
\details{The entire sample is split into bins of about equal number of observations and all observations regarding x and y.
Within the bin the values are aggregated. By default the mean is calculated.

There are several plotting methods related to twDEMC run. \itemize{
\item{ TODO: link methods  } 
\item{ the Gelman criterion: this method  } 
\item{ the theorectical minimum logLik-Value for significant model difference : \code{\link{getRLogLikQuantile}}  } 
}}


\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{twDEMCInt}}}
\examples{
data(twdemcEx1)
sample <- stackChains(twdemcEx1)
sample0 <- sample[ sample[,1] >= quantile(sample[,1],0.05), ]
#simple plot
#mtrace(plotMarginal2D)
plotMarginal2D( sample0, "a", "b", minN=1 )	# actually not a marginal
}
