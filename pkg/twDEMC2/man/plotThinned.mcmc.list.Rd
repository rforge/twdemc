\name{plotThinned.mcmc.list}
\alias{plotThinned.mcmc.list}
\title{plotThinned mcmc list}
\description{Plot of thinned coda's \code{mcmc.list}.}
\usage{
\method{plotThinned}{mcmc.list}(x, maxVars = 3, vars = 1:min(maxVars, ncol(x[[1]])), maxN = 500, suppressPlotWarnings = TRUE, ...)
}
\arguments{
  \item{x}{
}
  \item{maxVars}{maximum number of variable to put in one plot window}
  \item{vars}{index of variables to plot}
  \item{maxN}{maximum number of cases to plot, (thin before)}
  \item{suppressPlotWarnings}{if true, plot command is executed with suppressWarnings}
  \item{\dots}{
}
}
\details{default coda plot for mcmc.list squeezes all variables in one window}


\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{plotMarginal2D}}
\code{\link{twDEMCInt}}}

