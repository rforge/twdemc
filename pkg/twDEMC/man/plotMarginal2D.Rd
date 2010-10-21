\name{plotMarginal2D}
\alias{plotMarginal2D}
\title{plotMarginal2D}
\description{Plot the 2D marginal Likelihood of a sample}
\usage{
plotMarginal2D(smp, xCol = 2, yCol = 3, intCol = 1, ..., grains = 20, FUN = mean, col = rev(heat.colors(100)))
}
\arguments{
  \item{smp}{numeric matrix: first column log-Likelihood, other columns free parameters, see \code{\link{stackChains.twDEMC}}}
  \item{xCol}{index (column number or column name) of column for x ordinate}
  \item{yCol}{index (column number or column name) of column for y ordinate}
  \item{intCol}{index (column number or column name) of column of LogLikelihood values}
  \item{\dots}{additional arguments to FUN}
  \item{grains}{vector of length 1 or 2 giving the number of groups for x and y classes respectively}
  \item{FUN}{the function applied over log-Likelihoods in the classes}
  \item{col}{vector of colors}
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

