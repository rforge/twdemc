\name{constrainNStack}
\alias{constrainNStack}
\title{constrainNStack}
\description{Subsetting a sequence of parameter vectors. Keeps only the n cases in pss1 that are closest to thetaPrior.}
\usage{
constrainNStack(pss1, thetaPrior, n = nrow(pss1)\%/\%4, vars = names(thetaPrior), invCovarTheta = if (0 < length(thetaPrior)) solve(cov(pss1[, names(thetaPrior[vars]), drop = FALSE])) else numeric(0), returnAlpha = FALSE)
}
\arguments{
  \item{pss1}{numeric matrix: the stack to constrain, see details}
  \item{thetaPrior}{the point in parameter spcae for which to select the closest values}
  \item{n}{the number of rows in output, defaults to 1/4 of nrow(pss1)}
  \item{vars}{names or indices to constrain thetaPrior and invCovarTheta}
  \item{invCovarTheta}{the inverse of the covaraince matrix for thetaPrior, defaults to inverse of sample covariance}
  \item{returnAlpha}{switch to return also significance level}
}
\details{pss1 is a matrix with columns corresponding to variables 
and rows corresponding to cases. 
It is typicalla the result of \code{\link{stackChains.twDEMC}}}
\value{pss1[closestValues,]
if returnAlpha=TRUE then a list is returned with \describe{ 
\item{res}{ value as above }
\item{alpha}{ the significance level of the corresponding confidence ellipsis }} }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}
\code{\link{constrainCfStack}}}

