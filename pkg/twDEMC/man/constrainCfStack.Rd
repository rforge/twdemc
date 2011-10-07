\name{constrainCfStack}
\alias{constrainCfStack}
\title{constrainCfStack}
\description{Subsetting a sequence of parameter vectors. Keeps only the the cases in pss1 that are inside a confidence ellipsis around thetaPrior.}
\usage{
constrainCfStack(pss1, thetaPrior, alpha = 0.95, vars = names(thetaPrior), invCovarTheta = if (0 < length(thetaPrior)) solve(cov(pss1[, vars, drop = FALSE])) else numeric(0))
}
\arguments{
  \item{pss1}{numeric matrix: the stack to constrain, see details}
  \item{thetaPrior}{the point in parameter spcae for which to select the closest values}
  \item{alpha}{the conficence range of the ellipsis}
  \item{vars}{names or indices to constrain thetaPrior and invCovarTheta}
  \item{invCovarTheta}{the inverse of the covaraince matrix for thetaPrior, defaults to inverse of sample covariance}
}

\value{pss1[ withinConfidenceInterval, ]	}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}
\code{\link{constrainNStack}}}

