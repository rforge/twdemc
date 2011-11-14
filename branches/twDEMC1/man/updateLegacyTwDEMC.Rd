\name{updateLegacyTwDEMC}
\alias{updateLegacyTwDEMC}
\title{updateLegacyTwDEMC}
\description{update result of former version of twDEMC}
\usage{
updateLegacyTwDEMC(aTwDEMC)
}
\arguments{
  \item{aTwDEMC}{
}
}
\details{The function of the log of unnormalized density has been renamed from fLogLik to fLogDen
components "rLogLik","resFLogLikX" are renamed to "rLogDen","logDenCompX"

A twDEMC stores its call. Because argument names have been changed,
they are also changed within the \code{batchCall} entry: see \code{\link{updateLegacyTwDEMCArgs}}}


\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data( twdemcEx1)
aTwDEMC <- twdemcEx1
}
