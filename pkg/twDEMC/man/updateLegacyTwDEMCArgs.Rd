\name{updateLegacyTwDEMCArgs}
\alias{updateLegacyTwDEMCArgs}
\title{updateLegacyTwDEMCArgs}
\description{update list of arguments to of former version of twDEMC}
\usage{
updateLegacyTwDEMCArgs(argsF)
}
\arguments{
  \item{argsF}{a list of arguments to twDEMC}
}
\details{Arguments "fLogLik","argsFLogLik" are renamed to "fLogDen","argsFLogDenX"}


\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data( twdemcEx1)
argsF <- as.list(attributes(twdemcEx1)$batchCall)
updateLegacyTwDEMCArgs(argsF)
}
