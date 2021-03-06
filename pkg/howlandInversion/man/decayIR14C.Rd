\name{decayIR14C}
\alias{decayIR14C}
\title{decayIR14C}
\description{Calculate the atomic ratio for given years assuming only radioactive decay.}
\usage{
decayIR14C(yr, iR0 = c14Constants$iR14CStandard, yr0 = c14Constants$yr14CStandard, lambda = c14Constants$lambda)
}
\arguments{
  \item{yr}{numeric vector for years}
  \item{iR0}{isotopic ratio at yr0}
  \item{yr0}{year when iR0 is given}
  \item{lambda}{decay constant}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
yr=seq(-10000,2000,length.out=30)
plot( decayIR14C(yr)~yr )
}
