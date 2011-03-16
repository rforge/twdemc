\name{kinrespModel}
\alias{kinrespModel}
\title{Microbial respiration model}
\description{Respiration for given microbial parameters at given time.}
\usage{
kinrespModel(x0, r0, mumax, time, lambda = 0.9, YCO2 = 1.5)
}
\arguments{
  \item{x0}{initial microbial biomass (numeric scalar)}
  \item{r0}{initial microbial activity (numeric scalar)}
  \item{mumax}{maximum growth rate (numeric scalar)}
  \item{time}{time (numeric vector)}
  \item{lambda}{Ratio of growth associated (coupled) specific respiration to total specific respiration. Usually 0.9.}
  \item{YCO2}{Ratio of assimilated carbon per respired carbon.  Usually 1.5.}
}

\value{respiration at given time points (numeric vector)}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(dsKinrespTut)
ds <- subset(dsKinrespTut[order(dsKinrespTut$time),], replicate==1 & experiment=="9")
require(nlme)
gls1 <- gnls( resp ~ kinrespModel( x0, r0, mumax, time), ds, start=c(x0=140, r0=2e-3, mumax=0.24) )
plot( resp~time, data=ds)
lines(fitted(gls1)~ds$time)
}
