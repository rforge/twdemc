\name{meanInputFluctuating}
\alias{meanInputFluctuating}
\title{meanInputFluctuating}
\description{provide mean + normal year to year error}
\usage{
meanInputFluctuating(obsLitter, padj)
}
\arguments{
  \item{obsLitter}{list of datastream matrices leaf and root with first three columns time, obs, and sdObs}
  \item{padj}{parameters}
}
\details{generates a series of fluctuating input.
Mean and sd of both leaf and root litter are calculated from provied obs and sdObs assuming independent errors
Correlation between leaf and root input is taken from padj$corrLeafRootLitter. 
If this is missing, correlation of 0.8 is assumed. }


\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(Howland14C)
str(tmp <- meanInputFluctuating( Howland14C$litter ))
plot( obs ~ times, data=tmp$leaf)
plot( tmp$leaf[,2] ~ tmp$root[,2] )
}
