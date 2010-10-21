\name{meanInputFluctuating}
\alias{meanInputFluctuating}
\title{meanInputFluctuating}
\description{provide mean + normal year to year error}
\usage{
meanInputFluctuating(input, padj)
}
\arguments{
  \item{input}{list of datastream matrices with first two columns time and value}
  \item{padj}{parameters}
}



\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(Howland14C)
str(tmp <- meanInputFluctuating( Howland14C$litter ))
plot( obs ~ times, data=tmp$leaf)
plot( tmp$leaf[,2] ~ tmp$root[,2] )
}
