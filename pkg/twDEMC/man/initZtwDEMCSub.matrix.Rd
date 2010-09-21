\name{initZtwDEMCSub.matrix}
\alias{initZtwDEMCSub.matrix}
\title{initZtwDEMCSub matrix}
\description{generates an appropriate initial sample of parameter vectors for twDEMC from subsampling an array}
\usage{
initZtwDEMCSub.matrix(Zinit1, vars = colnames(Zinit1), nChains = 4, nPops = 1, m0 = calcM0twDEMC(length(unique(vars)), nPops, nChains), ...)
}
\arguments{
  \item{Zinit1}{the mcmc matrix to subsample (column variable, rows cases)}
  \item{vars}{which variables to keep}
  \item{nChains}{
}
  \item{nPops}{
}
  \item{m0}{number of required cases for initialization}
  \item{\dots}{
}
}

\value{an array of dimension suitable for Zinit for twDEMCInt}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}



\seealso{\code{\link{initZtwDEMCNormal}}}

