\name{twDEMC.array}
\alias{twDEMC.array}
\title{twDEMC array}
\description{Initialize \code{\link{twDEMCInt}} by array of initial population and remove those generations from results afterwards}
\usage{
\method{twDEMC}{array}(Zinit, ...)
}
\arguments{
  \item{Zinit}{initial population: a numeric array (d x M0 x Npop) see details in \code{\link{twDEMCInt}}}
  \item{\dots}{further arguments to \code{\link{twDEMCInt}}}
}

\value{result of \code{\link{twDEMCInt}} with initial M0-1 cases removed}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





