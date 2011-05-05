\name{calcSteadyK_ICBM1}
\alias{calcSteadyK_ICBM1}
\title{calcSteadyK ICBM1}
\description{calculte decay constants from assuming steady state and remaining parameters}
\usage{
calcSteadyK_ICBM1(Ctot, iY, parms)
}
\arguments{
  \item{Ctot}{SOM C-Stock}
  \item{iY}{steady state input}
  \item{parms}{list with entries kY,kO and cY and h}
}

\value{named numeric matrix with columns kY and kO, rows corresponding to longest input parameter (others are recycled)}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





