\name{calcSteadyK_ICBM1_mat}
\alias{calcSteadyK_ICBM1_mat}
\title{calcSteadyK ICBM1 mat}
\description{calculte decay constants from assuming steady state and remaining parameters}
\usage{
calcSteadyK_ICBM1_mat(Ctot, iY, parms)
}
\arguments{
  \item{Ctot}{SOM C-Stock}
  \item{iY}{steady state input}
  \item{parms}{named numeric vector with entries kY,kO and cY and h}
}

\value{named numeric matrix with columns kY and kO, rows corresponding to longest input parameter (others are recycled)}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





