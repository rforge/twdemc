\name{initStateSoilMod}
\alias{initStateSoilMod}
\title{initStateSoilMod}
\description{Creating initial state variables for Basic Colimitation submodel.}
\usage{
initStateSoilMod(xc12, cn, iR, modMeta)
}
\arguments{
  \item{xc12}{numeric vector: 12C mass for each pool \cr
If names are given these must comprise code{modMeta$rowNames}.}
  \item{cn}{numeric vector: carbon to nitrogen ratio for each pool \cr
If vector has length one, the same ratio is assumed for all pools.
If names are given, these must comprise code{modMeta$rowNames}}
  \item{iR}{numeric matrix (nPools, nIsotopes): atomic ratios \cr
For others be aware of the adjusted units to avoid bigger numerical errors.
If it is a vector (nIsotopes), the same atomic ratio is assumed for all pools   
If colnames/names are given they comprise code{modMeta$colNames}
Default 1 is added for columns c12 and n15.}
  \item{modMeta}{model meta information, see \code{\link{modMetaICBM1}}}
}

\value{Numeric matrix (nPool, nIsotopes) of state variable mass.}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





