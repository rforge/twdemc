\name{twModMeta-package}
\alias{twModMeta-package}
\alias{twModMeta}
\title{twModMeta}
\description{ODE-States consisting of several elements/isotopes.}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}


\details{
    This package supports the use of array based state variables with package deSolve.
    
    This often occurs when several isotopes of the same element need to be calculated. The matrix of elements in rows 
    and isotopes in columns, is converted to state vector with entries named: row_col
        
    An example \link{modMetaICBMDemo} with the simple ICBM carbon turnover model is provided using radiocarbon data.

\describe{\item{Basic meta-information structure and helper functions includes}{
\itemize{
\item setting up a structure about the models Metadata \code{\link{twCreateModMeta}}
\item adding auxiliary outputs to the Metadata, \code{\link{twSetModMetaAuxGroups}} 
\item initializing the state, \code{\link{initStateModMeta}} 
\item converting state vector to state matrix, \code{\link{twStateMatODERes}} 
\item constructing entry names \code{\link{twModElementNames}} 
} %itemize
}} %describe item

\describe{\item{Data}{
\itemize{
\item radiocarbon constants: \code{\link{c14Constants}}
\item radiocarbon atmospheric record of the bomb peak: \code{\link{delta14Catm}} and function \code{\link{fmAtmosphere}}.
} %itemize
}} %describe item

} %details


\keyword{ package }


