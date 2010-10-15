\name{twCreateModMeta}
\alias{twCreateModMeta}
\title{twCreateModMeta}
\description{Creating meta-information about rows (stateVars), isotopes (columns), and additonal outputs. }
\usage{
twCreateModMeta(rowNames, csts, iRUnits = 1, auxGroups = list())
}
\arguments{
  \item{rowNames}{character vector: names of the state variables}
  \item{csts}{list of character vectors: correponding to isotopes}
  \item{iRUnits}{numeric vector of units for different isotopes, ensuring equal magnitudes to avoid numerical errors (see details)\cr}
  \item{auxGroups}{namegroup (character) -> columnNames (character vector): column names of auxiliary outputs}
}
\details{\describe{\item{iRUnits}{ 
By default all units are 1.
If names are not given, is is assumed that order corresponds to colNames
If names are given, default 1 for c12 and n15 are added.
}}}
\value{a list with components \describe{
\item{rowNames}{character vector: names of the state variables}
\item{csts}{list of character vectors: correponding to isotopes}
\item{colNames}{character vector: names of the isotopes}
\item{elementNames}{character vector: names of the elements row_column }
\item{iRUnits}{numeric vector of units for different isotopes, ensuring equal magnitudes to avoid numerical errors (see details)\cr}	
\item{auxGroups}{namegroup (character) -> columnNames (character vector): column names of auxiliary outputs}
\item{auxOutputNames}{character vector: names of all auxiliary outputs}	
\item{auxOutputTemplate}{numeric vector, with names corresponding to auxOutputNames }	
\item{nRows}{number of rows}	
\item{nCols}{number of columns}	
\item{matrixTemplate}{numeric matrix: template for state variable matrix with corresponding row and column names}	
\item{nAux}{number of auxiliary outputs}	
}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





