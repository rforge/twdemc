\name{twModElementNames}
\alias{twModElementNames}
\title{twModElementNames}
\description{Construct element names for given rows and columns}
\usage{twModElementNames(rowNames, colNames)}
\arguments{
  \item{rowNames}{A character vector of row names.}
  \item{colNames}{A character vector of column names.}
}

\value{Character vector for each index (by row) rowName, colName}

\author{Thomas Wutzler}




\examples{
mmd <- modMetaICBMDemo()
twModElementNames("Y",mmd$csts$cis)
twModElementNames(c("Y","O"),mmd$csts$cis)
}
