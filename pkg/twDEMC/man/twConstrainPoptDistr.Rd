\name{twConstrainPoptDistr}
\alias{twConstrainPoptDistr}
\title{twConstrainPoptDistr}
\description{Constrain the information on parameters to selected parameters and add variance-covariance matrix and its inverse.}
\usage{
twConstrainPoptDistr(parNames,
parDistr,
corrMat = NULL)
}
\arguments{
  \item{parNames}{subsets of parameters: either character string, or indices, or boolean vector}
  \item{parDistr}{dataframe with columns trans, mu, and sigmaDiag, and variable names in rownames as returned by \code{link{twQuantiles2Coef}}}
  \item{corrMat}{correlations between parameters}
}



\author{Thomas Wutzler}





