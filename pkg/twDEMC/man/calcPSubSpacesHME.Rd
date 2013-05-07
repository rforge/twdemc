\name{calcPSubSpacesHME}
\alias{calcPSubSpacesHME}
\title{calcPSubSpacesHME}
\description{estimate of the integrated density of the subspace by harmonic mean estimate (HME)
subspace and population are used interchangeably here.}
\usage{calcPSubSpacesHME(ssLogDenL, qPop, iPopsSpace, nSamplePop = sapply(ssLogDenL, 
    length), qOmitLow = 0.05)}
\arguments{
  \item{ssLogDenL}{list (nPop) of vector: the (temperated) logDensity of each record in a}
  \item{qPop}{the previous quantile of samples within space}
  \item{iPopsSpace}{list (nSpace): of integer vectors: pops in space}
  \item{nSamplePop}{the number of samples in subspace}
  \item{qOmitLow}{quantile of low logDensities to omit in integral to stabilize HME estimate}
}



\author{Thomas Wutzler}





