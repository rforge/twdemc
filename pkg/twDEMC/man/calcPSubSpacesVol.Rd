\name{calcPSubSpacesVol}
\alias{calcPSubSpacesVol}
\title{calcPSubSpacesVol}
\description{estimate of the integrated density of the subspace by average density in volume
subspace and population are used interchangeably here.}
\usage{calcPSubSpacesVol(ssLogDenL, qPop, iPopsSpace, nSamplePop = sapply(ssLogDenL, 
    length))}
\arguments{
  \item{ssLogDenL}{list (nPop) of vector: the (temperated) logDensity of each record in a}
  \item{qPop}{the previous quantile of samples within space}
  \item{iPopsSpace}{list (nSpace): of integer vectors: pops in space}
  \item{nSamplePop}{the number of samples in subspace}
}



\author{Thomas Wutzler}





