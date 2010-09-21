\name{getAcceptedPos.twDEMCProps}
\alias{getAcceptedPos.twDEMCProps}
\title{getAcceptedPos twDEMCProps}
\description{Calculate the accepted state corresponding to the proposal.}
\usage{
getAcceptedPos.twDEMCProps(Y, acceptedRowName = "accepted")
}
\arguments{
  \item{Y}{matrix of proposals with row "accepted" and first step (column) initial state, and third dimension chains.}
  \item{acceptedRowName}{
}
}

\value{matrix of steps (steps x chains) whose state was the accepted state when evaluating proposal (row of Y)
the first step is NA because there is no previous accepted state}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





