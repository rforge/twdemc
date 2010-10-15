\name{twStateMatODERes}
\alias{twStateMatODERes}
\title{twStateMatODERes}
\description{returns the matrix of state variables for an ode output corresponding to modMeta}
\usage{
twStateMatODERes(out, iSteps, modMeta)
}
\arguments{
  \item{out}{result of ode to model of matrix statevariables}
  \item{iSteps}{integer vector index of outputs}
  \item{modMeta}{description of the model}
}

\value{array (nRow,nCol,length(iSteps)) of stacked state variable matrices for each step }

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}





