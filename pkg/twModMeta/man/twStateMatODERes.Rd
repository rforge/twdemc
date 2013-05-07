\name{twStateMatODERes}
\alias{twStateMatODERes}
\title{twStateMatODERes}
\description{converting a state in vector format to matrix (for deSolve output corresponding to modMeta)}
\usage{twStateMatODERes(out, iSteps, modMeta)}
\arguments{
  \item{out}{result of ode to model of matrix statevariables}
  \item{iSteps}{integer vector index of outputs}
  \item{modMeta}{description of the model}
}

\value{array (nRow,nCol,length(iSteps)) of stacked state variable matrices for each step }

\author{Thomas Wutzler}





