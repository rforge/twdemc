\name{twSetModMetaAuxGroups}
\alias{twSetModMetaAuxGroups}
\title{twSetModMetaAuxGroups}
\description{Setting auxiliare output items in modMeta.}
\usage{twSetModMetaAuxGroups(modMeta, auxGroupsNew, auxGroupsSolve = list())}
\arguments{
  \item{modMeta}{the data-structure (result of \code{\link{twCreateModMeta}}) to modify}
  \item{auxGroupsNew}{mapping namegroup (character) -> columnNames (character vector): column names of auxiliary outputs}
  \item{auxGroupsSolve}{same format as auxGroups, but pertaining to outputs calculated after integrating the model}
}

\value{modified modMeta, adjusted for auxGroups, auxOutputnames, auxOutputTemplate, and nAux.
In addition the following elements are added:
\itemize{
\item auxGroupsExt: auxGroups + auxGroupsSolve
\item auxGroupsExtNames: flat character vector of auxiliary output names
} }

\author{Thomas Wutzler}





