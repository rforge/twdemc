\name{combineTwDEMCPops}
\alias{combineTwDEMCPops}
\title{combineTwDEMCPops}
\description{combine given populations to one big chain with more cases}
\usage{
combineTwDEMCPops(pops,
popCases = integer(0),
mergeMethod = "random",
nInSlice = 4)
}
\arguments{
  \item{pops}{list of population objects}
  \item{popCases}{integer vector (sum(getNSamples(pops))): specifying for each case (row) from which population it is filled}
  \item{mergeMethod}{method of merging the pops, see \code{\link{twMergeSequences}}}
  \item{nInSlice}{sequence length from each population (only for mergeMethod \code{slice}}
}
\details{entry splits is discarded, as it is not generally determined.}
\value{
\item{pop}{the merged population}
\item{popCases}{integer vector (nSample): population that case is taken from
}}

\author{Thomas Wutzler}




\examples{
data(den2dCorEx)
runs0 <- subPops(den2dCorEx$mcSubspaces0, iSpaces=c(1))
getSpacesPop(runs0)	# all populatins belonging to space replicate 1
#pops <- runs0$pops
#mtrace(combineTwDEMCPops)
res <- combineTwDEMCPops(runs0$pops) 
#str(res)
head(res$pop$parms[,,1])

res <- combineTwDEMCPops(runs0$pops[3:4])
#str(res)

}
