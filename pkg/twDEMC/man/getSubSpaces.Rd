\name{getSubSpaces}
\alias{getSubSpaces}
\title{getSubSpaces}
\description{Recursively splitting the sample into subsamples and recording the upper and lower bounds}
\usage{
getSubSpaces(aSample, nSplit = 4, isBreakEarly = TRUE, argsFSplit = list(), pSub = 1, minPSub = 0.05)
}
\arguments{
  \item{aSample}{numeric matrix with parameters in columns: sample or subsample that is devided}
  \item{nSplit}{see \code{\link{findSplit}}, potentially modified due to argument \code{minPSub}}
  \item{isBreakEarly}{if set to FALSE all variable combinations are checked and an updated order of variables is returned}
  \item{argsFSplit}{further arguments to \code{\link{findSplit}}}
  \item{pSub}{the fraction that a subSample comprises within a bigger sample (used in recursive calls)}
  \item{minPSub}{minimum fraction a subSample, below which the sample is not split further}
}

\value{a list with entries
\item{spaces}{a list with an entry for each subspace. Each Entry is a list with
\item{sample}{numeric matrix: a subsample constrained to the subspace with col parameters}
\item{upperParBounds}{list with each entry numeric scalar: upper parameter bounds}
\item{lowerParBounds}{list with each entry numeric scalar: upper parameter bounds}
\item{pSub}{the proportion of the subSample to the overall Sample}
}}
\item{iVars}{integer vector: the indices of variables to check for splits. If \code{isBreakEarly=FALSE} its has been updated. NA if no split was found.}
\item{jVars}{integer vector: the indices of variables to check scales of variances given a split in iVar. If \code{isBreakEarly=FALSE} its has been updated. NA if no split was found.}
\item{iVars}{the index of variables to check for splits. If \code{isBreakEarly=FALSE} its has been updated.}
\item{jVars}{the index of variables to check scales of variances given a split in iVar. If \code{isBreakEarly=FALSE} its has been updated.
}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(den2dCorTwDEMC)
aSample <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
#mtrace(getSubSpaces)
subSpaces <- getSubSpaces(aSample, minPSub=0.4, isBreakEarly=FALSE )
str(subSpaces)
subSpaces <- getSubSpaces(aSample, minPSub=0.05)
(tmp <- sapply( subSpaces, function(subSpace){nrow(subSpace$sample)})/nrow(aSample))
#twUtestF(getSubSpaces)	# there are unit tests for this function
}
