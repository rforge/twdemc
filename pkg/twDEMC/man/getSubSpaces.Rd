\name{getSubSpaces}
\alias{getSubSpaces}
\title{getSubSpaces}
\description{Recursively splitting the sample into subsamples and recording the upper and lower bounds}
\usage{
getSubSpaces(aSample,
nSplit = 4,
isBreakEarly = TRUE,
isBreakEarlySubs = TRUE,
checkSlopesFirst = data.frame(ivar = 1,
j1AlphaSlope = 1,
jAlphaSlope = 2)[FALSE,
],
argsFSplit = list(),
pSub = 1,
minPSub = 0.05,
upperParBounds = numeric(0),
lowerParBounds = numeric(0),
splits = numeric(0),
maxNSample = 128,
verbose = FALSE)
}
\arguments{
  \item{aSample}{numeric matrix with parameters in columns: sample or subsample that is devided (do not include rLogLik in first column)}
  \item{nSplit}{see \code{\link{findSplit}}, potentially modified due to argument \code{minPSub}}
  \item{isBreakEarly}{if TRUE and argument \code{checkSlopesFirst} is given, then check slope angles for variables given in checkSlopesFirst first and if break is found, do not evaluate all the other slope angles}
  \item{isBreakEarlySubs}{same as argument \code{isBreakEarly} for recursive call from within}
  \item{checkSlopesFirst}{data.frame with entries ivar, j1AlphaSlope and j2AlphaSlope as returned in entry resD in result of ressplit}
  \item{argsFSplit}{further arguments to \code{\link{findSplit}}}
  \item{pSub}{the fraction that a subSample comprises within a bigger sample (used in recursive calls)}
  \item{minPSub}{minimum fraction a subSample, below which the sample is not split further}
  \item{upperParBounds}{named numeric vectors: giving upper parameter bounds: lowerBound < par <= upperBound
for exploring subspaces of the limiting distribution, see details}
  \item{lowerParBounds}{similar to upperParBounds: sample > bound}
  \item{splits}{named numeric numeric vector: history of splitting points}
  \item{maxNSample}{if given a value, then aSample is thinned before to given number of records (for efficiently calculating variances)}
  \item{verbose}{if TRUE then prints on every terminating subspace}
}

\value{a list with entries
\item{spaces}{a list with an entry for each subspace. Each Entry is a list with entries \itemize{
\item{sample}{numeric matrix: a subsample constrained to the subspace lb < val <= ub with col parameters}
\item{upperParBounds}{list with each entry numeric scalar: upper parameter bounds}
\item{lowerParBounds}{list with each entry numeric scalar: upper parameter bounds}
\item{pSub}{the proportion of the subSample to the overall Sample}
\item{splits}{named numeric vector of splitting points}
}}
\item{resD}{result details of \code{\link{findSplit}}
}}

\author{Thomas Wutzler}




\examples{
data(den2dCorEx)
aSample <- stackChains(thin(concatPops(den2dCorEx$mcBulk), start=300))[,-1]
#mtrace(getSubSpaces)
subSpaces <- getSubSpaces(aSample, minPSub=0.4 )
str(subSpaces)
subSpaces <- getSubSpaces(aSample, minPSub=0.05, argsFSplit=list(debugSequential=TRUE))
subSpaces <- getSubSpaces(aSample, minPSub=0.05)
(tmp <- sapply( subSpaces$spaces, function(subSpace){nrow(subSpace$sample)})/nrow(aSample)) # percentiles
(tmp <- sapply( subSpaces$spaces, "[[", "upperParBounds" )) # bounds
#visualize the splits
ss1 <- do.call(rbind, lapply( seq_along(subSpaces$spaces), function(i){cbind(spaceInd=i, subSpaces$spaces[[i]]$sample)}))
plot(b~a, as.data.frame(ss1), col=rainbow(length(subSpaces$spaces))[spaceInd] )
plot(b~a, as.data.frame(ss1), col=rainbow(length(subSpaces$spaces))[spaceInd], ylim=c(-5000,+5000) )
#twUtestF(getSubSpaces)	# there are unit tests for this function
}
