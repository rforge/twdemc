\name{getSubSpaces}
\alias{getSubSpaces}
\title{getSubSpaces}
\description{Recursively splitting the sample into subsamples and recording the upper and lower bounds}
\usage{
getSubSpaces(aSample, nSplit = 4, isBreakEarly = TRUE, isBreakEarlySubs = TRUE, checkSlopesFirst = data.frame(ivar = 1, j1AlphaSlope = 1, jAlphaSlope = 2)[FALSE, ], argsFSplit = list(), pSub = 1, minPSub = 0.05, splitHist = NA_real_[FALSE])
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
  \item{splitHist}{named numeric numeric vector: history of splitting points}
}

\value{a list with entries	
	\item{spaces}{a list with an entry for each subspace. Each Entry is a list with entries \itemize{	
		\item{sample}{numeric matrix: a subsample constrained to the subspace with col parameters}		
		\item{upperParBounds}{list with each entry numeric scalar: upper parameter bounds}		
		\item{lowerParBounds}{list with each entry numeric scalar: upper parameter bounds}		
		\item{pSub}{the proportion of the subSample to the overall Sample}		
		\item{splits}{named numeric vector of splitting points}		
		}}
	\item{resD}{result details of \code{\link{findSplit}}	
}}

\author{Thomas Wutzler <twutz@bgc-jena.mpg.de>}




\examples{
data(den2dCorTwDEMC)
aSample <- stackChains(thin(den2dCorTwDEMC, start=300))[,-1]
#mtrace(getSubSpaces)
subSpaces <- getSubSpaces(aSample, minPSub=0.4 )
str(subSpaces)
subSpaces <- getSubSpaces(aSample, minPSub=0.05, argsFSplit=list(debugSequential=TRUE))
subSpaces <- getSubSpaces(aSample, minPSub=0.05)
(tmp <- sapply( subSpaces$spaces, function(subSpace){nrow(subSpace$sample)})/nrow(aSample)) # percentiles
(tmp <- sapply( subSpaces$spaces, "[[", "upperParBounds" )) # bounds
#visualize the splits
ss1 <- do.call(rbind, lapply( seq_along(subSpaces$spaces), function(i){cbind(iSpace=i, subSpaces$spaces[[i]]$sample)}))
plot(b~a, as.data.frame(ss1), col=rainbow(length(subSpaces$spaces))[iSpace] )
plot(b~a, as.data.frame(ss1), col=rainbow(length(subSpaces$spaces))[iSpace], ylim=c(-5000,+5000) )
#twUtestF(getSubSpaces)	# there are unit tests for this function
}
