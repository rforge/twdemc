\name{updateBlockTwDEMC}
\alias{updateBlockTwDEMC}
\title{updateBlockTwDEMC}
\description{Perfrom one DEMC step, function to be called in remote process.}
\usage{updateBlockTwDEMC(xC, argsFUpdateBlock)}
\arguments{
  \item{xC}{numeric vector: current state that is used in density function of the block}
  \item{argsFUpdateBlock}{further argument provided for generating the update  \describe{
\item{compPosInDen}{ positions of the dimensions in x that are updated in this block}
\item{step}{proposed jump}
\item{rExtra}{extra portion in Metropolis decision because of selecting the jump}
\item{logDenCompC}{numeric vector: former result of call to same fLogDen}
\item{tempDenCompC}{numeric vector of length(logDenCompC): temperature for each density result component}
\item{fDiscrProp,argsFDiscrProp}{function and additional arguments applied to xProp, e.g. to round it to discrete values}
\item{argsFLogDen}{additional arguments to fLogDen and scalar factor applied to result of fLogDen}
\item{posLogDenInt}{the matching positions of intResCompNames within the the results components that are handled internally}
\item{ctrl$DRgamma}{ if !0 and >0 delayed Rejection (DR) (Haario06) is applied by jumping only DRgamma distance along the proposal }
\item{upperParBounds}{ named numeric vector, see \code{\link{twDEMCBlockInt}}  }
\item{lowerParBounds}{ named numeric vector, see \code{\link{twDEMCBlockInt}}  }
\item{fCalcComponentTemp}{ functiont to calculate temperature of result components, (way of transporting calcComponentTemp to remote process) }
\item{maxLogDen}{ integer vector (nResComp): maximum logDensity (overfitting control, usually -1/2 nObs) }
\item{iPop}{ just for debugging }
}}
}
\details{\describe{\item{internal Metropolis step}{
if posLogDenInt is given, then these components of result of fLogDen are handled
internally. Hence, for Metropolis step here operates only on other remaining components.
}}

\describe{\item{internal Metropolis step}{
if posLogDenInt is given, then these components of result of fLogDen are handled
internally. Hence, for Metropolis step here operates only on other remaining components.
}}}
\value{list with components
\item{accepted}{boolean scalar: if step was accepted}
\item{xC}{numeric vector: components of position in parameter space that are being updated (if accepted then the same as xP)}
\item{logDenCompC}{numeric vector: result components of fLogDen for current position (if accepted then the same as logDenCompP)}
\item{xP}{numeric vector: components of proposal that are being}
\item{logDenCompP}{numeric vector: result components of fLogDen for proposal, may inlcude attr "intermediate"
}}

\author{Thomas Wutzler}



\seealso{\code{\link{twDEMCBlockInt}}
\code{.updateBlocksTwDEMC}}

