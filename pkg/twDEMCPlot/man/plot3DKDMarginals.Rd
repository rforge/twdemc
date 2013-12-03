\name{plot3DKDMarginals}
\alias{plot3DKDMarginals}
\title{plot3DKDMarginals}
\description{Evaluates the Kernel-density regression on a grid and displays a rgl-figure.}
\usage{plot3DKDMarginals(form, ds, probs = c(0.2, 0.5, 0.75, 0.95), 
    ..., dims = 8, bandwidth = NULL, knotSpacing = "quantile", 
    nSampleBandwidth = 500, nSample = 1000, bwtype = "adaptive_nn", 
    tol = 0.1, ftol = 0.1)}
\arguments{
  \item{form}{formula of regression e.g. rLogDen~tvrY+tvrO+biasLitterLeaf}
  \item{ds}{dataframe or matrix with column names corresponding to variables in the formula}
  \item{probs}{the percentiles at which to plot contour surfaces for the response variable
the levels are calculated from the subsample }
  \item{\dots}{further argument to \code{\link{plot.twApply3DMesh}}}
  \item{dims}{vector of points across predictor space within each dimension: scales bad: n=prod(dims)~dims[1]^3}
  \item{bandwidth}{see adaptive bandwidth of npreg.
integer vector: bandwidth for each dimension
about Records to include for estimating bandwidth.
integer vector: bandwidth for each dimension
If only one value is supplied, then it is repeated for each predictor
integer vector: bandwidth for each dimension}
  \item{knotSpacing}{see argument in \code{\link{twApply3DMesh}}}
  \item{nSampleBandwidth}{the sample size to estimate the bandwidth}
  \item{nSample}{The larger the sample, the longer it takes
sample used for interpolation}
  \item{bwtype}{parameters for \code{\link{npregbw}}}
  \item{tol}{
}
  \item{ftol}{
}
}
\details{Uses package snowfall to split calculation to subprocesses.
Requires package np loaded in cluster (sfLibrary(np))}
\value{3D Array of Kernel density regression predictions over the mesh. See \code{\link{twApply3DMesh}}.}

\author{Thomas Wutzler}




\examples{
#Example: Nested contours of mixture of three tri-variate normal densities
# wrapped inside a function, because it takes so long.
.tmp.f <- function(){
    nmix3 <- function(x, y, z, m, s) {
        0.4 * dnorm(x, m, s) * dnorm(y, m, s) * dnorm(z, m, s) +
                0.3 * dnorm(x, -m, s) * dnorm(y, -m, s) * dnorm(z, -m, s) +
                0.3 * dnorm(x, m, s) * dnorm(y, -1.5 * m, s) * dnorm(z, m, s)
    }
    f <- function(x,y,z) nmix3(x,y,z,.5,.5)
    
    n <- 250
    x <- rnorm(n,.5,.5)
    y <- c(rnorm(n/2,.5,.5), rnorm(n/2,-.5,.5)) 
    zz <- rnorm(n,.5,.5)
    
    require("twMiscRgl")   # for plotting the original function, instead of the kernel-density re-estimate 
    plot(tmp <- twApply3DMesh(x,y,zz,f, nSample=1000, dims=10))	# just the points (nSample dominates over dims)
    plot(tmp, probs=seq(0.5, 0.95, len=4))
    tmpObs <- tmp$sample
    tmpObs$f2 <- tmpObs$f + rnorm(nrow(tmpObs),sd=0.025)
    plot(tmpObs$f2~tmpObs$f)
    plot3d( tmpObs, col=rev(heat.colors(20))[round(twRescale(tmpObs$f2,c(1,20)))]  )
    
    #reestimate density from disturbed sample
    # kernel density regression estimate
    #mtrace(plot3DKDMarginals)
    #sfInit(TRUE,cpus=2); sfLibrary(np)
    tmp2 <- plot3DKDMarginals( f2 ~ x+y+zz, tmpObs, dims=10, nSampleBandwidth=500  )
    tmp2b <- tmp2; tmp2b$sample<-NULL; plot(tmp2b)	# show the mesh instead of the sample
    plot(tmp2, probs=c(0.25,0.5,0.75,0.9), nDrawPoints=0)		#contours only
}
}
