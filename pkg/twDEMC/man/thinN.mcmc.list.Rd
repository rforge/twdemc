\name{thinN.mcmc.list}
\alias{thinN.mcmc.list}
\title{thinN mcmc list}
\description{Thin x so that each chain consists of about nThinnedRecords}
\usage{\method{thinN}{mcmc.list}(x, nThinnedRecords = 100, ...)}
\arguments{
  \item{x}{object to be thinned and plotted (mcmc.list)}
  \item{nThinnedRecords}{number of records in thinned object}
  \item{\dots}{
}
}
\details{Plotting mcmc.list with large n (>1000) takes long and involves much overplotting.
The same visual impression can be achieved using a thinned chain.}
\value{Thinned x}

\author{Thomas Wutzler}




\examples{FALSE}
