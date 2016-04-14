\name{twDEMC-package}
\alias{twDEMC-package}\alias{twDEMC}
\title{parallel DEMC}
\description{Parallel version of a differential evolution markov chain (DEMC) combined with simulated annealing (SA). DEMC is a variant of Monte Carlo Markov Chain (MCMC).}
\details{Method \code{\link{twDEMCSA}} automatically infers proper temperature decrease in simulated annealing (SA),based on mixing of different independent populations of DEMC chains.Method \code{\link{twDEMCBlock.twDEMCPops}} can be used to continue sampling with explicitly specifying starting and end temperature.Further functionality of the twDEMC package deals with \itemize{\item{ Generating an initial population for twDEMC: \code{\link{initZtwDEMCNormal}}  }\item{ Transforming the objects of type \code{twDEMCPops}: \code{\link{subset.twDEMCPops}}  }\item{ Transforming the objects of type \code{twDEMC} obtained by \code{\link{concatPops.twDEMCPops}}: \code{\link{subChains.twDEMC}}  }\item{ Transforming the parameter space: \code{\link{transOrigPopt.default}}  }\item{ Invoking fLogDen with proposal in a parallel load balanced way: \code{\link{twCalcLogDenPar}} and \code{\link{twCalcLogDensPar}}  }\item{ Calculating marginal aggregates: \code{\link{marginals1d}}  }\item{ Laplace approximation of the normalizing constant: \code{\link{normConstLaplace}} and \code{\link{normConstLaplaceBridge}}  }\item{ Plotting: \code{\link{plotThinned.mcmc.list}}  }\item{ Scaling logDensities by temperature: \code{\link{calcBaseTemp}}, \code{\link{calcTemperatedLogDen.matrix}}  }\item{ Ordering multiple log-Densities \code{\link{orderLogDen}}  }}%itemizeThe workhorse for parallelized sampling of different populations and chains is function \code{\link{twDEMCBlockInt}}, which is usually called by \code{\link{twDEMCBlock.array}} or \code{\link{twDEMCBlock.twDEMCPops}}}%details
\author{Thomas Wutzler}

\keyword{ package }


