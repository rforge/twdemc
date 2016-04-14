#require(testthat)
context("mvnormOpt")

test_that("mvnormOpt",{
            dimTheta <- 5L
            nObs <- 20L
            thetaMean <- rnorm(dimTheta)
            #.corG from GPSampling
            Sigma <- outer(1:dimTheta, 1:dimTheta, .corG, psi=dimTheta/2)
            cholSigma <- chol(Sigma)
            theta <- rmvnorm(nObs, thetaMean, Sigma)
            r1 <- dmvnorm(theta, thetaMean, Sigma, log=TRUE)
            r2 <- dmvnormOpt(theta, thetaMean, cholSigma=cholSigma, log=TRUE)
            expect_equal(r1,r2)
        })

