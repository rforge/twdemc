fCost <- function(theta){
    pred <- getUpdateFunction(predSpec)( theta, list(), xSparse, xRich )
    logDenComp <- logDenTheta( theta, list(pred=pred)
            , xSparse=xSparse, xRich=xRich, obs=obs, sdObs=sdObs )
    sum(logDenComp)
}
fCost(theta0)
res0 <- optim(theta0, fCost, control =list(fnscale=-1), hessian=TRUE ) # , method="SANN")
thetaBest <- res0$par
thetaCov <- -solve(res0$hessian)
sampleGradIgnore <- theta <- rmvnorm( 800, thetaBest, sigma=thetaCov )


# fix GP parameters: signal variance and correlation length 
sigma2d1 <- 1.5 * sdObs$y1^2
sigma2d2 <- 1.5 * sdObs$y2^2
psi1 <- diff(range(xSparse))/3 
psi2 <- diff(range(xRich))/3
# fix supporting locations
o1 <- seq(1, length(xSparse), length.out=4)
o2 <- seq(1, length(xRich), length.out=4)
# compute covariance matrix at supporing locations
fCor <- function(x1,x2,psi){  exp(-((x1-x2)/psi)^2) }
Kss1 <- sigma2d1 * outer( xSparse[o1], xSparse[o1], fCor, psi=psi1)
Kss2 <- sigma2d2 * outer( xRich[o2], xRich[o2], fCor, psi=psi2)
Kz1 <- Kss1 + diag(sdObs$y1^2, nrow=nrow(Kss1))
Kz2 <- Kss2 + diag(sdObs$y2^2, nrow=nrow(Kss2))
# compute covariance matrix at remaining locations
Krs1 <- sigma2d1 * outer( xSparse[-o1], xSparse[o1], fCor, psi=psi1)
Krs2 <- sigma2d2 * outer( xRich[-o2], xRich[o2], fCor, psi=psi2)

fCostGP <- function(theta){
    pred <- getUpdateFunction(predSpec)( theta, list(), xSparse, xRich )
    # compute expected value of model discrepancies
    z1 <- obs$y1 - pred$y1
    kZ_z1 <- solve(Kz1, z1[o1])
    delta1 <- numeric( length(pred$y1))
    delta1[o1] <- deltaS1 <- Kss1 %*% kZ_z1  
    delta1[-o1] <- deltaR1 <- Krs1 %*% kZ_z1
    z2 <- obs$y2 - pred$y2
    kZ_z2 <- solve(Kz2, z2[o2])
    delta2 <- numeric( length(pred$y2))
    delta2[o2] <- deltaS2 <- Kss2 %*% kZ_z2  
    delta2[-o2] <- deltaR2 <- Krs2 %*% kZ_z2
    # compute penalty term of model discrepancies
    logPDelta1 <- as.vector( -1/2* delta1[o1] %*% solve(Kss1, delta1[o1]) )    
    logPDelta2 <- as.vector( -1/2* delta2[o2] %*% solve(Kss2, delta2[o2]) )    
    # compute Likelihood of observations
    misfit1 <- obs$y1 - (pred$y1 + delta1)
    logPObs1 <- -1/2*sum((misfit1/sdObs$y1)^2)
    misfit2 <- obs$y2 - (pred$y2 + delta2)
    logPObs2 <- -1/2*sum((misfit2/sdObs$y2)^2)
    #
    ans <- logPObs1 + logPDelta1 + logPObs2 + logPDelta2 
    ans
}
fCostGP(theta0)
resGP <- optim(theta0, fCostGP, hessian=TRUE, control =list(fnscale=-1) ) 
thetaBestGP <- resGP$par
thetaCovGP <- -solve(resGP$hessian)
sampleGradGP <- theta <- rmvnorm( 800, thetaBestGP, sigma=thetaCovGP )



samplePostGrad <- rbind( doSamplePred(sampleGradIgnore,"gradIgnore") , doSamplePred(sampleGradGP,"gradGP") )
iStream <- "yRich"
g1 <- ggplot( subset(samplePostGrad, stream==iStream), aes(x=x ) ) +
        geom_point(aes(y=obs), col="grey50") +  
        facet_grid(~scenario, scales="free_x") +  
        geom_line(aes(y=median)) + 
        geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.4) +
        #ylab("respiration rate") + xlab("Temperature") +
        xlab("") + ylab("obs rich") + 
        theme_bw(base_size=baseSize)+
        theme(axis.title.x=element_blank()) +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()
                , axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
iStream <- "ySparse"
g2 <- ggplot( subset(samplePostGrad, stream==iStream), aes(x=x ) ) +
        geom_point(aes(y=obs), col="grey50") +  
        facet_grid(~scenario, scales="free_x") +  
        geom_line(aes(y=median)) + 
        geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.4) +
        #ylab("cumulated respiration") + xlab("Soil carbon stocks") +
        xlab("x") + ylab("obs sparse") + 
        theme_bw(base_size=baseSize)+
        theme(axis.text.x = element_blank(), axis.text.y = element_blank()
                , axis.ticks.x = element_blank(), axis.ticks.y = element_blank())  
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,1)))
print(g1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(g2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))

