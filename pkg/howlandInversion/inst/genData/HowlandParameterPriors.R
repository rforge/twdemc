

################################ The parameters ################################

# HowlandParameterPriors
.tmpf.biasDiffRespLitterfall <- function(){
	# calculating the sd of the bias of the difference between respiration and litterfall
	data(Howland14C)
	obs <- Howland14C$obsNutrientSite$respCum
	(resp <- mean(obs[,"obs"]))
	(sdResp <- sqrt(sum(obs[,"sdObs"]^2))/nrow(obs))		# variances sum up 
	obs <- Howland14C$litter$leaf
	(leafLitter <- mean(obs[,"obs"]))
	(sdLeafLitter <- sqrt(sum(obs[,"sdObs"]^2))/nrow(obs))		# variances sum up
	(sdDiff <- sqrt( sdResp^2 + sdLeafLitter^2 ))
}

# units are gC/m2 /yr
parms = parms0 = within( list(), {
		tLagLeaf=1	##<< less than 1 yr is recent-C, which is not accounted in the Howland study
		tLagRoot=5	##<< time lag between C-fixation and root turnover
		kY=1		##<< decay constant
		kO=1/30		
		h=0.01		##<< humification coefficient
		cY=0.1		##<< initial proportion of total C in young pool
		iROLayerCalcRelErr = 0.01 ##<< 1% error of reconstruction of O-Layer iR-ratio
		biasDiffRespLitterfall = 0  ##<< both litterfall and mean respiration may be biased,		
			##<< hence also root litter input - which is the difference betweent these two.
	})	
		
#-------------- a priori knowlege about the parameters ----------
# define distributions by mle and upper confidence bound,
upperBoundProb = 0.99
parmsBounds = list(
		tLagLeaf = parms$tLagLeaf * c(1,2)		
		,tLagRoot = parms$tLagRoot * c(1,3)
		,kY = parms$kY *c(1,10) 
		,kO = parms$kO *c(1,10)
		,h = c(parms$h, 0.977)
		,cY= c(parms$cY, 0.965)
		,biasDiffRespLitterfall= c(parms$biasDiffRespLitterfall, qnorm(upperBoundProb,sd=41.4))
	)
#which(sapply(parmsBounds,length)!=2)
		
#------- may have updated the mode in parms.bounds, transfer this value to parms and parms0
parms[ names(parmsBounds) ] <- sapply( parmsBounds, function(qList)qList[[1]] )
names(parms)[ which( !(names(parms) %in% names(parms0)))]
tmp.d <- unlist(parms[ names(parmsBounds) ]) - unlist(parms0[ names(parmsBounds) ])
tmp.d[ tmp.d != 0]
parms0 <- parms
		
#----------------- distributions of the variables
varDistr <- twVarDistrVec( names(parmsBounds) )
varDistr[] <- "lognorm"			#by default assume lognormal (0,Inf)
varDistr[c("h","cY")] <- "logitnorm" #logit-normal (0,1)
varDistr[c("biasDiffRespLitterfall")] <- "norm" #normal
# pStorageOnActivation is modeled lognorm, because increase in prior at 1 is not feasable

#----------------- calculate the standard mu and sd at normal scale from quantiles
#mtrace(twQuantiles2Coef)
#mtrace(
parDistr <- twQuantiles2Coef( parmsBounds, varDistr, upperBoundProb=upperBoundProb )

# do the flattest logitnormal distribution possible for cY and h
parDistr$mu[c("cY","h")] <-0
parDistr$sigmaDiag[c("cY","h")] <- twSigmaLogitnorm(unlist(parms0[c("cY","h")]))[,2]


HowlandParameterPriors <- list(
	parDistr=parDistr
	,parms0=parms0
)

save(HowlandParameterPriors, file="data/HowlandParameterPriors.RData")


#.tmp.f <- function(){
	HowlandParameterPriors <- NULL
	data(HowlandParameterPriors)
	str(HowlandParameterPriors)
	
	parDistr <- poptDistr <- HowlandParameterPriors$parDistr
	#poptDistr <- twConstrainPoptDistr(c("biasRespSum","F0","epsF","epsG"),parDistr)
	#ggplotDensity.poptDistr(poptDistr,parmsBounds=parmsBounds,doTransOrig=FALSE)
	ggplotDensity.poptDistr(poptDistr,parmsBounds=parmsBounds,pMin=0.025,plotUpperQuantile=FALSE)
}

#parms.tmp <- parms; parms.tmp$biasP14C=0.4; obsadj <- calcBiasedObs.hamer(obs, parms.tmp)



