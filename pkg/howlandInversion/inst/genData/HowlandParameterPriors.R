

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
		kY=kYNarrow=1	##<< decay constant
		kO=1/30		
		h=0.01		##<< humification coefficient
		cY=0.1		##<< initial proportion of total C in young pool
		iROLayerCalcRelErr = 0.01 ##<< 1% error of reconstruction of O-Layer iR-ratio
		biasDiffRespLitterfall = 0  ##<< both litterfall and mean respiration may be biased,		
			##<< hence also root litter input - which is the difference betweent these two.
		dO = 0		##<< first version change rate needs to be specified, second version rather varies h and calculates dO				
		biasLitterLeaf = 0  ##<< updated version of of.steady and also of.nonSteady account for bias in leaf and root litter separately		
		biasLitterRoot = 0  ##<<
	})	
		
#-------------- a priori knowlege about the parameters ----------
# define distributions by mle and upper confidence bound,
upperBoundProb = 0.99
parmsBounds = list(
		tLagLeaf = parms$tLagLeaf * c(1,2)		
		,tLagRoot = parms$tLagRoot * c(1,3)
		,kY = parms$kY *c(1,10) 
		,kYNarrow = parms$kY *c(1,1.2) 		# constrain kY very tightly simulating litter bag experiments
		,kO = parms$kO *c(1,10)
		,h = c(parms$h, 0.977)
		,cY= c(parms$cY, 0.965)
		,biasDiffRespLitterfall= c(parms$biasDiffRespLitterfall, qnorm(upperBoundProb,sd=41.4))
		#,dO=as.vector(c(parms$dO, Howland14C$obsNutrientSite$somStock[1,"obs"]*1/3 /50	)) 	# upper bound of increase rate is so that 1/3 of the stock can accumulate over 50yrs 
		,dO=as.vector(c(parms$dO, 50) )		##<< 20gC/m2/yr 
		,biasLitterLeaf= c(parms$biasLitterLeaf, qnorm(upperBoundProb,sd=41.4))
		,biasLitterRoot= c(parms$biasLitterRoot, qnorm(upperBoundProb,sd=80))
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
varDistr[c("biasDiffRespLitterfall","dO","biasLitterLeaf","biasLitterRoot")] <- "norm" #normal
# pStorageOnActivation is modeled lognorm, because increase in prior at 1 is not feasable

#----------------- calculate the standard mu and sd at normal scale from quantiles
#mtrace(twQuantiles2Coef)
#mtrace(
parDistr <- twQuantiles2Coef( parmsBounds, varDistr, upperBoundProb=upperBoundProb )

# do the flattest logitnormal distribution possible for cY and h
parDistr$mu[c("cY","h")] <-0
parDistr$sigmaDiag[c("cY","h")] <- 1.5 #twSigmaLogitnorm(unlist(parms0[c("cY","h")]))[,2]


HowlandParameterPriors <- HowlandParameterPriors0 <- list(
	parDistr=parDistr
	,parms0=parms0
)


.tmp.f <- function(){
	setInitialICBMPriors <- function(
		### modify priors to include or update distribution of dO (rate of accumulation) and Ctot0 (initial carbon stocks)
		priors=HowlandParameterPriors	##<< data structure to update
		, somStock=	Howland14C$obsNutrientSite$somStock[1,"obs"]	##<< observed stock
		, pIncreased=c(0,1/3)		##<< proportion of increase of stock within timePeriod: c(mode, upperBound)
		, timePeriod=50				##<< time perioed between t0 and observation of somStock
		, upperBoundProb=0.99
	){
		parmsBounds <- within( list(), {
				dO=somStock*pIncreased /timePeriod	
				Ctot0=somStock*(1+pIncreased)	
			})	
		# assuming no transformations (normal)
		varDistr <- twVarDistrVec( names(parmsBounds) )
		parDistr <- twQuantiles2Coef( parmsBounds, varDistr, upperBoundProb=upperBoundProb )
		for( pName in  names(parDistr$mu)){
			priors$parDistr$trans[[pName]] <- parDistr$trans[[pName]]
			priors$parDistr$mu[[pName]] <- priors$parms0[[pName]] <- parDistr$mu[[pName]]
			priors$parDistr$sigmaDiag[[pName]] <- parDistr$sigmaDiag[[pName]]
		}
		priors
	}
	
	
	HowlandParameterPriors <- setInitialICBMPriors(HowlandParameterPriors0, somStock=Howland14C$obsNutrientSite$somStock[1,"obs"])
}

save(HowlandParameterPriors, file="data/HowlandParameterPriors.RData")


.tmp.f <- function(){
	HowlandParameterPriors <- NULL
	data(HowlandParameterPriors)
	str(HowlandParameterPriors)
	
	parDistr <- poptDistr <- HowlandParameterPriors$parDistr; 
	#poptDistr <- twConstrainPoptDistr(c("biasRespSum","F0","epsF","epsG"),parDistr)
	#ggplotDensity.poptDistr(poptDistr,parmsBounds=parmsBounds,doTransOrig=FALSE)
	#options(error=recover)
	#options(error=NULL)
	ggplotDensity.poptDistr(poptDistr,pMin=0.025,plotUpperQuantile=FALSE)
}

#parms.tmp <- parms; parms.tmp$biasP14C=0.4; obsadj <- calcBiasedObs.hamer(obs, parms.tmp)



