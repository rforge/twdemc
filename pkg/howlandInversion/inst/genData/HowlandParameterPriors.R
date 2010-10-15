

################################ The parameters ################################

# HowlandParameterPriors

# units are gC/m2 /yr
parms = parms0 = within( list(), {
		tLagLeaf=1	##<< less than 1 yr is recent-C, which is not accounted in the Howland study
		tLagRoot=5	##<< time lag between C-fixation and root turnover
		kY=2.4		##<< decay constant
		kO=0.1		
		h=0.2		##<< humification coefficient
		cY=0.2		##<< initial proportion of total C in young pool
	})	
		
#-------------- a priori knowlege about the parameters ----------
# define distributions by mle and upper confidence bound,
upperBoundProb = 0.99
parmsBounds = list(
		tLagLeaf = parms$tLagLeaf * c(1,2)		
		,tLagRoot = parms$tLagRoot * c(1,3)
		,kY = parms$kY *c(1,1000) 
		,kO = parms$kO *c(1,1000)
		,h = c(parms$h, 0.8)
		,cY= c(parms$cY, 0.8)
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
#varDistr[c("biasRespSum","biasRespC14Sum")] <- "norm" 
# pStorageOnActivation is modeled lognorm, because increase in prior at 1 is not feasable

#----------------- calculate the standard mu and sd at normal scale from quantiles
#mtrace(twQuantiles2Coef)
parDistr <- twQuantiles2Coef( parmsBounds, varDistr, upperBoundProb=upperBoundProb )

HowlandParameterPriors <- list(
	parDistr=parDistr
	,parms0=parms0
)

save(HowlandParameterPriors, file="data/HowlandParameterPriors.RData")


.tmp.f <- function(){
	HowlandParameterPriors <- NULL
	data(HowlandParameterPriors)
	str(HowlandParameterPriors)
	
	parDistr <- poptDistr <- HowlandParameterPriors$parDistr
	#poptDistr <- twConstrainPoptDistr(c("biasRespSum","F0","epsF","epsG"),parDistr)
	#ggplotDensity.poptDistr(poptDistr,parmsBounds=parmsBounds,doTransOrig=FALSE)
	ggplotDensity.poptDistr(poptDistr,parmsBounds=parmsBounds,pMin=0.025,plotUpperQuantile=FALSE)
}

#parms.tmp <- parms; parms.tmp$biasP14C=0.4; obsadj <- calcBiasedObs.hamer(obs, parms.tmp)



