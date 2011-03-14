data(respWutzler10)
library(twKinresp)

# look for a dataset with high mumax range
kinrespRes <- kinrespGrowthphase(respWutzler10, weights=varPower(fixed=0.5))
rd.e <- getUnlimitedGrowthData(kinrespRes)
idExp <- "9"
serIds <- getSERId(rd.e)
serIdsUnique <- unique(serIds)
serId <- serIds[1]
serId <- "Fal_25_2"
mumaxRange <- lapply( serIdsUnique, function(serId){
		print(serId)
		fit1 <- try({
			rder <- subset(rd.e, serIds==serId)
			repFits <- kinrespRes$resRep[[serId]]
		 	fit1 <- fitKinrespReplicate(rder,start=coefKinresp(coef(repFits$fit)))
			diff(confintKinresp(confint(fit1))["mumax",])
		})
		if( inherits(fit1,"try-error") ) NA  else fit1 
	})
tmp <- unlist(mumaxRange); names(tmp) <- serIdsUnique
sort(tmp)

dsKinrespTut <-subset(rd.e, experiment %in% c("9","17","24"))

#rde <- subset(respWutzler10, experiment ==9)
#res4 <- kinrespGrowthphaseExperiment(rde, weights=varPower(fixed=0.5))
#dsKinrespTut <-getUnlimitedGrowthData(res4)
str(dsKinrespTut)
save(dsKinrespTut,file=file.path("data","dsKinrespTut.RData"))

names(res4$resRep)
coefKinresp(coef(res4$resRep[[1]]$fit))