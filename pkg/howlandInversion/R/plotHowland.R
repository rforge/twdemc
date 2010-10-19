plotHowlandFM <- function(
	### Plot the fraction modern result of a model-data fit for the Howland data.
	res		##<< result of of.HowlandSteady
	, obs	##<< observations
){
	matplot(res[,"time"], res[,c("F14C_Y","F14C_O","F14CT","respF14CT")], type="l", xlab="Nutrient Site, Time (yr)", ylab="" )
	mtext( "14C (Fraction modern)", side=2, line=2.5, las=0)
	#matplot(res[,"time"], res[,c("F14C_Y","F14C_O","respF14CT")], type="l" )
	lines( res[,"time"], attr(resOf,"iROLayer"), col="orange")
	#lines( res[,"time"], res[,"F14C_Y"], col="orange")
	#lines( res[,"time"], res[,"F14T"], col="orange")
	lines( delta2iR14C(delta14C)/model$modMeta$iR14CStandard ~ yr, data=delta14Catm, col="gray")
	points( obs$somOFM[,"times"],obs$somOFM[,"obs"], col="orange")
	points( obs$respFM[,"times"],obs$respFM[,"obs"], col="blue")
	legend( "topright", inset=0.05, xjust=1, yjust=1, bg="white"
		,legend=c( "RespSOM","O-Layer","Young","Old","SOM","Atmosphere")
		,lty=(c(4,1,1,2,3,1))
		,pch=(c(1,1,NA,NA,NA,NA))
		,col=(c("blue","orange","black","red","green","gray"))
	)
}

