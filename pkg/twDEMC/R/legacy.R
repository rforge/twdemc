# adapting result objects from former versions

updateLegacyTwDEMCArgs <- function(
	### update list of arguments to of former version of twDEMC
	argsF	##<< a list of arguments to twDEMC
){
	##details<< 
	## Arguments "fLogLik","argsFLogLik" are renamed to "fLogDen","argsFLogDenX"
	names(argsF)[ names(argsF) %in% c("fLogLik","argsFLogLik")] <- c("fLogDen","argsFLogDenX")
	argsF
}
attr(updateLegacyTwDEMCArgs,"ex") <- function(){
	data( twdemcEx1)
	argsF <- as.list(attributes(twdemcEx1)$batchCall)
	updateLegacyTwDEMCArgs(argsF)
}



updateLegacyTwDEMC <- function(
	### update result of former version of twDEMC
	aTwDEMC
){
	##details<<
	## The function of the log of unnormalized density has been renamed from fLogLik to fLogDen
	## components "rLogLik","resFLogLikX" are renamed to "rLogDen","resFLogDenX"
	names(aTwDEMC)[ names(aTwDEMC) %in% c("rLogLik","resFLogLikX")] <- c("rLogDen","resFLogDenX")
	
	##details<<
	## A twDEMC stores its call. Because argument names have been changed,
	## they are also changed within the \code{batchCall} entry: see \code{\link{updateLegacyTwDEMCArgs}}
	if( !is.null(bcall <- attributes(aTwDEMC)$batchCall) ){
		attributes(aTwDEMC)$batchCall <- as.call( updateLegacyTwDEMCArgs(as.list(bcall)))		
	}
	aTwDEMC
}
attr(updateLegacyTwDEMC,"ex") <- function(){
	data( twdemcEx1)
	aTwDEMC <- twdemcEx1
}





