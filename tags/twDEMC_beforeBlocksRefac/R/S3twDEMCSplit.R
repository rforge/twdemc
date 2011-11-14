setMethodS3("devide","twDEMC", function( 
		### Devide a twDEMC result into to two subPopulations
		x
		,split	##<< named numeric scalar: value and name of the parameter used to devide the sample
		,...
	){
		#chain <- x$parms[,,1]
		boLeft <- apply(x$parms,3,function(chain){ chain[names(split),] <= split })
		right <- left <- x
		
		#devide.twDEMC
		##seealso<<   
		## \code{\link{stackChains.twDEMC}}
		list( right=right, left=left )
	})
attr(devide.twDEMC,"ex") <- function(){
	data(twdemcEx1)
	#x <- twdemcEx1
	res <- devide(twdemcEx1,split=c(a=10.8))
	str(res$right)
}

