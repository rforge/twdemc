setupClusterHowlandDev <- function(
	### cluster setup for Howland with local development installation
	pkgDir = "../eclipse_R/howlandInversion"		#relative path to asom execution directory
){
	
	sfLibrary(twMisc)
	sfLibrary(twDEMC)
	sfLibrary(deSolve)
	#library(debug)
	#library(logitnorm)
	
	dynFilename <- file.path(pkgDir,"src",paste("howlandInversion", .Platform$dynlib.ext, sep = ""))
	sfExport("dynFilename")
	sfClusterEval( dyn.load(dynFilename) )
	dyn.load(dynFilename)
	#dyn.unload(dynFilename)

	#dsFileNames <- Sys.glob(file.path(pkgDir,"data","*.RData"))
	#for( dsFilename in dsFileNames) load(dsFilename)
	#dsNames <- twStripFileExt(basename(dsFileNames))
	#data( list=dsNames,  )
	#sfExport(list=dsNames)

	sfClusterCall( data, c("Howland14C"))	
	
	#tmp <- sapply(Sys.glob(file.path(pkgDir,"R","*.R")), source)
	tmp <- sapply(Sys.glob(file.path(pkgDir,"R","ICBM*.R")), sfSource)
	sfSource(file.path(pkgDir,"R","ModMeta.R"))	#init soil mod
	sfSource(file.path(pkgDir,"R","c14Utils.R"))	#init soil mod
	#sfSource(file.path(pkgDir,"R","of_Hamer_Bayes_FS.R"))	#further arguments to ofb.hamer and subfunctions
	
}

setupClusterMDIHamer <- function(
	### cluster setup for MDIHamer with installed package 
){ 
	sfLibrary(twMDIHamer) 
}
