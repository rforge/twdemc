#------------ generating package tw.DEMC
# inlcuding generation of Rd files from inline-docs
# install.packages("inlinedocs")

.tmp.setup <- function(){
	library(snowfall)
	sfInit(parallel=TRUE, cpus=4)
	sfLibrary(twMisc)
	sfLibrary(twDEMC)
	sfLibrary(deSolve)
	library(debug)
	library(logitnorm)
	library(ggplot2)
	
	dynFilename <- file.path("src",paste("mcmcMarginals", .Platform$dynlib.ext, sep = ""))
	sfExport("dynFilename")
	dyn.load("src/mcmcMarginals.dll")
	sfClusterEval( dyn.load(dynFilename) )
	#dyn.unload(dynFilename)
	
	tmp <- sapply(Sys.glob(file.path("R","*.R")), source)
	tmp <- sapply(Sys.glob(file.path("R","SoilMod*.R")), sfSource)
	sfSource(file.path("R","ModMeta.R"))	#init soil mod
	sfSource(file.path("R","of_Hamer_Bayes_FS.R"))	#further arguments to ofb.hamer
	dsNames <- twStripFileExt(basename(Sys.glob(file.path("data","*.RData"))))
	data( list=dsNames)
	
	
	#data(HamerLigninGlucoseRespUnc)
	#data(HamerKinRespPars)
	#data(HamerParameterPriors)
	
	
	#twUtestF(twDEMCBatch,"test.goodStart")
	#twUtestF(twDEMC)
	#twUtestF(constrainNStack)
	#twUtestF(ofb.hamer)	# takes long (refitting experiment)
	twUtestF()
}

.tmp.f <- function(){
	twUtestF("SoilMod_FS")
}

.tmp.compile <- function(){
	# compile and test dll on windows
	# only required when providing C or Fortran level code
	pkg <- "mcmcMarginals"
	dynFilenameLocal <- file.path("src",paste(pkg, .Platform$dynlib.ext, sep = ""))
	sfExport("dynFilenameLocal")
	sfClusterEval( dyn.unload(dynFilenameLocal) )
	dyn.unload(dynFilenameLocal)
	system(paste("R CMD SHLIB -o",dynFilenameLocal,"src/soilmod_fgsd.c src/soilmod_fs.c src/mic_c_part.c src/rc_helpers.c"))

	dyn.load(dynFilenameLocal)
	sfExport("dynFilenameLocal")
	sfClusterEval( dyn.load(dynFilenameLocal) )
	
	#needs to be done on host linux machine (pc026):
	system("rm src/*.o")
	system("R CMD SHLIB -o src/howlandInversion.so  src/icbm1.c src/soilmod_fgsd.c src/rc_helpers.c src/mic_c_part.c src/dllinit.c")

	#test compilation of soilmod_fgsd.c
	system("R CMD SHLIB src/soilmod_fgsd.c src/rc_helpers.c src/mic_c_part.c src/dllinit.c")
}

.tmp.genDoc <- function(){
	pkg <- "mcmcMarginals"
	# generate documentation	
	library(inlinedocs)
	unlink( "man/*.Rd", recursive=TRUE)	# take care, entire man directory deleted
	package.skeleton.dx(".")
	file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man" )	# copy descriptions of data

	
	# generate the HTML  files
	prevWd <- setwd("..")
	system(	paste("R CMD INSTALL --html ",pkg, sep="") )
	setwd(prevWd)
	
	# show in Browser
	htmlRoot <- file.path( system.file(package = pkg), "html" )
	html_viewer <- function(path) {
		browser <- getOption("browser")
		if(is.null(browser) && .Platform$OS.type == "windows")
			shell.exec(chartr("/", "\\", path))
		else browseURL(paste("file://", URLencode(path), sep=""))
	}
	html_viewer(file.path(htmlRoot,"00Index.html"))
	
}

system("rm src/*.o")
#R CMD check --no-vignettes --no-latex --no-install mcmcMarginals
#R CMD check --no-vignettes --no-latex --no-codoc mcmcMarginals
#R CMD INSTALL --html mcmcMarginals


.tmp.f <- function(){
	library(sos)
	fres1 <- findFn("mcmc") #also find MCMC
	fres2 <- findFn("Monte Carlo Markov Chain")
	fres3 <- findFn("demc")
	fres <- fres1 | fres2 | fres3;
	sapply( fres, length)
} 

.tmp.f <- function(){
	# get the newest own package versions
	install.packages(c("twMisc","logitnorm","twSnowfall"),repos="http://R-Forge.R-project.org")
}
