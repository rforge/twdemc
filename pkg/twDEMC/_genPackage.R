#------------ generating package twDEMC
# inlcuding generation of Rd files from inline-docs
# install.packages("inlinedocs")

.tmp.f <- function(){
	library(twMisc)
	library(twSnowfall)
	library(debug)
	library(abind)
	library(coda)
	library(mvtnorm)
	library(snowfall)
	library(logitnorm)
	library(ggplot2)
	sfInit(parallel=TRUE, cpus=4)
	
	tmp <- sapply(Sys.glob(file.path("R","*.R")), source)
	tmp <- sapply(Sys.glob(file.path("R","multiTemp.R")), sfSource)
	data( list=twStripFileExt(basename(Sys.glob(file.path("data","*.RData")))))
	
	twUtestF()
	
	twUtestF(transOrigPopt)
	twUtestF(twDEMCBatch)
	twUtestF(twDEMCBatch,"test.goodStart")
	twUtestF(twDEMCBatch,"test.badStart")
	twUtestF(twDEMCBatch,"test.saveAndRestart")
	twUtestF(twDEMC)
	twUtestF(constrainNStack)
}

.tmp.f <- function(){
	# generate documentation	
	library(inlinedocs)
	unlink( "man", recursive=TRUE)	# take care, entire man directory deleted
	package.skeleton.dx(".")
	file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man" )	# copy descriptions of data
	unlink(file.path("man","twDEMC.Rd"))   # else overwrites alias twDEMC to twDEMCInt 
}

#R CMD check --no-vignettes --no-latex --no-install twDEMC
#R CMD check --no-vignettes --no-latex --no-codoc twDEMC
#R CMD INSTALL --html twDEMC
