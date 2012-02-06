#------------ generating package twDEMC
# inlcuding generation of Rd files from inline-docs
# install.packages("inlinedocs")


.tmp.f <- function(){
	pkg<-"twDEMC"
	#library(twMiscRgl)
	library(twSnowfall)
	library(coda)
	library(mvtnorm)
	library(logitnorm)
	#library(ggplot2)
	#library(np)
	library(MASS)	#cov.rob in normConstLaplace
	
	library(debug)
	library(lattice)
	
	sfInit(parallel=TRUE, cpus=4)
	#sfInit(parallel=TRUE, cpus=2)
	
	for(filename in Sys.glob(file.path("R","*.R")) ){
		print(filename)
		source(filename)
	}
	#tmp <- sapply(Sys.glob(file.path("R","*.R")), source)
	tmp <- sapply(Sys.glob(file.path("R","multiTemp.R")), sfSource)
	data( list=twStripFileExt(basename(Sys.glob(file.path("data","*.RData")))))
	
	twUtestF()
	
	twUtestF(twCalcLogDenPar)
	twUtestF(transOrigPopt)
	twUtestF("constrainNStack")
	twUtestF("initZ")
	#twUtestF("S3twDEMCPops")
	twUtestF("S3twDEMC")
	twUtestF("twDEMC")
	twUtestF("twDEMCPops")
	twUtestF("getSubSpaces")
	twUtestF("findSplit")
	twUtestF("divideTwDEMC")

	twUtestF("twDEMC","test.badStartSeqData1D")
	twUtestF("findSplit","test.findSplit")
	twUtestF("divideTwDEMC","divideTwDEMCStep")
	twUtestF("divideTwDEMC","divideTwDEMC")
	
	
	twUtestF(twDEMCBatch,"test.goodStart")
	twUtestF(twDEMCBatch,"test.badStart")
	twUtestF(twDEMCBatch,"test.saveAndRestart")
	twUtestF(twDEMCBatch)
	twUtestF("runTwDEMC")
	twUtestF(constrainNStack)
	
	twUtestF(findSplit)
	twUtestF(getSubSpaces)
	twUtestF(findSplit)
}

.tmp.inlinedocs <- function(){
	# generate documentation
	
	# generate RD Files
	library(inlinedocs)
	unlink( file.path("man","*.Rd") )	
	package.skeleton.dx(".")
	unlink(file.path("man","twDEMC-package.Rd"))   # will be copied from genData 
	try(file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man" ), silent=TRUE)	# copy descriptions of data
	unlink(file.path("man","twDEMC.Rd"))   # else overwrites alias twDEMC to twDEMCBlockInt 
	
	# generate the HTML  files
	pkg<-"twDEMC"
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
	
	# copy to the generated html into working directory
	#file.copy( htmlRoot, ".", recursive=TRUE)
	
	require(twMisc)
	updateVersionAndDate()
}


.tmp.f <- function(){
	# generate documentation	
	library(inlinedocs)
	unlink( file.path("man","*.Rd"))	
	package.skeleton.dx(".")
	file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man" )	# copy descriptions of data
}

#R CMD check --no-vignettes --no-manual --no-install twDEMC
#R CMD check --no-vignettes --no-manual --no-codoc twDEMC
#R CMD INSTALL --html twDEMC
