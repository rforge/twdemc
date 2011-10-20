#------------ generating package twDEMC
# inlcuding generation of Rd files from inline-docs
# install.packages("inlinedocs")

pkg<-"twDEMC"

.tmp.f <- function(){
	library(twMiscRgl)
	library(twSnowfall)
	library(coda)
	library(mvtnorm)
	library(logitnorm)
	library(ggplot2)
	library(np)
	library(MASS)	#cov.rob
	
	library(debug)
	library(lattice)
	
	#sfInit(parallel=TRUE, cpus=4)
	sfInit(parallel=TRUE, cpus=2)
	
	tmp <- sapply(Sys.glob(file.path("R","*.R")), source)
	tmp <- sapply(Sys.glob(file.path("R","multiTemp.R")), sfSource)
	data( list=twStripFileExt(basename(Sys.glob(file.path("data","*.RData")))))
	
	twUtestF()
	
	twUtestF(transOrigPopt)
	twUtestF(twDEMCBatch)
	twUtestF("runTwDEMC")
	twUtestF(twDEMCBatch,"test.goodStart")
	twUtestF(twDEMCBatch,"test.badStart")
	twUtestF(twDEMCBatch,"test.saveAndRestart")
	twUtestF(twDEMC)
	twUtestF(constrainNStack)
	twUtestF("initZ")
	twUtestF("S3twDEMC")
}

.tmp.inlinedocs <- function(){
	# generate documentation
	
	# generate RD Files
	pkg<-"twDEMC"
	library(inlinedocs)
	unlink( file.path("man","*.Rd") )	
	package.skeleton.dx(".")
	try(file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man" ), silent=TRUE)	# copy descriptions of data
	unlink(file.path("man","twDEMC.Rd"))   # else overwrites alias twDEMC to twDEMCInt 
	
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
	
	# copy to the generated html into working directory
	file.copy( htmlRoot, ".", recursive=TRUE)
}


.tmp.f <- function(){
	# generate documentation	
	library(inlinedocs)
	unlink( file.path("man","*.Rd"))	
	package.skeleton.dx(".")
	file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man" )	# copy descriptions of data
}

#R CMD check --no-vignettes --no-latex --no-install twDEMC
#R CMD check --no-vignettes --no-latex --no-codoc twDEMC
#R CMD INSTALL --html twDEMC
