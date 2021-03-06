#------------ generating package tw.DEMC
# inlcuding generation of Rd files from inline-docs
# install.packages("inlinedocs")

.tmp.f <- function(){
	require(deSolve)
	
	require(debug)
	require(twMisc) # twStripFileExt
	#library(inlinedocs) #with twMisc

	tmp <- sapply(Sys.glob(file.path("R","*.R")), source)
	data( list=twStripFileExt(basename(Sys.glob(file.path("data","*.RData")))))
	#source(file.path("inst","cluster","setupCluster.R"))
	#mtrace(setupClusterHowlandDev)
	#setupClusterHowlandDev(pkgDir = ".")
	
	dynFilenameLocal <- file.path("src",paste("twModMeta", .Platform$dynlib.ext, sep = ""))
	dyn.load(dynFilenameLocal)
	
	twUtestF()
	
	twUtestF("ICBMDemo", "test.steadyStateHowland" )
	twUtestF("ICBMDemo")
}


.tmp.inlinedocs <- function(){
	# generate documentation
	data(delta14Catm)
	pkg <- "twModMeta"
	library(inlinedocs)
	htmlRoot <- file.path( system.file(package = pkg), "html" )	# may work only after the first R CMD INSTALL
	html_viewer <- function(path) {
		browser <- getOption("browser")
		if(is.null(browser) && .Platform$OS.type == "windows")
			shell.exec(chartr("/", "\\", path))
		else browseURL(paste("file://", URLencode(path), sep=""))
	}
	
	# generate RD Files
	unlink( file.path("man","*.Rd") )	# get rid of old documentation
	package.skeleton.dx(".")			# parse source code in R dir
	unlink( file.path("man","*-package.Rd") )	# get rid of old documentation
	try(file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man" ), silent=TRUE)	# copy descriptions of data
	
	# generate the HTML  files
	prevWd <- setwd("..")
	system(	paste("R CMD INSTALL --html ",pkg, sep="") )
	setwd(prevWd)
	
	# show in Browser
	html_viewer(file.path(htmlRoot,"00Index.html"))
	
}

.tmp.check <- function(){
	prevWd <- setwd("..")
	system(	paste("R CMD check --no-manual ",pkg, sep="") )
	setwd(prevWd)
}
	
.tmp.compile <- function(){
	# compile and test dll on windows
	dynFilenameLocal <- file.path("src",paste("twModMeta", .Platform$dynlib.ext, sep = ""))
	#sfExport("dynFilenameLocal")

	#sfClusterEval( dyn.unload(dynFilenameLocal) )
	dyn.unload(dynFilenameLocal)
	system("R CMD SHLIB -o src/twModMeta.dll src/icbmDemo.c src/rc_helpers.c")
	unlink(file.path("src","*.o"))

	#needs to be done on host linux machine (pc026):
	system("rm src/*.o")
	system("R CMD SHLIB -o src/twModMeta.so  src/icbmDemo.c src/soilmod_fgsd.c src/rc_helpers.c src/mic_c_part.c src/dllinit.c")
	
	dyn.load(dynFilenameLocal)
	#sfExport("dynFilenameLocal")
	#sfClusterEval( dyn.load(dynFilenameLocal) )
	
	
	#test compilation of icbm1.c
	#system("R CMD SHLIB src/icbm1.c src/rc_helpers.c")
}

.tmp.installCluster <- function(){
	# install from repository
	install.packages("twDEMC", repos="http://R-Forge.R-project.org")
}




