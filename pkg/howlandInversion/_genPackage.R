#------------ generating package tw.DEMC
# inlcuding generation of Rd files from inline-docs
# install.packages("inlinedocs")

.tmp.f <- function(){
	#library(twMisc)
	#library(snowfall)
	library(twDEMC)
	library(deSolve)
	sfInit(parallel=TRUE,cpus=4)
	#library(inlinedocs) #with twMisc
	
	tmp <- sapply(Sys.glob(file.path("R","*.R")), source)
	data( list=twStripFileExt(basename(Sys.glob(file.path("data","*.RData")))))
	source(file.path("inst","cluster","setupCluster.R"))
	#mtrace(setupClusterHowlandDev)
	setupClusterHowlandDev(pkgDir = ".")

	twUtestF()
	sfInit(parallel=FALSE)
	twUtestF()
	
	twUtestF("applyLB",test="test.sfFArgsApplyLB")
	twUtestF("applyLB")
}


.tmp.inlinedocs <- function(){
	# generate documentation	
	library(inlinedocs)
	unlink( file.path("man","*.Rd") )	
	package.skeleton.dx(".")
	file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man" )	# copy descriptions of data
}

.tmp.compile <- function(){
	# compile and test dll on windows
	dynFilenameLocal <- file.path("src",paste("howlandInversion", .Platform$dynlib.ext, sep = ""))
	sfClusterEval( dyn.unload(dynFilenameLocal) )
	dyn.unload(dynFilenameLocal)
	system("R CMD SHLIB -o src/howlandInversion.dll src/icbm1.c src/rc_helpers.c")
	dyn.load(dynFilenameLocal)
	sfExport("dynFilenameLocal")
	sfClusterEval( dyn.load(dynFilenameLocal) )
	
	#needs to be done on host linux machine (pc026):
	system("rm src/*.o")
	system("R CMD SHLIB -o src/howlandInversion.so  src/icbm1.c src/soilmod_fgsd.c src/rc_helpers.c src/mic_c_part.c src/dllinit.c")
	
	#test compilation of icbm1.c
	#system("R CMD SHLIB src/icbm1.c src/rc_helpers.c")
}




