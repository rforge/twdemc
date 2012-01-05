#------------ code to help to maintain the package code
# usage: paste single code-lines from inside the functions to the RGui or RTerm
# workspace must be the root of the package, i.e. above R and man

#-------- important adjustments specific for the user and working environment
# adjust the entries (especially Package and Version) in file DESCRIPTION before

# mdiBaseDir <- "/Net/Groups/C-Side/MDI"   # on unix
 mdiBaseDir <- "z:"
mdiBaseDir <- "m:"                         # on windows for most users
installDir <- if( file.access(.libPaths()[1],2)==0 ) .libPaths()[1] else  
  "c:/tmp"          # if standard libDir is not writeable, this gives the directory where the package will be installed

desc <- readLines("DESCRIPTION")
pkg <- gsub("(^Package:\\s+)(.+)$", "\\2", desc[grep("^Package:",desc)], perl=TRUE)
pkgVersion <- gsub("(^Version:\\s+)(.+)$", "\\2", desc[grep("^Version:",desc)], perl=TRUE)
mdiCodeDir <- file.path(mdiBaseDir,"_code/R") 					
mdiCodeHtmlDir <- file.path(mdiBaseDir,"public_html/code_doc")	
htmlRoot <- file.path( installDir, pkg, "html" )  # may work only after the first R CMD INSTALL



#.tmp.loadPackages <- function(){
	# loading libraries, sourcing the code and loading data
	# usually done on startup library(MyPackage)
	
	# this code uses several packages that need to be installed once
	# install.packages(c("RUnit","inlinedocs","R.methodsS3", "abind","twMisc"), repos = c("http://R-Forge.R-project.org","http://cran.rakanu.com/"), dep = TRUE)
	# in case that you use not the current R-version, you need to download the sources tarball from https://r-forge.r-project.org/R/?group_id=887
	#   upack the sources, and issue from a shell from the folder above extracted folder twMisc "R CMD INSTALL twMisc"
	
	#library(snowfall)
	#sfInit(parallel=TRUE,cpus=4)		# for parallel execution on 4 processors see (library snowfall)
	#library(debug)
	
	#function twStripFileExt resides in package twMisc, this needs to be installed once
	#install.packages(c("twMisc"), repos=c("http://cran.rakanu.com","http://R-Forge.R-project.org"))
	library(twMisc) 	# for twStipFileExt
	library(twDEMC)
	library(ggplot2)
	tmp <- sapply(Sys.glob(file.path("R","*.R")), source)		# load all the R code in file in dir R
	data( list=twStripFileExt(basename(Sys.glob(file.path("data","*.RData")))))	# load all the .RData files in dir data
#}





#######################################   build package     ####################################
.tmp.inlinedocs <- function(){
	# generate RD and HTML documentation
	require(inlinedocs)
	html_viewer <- function(path) {
		browser <- getOption("browser")
		if( .Platform$OS.type == "windows")
			shell.exec(chartr("/", "\\", path))
		else browseURL(paste("file://", path, sep=""))
	}
	
	# generate RD Files
	unlink( file.path("man","*.Rd") )	# get rid of old documentation
	package.skeleton.dx(".")			# parse source code in R dir
	try(file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man" ), silent=TRUE)	# copy descriptions of data

	# generate the HTML  files
	prevWd <- setwd("..")
	system(	paste('R CMD INSTALL --html -l "',installDir,'" ',pkg, sep="") )
	setwd(prevWd)
	
	# show in Browser
	html_viewer(file.path(htmlRoot,"00Index.html"))
  
	require(twMisc)
	updateVersionAndDate()
}

#######################################   build package     ####################################
.tmp.build.package <- function(){
  # create the zip file of the package, that can be installed from R-GUI and put it to mdi code directory
  # run from windows it will be a zip file
  # run from linux it will be a tgz file relating to the architecture
	prevWd <- setwd("..")
	system(	paste("R CMD INSTALL --build ",pkg, sep="") )
  setwd(prevWd)
}

.tmp.put2mdiCodeDir <- function(){
  # When done with local development, archive the code on the MDI code directory
  # All files need to be saved.
  # take care - the existing version will be overwritten

  destDir <- file.path(mdiCodeDir,pkg)  	# might need to change mounted path e.g. m:
  unlink(destDir, recursive=TRUE) 
  file.copy( file.path("..", pkg), mdiCodeDir, overwrite=TRUE, recursive=TRUE)
  
  #   copy also the buildt zip and/or tgz files generated above
  file.copy( Sys.glob( file.path("..",paste(pkg,'_',pkgVersion,'*',sep=""))), mdiCodeDir )

  # on unix it is possible to set permissions so that nobody other overwrites changes
  prevWd <- setwd(mdiCodeDir)
  system(  paste("chmod -R 755 ",pkg, sep="") )
  setwd(prevWd)
}

.tmp.putDesc2Web <- function(){
  # make the package accessible by wwww
  
	# mdi web-page documentation directory
	# access as http://www.bgc-jena.mpg.de/bgc-mdi/code_doc/<pkg>/html/00Index.html
	destDir <- file.path(mdiCodeHtmlDir,pkg)		# might need to change mounted path e.g. m:
  #unlink(destDir, recursive=TRUE) 
	unlink(Sys.glob(file.path(destDir,"*")), recursive=TRUE)
	try(dir.create(destDir, mode = "0775"))		# mode ignored on windows, need to set by hand, best on unix eg. dialog chmod 775 <dir>
  
  # copy the hmtl description generated above
  file.copy( htmlRoot, destDir, recursive=TRUE)
  
  # copy the pdf description generated above
  file.copy( file.path("..",paste(pkg,".pdf",sep="")), destDir )
	
  # copy the buildt zip and/or tgz files generated above
  file.copy( Sys.glob( file.path("..",paste(pkg,'_',pkgVersion,'*',sep=""))), destDir )
}

.tmp.unitTests <- function(){
	# executing Unit Tests before actual installation
	
	library(twMisc) 	# for twUTestF
	twUtestF()											# all unit tests		(displaying only summary output)
	twUtestF(respTempLlyodAndTaylor94)					# single unit test
	twUtestF(respTempLlyodAndTaylor94,"test.temp10")		# single test function  (displaying all the output)
	
	# let R check package consistency
	prevWd <- setwd("..")
	#system(	cmd <- paste("R CMD check --no-latex ",pkg, sep="") )
	system(	cmd <- paste("R CMD check --no-manual ",pkg, sep="") )
	setwd(prevWd)
	copy2clip(cmd)	# for allowing easy pasting into a shell
	
}

# for adding low-level C or Fortran Code to the package ask Thomas (twutz)


.tmp.compile <- function(){
	# compile and test dll on windows
	# only required when providing C or Fortran level code
	# replace howlandInversion with the basename of your dll
	
	dynFilenameLocal <- file.path("src",paste("howlandInversion", .Platform$dynlib.ext, sep = ""))
	sfExport("dynFilenameLocal")
	sfClusterEval( dyn.unload(dynFilenameLocal) )
	dyn.unload(dynFilenameLocal)
	system("R CMD SHLIB -o src/howlandInversion.dll src/icbm1.c src/rc_helpers.c")
	dyn.load(dynFilenameLocal)
	sfExport("dynFilenameLocal")
	sfClusterEval( dyn.load(dynFilenameLocal) )
	
	#needs to be done on host linux machine (pc026):
	system("rm src/*.o")
	system("R CMD SHLIB -o src/howlandInversion.so  src/icbm1.c src/soilmod_fgsd.c src/rc_helpers.c src/mic_c_part.c src/dllinit.c")
}





