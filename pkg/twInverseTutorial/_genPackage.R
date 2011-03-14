#------------ code to help to maintain the package code
# usage: paste single code-lines from inside the functions to the RGui or RTerm
# workspace must be the root of the package, i.e. above R and man

pkg <- "twInverseTutorial" 	# adjust to the name of your package, must correspond to the directory name and the DESCRIPTION file
mdiCodeDir <- file.path("z:","_code","R") 					# adjust z: to mounted path of mdi drive
mdiCodeHtmlDir <- file.path("z:","public_html","code_doc")	# adjust z: to mounted path of mdi drive
pkgVersion <- "1.0"		# must correspond to the version in DESCRIPTION file

.tmp.loadPackages <- function(){
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
	tmp <- sapply(Sys.glob(file.path("R","*.R")), source)		# load all the R code in file in dir R
	data( list=twStripFileExt(basename(Sys.glob(file.path("data","*.RData")))))	# load all the .RData files in dir data
}


.tmp.inlinedocs <- function(){
	# generate documentation
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
	try(file.copy( Sys.glob(file.path("inst","genData","*.Rd")), "man" ), silent=TRUE)	# copy descriptions of data

	# generate the HTML  files
	prevWd <- setwd("..")
	system(	paste("R CMD INSTALL --html ",pkg, sep="") )
	setwd(prevWd)
	
	# show in Browser
	html_viewer(file.path(htmlRoot,"00Index.html"))
}

.tmp.genPdfDoc <- function(){
	# generate pdf file
	prevWd <- setwd("..")
	system(	paste("R CMD Rd2pdf ",pkg, sep="") )
	setwd(prevWd)
}
	
.tmp.put2mdiCodeDir <- function(){
	# copy to mdi web-page documentation directory
	# access as http://www.bgc-jena.mpg.de/bgc-mdi/code_doc/<pkg>/html/00Index.html
	destDir <- file.path(mdiCodeHtmlDir,pkg)		# might need to change mounted path e.g. m:
	try(file.remove(destDir))
	dir.create(destDir, mode = "0775")		# mode ignored on windows, need to set by hand, best on unix eg. dialog
	file.copy( htmlRoot, destDir, recursive=TRUE)
	
	# create the zip file of the package, that can be installed from R-GUI and put it to mdi code directory
	prevWd <- setwd("..")
	system(	paste("R CMD INSTALL --build ",pkg, sep="") )
	file.copy( paste(pkg,"_",pkgVersion,".zip",sep=""), mdiCodeDir )	# adjust the version according to the DESCRIPTION file
	setwd(prevWd)
	
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






