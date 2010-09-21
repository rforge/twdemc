## Test unit 'calcDEMCTemp'

.setUp <-function () {}

.tearDown <- function () {}

test.tempBounds <- function (){
	T0 <- 10
	Tend <- 2
	temp <- calcDEMCTemp(T0,Tend,1,c(0,1))
	msg <- paste("temp=[",paste(temp,collapse=","),"]",sep="")
	checkEquals( length(temp), 2, msg )
	checkEquals( temp[1], T0, msg )
	checkEquals( temp[2], Tend, msg )
} 

test.tempMonotonous <- function (){
	nGen <- 10
	T0 <- 10
	Tend <- 2
	temp <- calcDEMCTemp(T0,Tend,nGen)
	msg <- paste("temp=[",paste(temp,collapse=","),"]",sep="")
	checkEquals( length(temp), nGen, msg )
	checkEquals( temp[nGen], Tend, msg )
	checkTrue( temp[1]<T0, msg )
	checkTrue( all(diff(temp)<0), msg )
	checkTrue( all(diff(diff(temp))>0), msg )
} 


tmp.f <- function(){
	#library(debug)
	currentPackage="tw.DEMC"
	mtrace(calcDEMCTemp)
	mtrace(twUtest)
	twUtestF(calcDEMCTemp,"test.tempBounds")
	twUtestF(calcDEMCTemp)
	twUtestF()
}

