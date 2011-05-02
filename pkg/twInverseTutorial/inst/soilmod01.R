#install.packages("FME")
require(FME)
#?FME

windows(width=3.4,height=2.8, pointsize=10, record=TRUE)
par( las=1 )					#also y axis labels horizontal
par(mar=c(2.2,2.8,1,0)+0.3 )  #margins
par(tck=0.02 )				#axe-tick length inside plots             
par(mgp=c(1.4,0.2,0) )  #positioning of axis title, axis labels, axis


## =======================================================================
## Jenny 49
## =======================================================================

Ne <- with( as.list(pars), k2/k1)	# equilibrium level of N
yini    <- c(N = 2*Ne)

pars    <- c(
	k1  = 0.0608	# 1/yr decay rate Nitrogen 
	,k2  = 38.0		# kg N / ha /yr, rate of N fixation
)

times   <- seq(0, 100, by = 1)	# years

### The Jenny model as derivatives
MJenny49 <- function(Time, State, Pars) {
	with(as.list(c(State, Pars)), {
			loss   <- k1*N
			fixation  <- k2
			dN     <- +fixation - loss
			return(list(c(dN)))
		})
}

# Solving the model
out     <- ode(func = MJenny49, y = yini, parms = pars, times = times)

# Plotting the results
matplot(out[,1], out[,2], type = "l", xlab = "time (yr)", ylab = "",	main = "Jenny 1949", lwd = 2)
mtext("N (kg/ha/yr)",2,2.2, las=0)
legend("topright", c("N stock"), col = 1:2, lty = 1:2, inset=c(0.05,0.05))




## =======================================================================
## ICBM
## =======================================================================

yini    <- c(
	Y = 1		##<< young carbon kg/ha
	,O=3		##<< old carbon kg/ha
)	 

pars    <- c(
	litter = 0.6	##<< kgC/ha/yr litter input
	,kY  = 1		##<< 1/yr decay rate young pool tvr of 1 year 
	,kO  = 1/40		##<< 1/yr decay rate of old pool
	,h = 0.15		##<< humification coefficient
)

yEq <- with(as.list(pars), litter/kY)
oEq <- with(as.list(pars), (h*yEq)/kO)


times   <- seq(0, 50, by = 1)	# years

# fluctuating inputs
inputs <- rnorm( length(times), pars["litter"], pars["litter"]*0.3 )
fInput <- approxfun(times, inputs)
plot( inputs ~ times )

MICBM <- function(Time, State, Pars, fInput) {
	with(as.list(c(State, Pars)), {
			decY  <- kY*Y
			decO  <- kO*O
			dY    <- +fInput(Time) -decY
			#dY    <- +litter -decY
			dO    <- +h*decY - decO
			resp  <- (1-h)*decY + decO
			return(list(c(dY,dO),c(decY=decY,decO=decO,resp=resp)))
		})
}

# Solving the model
out     <- ode(func = MICBM, y = yini, parms = pars, times = times, fInput=fInput)
tail(out)

# Plotting the results
matplot(out[,1], out[,2:3], type = "l", xlab = "time (yr)", ylab = "",	main = "ICBM", lwd = 2)
mtext("C (kg/ha/yr)",2,2.2, las=0)
legend("topright", c("Young","Old"), col = 1:2, lty = 1:2, inset=c(0.05,0.15))

## =======================================================================
## ICBM experiment of labelled litter
## =======================================================================
# assume we can mark litter inputs 
# and trace ratio of new Carbon in respiration and stocks

yini    <- c(
	Y = 1		##<< young carbon kg/ha
	,O=3		##<< old carbon kg/ha
	,YNew = 0
	,ONew = 0
)	 
