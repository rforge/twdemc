.tmp.f <- function(){
	install.packages("quantreg")
	library(quantreg)
	
	data(engel)
	z <- rq(foodexp ~ income, tau = -1)
	
	x.poor <- quantile(income, 0.1)
	x.rich <- quantile(income, 0.9)
	ps <- z$sol[1, ]
	qs.poor <- c(c(1, x.poor) %*% z$sol[4:5, ])
	qs.rich <- c(c(1, x.rich) %*% z$sol[4:5, ])
	par(mfrow = c(1, 2))
	plot(c(ps, ps), c(qs.poor, qs.rich), type = "n",
		xlab = expression(tau), ylab = "quantile")
	plot(stepfun(ps, c(qs.poor[1], qs.poor)), do.points = FALSE,
		add = TRUE)
	plot(stepfun(ps, c(qs.poor[1], qs.rich)), do.points = FALSE,
		add = TRUE, col.hor = "gray", col.vert = "gray")
	ps.wts <- (c(0, diff(ps)) + c(diff(ps), 0))/2
	ap <- akj(qs.poor, z = qs.poor, p = ps.wts)
	ar <- akj(qs.rich, z = qs.rich, p = ps.wts)
	plot(c(qs.poor, qs.rich), c(ap$dens, ar$dens),
		type = "n", xlab = "Food Expenditure", ylab = "Density")
	lines(qs.rich, ar$dens, col = "gray")
	lines(qs.poor, ap$dens, col = "black")
	legend("topright", c("poor", "rich"), lty = c(1,1), col = c("black", "gray"))
}





#--------------- conditional fitting
n <- 200
df <- 8
delta <- 8
set.seed(4003)
x <- sort(rt(n, df))
u <- runif(n)
v <- -log(1 - (1 - exp(-delta))/(1 + exp(-delta *
					pt(x, df)) * ((1/u) - 1)))/delta
y <- qt(v, df)


plot(x, y, col = "blue", cex = 0.25)
us <- c(0.25, 0.5, 0.75)
for (i in 1:length(us)) {
	u <- us[i]
	v <- -log(1 - (1 - exp(-delta))/(1 + exp(-delta *
						pt(x, df)) * ((1/u) - 1)))/delta
	lines(x, qt(v, df))
}
Dat <- NULL
Dat$x <- x
Dat$y <- y
deltas <- matrix(0, 3, length(us))
FrankModel <- function(x, delta, mu, sigma, df,
	tau) {
	z <- qt(-log(1 - (1 - exp(-delta))/(1 + exp(-delta *
							pt(x, df)) * ((1/tau) - 1)))/delta, df)
	mu + sigma * z
}
for (i in 1:length(us)) {
	tau = us[i]
	#mtrace(nlrq)
	fit <- nlrq(y ~ FrankModel(x, delta, mu, sigma,
			df = 8, tau = tau), data = Dat, tau = tau,
		start = list(delta = 5, mu = 0, sigma = 1),
		trace = TRUE)
	lines(x, predict(fit, newdata = x), lty = 2,
		col = "green")
	deltas[i, ] <- coef(fit)
}


#--------- exploring h
"lprq" <- function(x, y, h, m = 50, tau = 0.5) {
	xx <- seq(min(x), max(x), length = m)
	fv <- xx
	dv <- xx
	for (i in 1:length(xx)) {
		z <- x - xx[i]
		wx <- dnorm(z/h)
		r <- rq(y ~ z, weights = wx, tau = tau,
			ci = FALSE)
		fv[i] <- r$coef[1]
		dv[i] <- r$coef[2]
	}
	list(xx = xx, fv = fv, dv = dv)
}

library(MASS)
data(mcycle)
attach(mcycle)
plot(times, accel, xlab = "milliseconds", ylab = "acceleration")
hs <- c(1, 2, 3, 4)
for (i in hs) {
	h = hs[i]
	#mtrace(lprq)
	fit <- lprq(times, accel, h = h, tau = 0.5)
	lines(fit$xx, fit$fv, lty = i)
}
legend(45, -70, c("h=1", "h=2", "h=3", "h=4"),
	lty = 1:length(hs))

#------- use splines
library(splines)
plot(times, accel, xlab = "milliseconds", ylab = "acceleration",
	type = "n")
points(times, accel, cex = 0.75)
X <- model.matrix(accel ~ bs(times, df = 15))
taus <- 1:3/4
for (tau in taus) {
	fit <- rq(accel ~ bs(times, df = 15), tau = tau,
		data = mcycle)
	accel.fit <- X %*% fit$coef
	lines(times, accel.fit, col=match(tau, taus))
}


#--------- nonparametric fit
data(Mammals)
attach(Mammals)

x <- log(weight)
y <- log(speed)
plot(x, y, xlab = "Weight in log(Kg)", ylab = "Speed in log(Km/hour)",
	type = "n")
points(x[hoppers], y[hoppers], pch = "h", col = "red")
points(x[specials], y[specials], pch = "s", col = "blue")
others <- (!hoppers & !specials)
points(x[others], y[others], col = "black", cex = 0.75)
fit <- rqss(y ~ qss(x, lambda = 1), tau = 0.9)
plot(fit, add=TRUE)


## A bivariate example
#install.packages("tripack")
data(CobarOre)
attach(CobarOre)
fCO <- rqss(z ~ qss(cbind(x,y), lambda= .08), data=CobarOre)
plot(fCO)
plot(fCO, render = "rgl")

#demo(cobar)		#generating a movie with rgl


#-------------- additive models




