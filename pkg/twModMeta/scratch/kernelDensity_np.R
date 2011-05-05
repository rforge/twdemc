install.packages("np")
library("np")

?npudens
example(npudens)

#------- two dimensional example
# bimodal
new.xy <- rbind( MASS::mvrnorm(65, c(3,10), matrix( c(1,.6,.6,1),2) ),
	MASS::mvrnorm(35, c(6, 7), matrix( c(1,.6,.6,1), 2) ) )
#colnames(new.xy) <- c("x","y")
bw <- npudensbw(dat=new.xy)
npplot(bws=bw)

#--------- seven dimensions
library(mvtnorm)
help(package="mvtnorm", help_type = "html")
?dmvnorm
n=7
new.dat <- rbind( rmvnorm(400, rep(0,n), diag(0,n) ),
	MASS::mvrnorm(200, rep(2,n), diag(0.5,n)) )
#takes about one hour
bw2 <- npudensbw(dat=new.dat)
npplot(bws=bw2)






