sfClusterApplyLB <-
function (x, fun, ...) 
{
	sfCheck()
	checkFunction(fun)
	if (sfParallel()) 
		return(clusterApplyLB(sfGetCluster(), x, fun, ...))
	else return(lapply(x, fun, ...))
}
#<environment: namespace:snowfall>
#returns clusterApplyLB

clusterApplyLB <- 
function (cl, x, fun, ...) 
{
	argfun <- function(i) c(list(x[[i]]), list(...))
	dynamicClusterApply(cl, fun, length(x), argfun)
}

dynamicClusterApply <- 
function (cl, fun, n, argfun) 
{
	checkCluster(cl)
	p <- length(cl)
	if (n > 0 && p > 0) {
		submit <- function(node, job) sendCall(cl[[node]], fun, 
				argfun(job), tag = job)
		for (i in 1:min(n, p)) submit(i, i)
		val <- vector("list", n)
		for (i in 1:n) {
			d <- recvOneResult(cl)
			j <- i + min(n, p)
			if (j <= n) 
				submit(d$node, j)
			val[d$tag] <- list(d$value)
		}
		checkForRemoteErrors(val)
	}
}
