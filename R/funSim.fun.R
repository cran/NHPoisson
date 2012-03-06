funSim.fun <-
function(i,lambda,fun.name,fun.args=NULL,n=100)
{
	posNH<-simNHP.fun(lambda=lambda, n=n)$posNH
  	fun.args <- c(list(posNH), fun.args)
	aux<-do.call(fun.name,fun.args)
	return(aux)
}
