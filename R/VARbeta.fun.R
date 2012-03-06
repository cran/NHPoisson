VARbeta.fun <-
function(covariates, lambdafit)
{

	lambdafit<-lambdafit
	covariates<-covariates

	K<-dim(covariates)[2]
	InforM<-matrix(rep(0,K*K), ncol=K)

	calcVAR.fun<-function(j,i)
	{
		InforM[i,j]<<-sum(lambdafit*covariates[,i]*covariates[,j])
		InforM[j,i]<<-sum(lambdafit*covariates[,i]*covariates[,j])
	}

	calcVAR2.fun<-function(i)
	{
		aux<-sapply(c(1:i), FUN=calcVAR.fun, i)
	}

	aux2<-sapply(c(1:K), FUN=calcVAR2.fun)

	VARbeta<-solve(InforM)
	return(VARbeta)

}
