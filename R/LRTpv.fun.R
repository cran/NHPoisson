LRTpv.fun <-
function(fitPPobj)
{
	covariates<-fitPPobj$covariates
	ncov<-dim(covariates)[2]
	LRTpv<-NULL
	ncov<-dim(covariates)[2]
	if (fitPPobj$tind==TRUE)
	{
		covariates<-covariates[,2:ncov]
		ncov<-ncov-1#1 is substracted for the intercept
		for (j in c(1:ncov))
		{

			    llikR<-fitPP.fun(covariates=covariates[,-j], beta=fitPPobj$beta[-(j+1)],posE=fitPPobj$posE, inddat=fitPPobj$inddat,
					modCI=FALSE,modSim=TRUE, dplot=FALSE)$llik
			    difdev <- 2*(fitPPobj$llik - llikR)
     			    LRTpv[j]<-1-pchisq(difdev, 1)
		}
		
	}
	else
	{
		if (ncov>1)
		{
		for (j in c(1:ncov))
		{

			    llikR<-fitPP.fun(covariates=covariates[,-j], beta=fitPPobj$beta[-(j)],posE=fitPPobj$posE, inddat=fitPPobj$inddat,
					modCI=FALSE,modSim=TRUE, dplot=FALSE,tind=FALSE)$llik
			    difdev <- 2*(fitPPobj$llik - llikR)
     			    LRTpv[j]<-1-pchisq(difdev, 1)
		}
		}
		else stop('No test can be carried out since there is only one covariate and no intercept')

	}
	namcovariates<-fitPPobj$namcovariates	
	if (is.null(namcovariates)) namcovariates<-paste('Covariate',c(1:ncov))
	LRTpv<-matrix(LRTpv,ncol=1,dimnames=list(namcovariates,'p-values'))

	print(round(LRTpv,3))
	return(LRTpv)

}
