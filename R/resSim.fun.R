resSim.fun <-
function(i,lambda,covariates, beta,lint,t=NULL,
 tind='TRUE',typeI='Disjoint', typeRes='Pearson',h=NULL,n=100)
{
posNH<-simNHP.fun(lambda=lambda,n=n)$posNH
if (is.null(t)) t<-c(1:length(lambda))


mod<-fitPP.fun(covariates=covariates, beta=beta, posE=posNH,  tind=tind,
tim=t,modCI='FALSE', tit=NULL,modSim=TRUE,dplot='FALSE')


if (typeI=='Disjoint')
{
ResAux<-CalcResD.fun(obFPP=mod, h=h, nint=NULL,lint=lint, 
typeRes=typeRes,modSim='TRUE')

}
elseResAux<-CalcRes.fun(obFPP=mod,  h=h, lint=lint,typeRes=typeRes)

if (typeRes=='Raw') res<-ResAux$RawRes
else res<-ResAux$ScaRes$ScaRes
return(res)
}
