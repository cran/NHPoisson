globalval.fun <-
function(obFPP,lint=NULL,nint=NULL, covariates,
Xvar=NULL,namXvar=NULL, Xvart=NULL,namXvart=NULL, 
h=NULL, typeRes=NULL, typeResLV='Pearson',typeI='Disjoint', nsim=100,
clevel=0.95,resqqplot=FALSE, nintLP=100, tit='',
flow=0.5, addlow=FALSE,indgraph=FALSE,scax=NULL, scay=NULL,
legcex=0.5,ngen=100, cores=1, xlegend='topleft')
{
if ((is.null(lint))&(typeI=='Overlapping')) stop('Argument lint must be specified for Overlapping intervals')
if ((is.null(lint))&(is.null(nint))&(typeI=='Overlapping')) stop('one of arguments lint or nint must be specified for Disjoint intervals')
namcovariates<-obFPP$namcovariates
if(is.null(Xvar)) 
{
Xvar<-covariates
namXvar<-namcovariates
}
t<-obFPP$t
if (is.null(Xvart)==FALSE) Xvart<-as.matrix(Xvart)
if (is.null(Xvar)==FALSE) Xvar<-as.matrix(Xvar)
posE<-obFPP$posE
if (is.null(tit)) tit<-obFPP$tit
inddat<-obFPP$inddat
posEH<-transfH.fun(obFPP)$posEH
res<-unifres.fun(posEH)
graphresU.fun(unires=res$unires, posE=posE, Xvariables=Xvar,
namXv=namXvar,tit=tit, flow=flow, addlow=addlow,
indgraph=indgraph)

if (typeI=='Disjoint')
RRes<-CalcResD.fun(obFPP=obFPP, h=h, nint=nint,lint=lint, typeRes=typeRes)
else
RRes<-CalcRes.fun(obFPP=obFPP, h=h, typeRes=typeRes,lint=lint)

graphrate.fun(objres=RRes, tit=tit,scax=scax,scay=scay, xlegend=xlegend)

graphres.fun(objres=RRes,Xvariables=Xvart,typeRes=typeRes,
namXv=namXvart,indgraph=indgraph,addlow=addlow,tit=tit,flow=flow)

graphResCov.fun( Xvar=Xvar,obFPP=obFPP, h=h, nint=nintLP, tit=tit, 
typeRes=typeResLV,namX=namXvar,indgraph=indgraph)


if (resqqplot==TRUE) resQQplot.fun(nsim=nsim,objres=RRes, 
covariates=covariates, clevel=clevel,tit=tit,cores=cores, n=ngen)

return(RRes)
}
