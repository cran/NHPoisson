fitPP.fun <-
function(covariates=NULL,beta,posE=NULL,inddat=NULL,POTob=NULL,namcovariates=NULL,n=NULL,
tind='TRUE', tim=NULL, maxty='nlminb',method="Nelder-Mead", modCI='TRUE',CIty='Transf', clevel=0.95,
tit='', modSim='FALSE',dplot=TRUE,xlegend='topleft',lambdaxlim=NULL,lambdaylim=NULL)
{

if ( is.null(posE)&is.null(POTob) )
{
stop('Error: At least one of the arguments, posE or POTob must be specified')
}


if (is.null(POTob)==FALSE) 
{
T<-POTob$T
thres<-POTob$thres
n<-length(T)


exc<-as.numeric(T>thres)
inrachtx<-(c(0,diff(exc))==1)
numerachtx<-cumsum(inrachtx)[exc==1]
intentx<-T[exc==1]
posicintx<-c(1:n)[inrachtx==1]
posE<-posicintx+tapply(intentx,INDEX=numerachtx, FUN=which.max)-1

inddat<-1-exc
inddat[posE]<-1
}

if ( is.null(n)&is.null(covariates)&is.null(inddat) )
{
stop('Error: At least one of the arguments, n, covariates,inddat or POTob must be specified')
}
if (is.null(n)& is.null(inddat)) n<-dim(covariates)[1]
if (is.null(n)& is.null(covariates)) n<-length(inddat)

if (is.null(inddat)) inddat<-rep(1,n)
n<-length(inddat)
control<-sum(inddat)

if ((tind==FALSE)&(is.null(covariates)==TRUE))
{
stop('Error: Model without covariates and without constant term')
}

if (is.vector(covariates)) covariates<-matrix(covariates)

if (is.null(covariates)==FALSE)
{
if (n!=dim(covariates)[1])stop('Error:  Number of covariate observations  is not equal to n')

if ((tind==FALSE)&(length(beta)!=dim(covariates)[2]))
{
stop('Error:  Number of initial beta parameters  is not equal to the number of covariates')
}

if ((tind==TRUE)&(length(beta)!=dim(covariates)[2]+1))
{
stop('Error:  Number of initial beta parameters  is not equal to the number of covariates+constant term')
}
}


if (is.null(tim)==TRUE) tim<-c(1:n)

if (tind==TRUE)
{
covariates<-cbind(rep(1,n),covariates)
}
covariates<-as.matrix(covariates)
Bcovariates<-as.matrix(covariates[inddat==1,])
covariatest<-as.matrix(covariates[posE,])

likfun.fun<-function(beta, Bcovariates=Bcovariates, covariatest=covariatest)
{
mllikpois <- -sum(covariatest%*%beta) + sum(exp(Bcovariates%*%beta ))
return(as.double(mllikpois))
}

if (maxty=='nlminb')
{minpois<- nlminb(objective=likfun.fun, start=beta, Bcovariates=Bcovariates, 
covariatest=covariatest)
llik <--minpois$objective
}

if (maxty=='optim')
{minpois<- optim(par=beta, fn=likfun.fun, gr = NULL, Bcovariates=Bcovariates, 
covariatest=covariatest, method = method)
llik <--minpois$value
}
fbeta<-minpois$par

if (modSim==FALSE)
{
cat(fill = TRUE)
cat('Number of observations  not used in the estimation process: ', (n-control),fill=TRUE)
cat('Total number of time observations: ',n,fill=TRUE)
cat("Number of events: ", length(posE), fill = TRUE)
cat(fill = TRUE)
cat('Convergence code: ',minpois$convergence, fill=TRUE)
if (minpois$convergence==0) cat('Convergence attained', fill=TRUE)
if (minpois$convergence!=0) 
{
cat(minpois$message, fill=TRUE)
stop('Convergence  not attained')
}
cat('Loglikelihood: ',round(llik,3), fill=TRUE)
cat(fill = TRUE)
cat('Estimated beta: ',round(fbeta,3), fill=TRUE)
cat(fill = TRUE)
}

npar <- length(fbeta)
lambdafit <- exp(covariates%*%fbeta )

auxoutput<-list(llik=llik, npar=npar, beta =fbeta, inddat=inddat, 
lambdafit=lambdafit, posE=posE,covariates=covariates, namcovariates=namcovariates,
tit=tit,tind=tind, t=tim)

if (modCI==TRUE)
{

VARbeta<-VARbeta.fun(covariates=covariates, lambdafit=lambdafit)

if (CIty=='Transf')
{
CIlambda<-CItran.fun(VARbeta=VARbeta, lambdafit=lambdafit, covariates=covariates,
 clevel=clevel)
}

if (CIty=='Delta')
{
CIlambda<-CIdelta.fun(VARbeta=VARbeta, lambdafit=lambdafit, covariates=covariates,
 clevel=clevel)
}
UIlambda<-CIlambda$UIlambda
LIlambda<-CIlambda$LIlambda
auxoutput<-list(llik=llik, npar=npar, beta =fbeta, inddat=inddat, VARbeta=VARbeta,
lambdafit=lambdafit, LIlambda=LIlambda, UIlambda=UIlambda, posE=posE, covariates=covariates,
namcovariates=namcovariates,tit=tit,tind=tind, t=tim)
}

if (dplot==TRUE)
{
dev.new()
if ((modCI==TRUE)&(is.null(lambdaylim))) lambdaylim<-c(min(LIlambda, na.rm=TRUE), max(UIlambda, na.rm=TRUE))
plot(tim,lambdafit, ty='n',ylab='intensity', xlab='time', ylim=lambdaylim, xlim=lambdaxlim)
if (modCI==TRUE)
{
lines(tim,UIlambda, col='blue', lty=3)
lines(tim,LIlambda, col='red', lty=3)
legend(xlegend, legend=c('Fitted intensity', paste('Upper CI',CIty,sep=' '), paste('Lower CI',CIty,sep=' ')), col=c('black', 'blue', 'red'), lty=c(1,3,3), cex=0.8)
}
lines(tim,lambdafit, lwd=2)
mtext(paste(tit), outer = TRUE, line = -2,cex=1)
}


return(auxoutput)

}
