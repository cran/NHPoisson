\name{LRTpv.fun}
\alias{LRTpv.fun}

\title{Likelihood ratio test  of the covariates in a model}
 
\description{This function carries out, for each covariate in the model (except the intercept),  
a likelihood  ratio test comparing  the  original fitted NHPP with the model excluding 
that covariate from the linear predictor, and  provides the corresponding p-values.}

\usage{
LRTpv.fun(mlePP)
}

\arguments{
   \item{mlePP}{An object of class \code{\link{mlePP-class}}; usually,
the output  from \code{\link{fitPP.fun}}. The fitted model cannot include fixed parameters.}
}

\details{A LRT is  carried for all the covariates  in the linear predictor
except the intercept. If the model has not an intercept and there is only one covariate, no test can be carried out.
}

\value{A matrix with one column, which contains the LRT p-values for all the covariates in the model (except the intercept)
}





\seealso{ \code{\link{fitPP.fun}}, \code{\link{testlik.fun}}, \code{\link{dropAIC.fun}}, \code{\link{addAIC.fun}}}

\examples{
data(BarTxTn)
covB<-cbind(cos(2*pi*BarTxTn$dia/365), sin(2*pi*BarTxTn$dia/365), 
	BarTxTn$TTx,BarTxTn$Txm31,BarTxTn$Txm31**2)
BarEv<-POTevents.fun(T=BarTxTn$Tx,thres=318, 
	date=cbind(BarTxTn$ano,BarTxTn$mes,BarTxTn$dia))

mod1B<-fitPP.fun(tind=TRUE,covariates=covB, 
	posE=BarEv$Px, inddat=BarEv$inddat,
	tit="BAR Tx; cos, sin, TTx, Txm31, Txm31**2", 
	start=list(b0=-100,b1=1,b2=10,b3=0,b4=0,b5=0),dplot=FALSE, modCI=FALSE)

LRTpv.fun(mod1B)
}

