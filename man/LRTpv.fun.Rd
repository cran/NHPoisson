\name{LRTpv.fun}
\alias{LRTpv.fun}

\title{Calculates the p-value of a likelihood ratio test for each covariate in the NHPP}
 
\description{This function calculates, for each covariate included in the NHPP,  the p-value 
of a likelihood  ratio test comparing  the  original fitted NHPP with the model excluding 
a covariate from the linear predictor.}

\usage{
LRTpv.fun(fitPPobj)
}

\arguments{
   \item{fitPPobj}{The model  to be analyzed. Usually the output list from the function \code{\link{fitPP.fun}};  
it must have at least  five elements called:  covariates, tind, posE, inddat and beta }
}

\details{If the model has an intercept, a LRT is  carried for all the covariates  in the linear predictor
except the intercept. If the model has not an intercept and there is only one covariate, no test can be carried out.
}

\value{A matrix with one column, which contains the LRT p-values for all the covariates in the model (except the intercept)
}





\seealso{ \code{\link{fitPP.fun}}, \code{\link{testlik.fun}}}

\examples{


data(BarTxTn)
covB<-cbind(cos(2*pi*BarTxTn$dia/365), sin(2*pi*BarTxTn$dia/365), 
	BarTxTn$TTx,BarTxTn$Txm31,BarTxTn$Txm31**2)
BarEv<-POTevents.fun(T=BarTxTn$Tx,thres=318, 
	date=cbind(BarTxTn$ano,BarTxTn$mes,BarTxTn$dia))



mod1B<-fitPP.fun(tind="TRUE",covariates=covB, 
	posE=BarEv$Px, inddat=BarEv$inddat,
	tit="BAR Tx; cos, sin, TTx, Txm31, Txm31**2", 
	beta=c(-100,1,10,0,0,0),dplot=FALSE, modCI=FALSE)

LRTpv.fun(mod1B)

}

