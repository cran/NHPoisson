\name{summaryfit.fun}
\alias{summaryfit.fun}


\title{Summarizes the beta coefficients of the NHPP and their standar errors}

\description{This function  shows the estimated beta coefficients of the NHPP and their standar errors.}


\details{The  standar errors of the coefficients are calculated from the covariance matrix (VARbeta) 
which is calculated by \code{\link{fitPP.fun}}.}


\usage{summaryfit.fun(fitPPobj)}

\arguments{
  \item{fitPPobj}{The model  to be summarized. Usually the output list from the function \code{\link{fitPP.fun}};  
it must have at least two elements called, beta and VARbeta }
}



\value{A matrix with two columns, the first contains the estimated beta coefficients  of the NHPP and the second, their standard errors.}



\seealso{ \code{\link{fitPP.fun}}}

\examples{

data(BarTxTn)
covB<-cbind(cos(2*pi*BarTxTn$dia/365), sin(2*pi*BarTxTn$dia/365), 
	BarTxTn$TTx,BarTxTn$Txm31,BarTxTn$Txm31**2)
BarEv<-POTevents.fun(T=BarTxTn$Tx,thres=318, 
	date=cbind(BarTxTn$ano,BarTxTn$mes,BarTxTn$dia))


mod1B<-fitPP.fun(tind="TRUE",covariates=covB, 
	posE=BarEv$Px, inddat=BarEv$inddat,
	tit="BAR Tx; cos, sin, TTx, Txm31, Txm31**2", 
	beta=c(-100,1,10,0,0,0),dplot=FALSE)

summaryfit.fun(mod1B)

}

