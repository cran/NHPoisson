\name{testlik.fun}
\alias{testlik.fun}

\title{Likelihood ratio test to compare two nested models}

\description{This function performs a likelihood ratio test, a test to compare 
the fit of two models, where the first one (the null model ModR) is a particular case of the 
other (the alternative model ModG).}

\usage{ testlik.fun(ModG, ModR) }

\arguments{
  \item{ModG}{An object of class \code{\link{mlePP-class}}; usually,
the output  from \code{\link{fitPP.fun}}.}
  \item{ModR}{An object of class \code{\link{mlePP-class}}; usually,
the output  from \code{\link{fitPP.fun}}.}}

\details{The test statistic is twice the difference in the log-likelihoods 
of the models.
Under the null, the  statistic follows a \eqn{\chi^2} distribution with degrees 
of freedom df2-df1,the number of parameters of modG and modR respectively.}

\value{ A list with elements
\item{pv}{P-value of the likelihood ratio test.}
\item{ModG}{Input argument.}
\item{ModR}{Input argument.} }


\seealso{ \code{\link{fitPP.fun}},\code{\link{LRTpv.fun}}}

\references{

Casella, G. and Berger, R.L., (2002). \emph{Statistical inference.} Brooks/Cole.

}

\examples{

##The alternative model modB is specified  by the output fitPP.fun
##The null model modBR is specified  by a list with elments llik and npar


data(BarTxTn)

covB<-cbind(cos(2*pi*BarTxTn$dia/365), sin(2*pi*BarTxTn$dia/365), 
	BarTxTn$TTx,BarTxTn$Txm31,BarTxTn$Txm31**2)


modB<-fitPP.fun(tind=TRUE,covariates=covB, POTob=list(T=BarTxTn$Tx, thres=318),
	tim=c(1:8415), tit="BAR Tx; cos, sin, TTx, Txm31, Txm31**2", 
	start=list(b0=-100,b1=1,b2=10,b3=0,b4=0,b5=0),dplot=FALSE,modCI=TRUE,
	modSim=TRUE)

modBR<-fitPP.fun(tind=TRUE,covariates=covB[,1:4], 
	POTob=list(T=BarTxTn$Tx, thres=318),
	tim=c(1:8415), tit="BAR Tx; cos, sin, TTx, Txm31", 
	start=list(b0=-100,b1=1,b2=10,b3=0,b4=0),dplot=FALSE,modCI=TRUE,	modSim=TRUE)


aux<-testlik.fun(ModG=modB,ModR=modBR)


}

