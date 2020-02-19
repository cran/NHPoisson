\name{confintAsin.fun}
\alias{confintAsin.fun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute confidence intervals for  the  \eqn{\beta} parameters}

\description{ This function computes confidence intervals for  the  \eqn{\beta} parameters.}

\usage{confintAsin.fun(mlePP, level = 0.95)}


\arguments{
  \item{mlePP}{A \code{"\linkS4class{mlePP}"}-class object; usually the output from \code{\link{fitPP.fun}}.}

  \item{level}{The confidence level required for the intervals.}
}


\details{ The confidence intervals  calculated by this function are based on the asymptotic normal 
approximation of th MLE of the \eqn{\beta} parameters, that is
\eqn{( \hat \beta -z_{(1-\alpha/2)}s.e.(\hat \beta ),  \hat \beta +z_{(1-\alpha/2)} s.e.(\hat \beta ) )}
with \eqn{\alpha=1-level}
}



\value{A matrix with two columns, the first contains the lower limits of the confidence intervals 
of all the parameters and the second the upper limits.}


\references{
Casella, G. and Berger, R.L., (2002). \emph{Statistical inference.} Brooks/Cole.

Cebrian, A.C., Abaurrea, J. and Asin, J. (2015). {NHPoisson}: An {R} Package for 
Fitting and Validating Nonhomogeneous Poisson Processes. 
\emph{Journal of Statistical Software}, 64(6), 1-24.

}


\seealso{\code{\link{confint}}, \code{\link{VARbeta.fun}}
}

\examples{
data(BarTxTn)

covB<-cbind(cos(2*pi*BarTxTn$dia/365), sin(2*pi*BarTxTn$dia/365), 
	BarTxTn$TTx,BarTxTn$Txm31,BarTxTn$Txm31**2)

BarEv<-POTevents.fun(T=BarTxTn$Tx,thres=318, 
	date=cbind(BarTxTn$ano,BarTxTn$mes,BarTxTn$dia))

mod1B<-fitPP.fun(covariates=covB, 
	posE=BarEv$Px, inddat=BarEv$inddat,
	tit="BAR Tx; cos, sin, TTx, Txm31, Txm31**2", 
	start=list(b0=-100,b1=1,b2=-1,b3=0,b4=0,b5=0))

confintAsin.fun(mod1B)

}
