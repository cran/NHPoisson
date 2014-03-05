\name{CItran.fun}
\alias{CItran.fun}


\title{Confidence intervals for  \eqn{\lambda(t)} based on transformation }

\description{
Given the  \eqn{\beta} variance matrix (or its estimation),  an approximate
confidence interval for each \eqn{\lambda(t)} is calculated using a transformation of 
the confidence interval for the linear
predictor \eqn{\nu(t)=\textbf{X(t)} \beta}. The transformation is \eqn{\exp(I_i)}, 
where \eqn{I_i} are the extremes of the intervals for \eqn{\nu(t)}.}

\usage{CItran.fun(VARbeta, lambdafit, covariates, clevel = 0.95)}


\arguments{
  \item{VARbeta}{(Estimated) Variance matrix of the \eqn{\beta}  parameter vector of
the PP intensity.}

  \item{lambdafit}{Numeric vector of fitted values of the PP intensity  
\eqn{\hat \lambda(t)}.}

  \item{covariates}{Matrix of covariates to estimate the  PP intensity.}

  \item{clevel}{Confidence level of  the confidence intervals. A value in  the interval 
(0,1).}
}



\value{A list with  elements
  \item{LIlambda}{Numeric vector of the lower values of the intervals.}
  \item{UIlambda}{Numeric vector of the upper values of the intervals.}
  \item{lambdafit}{Input argument.}
}


\references{
Casella, G. and Berger, R.L., (2002). \emph{Statistical inference.} Brooks/Cole.
}

\note{\code{fitPP.fun} calls \code{CItran.fun}  when  the argument is \emph{CIty='Transf'}.}


\seealso{\code{\link{CIdelta.fun}}, \code{\link{fitPP.fun}}, \code{\link{VARbeta.fun}}
}

\examples{
aux<-CItran.fun(VARbeta=0.01, lambdafit=exp(rnorm(100)), covariates=matrix(rep(1,100)),
	 clevel=0.95)
}