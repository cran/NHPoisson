\name{fitPP.fun}
\alias{fitPP.fun}
\title{Fit  a non homogeneous Poisson Process}

\description{
This function fits by maximum likelihood a NHPP where the  intensity \eqn{\lambda(t)}  
is formulated  as a function of covariates.  It also calculates and plots 
approximate confidence intervals for \eqn{\lambda(t)}.}

\usage{
fitPP.fun(covariates = NULL, beta, posE = NULL, inddat = NULL, POTob = NULL, 
namcovariates = NULL, n = NULL, tind = "TRUE", tim = NULL, maxty='nlminb',
method="Nelder-Mead", modCI = "TRUE", CIty = "Transf", clevel = 0.95,
 tit = '', modSim = "FALSE", dplot = TRUE, xlegend = "topleft",
lambdaxlim=NULL,lambdaylim=NULL)
}

\arguments{

  \item{covariates}{Matrix of the covariates to be included in the 
linear predictor of the PP intensity (each column is a covariate).}

  \item{beta}{Numeric vector of the initial values  for the estimation of 
the \eqn{\beta} parameters. }

  \item{posE}{Optional (see Details section). Numeric vector of the position of the PP occurrence points.}

  \item{inddat}{Optional (see Details section). Index vector equal to 1 for the observations  used in the estimation process
By default, all the observations are considered.}

  \item{POTob}{Optional (see Details section). List with  elements T and thres 
that defines the PP resulting from a POT approach; 
see \code{\link{POTevents.fun}} for more details.}

  \item{namcovariates}{Optional. Vector of  the names of the variables in covariates.}

  
\item{n}{Optional. Number of observations  in the observation period; 
it is only neccessary if POTob, inddat and covariates are NULL.}

  \item{tind}{Logical flag. If it is TRUE, an independent term is fitted in the 
linear predictor.}

  \item{tim}{Optional. Time vector  of  the observation period.
By default, a vector 1,...n is considered.}


   \item{maxty}{Label indicating  the function used to maximize the loglikelihood function.
 If it is 'nlminb'  the function \code{\link{nlminb}} is used (the default option). If it is 'optim'
the function \code{\link{optim}} is used; in this case, the method of optimization 
can be chosen with the  following argument.}

   \item{method}{Label indicating  the algorithm of optimization if 
argument \code{maxty} is 'optim'. See \code{\link{optim}} for details. The default option
is 'Nelder-Mead'.}

  \item{modCI}{Logical flag. If it is TRUE, confidence intervals  
for \eqn{\lambda(t)} values are calculated.}

  \item{CIty}{Label indicating  the method to calculate the approximate
confidence intervals  for \eqn{\lambda(t)}. It can be 'Delta' for the delta method
or 'Transf' for transformed asymptotic intervals (default); see \code{\link{CIdelta.fun}} and
\code{\link{CItran.fun}} for details.}

  \item{clevel}{Confidence level of the  confidence intervals.}

  \item{tit}{Character string. A title for the plot.}

  \item{modSim}{Logical flag. If it is FALSE, information on the  
estimation process is shown on the screen.
 For simulation process, the option TRUE should be prefered.}

  \item{dplot}{Logical flag. If it is TRUE,  the fitted intensity is plotted.}

  \item{xlegend}{Label indicating the position where the legend on 
the graph will be located. }

   \item{lambdaxlim}{Optional. Numeric vector of length 2, giving the  lowest and highest values which determine the x  range.}

    \item{lambdaylim}{Optional. Numeric vector of length 2, giving the  lowest and highest values which determine the y  range.}

}

\details{
A Poisson process (PP)  is usually specified by a vector containing the  occurrence 
points of the process \eqn{(t_i)_{i=1}^k},  (argument posE). 
Since PP are often used in the framework of POT models, \code{fitPP.fun} also 
provides the possibility of 
using as input the  series of the observed values  in a POT model 
\eqn{(x_i)_{i=1}^n} and the threshold used to define the  extreme events 
(argument POTob). 


In the case of PP defined by a POT approach, 
the observations of the extreme events which are 
not defined as the occurrence point are not considered in the estimation. This is done
 through the argument inddat, see \code{\link{POTevents.fun}}. If  the input is provided via argument POTob, index inddat
is calculated automatically. See \cite{Coles (2001)} for more details on the POT approach.


The estimation of the \eqn{\beta} covariance matrix is based on the 
asymptotic distribution of the MLE \eqn{\hat \beta}, and calculated as the inverse of the hessian.
Confidence intervals for \eqn{\lambda(t)} can be calculated using  two approaches, 
the  delta method or a transformation of the confidence interval for the 
linear predictor \eqn{\nu(t)=\textbf{X(t)} \beta}.
The interval for \eqn{\nu(t)} is also based on the asymptotic properties of the 
MLE  \eqn{\hat \nu(t)}. See \cite{Casella (2002)} for more details on ML theory and delta method.}



\value{
   A list with elements
  \item{llik }{Value of the loglikelihood function.}
  \item{npar }{Number of estimated parameters.}
  \item{beta }{Vector of the MLE \eqn{\hat \beta}.}
  \item{inddat }{Input argument.}
  \item{VARbeta }{Covariance matrix of the \eqn{\beta} parameters}
  \item{lambdafit }{Vector of the fitted intensity \eqn{\hat \lambda(t)}.}
  \item{LIlambda }{Vector of lower extremes of the CI}
  \item{LUlambda }{Vector of upper extremes of the CI.}
  \item{posE }{Input argument.}
  \item{namcovariates }{Input argument.}
  \item{tit }{Input argument.}
  \item{tind }{Input argument.}
}

\note{A homogeneous Poisson process (HPP) can be fitted as a particular case,  
using  an intensity defined by  only an independent term  and no covariates.}


\seealso{ \code{\link{POTevents.fun}}, \code{\link{globalval.fun}},
\code{\link{VARbeta.fun}}, \code{\link{CItran.fun}}, \code{\link{CIdelta.fun}} }

\references{

Coles, S. (2001). \emph{An introduction to statistical modelling of extreme
values.} Springer.

Casella, G. and Berger, R.L., (2002). \emph{Statistical inference.} Brooks/Cole.

}

\examples{
#model fitted  using as input posE and inddat and  no confidence intervals 

data(BarTxTn)
covB<-cbind(cos(2*pi*BarTxTn$dia/365), sin(2*pi*BarTxTn$dia/365), 
	BarTxTn$TTx,BarTxTn$Txm31,BarTxTn$Txm31**2)
BarEv<-POTevents.fun(T=BarTxTn$Tx,thres=318, 
	date=cbind(BarTxTn$ano,BarTxTn$mes,BarTxTn$dia))


mod1B<-fitPP.fun(tind='TRUE',covariates=covB, 
	posE=BarEv$Px, inddat=BarEv$inddat,
	tit='BAR Tx; cos, sin, TTx, Txm31, Txm31**2', 
	beta=c(-100,1,10,0,0,0))

#model fitted  using as input  a list from POTevents.fun and with  confidence intervals 
tiempoB<-BarTxTn$ano+rep(c(0:152)/153,55)

mod2B<-fitPP.fun(tind='TRUE',covariates=covB, 
	POTob=list(T=BarTxTn$Tx, thres=318),
	tim=tiempoB, tit='BAR Tx; cos, sin, TTx, Txm31, Txm31**2', 
	beta=c(-100,1,10,0,0,0),CIty='Delta',modCI=TRUE,
	modSim=TRUE)
}