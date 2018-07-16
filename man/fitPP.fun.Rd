\name{fitPP.fun}
\alias{fitPP.fun}
\title{Fit of a non homogeneous Poisson Process}

\description{
This function fits by maximum likelihood a NHPP where the  intensity \eqn{\lambda(t)}  
is formulated  as a function of covariates.  It also calculates and plots 
approximate confidence intervals for \eqn{\lambda(t)}.}

\usage{
fitPP.fun(covariates = NULL, start, fixed=list(), posE = NULL, inddat = NULL, 
POTob = NULL,  nobs = NULL, tind = TRUE, tim = NULL, minfun="nlminb",
 modCI = "TRUE", CIty = "Transf", clevel = 0.95,
 tit = "", modSim = "FALSE", dplot = TRUE, xlegend = "topleft",
lambdaxlim=NULL,lambdaylim=NULL,...)
}

\arguments{

  \item{covariates}{Matrix of the covariates to be included in the 
linear predictor of the PP intensity (each column is a covariate). It is advisable to give 
names to the columns of this matrix
 (using \code{\link{dimnames}}), since they  will be used in the output.  
Otherwise the default names 'Covariate i' are used.
The offset covariates must be included in this matrix.  A maximum of 50 covariates are allowed.}

  \item{start}{Named list of the initial values  for the estimation of 
the \eqn{\beta} parameters (including fixed parameters).
 The names of the list  must be (compulsory): b0 (for the intercept), b1 (for the first column in \code{covariates}), b2 (for the second column), 
 b3 (for the third column), etc.}


  \item{fixed}{Named list of the fixed \eqn{\beta} parameters. The elements of this list  must be elements 
of the list \code{start}. }


  \item{posE}{Optional (see Details section). Numeric vector of the position of the PP occurrence points.}

  \item{inddat}{Optional (see Details section). Index vector equal to 1 for the observations  used in the estimation process
By default, all the observations are considered.}

  \item{POTob}{Optional (see Details section). List with  elements T and thres 
that defines the PP resulting from a POT approach; 
see \code{\link{POTevents.fun}} for more details.}


  
\item{nobs}{Optional. Number of observations  in the observation period; 
it is only neccessary if POTob, inddat and covariates are NULL.}

  \item{tind}{Logical flag. If it is TRUE, an independent term is fitted in the 
linear predictor. It cannot be a character string, so TRUE and not'TRUE' should be used.}

  \item{tim}{Optional. Time vector  of  the observation period.
By default, a vector 1,...n is considered.}


   \item{minfun}{Label indicating  the function to minimize the negative of the loglikelihood function.
There are two possible values: "nlminb"  (the default option) and  "optim". In the last case, 
the method of optimization can be chosen with an additional method argument.}



  \item{modCI}{Logical flag. If it is TRUE, confidence intervals  
for \eqn{\lambda(t)} values are calculated.}

  \item{CIty}{Label indicating  the method to calculate the approximate
confidence intervals  for \eqn{\lambda(t)}. It can be "Transf" for transformed asymptotic intervals (default) or 
"Delta" for the delta method; see \code{\link{CItran.fun}}  and \code{\link{CIdelta.fun}} for details.}

  \item{clevel}{Confidence level of the  confidence intervals.}

  \item{tit}{Character string. A title for the plot.}

  \item{modSim}{Logical flag. If it is FALSE, information on the  
estimation process is shown on the screen.
 For simulation process, the option TRUE should be preferred.}

  \item{dplot}{Logical flag. If it is TRUE,  the fitted intensity is plotted.}

  \item{xlegend}{Label indicating the position where the legend on 
the graph will be located. }

   \item{lambdaxlim}{Optional. Numeric vector of length 2, giving the  lowest and highest values which determine the x  range.}

    \item{lambdaylim}{Optional. Numeric vector of length 2, giving the  lowest and highest values which determine the y  range.}
 
    \item{...}{Further arguments to pass to \code{\link{optim}} or to \code{\link{nlminb}} (depending on the value of  the 
\code{minfun} argument).}

}

\details{
A Poisson process (PP)  is usually specified by a vector containing the  occurrence 
points of the process \eqn{(t_i)_{i=1}^k},  (argument posE). 
Since PP are often used in the framework of POT models, \code{\link{fitPP.fun}} also 
provides the possibility of 
using as input the  series of the observed values  in a POT model 
\eqn{(x_i)_{i=1}^n} and the threshold used to define the  extreme events 
(argument POTob). 


In the case of PP defined by a POT approach, 
the observations of the extreme events which are 
not defined as the occurrence point are not considered in the estimation. This is done
 through the argument inddat, see \code{\link{POTevents.fun}}. If  the input is provided via argument POTob, index inddat
is calculated automatically. See \cite{Coles (2001)} for more details on the POT approach.

The maximization of the loglikelihood function can be done using  two different optimization  routines, 
\code{\link{optim}} or \code{\link{nlminb}}, selected in the argument \code{minfun}. Depending on 
the covariates included in the function, one routine can  succeed to converge when the other fails.

 This function allows  us to  keep fixed some \eqn{\beta} parameters  (offset terms).  This can be 
used to specify an a priori known component to be included in the linear predictor during fitting. The fixed parameters
must be specified in the \code{fixed}  argument (and also   in \code{start}); 
the fixed covariates must be included as columns of \code{covariates}.

 

The estimation of the \eqn{\hat \beta} covariance matrix is based on the 
asymptotic distribution of the MLE \eqn{\hat \beta}, and calculated as the inverse of the negative of the  hessian matrix.
Confidence intervals for \eqn{\lambda(t)} can be calculated using  two approaches
specified in the argument \code{CIty}. See \cite{Casella (2002)} for more details on ML theory and delta method.}



\value{
An object of class  \code{\linkS4class{mlePP}}, which is a subclass of \code{\linkS4class{mle}}. 
Consequently, many of the generic functions with \code{mle} methods, such as 
\code{\link{logLik}} or \code{\link{summary}}, can be applied to the output of this function.  Some other generic 
functions related to fitted models, such as \code{AIC} or \code{BIC}, can  also  be applied to \code{mlePP} objects. 
}

\note{A homogeneous Poisson process (HPP) can be fitted as a particular case,  
using  an intensity defined by  only an intercept  and no covariate.}


\seealso{ \code{\link{POTevents.fun}}, \code{\link{globalval.fun}},
\code{\link{VARbeta.fun}}, \code{\link{CItran.fun}}, \code{\link{CIdelta.fun}} }

\references{

Cebrian, A.C., Abaurrea, J. and Asin, J. (2015). {NHPoisson}: An {R} Package for 
Fitting and Validating Nonhomogeneous Poisson Processes. 
\emph{Journal of Statistical Software}, 64(6), 1-25.


Coles, S. (2001). \emph{An introduction to statistical modelling of extreme
values.} Springer.

Casella, G. and Berger, R.L., (2002). \emph{Statistical inference.} Brooks/Cole.

Kutoyants Y.A. (1998).\emph{Statistical inference for spatial Poisson processes.} 
Lecture notes in Statistics 134. Springer.

}

\examples{
#model fitted  using as input posE and inddat and  no confidence intervals 

data(BarTxTn)
covB<-cbind(cos(2*pi*BarTxTn$dia/365), sin(2*pi*BarTxTn$dia/365), 
	BarTxTn$TTx,BarTxTn$Txm31,BarTxTn$Txm31**2)
BarEv<-POTevents.fun(T=BarTxTn$Tx,thres=318, 
	date=cbind(BarTxTn$ano,BarTxTn$mes,BarTxTn$dia))


mod1B<-fitPP.fun(covariates=covB, 
	posE=BarEv$Px, inddat=BarEv$inddat,
	tit="BAR Tx; cos, sin, TTx, Txm31, Txm31**2", 
	start=list(b0=-100,b1=1,b2=-1,b3=0,b4=0,b5=0))

#model fitted  using as input  a list from POTevents.fun and with  confidence intervals 

tiempoB<-BarTxTn$ano+rep(c(0:152)/153,55)

mod2B<-fitPP.fun(covariates=covB, 
	POTob=list(T=BarTxTn$Tx, thres=318),
	tim=tiempoB, tit="BAR Tx; cos, sin, TTx, Txm31, Txm31**2", 
	start=list(b0=-100,b1=1,b2=-1,b3=0,b4=0,b5=0),CIty="Delta",modCI=TRUE,
	modSim=TRUE)

#model  with a fixed parameter (b0)

mod1BF<-fitPP.fun(covariates=covB, 
	posE=BarEv$Px, inddat=BarEv$inddat,
	tit="BAR Tx; cos, sin, TTx, Txm31, Txm31**2", 
	start=list(b0=-89,b1=1,b2=10,b3=0,b4=0,b5=0), 
	fixed=list(b0=-100))
}