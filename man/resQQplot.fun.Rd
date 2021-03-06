\name{resQQplot.fun}
\alias{resQQplot.fun}
\alias{resSim.fun}



\title{Perform a qqplot for the residuals of a NHPP}

\description{This function  performs a qqplot comparing the empirical quantiles of the 
residuals with the expected quantiles under the fitted NHPP,
calculated by a Monte Carlo approach.

It calls the auxiliary function \code{resSim.fun} 
(not intended for the users), see Details section.}


\usage{
resQQplot.fun(nsim, objres, covariates, clevel = 0.95, cores = 1, 
tit ="", fixed.seed=NULL, histWgraph=TRUE)}

\arguments{
  \item{nsim}{Number of simulations for the calculations.}

  \item{objres}{A list with the same elements of the output list from the
function \code{\link{CalcRes.fun}} or \code{\link{CalcResD.fun}}.}

  \item{covariates}{Matrix of covariates to fit the NHPP (each column is a covariate).}

  \item{clevel}{Confidence level of the  residual envelope.}

  \item{cores}{Optional. Number of cores of the computer to be used in the 
calculations. Default: one core is used.}

  \item{tit}{Character string. A title for the plot.}

  \item{fixed.seed}{An integer or NULL. If it is an integer, that is the value  used to set the seed 
in  random  generation processes. It it is NULL, a random seed is used.}

  \item{histWgraph}{ Logical flag.  Only used in Windows platforms. If it is TRUE,  a new graphical device is opened
 with the option \code{record=TRUE}.}

}

\details{
The expected quantiles are calculated as the median values  of the simulated samples. 
Confidence intervals for each quantile \eqn{r_{(i)}} with pointwise  significance level 
\code{clevel}  are calculated as quantiles of probability 1-clevel /2 and clevel/2 
of the simulated sample for each residual. 

All type of residuals (disjoint or overlapping and Pearson or raw residuals) are 
supported by this function. However, the  qqplot for overlapping residuals can  be
a high time consuming process. So,  disjoint
residuals should be prefered in this function.

The  auxiliary function \code{resSim.fun} generates a NHPP with intensity \eqn{\lambda(t)},
 fits the model using the covariate matrix  and calculates the residuals.}


\value{
A list with elements
  \item{resmed}{Numeric vector containing the mean of the simulated residuals 
in each point.}
  \item{ressup}{Numeric vector of the upper values of the simulated envelopes.}
  \item{resinf}{Numeric vector of the lower values of the simulated envelopes.}
  \item{objres}{Input argument.}
  \item{nsim}{Input argument.}
  \item{fixed.seed}{Input argument.}
}



\references{

Cebrian, A.C., Abaurrea, J. and Asin, J. (2015). {NHPoisson}: An {R} Package for 
Fitting and Validating Nonhomogeneous Poisson Processes. 
\emph{Journal of Statistical Software}, 64(6), 1-24.

}


\seealso{\code{\link{simNHP.fun}}, \code{\link{GenEnv.fun}}}

\examples{

X1<-rnorm(500)
X2<-rnorm(500)

aux<-fitPP.fun(tind=TRUE,covariates=cbind(X1,X2), 
	posE=round(runif(40,1,500)), inddat=rep(1,500),
	tim=c(1:500), tit="Simulated example", start=list(b0=1,b1=0,b2=0),dplot=FALSE)

auxRes<-CalcResD.fun(mlePP=aux,lint=50)


#if we want reproducible results, we can fixed the seed in the generation process
#(the number of cores used in the calculations must also be the same to reproduce
# the result)

auxqq<-resQQplot.fun(nsim=50,objres=auxRes, covariates=cbind(X1,X2), fixed.seed=123)

}
