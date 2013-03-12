\name{CalcResD.fun}
\alias{CalcResD.fun}

\title{Calculate NHPP residuals on disjoint intervals}

\description{This function calculates raw and scaled  residuals of a NHPP based on 
disjoint intervals. The scaled residuals can be Pearson or any other type of scaled residuals
defined by the function \eqn{h(t)}.}


\usage{
CalcResD.fun(obFPP, h = NULL, nint = NULL, lint = NULL, typeRes = NULL,
 modSim = "FALSE")
}



\arguments{
  \item{obFPP}{A list with (at least) elements posE, inddat and lambdafit; for example,
the output list from \code{fitPP.fun}, see that function for the details.}

  \item{lint}{Optional. Length of the  intervals to calculate 
the residuals.}

  \item{h}{Optional. Weight function to calculate the scaled residuals. By default,  
Pearson residuals with \eqn{h(t)=1/\sqrt{\hat \lambda(t)}} are calculated.}

  \item{typeRes}{Optional. Label indicating the type of scaled residuals. By default,  
Pearson residuals  are calculated and label is 'Pearson'.}

  \item{modSim}{Logical flag. If it is FALSE, some 
information on the intervals is shown on the screen.}

  \item{nint}{Number of intervals used to calculate the residuals. Intervals with 
the same length are considered. Only one of 
lint or nint must be specified.}
}


\details{

The  intervals used to calculate the residuals can be specified either
by  nint or lint; only one of the arguments must be provided. If nint is specified,
 intervals of equal length are calculated.

The raw residuals are  based on the increments of 
the raw process \eqn{R(t)=N_t-\int_0^t\hat\lambda(u)du} 
in   disjoint intervals  \eqn{(l_1, l_2)} centered on t:

\deqn{r'(l_1, l_2)=R(l_2)-R(l_1)=\sum_{t_i \in (l_1,l_2)}I_{t_i}-\int_{l_1}^{l_2} \hat \lambda(u)du.}


Residuals \eqn{r'(l_1, l_2)} are made 'instantaneous' dividing by the
length of the intervals (specified by the argument lint),
\eqn{r(l_1, l_2)=r'(l_1,l_2)/(l_2-l_1)}.


The function also calculates  the residuals scaled  with  the function \eqn{h(t)}

\deqn{r_{sca}(l_1, l_2)=\sum_{t_i \in (l_1,l_2)}h_{t_i}-\int_{l_1}^{l_2} h(u) \hat \lambda(u)du.}

By default,  Pearson residuals with \eqn{h(t)=1/\sqrt{\hat \lambda(t)}} are calculated.}

\value{A list with elements

  \item{RawRes}{Numeric vector of the raw residuals.}

  \item{ScaRes}{A list with elements ScaRes (vector of the scaled residuals) 
	and typeRes (name of the type of scaled residuals).}

  \item{emplambda}{Numeric vector of the empirical estimator of the PP intensity   on the 
considered intervals.}

  \item{fittedlambda}{Numeric vector of the   sum  of the intensities 
\eqn{\hat \lambda(t)}  on the considered intervals, divided by the length of the interval.}


  \item{lintV}{Numeric vector of the  exact length of each interval. 
The exact length is defined as the number of observations  in each interval used in the 
estimation (observations with inddat=1).}

  \item{lint}{Input argument.}

  \item{nint}{Input argument.}

  \item{pm}{Numeric vector of the mean point of the intervals.}

  \item{typeI}{Label indicating  the type  of intervals
	used to calculate the residuals, 'Disjoint' .}

  \item{h}{Input argument.}

  \item{obFPP}{Input argument.}
}



\references{

Abaurrea, J., Asin, J., Cebrian, A.C. and Centelles, A.
(2007). Modeling and forecasting extreme heat events in the central
Ebro valley, a continental-Mediterranean area. \emph{Global and
Planetary Change}, 57(1-2), 43-58.

Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005).
Residual analysis for spatial point processes.
\emph{Journal of the Royal Statistical Society}, Series B 67,617-666.



Brillinger, D. (1994). Time series, point processes and hybrids. \emph{Can. J. Statist.},
22, 177-206.

Lewis, P. (1972). Recent  results in the statistical analysis of univariate point
processes. In \emph{Stochastic point processes} (Ed. P. Lewis), 1-54. Wiley.


}




\seealso{\code{\link{CalcRes.fun}},  \code{\link{unifres.fun}}, 
\code{\link{graphres.fun}}}

\examples{

X1<-rnorm(1000)
X2<-rnorm(1000)


modE<-fitPP.fun(tind='TRUE',covariates=cbind(X1,X2), 
	posE=round(runif(40,1,1000)), inddat=rep(1,1000),
	tim=c(1:1000), tit='Simulated example',beta=c(1,0,0),
	dplot=FALSE,modCI=FALSE,modSim=TRUE)

#Residuals, based on 20 disjoint intervals of length 50, from the fitted NHPP modE  

ResDE<-CalcResD.fun(obFPP=modE,lint=50)



}

