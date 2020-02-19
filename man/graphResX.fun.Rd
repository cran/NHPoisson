\name{graphResX.fun}

\alias{graphResX.fun}


\title{Perform a lurking  variable plot}

\description{This function performs a lurking  variable plot to analyze the residuals  in terms
of different levels of the variable. }


\usage{graphResX.fun(X, nint, mlePP, typeRes = "Pearson", h = NULL, namX = NULL)}


\arguments{
  \item{X}{Numeric vector, the variable for the lurking variable plot.}
  \item{nint}{Number of intervals or levels the variable  is divided into.}

  \item{mlePP}{An object of class \code{\link{mlePP-class}}; usually,
the output  from \code{\link{fitPP.fun}}.}

  \item{typeRes}{Label indicating the type of residuals ('Raw' or any type of 
scaled residuals such as 'Pearson').}

  \item{h}{Optional. Weight function  used to calculate the scaled residuals (if
typeRes is not equal to 'Raw'). By default,  Pearson residuals with 
\eqn{h(t)=1/\sqrt{\hat \lambda(t)}} are calculated. \eqn{\hat \lambda(t)} is provided
by  the lambdafit  slot in mlePP.}

  \item{namX}{Optional. Name of variable X.}

}

\details{The residuals for different levels of the variable are analyzed.
For a variable \eqn{X(t)}, the considered levels are
\deqn{W(P_{X,j}, P_{X,j+1})=\{ t:  P_{X,j} \le X(t) < P_{X,j+1} \}}
where  \eqn{P_{X,i}} is the sample j-percentile of  X. This type of plot is 
specially useful for variables which are not a  monotonous function of time.


In the case typeRes='Raw'  or typeRes='Pearson',  envelopes for the residuals 
are also plotted.
The envelopes are based on an  approach analogous to the one  in
Baddeley et al. (2005) for spatial  Poisson processes.   
The envelopes for raw residuals are
 
\deqn{\pm {2 \over l_W} \sqrt{\sum_i \hat \lambda(i)}}

 where index i  runs over  the integers in the level \eqn{W(P_{X,j}, P_{X,j+1})},
and \eqn{l_W} is its length (number of observations in \eqn{W}). 
The envelopes for the Pearson residuals are,
\deqn{\pm 2/\sqrt{l_W}.}
}


\value{A list with elements

  \item{Xres}{Vector of  residuals. }
  \item{xm}{Vector of the  mean value  of the variable in each interval.}
  \item{pc}{Vector of the quantiles that define the levels of the variable.}
  \item{typeRes}{Input argument.}
  \item{namX}{Input argument.}
  \item{lambdafit}{Input argument.}
  \item{posE}{Input argument.}
}


\references{

Atkinson, A. (1985). \emph{Plots, transformations and regression.} Oxford 
University Press.

Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005).
Residual analysis for spatial point processes.
\emph{Journal of the Royal Statistical Society}, Series B 67, 617-666.

Cebrian, A.C., Abaurrea, J. and Asin, J. (2015). {NHPoisson}: An {R} Package for 
Fitting and Validating Nonhomogeneous Poisson Processes. 
\emph{Journal of Statistical Software}, 64(6), 1-24.

}



\seealso{\code{\link{graphResCov.fun}}, \code{\link{graphres.fun}}  }

\examples{

##Simulated process not related to variable X
##Plots dividing the  variable into  50 levels

X1<-rnorm(500)
X2<-rnorm(500)
auxmlePP<-fitPP.fun(posE=round(runif(50,1,500)), inddat=rep(1,500),
	covariates=cbind(X1,X2),start=list(b0=1,b1=0,b2=0))



##Raw residuals
res<-graphResX.fun(X=rnorm(500),nint=50,mlePP=auxmlePP,typeRes="Raw")

##Pearson residuals
res<-graphResX.fun(X=rnorm(500),nint=50,mlePP=auxmlePP,typeRes="Pearson")


}
