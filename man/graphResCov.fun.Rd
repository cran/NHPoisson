\name{graphResCov.fun}

\alias{graphResCov.fun}

\title{Perform lurking variable plots for a set of variables}

\description{This function performs  lurking variable plots 
for a set of variables. The function 
\code{\link{graphResX.fun}}  performs the lurking  variable plot for one variable and
\code{graphResCov.fun}  calls this function for a set of variables;
see \code{\link{graphResX.fun}} for details.}

\usage{
graphResCov.fun(Xvar, nint, mlePP, h = NULL, typeRes = "Pearson", namX = NULL, 
 histWgraph=TRUE, plotDisp=c(2,2), tit = "")}

\arguments{
  \item{Xvar}{Matrix of variables (each column is a variable).}

  \item{nint}{Number of intervals each  covariate  is divided into to perform
 the lurking variable plot.}

  \item{mlePP}{An object of class \code{\link{mlePP-class}}; usually,
the output  from \code{\link{fitPP.fun}}.}

  \item{typeRes}{Label indicating the type of residuals ("Raw" or any type of 
scaled residuals such as "Pearson") used in the plots.}

  \item{h}{Optional. Weight function  used to calculate the scaled residuals (if
typeRes is not equal to "Raw"). By default,  Pearson residuals with 
\eqn{h(t)=1/\sqrt{\hat \lambda(t)}} are calculated. \eqn{\hat \lambda(t)} is provided
by element lambdafit in mlePP.}

  \item{namX}{Optional. Vector of  the names of the variables in Xvar.}

  \item{histWgraph}{ Logical flag.   If it is TRUE,  a new graphical device is opened
 with the option \code{record=TRUE}, so that the history of all plots is recorded in the new device. 
 This option may not work on some platforms; for example, RStudio does not allow the user to open
 new  graphical devices.}


  \item{plotDisp}{ A vector of the form \code{c(nr, nc)}.  The lurking variable plots will be drawn in 
a nr\eqn{\times}nc array. It is used as argument \code{mfrow} in \code{\link{par}}. By default,  
a 2 \eqn{\times} 2 window is used.}

  \item{tit}{Character string. A title for the plot.}
}



\value{A list with elements
  \item{mXres}{Matrix of  residuals (each column contains the  residuals of a variable).}
  \item{mXm}{Matrix of   mean values (each column contains
the  mean values of a variable in each interval).}
  \item{mXpc}{Matrix of the quantiles that define the intervals of each variable
(each column contains the  quantiles of one variable).}
  \item{nint}{Input argument.}
  \item{mlePP}{Input argument.}
}



\references{

Atkinson, A. (1985). \emph{Plots, transformations and regression.} Oxford 
University Press.

Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005).
Residual analysis for spatial point processes.
\emph{Journal of the Royal Statistical Society}, Series B 67,617-666.

Cebrian, A.C., Abaurrea, J. and Asin, J. (2015). {NHPoisson}: An {R} Package for 
Fitting and Validating Nonhomogeneous Poisson Processes. 
\emph{Journal of Statistical Software}, 64(6), 1-24.

}



\seealso{ \code{\link{graphResX.fun}}, \code{\link{graphres.fun}}  }

\examples{



#Simulated process without any relationship with variables Y1 and Y2
#The plots are performed dividing the  variables into  50 intervals
#Raw residuals. 

X1<-rnorm(500)
X2<-rnorm(500)
auxmlePP<-fitPP.fun(posE=round(runif(50,1,500)), inddat=rep(1,500),
	covariates=cbind(X1,X2),start=list(b0=1,b1=0,b2=0))

Y1<-rnorm(500)
Y2<-rnorm(500)
res<-graphResCov.fun(mlePP=auxmlePP, Xvar=cbind(Y1,Y2), nint=50,  
	typeRes="Raw",namX=c("Y1","Y2"),plotDisp=c(2,1))

#If more variables were specified in the argument Xvar, with
#the same 2X1 layout specified in plotDisp, the resulting plots could be 
#scrolled up and down with the "Page Up" and "Page Down" keys.

}