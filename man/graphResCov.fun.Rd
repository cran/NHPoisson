\name{graphResCov.fun}

\alias{graphResCov.fun}

\title{Perform lurking variable plots for a set of variables}

\description{This function performs  lurking variable plots 
for a set of variables. The function 
\code{\link{graphResX.fun}}  performs the lurking  variable plot for one variable and
\code{graphResCov.fun}  calls this function for a set of variables;
see \code{\link{graphResX.fun}} for details.}

\usage{
graphResCov.fun(Xvar, nint, obFPP, h = NULL, typeRes = "Pearson", namX = NULL, 
indgraph = "FALSE", tit = "")}

\arguments{
  \item{Xvar}{Matrix of variables (each column is a variable).}

  \item{nint}{Number of intervals each  covariate  is divided into to perform
 the lurking variable plot.}

  \item{obFPP}{A list with (at least) elements posE, inddat and lambdafit; for example,
the output list from \code{\link{fitPP.fun}}, see that function for details.}

  \item{typeRes}{Label indicating the type of residuals ('Raw' or any type of 
scaled residuals such as 'Pearson') used in the plots.}

  \item{h}{Optional. Weight function  used to calculate the scaled residuals (if
typeRes is not equal to 'Raw'). By default,  Pearson residuals with 
\eqn{h(t)=1/\sqrt{\hat \lambda(t)}} are calculated. \eqn{\hat \lambda(t)} is provided
by element lambdafit in obFPP.}

  \item{namX}{Optional. Vector of  the names of the variables in Xvar.}

  \item{indgraph}{Logical flag. If it is TRUE, each plot is carried out in an
individual window. By default, windows with  2 \eqn{\times} 2 plots are used.}

  \item{tit}{Character string. A title for the plot.}
}



\value{A list with elements
  \item{mXres}{Matrix of  residuals (each column contains the  residuals of a variable).}
  \item{mXm}{Matrix of   mean values (each column contains
the  mean values of a variable in each interval).}
  \item{mXpc}{Matrix of the quantiles that define the intervals of each variable
(each column contains the  quantiles of one variable).}
  \item{nint}{Input argument.}
  \item{obFPP}{Input argument.}
}



\references{

Atkinson, A. (1985). \emph{Plots, transformations and regression.} Oxford 
University Press.

Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005).
Residual analysis for spatial point processes.
\emph{Journal of the Royal Statistical Society}, Series B 67,617-666.

}



\seealso{ \code{\link{graphResX.fun}}, \code{\link{graphres.fun}}  }

\examples{



#Simulated process without any relationship with variables X1 and X2
#The plots are performed dividing the  variables into  50 intervals
#Raw residuals. 

X1<-rnorm(500)
X2<-rnorm(500)

res<-graphResCov.fun(obFPP=list(posE=round(runif(50,1,500)), inddat=rep(1,500),
 lambdafit=runif(500,0,0.2)), X=cbind(X1,X2), nint=50,  
	typeRes="Raw",namX=c("X1","X2"),indgraph="FALSE")

}