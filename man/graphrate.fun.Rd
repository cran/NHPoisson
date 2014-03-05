\name{graphrate.fun}
\alias{graphrate.fun}

\title{Plot fitted and  empirical  PP occurrence rates}

\description{This function calculates the empirical  and the cumulative fitted occurrence rate
 of a PP  on overlapping or disjoint intervals and plot them versus time.}



\usage{
graphrate.fun(objres = NULL, fittedlambda = NULL, emplambda = NULL, t = NULL, 
lint = NULL, typeI = "Disjoint", tit = "", scax = NULL, scay = NULL, 
xlegend = "topleft")}

\details{Either the argument objres or the set of arguments (fittedlambda, emplambda, t) 
must be specified. If objres is provided, fittedlambda, emplambda, t,lint and typeI  
are ignored.

In order to make comparable the empirical and the fitted occurrence rates,
a cumulative fitted rate must be used. That means that argument fittedlambda 
must be the  sum of the
intensities fitted by the model  over  the same interval where the 
empirical rates have been calculated.}



\arguments{

 \item{objres}{Optional. A list with (at least) elements fittedlambda, emplambda, t, 
and typeI. 
For example, the output from \code{\link{CalcRes.fun}} or \code{\link{CalcResD.fun}}; 
see those functions for details.}

  \item{fittedlambda}{Optional. Numeric vector of the  cumulative  
fitted intensities \eqn{\hat \lambda(t)} over the considered intervals
(and usually divided by the length of the interval).}

  \item{emplambda}{Optional. Numeric vector of the  empirical PP 
occurrence rates estimated  over the considered intervals (usually divided by the
length of the interval).}

  \item{t}{Optional. Time vector of the  PP observation period.}

  \item{lint}{Optional. Length of the  intervals used to calculate 
the empirical and the (cumulative) fitted occurrence intensities.}

  \item{typeI}{Label indicating  the type ('Overlapping' or 'Disjoint') of  
the intervals.}

  \item{tit}{Character string. A title for the plot.}

  \item{scax}{Optional. Vector of two values giving  the range of values for 
the x-axis. An adequate range is selected by default.}


  \item{scay}{Optional. Vector of two values giving  the range of values for the 
y-axis. An adequate range is selected by default.}

  \item{xlegend}{Label indicating the position where the legend on 
the graph will be located.}
}


\seealso{ \code{\link{CalcRes.fun}}, \code{\link{CalcResD.fun}}}

\examples{
##plot of rates based on overlapping intervals
graphrate.fun(emplambda=runif(500,0,1), fittedlambda=runif(500,0,1), 
	t=c(1:500), lint=100, tit="Example", typeI="Overlapping")

#plot of rates based on disjoint intervals
graphrate.fun(emplambda=runif(50,0,1), fittedlambda=runif(50,0,1), 
	t=c(1:50), lint=10, tit="Example", typeI="Disjoint")

#Example using objres as input. In this example X1 has no influence on the rate;
#consequently the fitted rate is almost a constant.

X1<-rnorm(1000)

modE<-fitPP.fun(tind="TRUE",covariates=cbind(X1), 
	posE=round(runif(40,1,1000)), inddat=rep(1,1000),
	tim=c(1:1000), tit="Simulated example", beta=c(1,0),
	modCI=FALSE,modSim=TRUE,dplot=FALSE)

ResDE<-CalcResD.fun(obFPP=modE,lint=50)

graphrate.fun(ResDE, tit="Example")



}