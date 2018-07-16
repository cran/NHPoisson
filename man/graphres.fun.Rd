\name{graphres.fun}
\alias{graphres.fun}

\title{Plot of the NHPP residuals}

\description{This function plots the residuals of a NHPP (raw or scaled, 
overlapping or disjoint) 
versus time or other variables which are monotonous functions of time.}

\usage{
graphres.fun(objres = NULL, typeRes = "Raw", t = NULL, res = NULL, lint = NULL, 
posE = NULL, fittedlambda = NULL, typeI = "Disjoint", Xvariables = NULL, 
namXv = NULL, histWgraph=TRUE, plotDisp=c(2,2), addlow = FALSE, lwd = 2, 
tit = "", flow = 0.5, xlegend = "topleft", legcex = 0.5)}

\arguments{
 \item{objres}{Optional. A list with (at least) elements t, typeI and Rawres  and/or ScaRes, depending on 
the value of typeRes. For example, the output list from the functions 
\code{\link{CalcRes.fun}} or \code{\link{CalcResD.fun}}; see those functions for details.}

  \item{typeRes}{Label indicating the type of residuals ("Raw" or any type of scaled 
residuals such as "Pearson").}

  \item{t}{Optional. Time vector of the PP observation period.}

  \item{res}{Optional. Vector of residuals.}

  \item{lint}{Optional. Length of the  intervals used to calculate the residuals.}

  \item{posE}{Optional. Numeric vector of the PP occurrence times. 
Only used when typeI = "Overlapping".}

 \item{fittedlambda}{Optional. Vector of the cumulative  fitted PP intensity 
over the intervals. Used  to calculate the envelopes when typeRes="Raw".}

  \item{typeI}{Label indicating  the type ("Overlapping" or "Disjoint") of intervals.}

  \item{Xvariables}{Optional. Matrix  of the variables for
 the residual plots (each column is a variable).}

  \item{namXv}{Optional. Vector of the  names of the variables in Xvariables.}

  \item{histWgraph}{ Logical flag.   If it is TRUE,  a new graphical device is opened
 with the option \code{record=TRUE}, so that the history of all plots is recorded in the new device. 
 This option may not work on some platforms; for example, RStudio does not allow the user to open
 new  graphical devices.}


  \item{plotDisp}{ A vector of the form \code{c(nr, nc)}.  The residual plots will be drawn in 
a nr\eqn{\times}nc array. It is used as argument \code{mfrow} in \code{\link{par}}. By default,  
a 2 \eqn{\times} 2 window is used.}

  \item{tit}{Character string. A title for the plots.}

  \item{addlow}{Logical flag. If it is TRUE, a lowess is added to the residual plots.}

  \item{lwd}{Argument lwd   for plotting the lowess lines, see \code{\link{par}}
for details.}


  \item{flow}{Argument f for the lowess, see \code{\link{lowess}} for details.}

  \item{xlegend}{Label giving the position of the graph where the legend will 
be located.}

  \item{legcex}{Argument cex for the legend, see \code{\link{par}} for details.}
}

\details{Either argument objres or pair of arguments (t,res) must be specified. 
If objres is provided, arguments t,res, typeRes, typeI, posE and  fittedlambda
are ignored.  

A residual plot versus time is always performed. 
These plots are intended for time or variables which are monotonous functions,
since residuals are  calculated over a given time interval and plotted versus the
value of the variables in the mean point of the interval.

A smoother (lowess) of the residuals can be optionally added to the plots. 
In the case of overlapping intervals,
the residuals of the occurrence points   are  marked differently from the rest.
In the case typeRes="Raw" (if argument fittedlambda is available) or 
typeRes="Pearson",  envelopes for the residuals are also plotted.
The envelopes are based on an  approach analogous to the one shown in
Baddeley et al. (2005) for spatial  Poisson processes.  
The envelopes for raw residuals are,
 
\deqn{\pm {2 \over l_2-l_1} \sqrt{\sum_{ i \in (l_1,l_2)} \hat \lambda(i)} }

 where index i  runs over  the integers in the interval \eqn{(l_1,l_2)}. 
The envelopes for the Pearson residuals are,
\deqn{\pm 2/\sqrt{l_2-l_1}.}  

These plots  allow us to analyze the effect on the intensity, of the covariates included
in the model or other potentially influent variables. They 
show  if the mean or the dispersion of the residuals vary sistematically, see  for 
example  residual analysis in Atkinson (1985) or Collett (1994).}



\references{

Atkinson, A. (1985). \emph{Plots, transformations and regression.} Oxford 
University Press.

Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005).
Residual analysis for spatial point processes.
\emph{Journal of the Royal Statistical Society, Series B}, 67, 617-666.

Cebrian, A.C., Abaurrea, J. and Asin, J. (2015). {NHPoisson}: An {R} Package for 
Fitting and Validating Nonhomogeneous Poisson Processes. 
\emph{Journal of Statistical Software}, 64(6), 1-25.


Collett, D. (1994). Modelling survival data in medical research. Chapman & Hall.

}

\seealso{\code{\link{graphrate.fun}}}

\examples{

#Example using objres as input

X1<-c(1:1000)**0.5

modE<-fitPP.fun(tind=TRUE,covariates=cbind(X1), 
	posE=round(runif(40,1,1000)), inddat=rep(1,1000),
	tim=c(1:1000), tit="Simulated example", start=list(b0=1,b1=0),
	modSim = TRUE, dplot = FALSE)

ResDE<-CalcResD.fun(mlePP=modE,lint=50)
graphres.fun(objres=ResDE, typeRes="Raw", Xvariables=cbind(X1),
	namXv=c("X1"), plotDisp=c(2,1), addlow=TRUE,tit="Example")


#Example using the set of arguments res, t and fittedlambda as input
#In this case, with typeI="Disjoint", only values of t, fittedlambda and Xvariables 
#in the midpoint of the intervals must be provided.

#Since   a 1X1 layout is  specified in plotDisp and only one  
#graphical device is opened by default, the two  resulting plots can be scrolled  
#up and down  with the "Page Up" and "Page Down" keys.

X1<-c(1:500)**0.5
graphres.fun(res=rnorm(50),posE=round(runif(50,1,500)),
	fittedlambda=runif(500,0,1)[seq(5,495,10)],
	t=seq(5,495,10), typeRes="Raw", typeI="Disjoint",Xvariables=X1[seq(5,495,10)],
	namXv=c("X1"), plotDisp=c(1,1), tit="Example 2",lint=10)
}