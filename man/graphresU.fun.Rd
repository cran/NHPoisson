\name{graphresU.fun}
\alias{graphresU.fun}

\title{Validation analysis of PP uniform (generalized) residuals}

\description{This function checks the properties that must
be fulfilled  by the uniform (generalized) residuals of a PP: uniform
character and uncorrelation. Optionally, the existence of patterns
versus covariates or potentially influent variables can be graphically analyzed.}



\usage{graphresU.fun(unires, posE,  Xvariables = NULL, namXv = NULL, flow = 0.5,
 tit = "", addlow = TRUE, histWgraph=TRUE, plotDisp=c(2,2), indgraph = FALSE)}

\arguments{
  \item{unires}{Numeric vector of the uniform residuals.}

  \item{posE}{Numeric vector of the occurrence times of the PP.}

  \item{Xvariables}{Matrix of variables to perform the residual plots (each column
is a variable).}

  \item{namXv}{Optional. Vector of names of the variables in Xvariables.}

  \item{tit}{Character string. A title for the plot.}

  \item{addlow}{Logical flag. If it is TRUE, a lowess is added to the plots.}

  \item{flow}{Argument f for the lowess smoother; see \code{\link{lowess}} 
for details.}


  \item{histWgraph}{ Logical flag.   If it is TRUE,  a new graphical device is opened
 with the option \code{record=TRUE}, so that the history of all plots is recorded in the new device. 
 This option may not work on some platforms; for example, RStudio does not allow the user to open
 new  graphical devices.}

  \item{plotDisp}{ A vector of the form \code{c(nr, nc)}.  The residual versus variables plots will be drawn in 
a nr\eqn{\times}nc array. It is used as argument \code{mfrow} in \code{\link{par}}. By default,  
a 2 \eqn{\times} 2 layout is used.}

  \item{indgraph}{Logical flag. If it is TRUE, the validation plots (except the residuals versus variables
 plots) are carried out in four1 \eqn{\times} 1 layouts. By default,  a 2 \eqn{\times} 2 layout is used.}
}

\details{The validation analysis  of the uniform
character consists in a  uniform Kolmogorov-Smirnov test and 
a qqplot with a 95\%  confidence band based on a beta distribution. 
The analysis of the serial correlation is based on the Pearson correlation 
coefficient, Ljung-Box  tests and a lagged  serial correlation plot. 
An index plot of the residuals and residual plots  versus the variables 
in argument Xvariables are performed to analyze the effect of
covariates   or  other potentially influent variables. These  plots
will show  if the mean or dispersion of the residuals vary sistematically, see model
diagnostic of Cox-Snell residuals in Collett (1994) for more details.}

\references{
Abaurrea, J., Asin, J., Cebrian, A.C. and Centelles, A.
(2007). Modeling and forecasting extreme heat events in the central
Ebro valley, a continental-Mediterranean area. \emph{Global and
Planetary Change}, 57(1-2), 43-58.

Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005).
Residual analysis for spatial point processes.
\emph{Journal of the Royal Statistical Society, Series B}, 67, 617-666.

Cebrian, A.C., Abaurrea, J. and Asin, J. (2015). {NHPoisson}: An {R} Package for 
Fitting and Validating Nonhomogeneous Poisson Processes. 
\emph{Journal of Statistical Software}, 64(6), 1-25.

Collett, D. (1994). Modelling survival data in medical research. Chapman \& Hall.

Ogata, Y. (1988). Statistical models for earthquake occurrences and residual
analysis for point processes. \emph{Journal of  the American Statistical Association},
83(401), 9-27.
}



\seealso{
 \code{\link{unifres.fun}}, \code{\link{transfH.fun}}
}

\examples{

#Since  only one graphical device is opened  and the argument histWgraph 
#is TRUE by default, the resulting residual plots  (three pages with the 
#considered 1X2 layout for the residual versus variables plot)   
#can be scrolled up and down with the "Page Up" and "Page Down" keys.

X1<-rnorm(500)
X2<-rnorm(500)

graphresU.fun(unires=runif(30,0,1),posE=round(runif(30,0,500)), 
  Xvariables=cbind(X1,X2), namXv=c("X1","X2"),tit="Example",flow=0.7,plotDisp=c(1,2))

}