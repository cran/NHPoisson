\name{graphresU.fun}
\alias{graphresU.fun}

\title{Validation analysis of PP uniform (generalized) residuals}

\description{This function checks the properties that must
be fulfilled  by the uniform (generalized) residuals of a PP: uniform
character and uncorrelation. Optionally, the existence of patterns
versus covariates or potentially influent variables can be graphically analyzed.}



\usage{graphresU.fun(unires, posE,  Xvariables = NULL, namXv = NULL, flow = 0.5,
 tit = "", addlow = TRUE, indgraph = FALSE)}

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

  \item{indgraph}{Logical flag. If it is TRUE, plots are carried out in 
individual windows. By default,  2 \eqn{\times} 2 windows are used.}
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

Collett, D. (1994). Modelling survival data in medical research. Chapman \& Hall.

Ogata, Y. (1988). Statistical models for earthquake occurrences and residual
analysis for point processes. \emph{Journal of  the American Statistical Association},
83(401), 9-27.
}



\seealso{
 \code{\link{unifres.fun}}, \code{\link{transfH.fun}}
}

\examples{

X1<-rnorm(500)
X2<-rnorm(500)

graphresU.fun(unires=runif(30,0,1),posE=round(runif(30,0,500)), 
	Xvariables=cbind(X1,X2), namXv=c("X1","X2"),tit="Example",flow=0.7)

}