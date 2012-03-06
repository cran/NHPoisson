\name{globalval.fun}
\alias{globalval.fun}

\title{Perform a global validation analysis for a NHPP}

\description{This function performs a thorough validation analysis for a fitted NHPP. 
It calculates
the (generalized) uniform  and the raw (or scaled) residuals, performs residual plots
for the uniform  residuals, and time residual and lurking variable plots for the 
raw  or scaled residuals.
It also  plots the fitted and  empirical estimations of the NHPP intensity.
Optionally, it also performs a  residual QQplot.}

\details{If typeI='Overlapping', argument lint is compulsory. If 
typeI='Disjoint', only one of the arguments lint or  nlint must be specified.

Variables in argument Xvar are used for the plots of the
(generalized) uniform residuals and for the lurking variable plots.
If Xvar is NULL, the variables in argument covariates are used.}

\usage{
globalval.fun(obFPP, lint = NULL, nint = NULL, covariates, Xvar = NULL, 
namXvar = NULL, Xvart = NULL, namXvart = NULL,  h = NULL, typeRes = NULL,
typeResLV='Pearson',typeI = "Disjoint", nsim = 100, clevel = 0.95, 
resqqplot = FALSE, nintLP = 100, tit = '', flow = 0.5, addlow = FALSE, 
indgraph = FALSE, scax = NULL, scay = NULL, legcex = 0.5, ngen = 100, 
cores = NULL, xlegend = "topleft")
}

\arguments{


  \item{obFPP}{A list with  the elements  of the output from \code{\link{fitPP.fun}}, 
see that function for details.}


  \item{lint}{Length of the  intervals  used to calculate the residuals.}

  \item{nint}{Number of intervals  used to calculate the residuals. 
Intervals of equal length are considered.  Only used if typeI='Disjoint'.
In that case, only one of  the arguments lint or nint must be specified.}

  \item{covariates}{Matrix of covariates used to fit the NHPP specified in objFPP
(each column is a covariate).}

   \item{Xvar}{Optional. Matrix of the lurking variables (each column is a variable). 
They are also used  for the plots of the (generalized) uniform residuals.}

  \item{namXvar}{Optional. Vector of names of the variables in Xvar.}


  \item{Xvart}{Optional. Matrix of the variables for the 
residual plots  for monotonous variables (each column is a variable). 
A time  plot is performed in all the cases.}

  \item{namXvart}{Optional. Vector of names of the variables in Xvart.}

  \item{h}{Optional. Weight function  to calculate the scaled residuals.
By default,  Pearson residuals with \deqn{h(t)=1/\sqrt{\hat \lambda(t)}} are calculated.
This function is used to calculate both the  scaled residuals and the residuals 
for the lurking variables (except if typeResLV='Raw').}


  \item{typeRes}{Optional. Label indicating the type of scaled residuals. By default,  
Pearson residuals  are calculated and label is 'Pearson'.}

  \item{typeResLV}{Label indicating the type of residuals ('Raw' or any type of 
scaled residuals such as 'Pearson') to calculate the residuals for the lurking variable
plots.}


 \item{typeI}{Label indicating  the type ('Overlapping' or 'Disjoint') of intervals
	used to calculate the residuals.}

  \item{clevel}{Confidence level of the  residual envelopes.} 

  \item{resqqplot}{Logical  flag. It is is TRUE, a residual qqplot is carried out.}

  \item{nsim}{Number of simulations  for the  residual qqplot.}


  \item{nintLP}{Number of levels considered in the lurking variables. It is used as argument
nint in the call of the function \code{\link{graphResCov.fun}}.}

  \item{tit}{Character string. A title for the plot.}

  \item{addlow}{Logical flag. If it is TRUE, a lowess is added in the residual plots.}
 
  \item{flow}{Argument f for the lowess smoother of the raw (or scaled) residual 
plots, see \code{\link{lowess}}.}

  \item{indgraph}{Logical flag. If it is TRUE, plots are carried out in 
individual windows. By default, windows with  2 \eqn{\times} 2 plots are used.}

  \item{scax}{Optional. Vector of two values indicating  the range of values for the x-axis 
in the fitted and empirical  rate plot. An adequate range is selected by default.}

  \item{scay}{Optional. Vector of two values indicating  the range of values for the x-axis
in the fitted and empirical  rate plot. An adequate range is selected by default.}

  \item{legcex}{cex argument  for the legend in the  residual time plots 
(see \code{\link{par}} for details).}

  \item{ngen}{Argument n  used in  the call of the function 
	\code{\link{resQQplot.fun}}; see that function for details.}


  \item{cores}{Optional. Number of cores of the computer
to be used in the calculations. Default: all cores available in the computer are used.}


  \item{xlegend}{Argument xlegend used in the call of the function 
\code{\link{graphrate.fun}}; see that function for details.}


}


\value{A list with the same elements that \code{\link{CalcRes.fun}} or 
\code{\link{CalcResD.fun}}
 (depending on the value of the argument typeI).}



\seealso{
\code{\link{graphres.fun}}, \code{\link{graphrate.fun}}, \code{\link{resQQplot.fun}},
\code{\link{graphResCov.fun}}, \code{\link{graphresU.fun}}}

\examples{

data(BarTxTn)

covB<-cbind(cos(2*pi*BarTxTn$dia/365), sin(2*pi*BarTxTn$dia/365), 
	BarTxTn$TTx,BarTxTn$Txm31,BarTxTn$Txm31**2)


modB<-fitPP.fun(tind='TRUE',covariates=covB, 
	POTob=list(T=BarTxTn$Tx, thres=318),
	tit='BAR Tx; cos, sin, TTx, Txm31, Txm31**2', 
	beta=c(-100,1,10,0,0,0),CIty='Transf',modCI=TRUE,
	modSim=TRUE,dplot=FALSE)

aux<-globalval.fun(obFPP=modB,lint=153,	covariates=covB,typeI='Disjoint', 
	typeRes='Raw',typeResLV='Raw',	resqqplot=FALSE)

##If typeRes and typeResLV are not specified, Pearson residuals are calculated
##by default.

aux<-globalval.fun(obFPP=modB,lint=153,	covariates=covB,typeI='Disjoint', 
	resqqplot=FALSE,Xvar=covB[,1])

}
