\name{emplambda.fun}
\alias{emplambda.fun}


\title{Empirical occurrence rates of a NHPP on overlapping intervals}

\description{This function calculates the empirical  occurrence  rates
of a point process on overlapping intervals.  The empirical rate centered in each
time of the observation period is calculated using intervals of a given length. A plot of the empirical 
rate over time can be performed optionally. 
}

\usage{
emplambda.fun(posE, t, lint, plotEmp = TRUE, inddat = NULL, tit = '', 
scax = NULL, scay = NULL)}



\arguments{
  \item{posE}{Numeric vector of the  position of the occurrence points of 
the NHPP (or any point process in time).}

  \item{t}{Time index of the observation period. The simplest option is  1,...,n with 
n the length of the period.}

  \item{lint}{Length of the intervals used to calculate the rates.}


  \item{plotEmp}{Logical flag. If it is TRUE, a plot of the empirical rate 
is carried out.}

  \item{inddat}{Optional. Index vector equal to 1 for the observations used in the estimation process
By default, all the observations are considered, see \code{\link{POTevents.fun}}.}

  \item{tit}{Character string. A  title  for the plot.}

  \item{scax}{Optional.  A two element vector indicating the x-scale for the plot.}

  \item{scay}{Optional.  A two element vector indicating the y-scale for the plot.}
}


\value{A list with elements
  \item{emplambda }{Vector of the empirical rates.}
  \item{lint }{Input argument.}
}


\seealso{\code{\link{emplambdaD.fun}}, \code{\link{fitPP.fun}}, \code{\link{POTevents.fun}}}

\examples{

data(BarTxTn)

BarEv<-POTevents.fun(T=BarTxTn$Tx,thres=318, 
	date=cbind(BarTxTn$ano,BarTxTn$mes,BarTxTn$dia))

# empirical rate based on overlapping intervals
emplambdaB<-emplambda.fun(posE=BarEv$Px,inddat=BarEv$inddat, t=c(1:8415), 
	lint=153,  tit='Barcelona')


}