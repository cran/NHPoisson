\name{emplambdaD.fun}
\alias{emplambdaD.fun}

\title{Empirical occurrence rates of a NHPP on disjoint intervals}

\description{This function calculates the empirical  occurrence  rates
of a point process using  disjoint intervals.  The rate is
assigned to the mean point of the  interval. A plot of the empirical 
rate over time can be performed optionally. }


\usage{emplambdaD.fun(posE, t, lint=NULL, nint = NULL, plotEmp = TRUE, inddat = NULL, 
tit = '', scax = NULL, scay = NULL)}


\details{The  intervals  can be specified either
by  nint or lint; only one of the arguments must be provided.}

\arguments{
  \item{posE}{Numeric vector of the  position of the occurrence points of 
the NHPP (or any point process in time).}

  \item{t}{Time index of the observation period. The simplest option is  1,...,n with 
n the length of the period.}

  \item{lint}{Optional (alternative argument to nint). Length of the intervals used to calculate 
the rates.}

  \item{nint}{Optional (alternative argument to lint). Number of intervals 
(of equal length) used to to calculate the rates.
It is an alternative way  to lint for identifying the	intervals. }

  \item{plotEmp}{Logical flag. If it is TRUE, a plot of the empirical rate 
is carried out.}

  \item{inddat}{Optional. Index vector equal to 1 for the observations used in the 
estimation process. By default, all the observations are considered, 
see \code{\link{POTevents.fun}}.}

  \item{tit}{Character string. A  title  for the plot.}

  \item{scax}{Optional.  A two element vector indicating the x-scale for the plot.}

  \item{scay}{Optional.  A two element vector indicating the y-scale for the plot.}
}


\value{A list with elements
  \item{emplambda }{Vector of the empirical rates.}
  \item{lint }{Input argument.}
  \item{nint }{Input  argument.}
}


\seealso{\code{\link{emplambda.fun}}, \code{\link{fitPP.fun}}, 
\code{\link{POTevents.fun}}}

\examples{

data(BarTxTn)

BarEv<-POTevents.fun(T=BarTxTn$Tx,thres=318, 
	date=cbind(BarTxTn$ano,BarTxTn$mes,BarTxTn$dia))


# empirical rate based on disjoint intervals using nint to specify the intervals
emplambdaDB<-emplambdaD.fun(posE=BarEv$Px,inddat=BarEv$inddat, t=c(1:8415), 
	nint=55)

# empirical rate based on disjoint intervals using lint to specify the intervals
emplambdaDB<-emplambdaD.fun(posE=BarEv$Px,inddat=BarEv$inddat, t=c(1:8415), 
	lint=153)
}