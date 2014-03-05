\name{POTevents.fun}
\alias{POTevents.fun}
\title{Calculate extreme events using a POT approach}

\description{This function calculates the characteristics of the extreme events
of a series  \eqn{(x_i)}  defined using a peak over threshold (POT) method
 with an extreme threshold. The initial  and the maximum intensity positions,  
the mean excess, the maximum excess  and the length of each event are calculated. 
}

\usage{POTevents.fun(T, thres, date = NULL)}

\arguments{
  \item{T}{Numeric vector, the series \eqn{(x_i)} to calculate the
extreme events.}

  \item{thres}{Threshold value  used to define the extreme events.}

  \item{date}{Optional. A vector or matrix indicating the date of each observation.}
}


\details{One of the elements of the output from this function is  a vector (inddat)
which marks the observations  that should be used in the estimation of a point  process,
resulting from a POT approach. The  observations  to be considered in the estimation
are marked with 1 and correspond to the non occurrence  observations and to
a single occurrence point per event. The occurence point is defined as the point 
where maximum intensity of the event  occurs.The 
observations in an extreme event which are not the occurrence point are marked with 0
and treated as non observed. }

\value{A list with components
  \item{Pi }{Vector of  the initial points of the extreme events.}
  \item{datePi }{Date of the initial points Pi.}
  \item{Px }{Vector of the points of  maximum excess of the extreme events.}
  \item{datePx }{Vector of the date of the maximum excess points Px.}
  \item{Im }{Vector of the mean excesses (over the threshold) of the extreme events.}
  \item{Ix }{Vector of the maximum excesses (over the threshold) of the extreme events.}
  \item{L }{Vector of the lengths of the extreme events.}
  \item{inddat }{Index equal to 1 in the observations used in the estimation process 
and  to 0  in the others.}
}



\seealso{\code{\link{fitPP.fun}}}

\examples{
data(BarTxTn)
dateB<-cbind(BarTxTn$ano,BarTxTn$mes,BarTxTn$diames)
BarEv<-POTevents.fun(T=BarTxTn$Tx,thres=318, date=dateB)
}