\name{simNHP.fun}
\alias{simNHP.fun}
\alias{buscar}

\title{Generate the occurrence points of a NHPP}

\description{This function generates the occurrence times of the points  
of a NHPP with  a given  time-varying intensity \eqn{\lambda(t)}, 
in a period (0, T). The length of argument lambda determines T,
the length of the observation period.

It calls the auxiliary function \code{buscar} (not intended
for the users), see Details section.}




\usage{simNHP.fun(lambda, fixed.seed=NULL)}



\arguments{
  \item{lambda}{Numeric vector, the time varying intensity \eqn{\lambda(t)}
to generate the NHPP.}

  \item{fixed.seed}{An integer or NULL. If it is an integer, that is the value  used to set the seed 
in  random  generation processes. It it is NULL, a random seed is used.}

}

\details{The generation of the NHPP points consists in two steps. 
First,  the points of a homogeneous PP  of intensity 1 are generated using 
independent exponentials. Then, the  homogeneous occurrence times are transformed into 
the points  of a  non homogeneous process with intensity \eqn{\lambda(t)}.
This transformation is performed by the auxiliary function \code{buscar} 
(not intended for the user).}

\value{A list with elements
  \item{posNH }{Numeric vector of the occurrences times of the NHPP generated in 
the observation period (0,T).}
  \item{lambda }{Input argument.}
  \item{fixed.seed}{Input argument.}
}



\seealso{ \code{\link{GenEnv.fun}}, \code{\link{resQQplot.fun}}}

\references{

Cebrian, A.C., Abaurrea, J. and Asin, J. (2015). {NHPoisson}: An {R} Package for 
Fitting and Validating Nonhomogeneous Poisson Processes. 
\emph{Journal of Statistical Software}, 64(6), 1-24.

Ross, S.M. (2006). \emph{Simulation.} Academic Press.

}

\examples{
#Generation of the occurrence times of a homogeneours PP with  constant intensity 
#0.01 in a period  of time of length 1000

aux<-simNHP.fun(lambda=rep(0.01,1000))
aux$posNH

#if we want reproducible results, we can fixed the seed in the generation process
aux<-simNHP.fun(lambda=rep(0.01,1000),fixed.seed=123)
aux$posNH

#and the result is:
#  [1]  85 143 275 279 284 316 347 362 634 637 738 786 814 852 870 955


#Generation of  the occurrence times of a NHPP with  time-varying intensity t in 
#a period  of time of length 500

t<-runif(500, 0.01,0.1)
aux<-simNHP.fun(lambda=t)
aux$posNH
}



