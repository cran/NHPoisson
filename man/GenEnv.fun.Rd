\name{GenEnv.fun}
\alias{GenEnv.fun}
\alias{funSim.fun}
\alias{posk.fun}


\title{Simulated envelopes}

\description{This function calculates a  point estimation and an envelope
for a  given statistic using a Monte Carlo approach. The statistic must be a function 
of the occurrence points of a NHPP.

It calls the auxiliary function \code{funSim.fun} (not intended
for the users), see Details section.}


\details{The auxiliary function \code{funSim.fun} 
generates a  simulated sample  of the occurrence points in a NHPP
and  calculates the corresponding statistic using the simulated points.}


\usage{
GenEnv.fun(nsim, lambda, fun.name, fun.args = NULL, clevel = 0.95, 
cores = 1, fixed.seed=NULL)}

\arguments{
  \item{nsim}{Number of simulations for the calculations.}

  \item{lambda}{Numeric vector of the intensity \eqn{\lambda(t)}
(or \eqn{\hat \lambda(t)}) of the NHPP.}

  \item{fun.name}{Name of the function defining the statistic to be estimated.}

  \item{fun.args}{Additional arguments for the function fun.name.}

  \item{clevel}{Confidence level of the envelope.}

  \item{cores}{Optional.  Number of cores of the computer
to be used in the calculations. Default: one core is used.}

  \item{fixed.seed}{An integer or NULL. If it is an integer, that is the value  used to set the seed 
in  random  generation processes. It it is NULL, a random seed is used.}

}



\value{ A list with elements
  \item{valmed}{Point estimation (mean value) of the statistic  to be calculated.}
  \item{valinf}{Lower value of the simulated CI.}
  \item{valsup}{Upper value of the simulated CI.}
  \item{lambda}{Input argument.}
  \item{nsim}{Input argument.}
  \item{nsimval}{Number of valid simulations (used in the calculation of the CI and the point estimation).}
  \item{fixed.seed}{Input argument.}	
}



\seealso{ \code{\link{simNHP.fun}}, \code{\link{resQQplot.fun}}}

\examples{
# Calculation of the point estimation and a 95% CI based on 100 simulations 
#for the second occurrence time of a NHPP with intensity lambdat.
#posk.fun(x, k) is a function that returns the value in the row k of vector x.

lambdat<-runif(1000,0.01,0.02)
aux<-GenEnv.fun(lambda=lambdat,fun.name="posk.fun",fun.args=2,nsim=100)

#if we want reproducible results, we can fixed the seed in the generation process
#(the number of cores used in the calculations must also be the same to reproduce 
#the result)

aux<-GenEnv.fun(lambda=lambdat,fun.name="posk.fun",fun.args=2,nsim=100,fixed.seed=123)

#the result (with 1 core): Lower interval:  25.55; Mean value:  136.06; Upper interval:  288

}

