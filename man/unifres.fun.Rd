\name{unifres.fun}

\alias{unifres.fun}

\title{Calculate exponential and uniform (generalized) residuals of a HPP}

\description{
This function calculates the exponential \eqn{d_i} and the uniform (generalized)
residuals \eqn{u_i} of a 
HPP, using the occurrence points \eqn{t_i}.}

\usage{unifres.fun(posEH)}

\arguments{
  \item{posEH}{Numeric vector, the occurrence points of a HPP.}
}

\details{
The exponential residuals of a HPP are defined as the inter-event distances
\eqn{d_i=t_i-t_{i-1}}, that are an i.i.d. exponential sample. The series
\eqn{d_i} is  an example of the generalized residuals proposed by Cox and Snell (1968).
The uniform residuals, defined as the function  \eqn{\exp(-d_i)}
of the exponential residuals, are an i.i.d. uniform sample, see Ogata (1988).}

\value{A list with elements
  \item{expres}{Numeric vector of the exponential residuals.}
  \item{unires}{Numeric vector of the uniform residuals.}
  \item{posEH}{Input argument.}
}

\references{
Abaurrea, J., Asin, J., Cebrian, A.C. and Centelles, A.
(2007). Modeling and forecasting extreme heat events in the central
Ebro valley, a continental-Mediterranean area. \emph{Global and
Planetary Change}, 57(1-2), 43-58.


Cox, D. R. and Snell, E. J. (1968). A general definition of residuals.
\emph{Journal of  the Royal Statistical Society, series B}, 30(2), 248-275.
83(401), 9-27.

Ogata, Y. (1988). Statistical models for earthquake occurrences and residual
analysis for point processes.\emph{Journal of  the American Statistical Association},
83(401), 9-27.

}




\seealso{\code{\link{transfH.fun}}, \code{\link{graphresU.fun}}}

\examples{
## generates the occurrence times of a homogeneours PP with  constant intensity 0.01 
## and calculates de residuals

aux<-simNHP.fun(lambda=rep(0.01,1000))

res<-unifres.fun(aux$posNH)
}