\name{VARbeta.fun}
\alias{VARbeta.fun}

\title{Calculate the variance matrix of the \eqn{\hat \beta} vector.}

\description{This function estimates  the variance matrix of the MLE of the vector of 
\eqn{\beta} parameters, using the asymptotic distribution of the MLE.}

\usage{VARbeta.fun(covariates, lambdafit)}

\arguments{
  \item{covariates}{Matrix of covariates (each column is a covariate).}
  \item{lambdafit}{Numeric vector, the fitted PP intensity \eqn{\hat \lambda(t)}.}
}

\details{The variance matrix is calculated as the inverse of the hessian. The inverse of the matrix
is calculated using the  solve function. If this function  leads to an error in the calculation, the
inverse is tried to be calculated via its Cholesky decomposition. If this option also fails, 
the variance matrix is not estimated.}

\value{
  \item{VARbeta }{Variance matrix of the \eqn{\hat \beta} vector.}
}

\note{The function \code{\link{fitPP.fun}} calls this function.}

\references{
Casella, G. and Berger, R.L., (2002). \emph{Statistical inference.} Brooks/Cole.

}


\seealso{ \code{\link{CItran.fun}}, \code{\link{CIdelta.fun}}}

\examples{

lambdafit<-runif(100,0,1)
X<-cbind(rep(1,100),rnorm(100),rnorm(100))

aux<-VARbeta.fun(covariates=X, lambdafit=lambdafit)


}