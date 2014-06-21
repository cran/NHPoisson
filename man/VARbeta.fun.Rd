\name{VARbeta.fun}
\alias{VARbeta.fun}

\title{Calculate the covariance matrix of the \eqn{\hat \beta} vector.}

\description{This function estimates  the covariance matrix of the ML estimators of the  
\eqn{\beta} parameters, using the asymptotic distribution and properties of the ML estimators.}

\usage{VARbeta.fun(covariates, lambdafit)}

\arguments{
  \item{covariates}{Matrix of covariates (each column is a covariate).}
  \item{lambdafit}{Numeric vector, the fitted PP intensity \eqn{\hat \lambda(t)}.}
}

\details{The covariance matrix is calculated as the inverse of the  negative of the hessian matrix. The inverse of the matrix
is calculated using the  solve function. If this function  leads to an error in the calculation, the
inverse is calculated via its Cholesky decomposition. If this option also fails, 
the covariance matrix is not estimated and a  matrix  of dimension \eqn{0 \times 0} is returned. }

\value{
  \item{VARbeta }{Coariance matrix of the \eqn{\hat \beta} vector.  It has
 an attribute, called 'CalMethod' which shows the method  used to calculate the  inverse of the matrix: 
'Solve function', 'Cholesky' or 'Not possible'.}
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