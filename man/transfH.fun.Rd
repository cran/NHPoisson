\name{transfH.fun}
\alias{transfH.fun}

\title{Transforming  a NHPP into a HPP}

\description{This function transforms the  points \eqn{t^{NH}_i} of a NHPP into 
the occurrence points  
\eqn{t^{H}_i} of a HPP  of rate 1.}

\usage{transfH.fun(mlePP)}

\arguments{
  \item{mlePP}{An object of class \code{\link{mlePP-class}}; usually,
the output  from \code{\link{fitPP.fun}}.}
}


\value{A list with elements

\item{posEH}{Numeric vector of the transformed occurrence times 
of the HPP.}
\item{posE}{Slot of the  input argument mlePP.}
\item{lambdafit}{Slot of the  input argument mlePP.}
\item{inddat}{Slot of the  input argument mlePP.}
}

\details{Transformation of the NHPP points \eqn{t^{NH}_i}  into 
the HPP points \eqn{t^{H}_i}  is based on
 the time scale transformation,  
\deqn{t^H_i=\int_0^{t^{NH}_i}\lambda(t)dt.}
(usually the estimated value \eqn{\hat \lambda(t)} is used in the transformation.)}


\references{

Cebrian, A.C., Abaurrea, J. and Asin, J. (2015). {NHPoisson}: An {R} Package for 
Fitting and Validating Nonhomogeneous Poisson Processes. 
\emph{Journal of Statistical Software}, 64(6), 1-25.

Cox, D.R., Isham, V., 1980. \emph{Point Processes}. Chapman and Hall.

Daley, D. and D. Vere-Jones (2003). 
\emph{An Introduction to the Theory of Point Processes}. Springer.


}


\seealso{\code{\link{simNHP.fun}}}

\examples{

X1<-rnorm(500)
X2<-rnorm(500)
auxmlePP<-fitPP.fun(posE=round(runif(50,1,500)), inddat=rep(1,500),
	covariates=cbind(X1,X2),start=list(b0=1,b1=0,b2=0))


posEH<-transfH.fun(auxmlePP)
}