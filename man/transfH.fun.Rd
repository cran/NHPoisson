\name{transfH.fun}
\alias{transfH.fun}

\title{Transform  a NHPP into a HPP}

\description{This function transforms the  points \eqn{t^{NH}_i} of a NHPP into 
the occurrence points  
\eqn{t^{H}_i} of a HPP  of rate 1.}

\usage{transfH.fun(obFPP)}

\arguments{
  \item{obFPP}{A list with (at least)  elements posE, lambdafit and inddat. 
For example, the output list from \code{\link{fitPP.fun}}; 
see this function for details.}
}


\value{A list with elements

\item{posEH}{Numeric vector of the transformed occurrence times 
of the HPP.}
\item{posE}{Element of the  input argument obFPP.}
\item{lambdafit}{Element of the  input argument obFPP.}
\item{inddat}{Element of the  input argument obFPP.}
}

\details{Transformation of the NHPP points \eqn{t^{NH}_i}  into 
the HPP points \eqn{t^{H}_i}  is based on
 the time scale transformation,  
\deqn{t^H_i=\int_0^{t^{NH}_i}\lambda(t)dt.}
(usually the estimated value \eqn{\hat \lambda(t)} is used in the transformation.)}


\references{
Daley, D. and D. Vere-Jones (2003). 
\emph{An Introduction to the Theory of Point Processes}. Springer.

Cox, D.R., Isham, V., 1980. \emph{Point Processes}. Chapman and Hall.

}


\seealso{\code{\link{simNHP.fun}}}

\examples{

dat<-list(posE=runif(50, 0,500), lambdafit=rep(0.01,500), inddat=rep(1,500))
posEH<-transfH.fun(dat)
}