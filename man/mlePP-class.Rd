\name{mlePP-class}
\Rdversion{1.1}
\docType{class}
\alias{mlePP-class}

\title{Class \code{"mlePP"} for results of maximum likelihood estimation of Poisson processes with covariates}

\description{This class encapsulates  the output from the  maximum likelihood  estimation of  a Poisson process where the 
intensity is modeled as a linear function of  covariates.}

\section{Objects from the Class}{
Objects can be created by calls of the form  \code{new("mlePP", ...)}, but most often as the 
result of a call to \code{\link{fitPP.fun}}.  }


\section{Slots}{
  \describe{
    \item{\code{call}:}{Object of class \code{"language"}.  The call to \code{\link{fitPP.fun}}.}
    \item{\code{coef}:}{Object of class \code{"numeric"}. The estimated coefficientes of the model. }
    \item{\code{fullcoef}:}{Object of class \code{"numeric"}. The  full coefficient vector, including the fixed 
	parameters of the model. It has an attribute, called 'TypeCoeff' which shows  the names of the fixed parameters.}
    \item{\code{vcov}:}{Object of class \code{"matrix"}. Approximate variance-covariance matrix of the estimated coefficients.
	 It has an attribute, called 'CalMethod' which shows the method  used to calcualte the  inverse of the information matrix: 
	'Solve function', 'Cholesky', 'Not possible' or 'Not required' if  \code{modCI=FALSE}.}
    \item{\code{min}:}{Object of class \code{"numeric"}. Minimum value of objective function, that is the 
	negative of the loglikelihood function. }
    \item{\code{details}:}{Object of class \code{"list"}. The output returned from \code{\link{optim}}. 
	If  \code{\link{nlminb}} is used to minimize the function, it is NULL.}
    \item{\code{minuslogl}:}{Object of class \code{"function"}.  The negative  of the loglikelihood function.}
    \item{\code{nobs}:}{Object of class \code{"integer"}.  The number of observations.}
    \item{\code{method}:}{Object of class \code{"character"}.  It is a bit different from the slot in  the extended 
	class 	\code{\linkS4class{mle}}: here,
	it is the input argument \code{minfun} of \code{\link{fitPP.fun}} instead of the  method used 
	in \code{\link{optim}} (this information already appears in \code{details}).}
    \item{\code{detailsb}:}{Object of class \code{"list"}.The output returned from \code{\link{nlminb}}. 
	If  \code{\link{optim}} is used to minimize the function, it is NULL.}
    \item{\code{npar}:}{Object of class \code{"integer"}. Number of estimated parameters.}
    \item{\code{inddat}:}{Object of class \code{"numeric"}.  Input argument  of \code{\link{fitPP.fun}}.}
    \item{\code{lambdafit}:}{Object of class \code{"numeric"}. Vector of the fitted intensity \eqn{\hat \lambda(t)}.}
    \item{\code{LIlambda}:}{Object of class \code{"numeric"}. Vector of lower limits of the CI.}
    \item{\code{UIlambda}:}{Object of class \code{"numeric"}. Vector of upper limits of the CI.}
    \item{\code{convergence}:}{Object of class \code{"integer"}. A code of convergence. 0 indicates successful convergence. }
    \item{\code{posE}:}{Object of class \code{"numeric"}. Input argument  of \code{\link{fitPP.fun}}.}
    \item{\code{covariates}:}{Object of class \code{"matrix"}. Input argument  of \code{\link{fitPP.fun}}.} 
    \item{\code{fixed}:}{Object of class \code{"list"}. Input argument  of \code{\link{fitPP.fun}}.} 
    \item{\code{tit}:}{Object of class \code{"character"}. Input argument of \code{\link{fitPP.fun}}.} 
    \item{\code{tind}:}{Object of class \code{"logical"}. Input argument of \code{\link{fitPP.fun}}.}
    \item{\code{t}:}{Object of class \code{"numeric"}. Input argument of \code{\link{fitPP.fun}}.} 

  }
}
\section{Extends}{
Class \code{"\linkS4class{mle}"}, directly.}

\section{Methods}{
Most of the  S4 methods in \pkg{stats4}  for the S4-class \code{mle} can be used.  Also  a  \code{mle} method 
for  the generic function \code{\link{extractAIC}} and a version of the \code{\link{profile}}  
\code{\linkS4class{mle}}  method adapted to the  \code{\linkS4class{mlePP}} objects are available:
 \describe{
    \item{coef}{\code{signature(object = "mle")}}
    \item{logLik}{\code{signature(object = "mle")}}
    \item{nobs}{\code{signature(object = "mle")}}
    \item{show}{\code{signature(object = "mle")}}
    \item{summary}{\code{signature(object = "mle")}}
    \item{update}{\code{signature(object = "mle")}}
    \item{vcov}{\code{signature(object = "mle")}}
    \item{confint}{\code{signature(object = "mle")}}
    \item{extractAIC}{\code{signature(object = "mle")}}
    \item{profile}{\code{signature(fitted = "mlePP")}}
}

Some other generic functions related to fitted models, such as \code{AIC} or \code{BIC}, can  also  
be applied to \code{mlePP} objects. 
}

\note{Let us remind that, as in all the S4-classes,  the symbol @  must be used instead of $ to name the slots: mlePP@covariates, 
mlepp@lambdafit, etc.}


\seealso{\code{\link{fitPP.fun}},  \code{\linkS4class{mle}} }

\examples{
showClass("mlePP")
}

