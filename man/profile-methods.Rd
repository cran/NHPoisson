\name{profile-methods}
\docType{methods}
\alias{profile-methods}
\alias{profile,mlePP-method}

\title{ Method \code{mlePP} for Function \code{profile}}

\description{ Method for  generic function \code{\link{profile}} for  objects of the S4-class 
\code{\linkS4class{mlePP}}. It is almost  identical to the method \code{\linkS4class{mle}} for this function in \pkg{stats4}, 
 but  small changes  have to be done due to the differences  in the arguments of  the functions
\code{mle}  and \code{fitPP.fun}.  In order to profile an \code{mlePP} object,  its \code{vcov} slot cannot be missing. 
That means that if the function \code{fitPP.fun} is used to create the object, the argument \code{modCI=TRUE} must be used.
}

\section{Methods}{
\describe{
\item{\code{signature(fitted = "mlePP")}}{}
}}


\keyword{methods}

