\name{stepAICmle.fun}
\alias{stepAICmle.fun}
\alias{checkdim}


\title{Choose  the best PP model by AIC in a stepwise algorithm}

\description{Performs stepwise model selection by AIC for Poisson proces models 
estimated by maximum likelihood. 

It calls the auxiliary function \code{checkdim} (not intended
for the users).}

\usage{stepAICmle.fun(ImlePP, covariatesAdd = NULL, startAdd = NULL, 
direction = "forward", ...)}


\arguments{

  \item{ImlePP}{A \code{\linkS4class{mlePP}}-class object; usually the output from \code{\link{fitPP.fun}}. 
It defines the initial model of the stepwise algorithm.  The fitted model cannot include fixed parameters.}

  \item{covariatesAdd}{Matrix   of the potential covariates to be added to the model; each column  must 
contain a covariate. In the 'forward'  and the 'both' directions, it is compulsory to assign a matrix to this argument.  
It is advisable to give names to the columns of 
this matrix (using \code{\link{dimnames}}) since,   they will be used  in the output. 
Otherwise the default names 'New Covariate i' are used.}

  \item{startAdd}{Optional. The vector of initial values  for the estimation  of  the coefficients 
of each potential covariate. If it  is NULL,  initial values  equal to 0 are used.}

  \item{direction}{Label indicating the direction of the algortihm: 'forward' (the default), 'backward' or 'both'.}

    \item{...}{Further arguments to pass to \code{\link{addAIC.fun}} and \code{\link{dropAIC.fun}}, 
for example the constant k for the AIC calculation}

}



\details{Three directions, forward, backward and both, are implemented. The initial model is given by 
 \code{ImlePP} and the algorithm stops 
 when none of the covariates eliminated from the model 
 or  added from the potential covariates set (argument \code{covariatesAdd} ) improves the model 
fitted in the  previous step, according to the AIC.  For the  'both'  and 'forward' directions, the argument \code{covariatesADD}
is compulsary, and the default NULL leads to an error.

In the 'both' direction,   'forward' and 'backward' steps are carried out alternatively. In the 'forward' direction, 
the initial model usually contains only the intercept.
}

\value{A \code{\linkS4class{mlePP}}-class object,  the fit of the final PP model selectecd by the algorithm.}

\references{

Casella, G. and Berger, R.L., (2002). \emph{Statistical inference.} Brooks/Cole.

Cebrian, A.C., Abaurrea, J. and Asin, J. (2015). {NHPoisson}: An {R} Package for 
Fitting and Validating Nonhomogeneous Poisson Processes. 
\emph{Journal of Statistical Software}, 64(6), 1-24.

Venables, W. N. and Ripley, B. D. (2002). \emph{Modern Applied Statistics with S.} Fourth edition. Springer.
}



\seealso{
\code{\link{addAIC.fun}}, \code{\link{dropAIC.fun}}, \code{\link{testlik.fun}} }

\examples{


data(BarTxTn)

BarEv<-POTevents.fun(T=BarTxTn$Tx,thres=318, 
	date=cbind(BarTxTn$ano,BarTxTn$mes,BarTxTn$dia))

#The initial model contains only the inercept
 mod1Bind<-fitPP.fun(covariates=NULL, posE=BarEv$Px, inddat=BarEv$inddat,
	tit='BAR  Intercept ', 	start=list(b0=1))
#the potential covariates
covB<-cbind(cos(2*pi*BarTxTn$dia/365), sin(2*pi*BarTxTn$dia/365), 
	BarTxTn$TTx,BarTxTn$Txm31,BarTxTn$Txm31**2)
dimnames(covB)<-list(NULL,c('cos','sin','TTx','Txm31', 'Txm31**2'))

bb<-stepAICmle.fun(ImlePP=mod1Bind, covariates=covB, startAdd=c(1,-1,0,0,0), 
direction='both')

}
