\name{dropAIC.fun}
\alias{dropAIC.fun}

\title{Calculate the AIC for all one-covariate deletions from the  current model}

\description{ This function fits all models  obtained from the current model by deleting one covariate 
(except the intercept),  and calculates their AIC value. 
It selects the best covariate to be deleted, according to the AIC value.}

\usage{dropAIC.fun(mlePP, modSim = FALSE,...)}

\arguments{
 \item{mlePP}{A \code{"\linkS4class{mlePP}"}-class object; usually the output from \code{\link{fitPP.fun}}. 
It defines the current model.  The fitted model cannot include fixed parameters.}


  \item{modSim}{Logical flag. If it is FALSE, information about the process is shown on the screen.
 For automatic selection processes, the option TRUE should be preferred.}

    \item{...}{Further arguments to pass to \code{\link{AIC}}, for example the constant k for the AIC calculation.}
}


\details{The definition of AIC uses constant k=2, but  a different value  k can be passed as
an additional argument. The best covariate to be deleted  is the one  whose deletion leads to  the model with the lowest AIC value
 and it  improves the current model if the new AIC  is lower than the current one. }


\value{ A list with the following components
  \item{AICadd }{Vector of the AIC values obtained from deleting each covariate  of the current model.}
  \item{posminAIC}{An integer indicating the number of the column  of the covariates matrix with the covariate leading to the minimum AIC.}
  \item{namecov}{Name of the covariate leading to the minimum AIC.}
  \item{AICcurrent}{AIC value of the current (initial) model.}
}

\references{

Casella, G. and Berger, R.L., (2002). \emph{Statistical inference.} Brooks/Cole.

Venables, W. N. and Ripley, B. D. (2002). \emph{Modern Applied Statistics with S.} Fourth edition. Springer.
}



\seealso{
\code{\link{addAIC.fun}}, \code{\link{stepAICmle.fun}}, \code{\link{LRTpv.fun}} }

\examples{

data(BarTxTn)

BarEv<-POTevents.fun(T=BarTxTn$Tx,thres=318, 
	date=cbind(BarTxTn$ano,BarTxTn$mes,BarTxTn$dia))

covB<-cbind(cos(2*pi*BarTxTn$dia/365), sin(2*pi*BarTxTn$dia/365), 
	BarTxTn$TTx,BarTxTn$Txm31,BarTxTn$Txm31**2)

dimnames(covB)<-list(NULL,c('cos','sin','TTx','Txm31', 'Txm31**2'))

mod1B<-fitPP.fun(covariates=covB, posE=BarEv$Px, inddat=BarEv$inddat,
	tit="BAR Tx; cos, sin, TTx, Txm31, Txm31**2", 
	start=list(b0=-100,b1=1,b2=10,b3=0,b4=0,b5=0))


aux<-dropAIC.fun(mod1B)
}
