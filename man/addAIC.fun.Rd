\name{addAIC.fun}
\alias{addAIC.fun}

\title{Calculate the AIC  for all one-covariate additions to  the current model}

\description{ This function fits all models that differ from the current model by adding a single covariate from 
those supplied,  and calculates their AIC value. It selects the best covariate to be added to the model, 
 according to the AIC.}

\usage{addAIC.fun(mlePP, covariatesAdd, startAdd = NULL, modSim = FALSE,...)}


\arguments{
  \item{mlePP}{A \code{"\linkS4class{mlePP}"}-class object; usually the output from \code{\link{fitPP.fun}}. 
	It defines the current model.  The fitted model cannot include fixed parameters.}

  \item{covariatesAdd}{Matrix   of the potential covariates to be added to the model; each column  must 
contain a covariate.}

  \item{startAdd}{Optional. The vector of initial values  for the estimation  algorithm of  the coefficients 
of each potential covariate. If it  is NULL,  initial values  equal to 0 are used. Remark
that in contrast to argument \code{start} of \code{\link{fitPP.fun}}, \code{startAdd} is a  numeric vector not a list.}

  \item{modSim}{Logical flag. If it is FALSE, information about the process is shown on the screen.
 For automatic selection processes, the option TRUE should be preferred.}

    \item{...}{Further arguments to pass to \code{\link{AIC}}, for example the constant k for the AIC calculation.}
}

\details{The definition of AIC uses constant k=2, but  a different value  k can be passed as
an additional argument. The best covariate to be added  is the one which leads to  the model with the lowest AIC value
 and it  improves the current model if the new AIC is lower than the current one. }


\value{ A list with the following components
  \item{AICadd }{Vector of the AIC values obtained from  adding to the current model each covariate in  \code{covariatesAdd}.}
  \item{posminAIC}{An integer indicating the number of the column  of covariatesAdd with the covariate leading to the minimum AIC.}
  \item{namecov}{Name of the covariate leading to the minimum AIC.}
  \item{AICcurrent}{AIC value of the current (initial) model.}
  \item{newCoef}{A (named) list with   the  initial value for the coefficient
of the best covariate to be added. It is used in \code{\link{stepAICmle.fun}}.}

}





\seealso{
\code{\link{dropAIC.fun}}, \code{\link{stepAICmle.fun}}, \code{\link{LRTpv.fun}} }

\examples{

data(BarTxTn)

BarEv<-POTevents.fun(T=BarTxTn$Tx,thres=318, 
	date=cbind(BarTxTn$ano,BarTxTn$mes,BarTxTn$dia))

#The initial model contains only the intercept
 mod1Bind<-fitPP.fun(covariates=NULL, posE=BarEv$Px, inddat=BarEv$inddat,
	tit='BAR  Intercept ', 	start=list(b0=1))
#the potential covariates
covB<-cbind(cos(2*pi*BarTxTn$dia/365), sin(2*pi*BarTxTn$dia/365), 
	BarTxTn$TTx,BarTxTn$Txm31,BarTxTn$Txm31**2)
dimnames(covB)<-list(NULL,c('cos','sin','TTx','Txm31', 'Txm31**2'))

aux<-addAIC.fun(mod1Bind, covariatesAdd=covB)
}
