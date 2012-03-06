transfH.fun <-
function(obFPP)
{
	inddat<-obFPP$inddat
	posE<-obFPP$posE
	lambdafit<-obFPP$lambdafit
	lambdafitc<-lambdafit*inddat

	Ilambda <- inddat*cumsum(lambdafitc)
	posEH <- Ilambda[posE]-lambdafitc[posE]

	return(list(posEH=posEH,posE=posE, lambdafit=lambdafit, inddat=inddat))
}
