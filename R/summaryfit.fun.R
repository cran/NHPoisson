summaryfit.fun <-
function(fitPPobj)
{

	if (is.null(fitPPobj$beta)) stop('The element beta of the object fitPPobj is missing')
	if (is.null(fitPPobj$VARbeta)) stop('The element VARbeta of the object fitPPobj is missing')

	ii<-c(1:length(beta))
	SD_beta<-diag(fitPPobj$VARbeta)**0.5
	beta<-fitPPobj$beta
	auxmat<-cbind(beta, SD_beta)
	print(round(auxmat,4))
	return(auxmat)
}
