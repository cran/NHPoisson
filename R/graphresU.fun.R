graphresU.fun <-
function(unires, posE, Xvariables=NULL, namXv=NULL, flow=0.5, tit='',
	addlow=TRUE,indgraph=FALSE)
{
	tt<-c(1:length(unires))
	if (is.null(Xvariables)==FALSE) Xvariables<-as.matrix(Xvariables)
	if (is.null(Xvariables)==TRUE) nXv<-0
	else nXv<-dim(Xvariables)[2]
	Xvariablest<-as.matrix(Xvariables[posE,])
	n<-length(unires)

	if (indgraph==FALSE)
	{
	dev.new()
	par(mfrow=c(2,2))
	plot(tt, unires, cex = 1, xlab = "index", 
			ylab = "uniform residuals", type = "n")
	points(tt, unires, cex = 0.3,pch=16)
	if (addlow==TRUE)
	{	aux<-lowess(tt,unires,f=flow)
		lines(aux$x,aux$y)
	}

	iXv<-1
	igraf<-1

	while (iXv<=nXv)
	{
		if ((igraf-4*floor(igraf/4))==0) #it check if the number of performed plots is multiple of 4
		{
			mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)
			dev.new()
			par(mfrow=c(2,2))
		}
		plot(Xvariablest[,iXv], unires, cex = 1, xlab = namXv[iXv], 
			ylab = "uniform residuals", type = "n")
		points(Xvariablest[,iXv], unires, cex = 0.3,pch=16)
		if (addlow==TRUE)
		{	aux<-lowess(Xvariablest[,iXv],unires,f=flow)
			lines(aux$x,aux$y)
		}
		igraf<-igraf+1
		iXv<-iXv+1
	}
	mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)

	res<-unires
	lagres<-c(NA, res[1:(length(res)-1)])
	dev.new()
	par(mfrow=c(2,2))
	plot(lagres,res,ylab='residuals (t)', xlab='residuals(t-1)', cex=0.5)
	obreg<-lm(res~lagres, na.action=na.exclude)
	lines(lagres,predict(obreg))
	corr<-cor.test(res,lagres,na.action='remove')
	title(sub=paste(' Pearson coef.: ',round(corr$estimate,2),'. P-value: ',round(corr$p.value,3),sep=''), cex=0.7)

	auxacf<-acf(res, plot=T, xlab='lag',main='')

	ret<-10*log10(n) 
	a<-NULL
	indret<-c(1:ret)
	for (i in indret)
	{
		a[i]<-Box.test(res,type='Ljung-Box',lag=i)$p.value
	}
	plot(indret, a, pch=16, xlab='lag', ylab='p-values', cex=0.5, ylim=c(0,1))
	abline(h=0.05,col='red')
	title(sub=paste('Ljung-Box p-values ',sep=''), cex=0.7)
	mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)

	qqPlot(unires,dist='unif')
	ksres<-ks.test(unires, y = "punif", min = 0, max = 1)
	title(sub=paste('KS (uniform) p-value: ',round(ksres$p.value,3),sep=''),cex=0.7)
	mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)
	}

	if(indgraph==TRUE)
	{
	dev.new()
	plot(tt, unires, cex = 1, xlab = "index", 
			ylab = "uniform residuals", type = "n")
	points(tt, unires, cex = 0.3,pch=16)
	if (addlow==TRUE)
	{	aux<-lowess(tt,unires,f=flow)
		lines(aux$x,aux$y)
	}
	mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)

	iXv<-1
	while (iXv<=nXv)
	{
		dev.new()
		plot(Xvariablest[,iXv], unires, cex = 1, xlab = namXv[iXv], 
			ylab = "uniform residuals", type = "n")
		points(Xvariablest[,iXv], unires, cex = 0.3,pch=16)
		if (addlow==TRUE)
		{	aux<-lowess(Xvariablest[,iXv],unires,f=flow)
			lines(aux$x,aux$y)
		}
		mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)

		iXv<-iXv+1
	}

	res<-unires
	lagres<-c(NA, res[1:(length(res)-1)])
	dev.new()
	plot(lagres,res,ylab='residuals (t)', xlab='residuals(t-1)', cex=0.5)
	obreg<-lm(res~lagres, na.action=na.exclude)
	lines(lagres,predict(obreg))
	corr<-cor.test(res,lagres,na.action='remove')
	mtext(paste(" Model: ", tit,'. Pearson coef.: ',round(corr$estimate,2),
			'. P-value: ',round(corr$p.value,3), sep=''), outer = TRUE, line = -2,cex=1)

	dev.new()
	auxacf<-acf(res,plot=T,main='',xlab='lag')
	mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)

	dev.new()
	ret<-10*log10(n) 
	a<-NULL
	indret<-c(1:ret)
	for (i in indret)
	{
		a[i]<-Box.test(res,type='Ljung-Box',lag=i)$p.value
	}
	plot(indret, a, pch=16, xlab='lag', ylab='p-values', cex=0.5, ylim=c(0,1))
	abline(h=0.05,col='red')
	mtext(paste('Ljung-Box p-values ',sep=''), outer = TRUE, line = -2,cex=1)

	dev.new()
	qqPlot(unires,dist='unif')
	ksres<-ks.test(unires, y = "punif", min = 0, max = 1)
	mtext(paste(" Model: ", tit,'. KS (uniform) p-value: ',round(ksres$p.value,3), sep=''), outer = TRUE, line = -2,cex=1)

	}
	return(NULL)
}
