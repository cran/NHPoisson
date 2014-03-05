graphResCov.fun <-
function(Xvar, nint,obFPP, h=NULL, typeRes='Pearson',
namX=NULL, indgraph='FALSE', tit='')
{
Xvar<-as.matrix(Xvar)
n<-dim(Xvar)[1]
if (is.null(tit)) tit<-obFPP$tit


mXres<-NULL
mXm<-NULL
mXpc<-NULL

nXv<-dim(Xvar)[2]
if (indgraph==TRUE)
{
iXv<-1
while (iXv<=nXv)
{
dev.new()
auxX<-graphResX.fun( X=Xvar[,iXv],nint=nint,obFPP=obFPP,
h=h,typeRes=typeRes, namX=namX[iXv])
mXres<-cbind(mXres, auxX$Xres)
mXm<-cbind(mXm, auxX$Xm)
mXpc<-cbind(mXpc, auxX$pc)
iXv<-iXv+1
}

}


if (indgraph==FALSE)
{
dev.new()
par(mfrow=c(2,2))

auxX<-graphResX.fun(X=Xvar[,1], nint=nint, obFPP=obFPP,
typeRes=typeRes, namX=namX[1])
mXres<-cbind(mXres, auxX$Xres)
mXm<-cbind(mXm, auxX$Xm)
mXpc<-cbind(mXpc, auxX$pc)

iXv<-2
igraf<-1
while (iXv<=nXv)
{
if ((igraf-4*floor(igraf/4))==0) #it check if the number of performed plots is multiple of 4
{
mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=0.8)
mtext(paste(typeRes, " residuals  ", sep=' '), outer = TRUE, line = -3,cex=0.7)
dev.new()
par(mfrow=c(2,2))
}

auxX<-graphResX.fun(X=Xvar[,iXv], nint=nint,obFPP=obFPP,
typeRes=typeRes, namX=namX[iXv])
mXres<-cbind(mXres, auxX$Xres)
mXm<-cbind(mXm, auxX$Xm)
mXpc<-cbind(mXpc, auxX$pc)

iXv<-iXv+1
igraf<-igraf+1
}
mtext(paste(" Model: ", tit, sep=''), outer = TRUE, line = -2,cex=1)
mtext(paste(typeRes, " residuals. Number of intervals:  ",nint, sep=' '), outer = TRUE, line = -3,cex=0.7)
}


return(list(mXres=mXres,mXm=mXm,mXpc=mXpc,nint=nint, obFPP=obFPP))

}
