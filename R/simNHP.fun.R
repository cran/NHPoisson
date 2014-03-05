simNHP.fun <-
function(lambda,n=100)
{

Tfinal<-length(lambda)
lambdacum<-cumsum(lambda)
distexp<-rexp(n,1)
posH<-cumsum(distexp)
posNHaux<-apply(as.matrix(posH),MARGIN=1,FUN=buscar, lambdacum)

nn<-n
while(posNHaux[nn]<Tfinal)
{
distexp<-rexp(n,1)
posH<-posH[n]+cumsum(distexp)
posNHaux<-c(posNHaux,apply(as.matrix(posH),MARGIN=1,FUN=buscar, lambdacum))
nn<-nn+n
}

posNH<-posNHaux[posNHaux<Tfinal]
return(list(posNH=posNH, lambda=lambda))
}
