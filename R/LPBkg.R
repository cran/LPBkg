
Legj <- function(u,m){
  m <- min( length(unique(u))-1, m   )
  leg_js <-  slegendre.polynomials(m,normalized=TRUE)
  X <- polynomial.values(leg_js, u)
  S <- do.call(cbind,X)[,-1]}



denoise <- function(LP,n,method){
  M<-length(LP)
  LP2s <- sort(LP^2,decreasing=TRUE,index=TRUE)$x
  criterion <- rep(0,length(LP2s))
  if(method=="AIC"){ penalty <- 2}
  if(method=="BIC"){ penalty <- log(n)}
  criterion[1] <- LP2s[1] - penalty/n
  if(criterion[1]< 0){ return(rep(0,M))  }
  for(k in 2:M){
    criterion[k] <- sum(LP2s[1:k]) - k*penalty/n
  }
  LPsel<-LP
  LPsel[LP^2<LP2s[which.max(criterion)]] <- 0
  return(LPsel)
}


c_alpha2<-function(M,IDs,alpha=0.05,c_interval=c(1,10)){#m=4
  sleg.f <-  slegendre.polynomials(M,normalized=TRUE)
  Der<- polynomial.derivatives( sleg.f)
  Derf<-0
  for(j in (IDs+1)){
    Derf<-paste( Derf,"+(",toString(Der[[j]]),")^2")
  }
  sqDerf<-paste("sqrt(",Derf,")")
  integrand<-function(x)eval(parse(text=sqDerf))
  integrand<-Vectorize(integrand,"x")
  k0<-integrate(integrand,lower=0,upper=1)$value
  whichC<-function(C,k0){2*(1-pnorm(C,0,1))+k0*exp(-C^2/2)/pi-alpha}
  Calp<-uniroot(whichC, interval=c_interval,k0=k0)$root
  return( Calp)
}



BestM<-function(data,g, Mmax=20,range=c(min(data),max(data))){
  G<-function(y)integrate(g,lower=range[1],upper=y)$value
  G<-Vectorize(G)
  u<-G(data)
  n<-length(data)
  p.val<-c()
  S<- as.matrix(Legj(u,Mmax) )
  LP <- apply(S,FUN="mean",2)
  for(m in 1:Mmax){
    deviance<-n*sum(LP[1:m]^2)
    p.val[m]<-pchisq( deviance,m,lower.tail = FALSE)}
  names(p.val)<-1:Mmax
  return(list(pvals=p.val,minp=min(p.val),Msel=which.min(p.val)))}


dhatL2<-function(data,g, M=6, Mmax=NULL,
                 smooth=TRUE,criterion="AIC",
                 hist.u=TRUE,breaks=20,
                 ylim=c(0,2.5),range=c(min(data),max(data)),
                 sigma=2){
  bluetrans <- rgb(0, 0, 250, 50, maxColorValue=300)
  pinktrans <- rgb(30, 0, 10, 30, maxColorValue=300)
  G<-function(y)integrate(g,lower=range[1],upper=y)$value
  G<-Vectorize(G)
  u<-G(data)
  xx<-seq(range[1],range[2],by=0.01)
  uu<-G(xx)
  n<-length(u)
  S<- as.matrix(Legj(u=u,m=M))
  LP <- apply(S,FUN="mean",2)
  if(smooth==TRUE)LP <- denoise(LP,n=n,criterion)
  IDS0<-which(LP==0)
  IDS1<-which(LP!=0)
  z.L2<- approxExtrap(u,1 + S%*%LP,xout=c(0,u,1))
  dhat<- approxfun(z.L2$x,z.L2$y, rule=2)
  covs<-cov(S)*(n-1)/n^2
  covs[IDS0,]<-0
  covs[,IDS0]<-0
  vec1<-rep(0,length(u))
  for(j in 1:M){
    for(k in 1:M){
      vec1<-vec1+S[,j]*S[,k]*covs[j,k]
    }}
  z.SE.L2<- approxExtrap(u,sqrt(vec1),xout=c(0,u,1))
  SEdhat<- approxfun(z.SE.L2$x,z.SE.L2$y, rule=2)
  sigmas0<-rep(1/n,M)
  sigmas0[IDS0]<-0
  z.SE0.L2<- approxExtrap(u,sqrt( apply(S^2%*%sigmas0,1,sum) ),xout=c(0,u,1))
  SEdhatH0<- approxfun(z.SE0.L2$x,z.SE0.L2$y, rule=2)
  alpha=pnorm(sigma,0,1,lower.tail=FALSE)
  if(Mmax>1&smooth==FALSE)alpha=alpha/Mmax
  if(Mmax>1&smooth==TRUE)alpha=2*alpha/(Mmax+M*(M-1))
  qa<-c_alpha2(M,IDS1,alpha=alpha,c_interval=c(1,20))
  LBf1<-function(u)1-qa*SEdhatH0(u)
  UBf1<-function(u)1+qa*SEdhatH0(u)
  SE<-function(u)SEdhat(u)
  deviance<-n*sum(LP^2)
  pval<-pchisq(deviance,sum(LP!=0),lower.tail = FALSE)
  adj_pval<-NA
  if(Mmax>1&smooth==FALSE)adj_pval=pval*Mmax
  if(smooth==TRUE)adj_pval=pval*(Mmax+M*(M-1))/2
  if(hist.u==TRUE){
    oldpar <- par(mfrow=c(1,1),mar=c(5,5,1,1))
    on.exit(par(oldpar))
    par(mfrow=c(1,1),mar=c(5,5,1,1))
         hist(u,breaks=breaks,prob=TRUE,main=" ",ylim=ylim,
         ylab="Comparison density",xlab=expression(G(x)),cex.axis=1,cex.lab=1.4,xaxt='n' )
    lines(uu,dhat(uu),col="dodgerblue1",lwd=2)
    polygon(c(uu,rev(uu)),c(LBf1(uu),rev(UBf1(uu))),col=pinktrans, border = FALSE)
    polygon(c(uu,rev(uu)),c(dhat(uu)-SE(uu),rev(dhat(uu)+SE(uu))),col=bluetrans, border = FALSE)
    abline(h=1,lwd=2,lty=2,col="tomato3")
    Gx<-quantile(data,probs=seq(0,1,by=0.1))
    Axis(side=1, at=seq(0,1,by=0.1),
         labels=paste("G(",round(Gx,4),")",sep=""),
         cex.axis=1,col="black",col.lab="black" )
  }
  dhat_num<-function(x)dhat(G(x))
  dhat.x<-function(x)dhat_num(x)
  f<-function(x)g(x)*dhat(G(x))
  return(list(Deviance=deviance,Dev_pvalue=pval,
              Dev_adj_pvalue=adj_pval,dhat=dhat,
              kstar=sum(LP!=0),
              dhat.x=dhat.x,SE=SE,LBf1=LBf1,
              UBf1=UBf1,f=f,
              u=u,LP=LP,G=G
  ))}
