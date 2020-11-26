
#library(devtools)
#install_github('seonjoo/smm')
rm(list=ls())
library(parallel) ## Multicore Compuataion
library(smm)
library(MASS)
library(corpcor)
 
runsim1<-function(seednum, N=100,
                  lambda1=exp(seq(0,-6, length=15)),
                  lambda2=exp(c(-1)),
                  alpha=c(0.5,0.75,0.95),
                  K=5,multicore=1,
                  max.iter=50,
                  tol=10^(-4),
                  covmat=diag(rep(1,10))){
  set.seed(seednum)
  print(seednum)
  
  v=ncol(covmat)-2
  tmpmat = mvrnorm(n=N, mu=rep(0,v+2),Sigma=covmat)
  x=as.matrix(tmpmat[,1])
  y=tmpmat[,2]
  m=as.matrix(tmpmat[,1:v + 2])
 
  xindx=which(p.adjust(apply(m,2,function(xx)cor.test(x,xx)$p.value), method='BH')<0.05)
  yindx=which(p.adjust(apply(m,2,function(xx)cor.test(y,xx)$p.value), method='BH')<0.05)
  xindx_adjust=xindx[xindx%in%yindx]
  
  xindx=which(apply(m,2,function(xx)cor.test(x,xx)$p.value)<0.05)
  yindx=which(apply(m,2,function(xx)cor.test(y,xx)$p.value)<0.05)
  xindx_p05=xindx[xindx%in%yindx]
 
  xindx=which(apply(m,2,function(xx)cor.test(x,xx)$p.value)<0.1)
  yindx=which(apply(m,2,function(xx)cor.test(y,xx)$p.value)<0.1)
  xindx_p1=xindx[xindx%in%yindx]
 
  xindx=which(apply(m,2,function(xx)cor.test(x,xx)$p.value)<0.2)
  yindx=which(apply(m,2,function(xx)cor.test(y,xx)$p.value)<0.2)
  xindx_p2=xindx[xindx%in%yindx]
  
#round(cor(x,m[,1:3]),2)
#  round(cor(y,m[,1:3]),2)
#  round(unlist(lapply(1:3,function(xx)pcor.test(m[,xx],y,cbind(m[,1:3][,-xx],x))$estimate)),3)
  time1<-NA;fit.cv<-c()

  time1<-system.time(
  fit.cv<-cv.sparse.mediation.sgrplasso(x,m,y,tol=tol,K=K,max.iter=max.iter,
                                       lambda1 = lambda1,lambda2=lambda2,
                                       alpha=alpha,multicore=multicore,
                                       non.zeros.stop=max(150,V/4)))
  time2<-system.time(fit<-sparse.mediation.grplasso(x,m,y,tol=tol,max.iter=max.iter,
                                                    lambda1 =fit.cv$cv.lambda1,
                                                    lambda2=fit.cv$cv.lambda2,
                                                    alpha=fit.cv$cv.alpha,
                                                    non.zeros.stop=max(150,V/4)))
  return(list(time=time1,time2=time2,cvfit=fit.cv,fit=fit,xindx_adjust=xindx_adjust,xindx_p05=xindx_p05,xindx_p1=xindx_p1, xindx_p2=xindx_p2))
}

B=100
inst_run<-function(N=200, V=50,nrepeat=4,K=5, n.cores=12){
	cmat=matrix(0,V+2,V+2);
	a=c(0.5,0.4,0.3)
	totef=0.5
	b=a*(totef+sqrt((1-totef^2)*(1-a^2)))
	cmat[1,2]<-totef
	cmat[1,1:3+2]<-a
	cmat[2,1:3+2]<-b
	cmat=cmat+t(cmat)
	diag(cmat)<-1
  
  z<-mclapply( 1:nrepeat +20200512, function(x){runsim1(seednum=x,N=N, K=K,covmat= cmat)},mc.cores= n.cores)
  return(z)
}

for (V in c(400,800)){
  for( n in c(100,300)){
    print(c(n,V))
    eval(parse(text=paste('z_n',n,'_V',V,'=inst_run(N=',n,',V=',V,',nrepeat=100)',sep='')))
    eval(parse(text=paste("save(z_n",n,"_V",V,",file='z_n",n,"_V",V,".Rdata')",sep="")))
      }
}

