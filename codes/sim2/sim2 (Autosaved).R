library(parallel) ## Multicore Compuataion
library(smm)
library(MASS)
library(corpcor)
runsim1<-function(seednum, N=100,
                  lambda1=exp(seq(0,-6, length=15)),
                  lambda2=exp(c(0)),
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
  time1<-NA;fit.cv<-c()

  time1<-system.time(
  fit.cv<-cv.sparse.mediation.sgrplasso(x,m,y,tol=tol,K=K,max.iter=max.iter,
                                       lambda1 = lambda1,lambda2=lambda2,
                                       alpha=alpha,multicore=multicore,
                                       non.zeros.stop=max(150,v/4)))
  time2<-system.time(fit<-sparse.mediation.grplasso(x,m,y,tol=tol,max.iter=max.iter,
                                                    lambda1 =fit.cv$cv.lambda1,
                                                    lambda2=fit.cv$cv.lambda2,
                                                    alpha=fit.cv$cv.alpha,
                                                    non.zeros.stop=max(150,v/4)))
  return(list(time=time1,time2=time2,cvfit=fit.cv,fit=fit))
}

inst_run<-function(N=200, V=50,nrepeat=4,K=5, n.cores=16){
  pcmat=matrix(0,5,5)
  pcmat[1,]<-c(1,0,0.2894310,0.12851795,0.0176913)
  pcmat[2,3:5]<-c(0.5,0.4,0.3)
  pcmat[3:5,3:5]<-0.05
  pcmat=pcmat+t(pcmat)
  #diag(pcmat)<- 1
  #eigen(pcmat)$values
  diag(pcmat)<- -1
  #cov2cor(pseudoinverse(-pcmat))
  cmat=matrix(0,V+2,V+2);
  cmat[1:5,1:5]<-cov2cor(pseudoinverse(-pcmat))
	for (k in 1:10){cmat[(1:10)+5+(k-1)*10,(1:10)+5+(k-1)*10]<-0.5}
  cmat[1:10+5,1]<-0.5
  cmat[1,1:10+5]<-0.5
	diag(cmat)<- 1
#eigen(cmat)$values
  z<-mclapply( (1:nrepeat)*2 +20200513, function(num){runsim1(seednum=num,N=N, K=K,covmat= cmat)},mc.cores= n.cores)
  return(z)
}

for (V in c(800)){
  for( n in c(100,300)){
    print(c(n,V))
    eval(parse(text=paste('z_n',n,'_V',V,'=inst_run(N=',n,',V=',V,',nrepeat=100)',sep='')))
    eval(parse(text=paste("save(z_n",n,"_V",V,",file='z_n",n,"_V",V,".Rdata')",sep="")))
    }
}

