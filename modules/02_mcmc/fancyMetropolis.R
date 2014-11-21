normalMetropolis<-function(y=NA,nsteps=1000,mu0=0,tau0=0.001,alpha=0.01,beta=0.01,pmu=0.1,ptau=0.2){
  ## matrix to store MCMC samples
  mcmcsam<-matrix(NA,nrow=nsteps,ncol=2)
  colnames(mcmcsam)<-c("mu","tau")
  ## intialize mu by sampling one value of y
  mcmcsam[1,1]<-sample(y,1)
  ## initialize tau by sampling a value similar to the empirical sd
  s<-sd(y)
  mcmcsam[1,2]<-1/(runif(1,s-(0.1*s),s+(0.1*s)))^2
  ## main loop for gibbs sample, start plot
  ybar<-mean(y)
  vary<-var(y)
  plot(c(ybar-0.23*ybar,ybar+0.23*ybar),c(0,1/vary+1.2*(1/vary)),xlab="mu",ylab="tau",type='n')
  for (i in 2:nsteps){
    ## update mu, random walk, uniform proposal
    m<-mcmcsam[i-1,1]
    mstar<-m+runif(1,-pmu,pmu)
    
    points(c(m,mstar),rep(mcmcsam[i-1,2],2),type='l',col='gray');Sys.sleep(0.3)
    
    pr<-sum(dnorm(x=y,mean=m,sd=sqrt(1/mcmcsam[i-1,2]),log=TRUE)) + dnorm(x=m,mean=mu0,sd=1/sqrt(tau0),log=TRUE)
    prstar<-sum(dnorm(x=y,mean=mstar,sd=sqrt(1/mcmcsam[i-1,2]),log=TRUE)) + dnorm(x=mstar,mean=mu0,sd=1/sqrt(tau0),log=TRUE)
    ## caclulate Metropolis ratio
    if ((prstar-pr) > log(runif(1,0,1))){
        mcmcsam[i,1]<-mstar
        points(c(m,mstar),rep(mcmcsam[i-1,2],2),type='l',col='black');Sys.sleep(0.3)
    }    
    else {
        mcmcsam[i,1]<-m
        points(mstar,mcmcsam[i-1,2],col="red",pch=4);Sys.sleep(0.2)
        points(c(m,mstar),rep(mcmcsam[i-1,2],2),type='l',col='white')
        points(mstar,mcmcsam[i-1,2],col="white",pch=4)
        
    }
    ## update tau, random walk, uniform proposal
    tau<-mcmcsam[i-1,2]
    taustar<-tau+runif(1,-ptau,ptau)
    if (taustar < 0){## cannot have negatuve precision
        taustar<-tau
    }
    
    points(rep(mcmcsam[i,1],2),c(tau,taustar),type='l',col='gray');Sys.sleep(0.3)

    
    pr<-sum(dnorm(x=y,mean=mcmcsam[i,1],sd=sqrt(1/tau),log=TRUE)) + dgamma(x=tau,shape=alpha,rate=beta,log=TRUE)
    prstar<-sum(dnorm(x=y,mean=mcmcsam[i,1],sd=sqrt(1/taustar),log=TRUE)) + dgamma(x=taustar,shape=alpha,rate=beta,log=TRUE)
    ## caclulate Metropolis ratio
    if ((prstar-pr) > log(runif(1,0,1))){
        mcmcsam[i,2]<-taustar
        points(rep(mcmcsam[i,1],2),c(tau,taustar),type='l',col='black');Sys.sleep(0.3)

    }    
    else {
        mcmcsam[i,2]<-tau
        points(mcmcsam[i,1],taustar,col="red",pch=4);Sys.sleep(0.2)
        points(rep(mcmcsam[i,1],2),c(tau,taustar),type='l',col='white')
        points(mcmcsam[i,1],taustar,col="white",pch=4)

    }
  }
  return(mcmcsam)
}
