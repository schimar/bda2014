normalGibbs<-function(y=NA,nsteps=1000,mu0=0,tau0=0.001,alpha=0.1,beta=0.001){
  ## calculate and store empirical mean and precision
  n<-length(y)
  ybar<-mean(y)
  v<-1/var(y)
  ## matrix to store MCMC samples
  mcmcsam<-matrix(NA,nrow=nsteps,ncol=2)
  colnames(mcmcsam)<-c("mu","tau")
  ## intialize mu by sampling one value of y
  mcmcsam[1,1]<-sample(y,1)
  ## gibbs sample to intialize tau
  varymu <- sum( (y - mcmcsam[1,1])^2)
  mcmcsam[1,2]<-rgamma(1,shape=n/2 + alpha,rate= 0.5 * varymu + beta)
  ## main loop for gibbs sample
  for (i in 2:nsteps){
    ## update mu | tau
    muy<-(tau0 * mu0 + n * ybar * mcmcsam[(i-1),2])/(tau0 + n * mcmcsam[(i-1),2])
    tauy<-tau0 + n * mcmcsam[(i-1),2]
    mcmcsam[i,1]<-rnorm(1,muy,sqrt(1/tauy))
    ## update tau | mu
    varymu <- sum( (y - mcmcsam[i,1])^2)
    mcmcsam[i,2]<-rgamma(1,shape=n/2 + alpha,rate= 0.5 * varymu + beta)
  }
  return(mcmcsam)
}
