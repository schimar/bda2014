## input gibbs function
source("normalGibbs.R")

## simulate data, which we will then analyze
y<-rnorm(n=25,mean=3,sd=0.5)

## run mcmc
out<-normalGibbs(y=y,nstep=5000,mu0=0,tau0=1e-5,alpha=0.01,beta=0.01)

## plot posterior distribution
plot(out[1:50,1],out[1:50,2],type='l',xlab="mu",ylab="tau")
plot(out[,1],out[,2],type='l',xlab="mu",ylab="tau")

plot(out[,1],out[,2],col="gray",xlab="mu",ylab="tau")
library(MASS)
d<-kde2d(out[,1],out[,2]) ## kernel density estimation
contour(d,add=T)
  
## marginal posterior distributions
hist(out[,1],xlab="mu")
hist(out[,2],xlab="tau")

## summaries
apply(out,2,quantile,probs=c(0.5,0.025,0.975))
apply(out,2,mean)
