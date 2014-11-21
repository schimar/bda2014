## lets set the true mean to 7 and s.d. to 1
mu<-7
sigma<-0.5
tau<-1/(sigma^2) ## precision 1/(sigma^2)

## now lets draw 20 random deviates
N<-20
y<-rnorm(n=N,mean=mu,sd=sigma)

## now lets visualize the posterior for mu (assume tau is known) with different normal priors
ybar<-mean(y)

mu0<-0
tau0<-1e-6
sd0<-sqrt(1/tau0)

## plot prior on two different scales
x<-seq(-25,25,0.1)
plot(x,dnorm(x=x,mean=mu0,sd=sd0),type='l')

x<-seq(4,9,0.01)
plot(x,dnorm(x=x,mean=mu0,sd=sd0),type='l')

## calculate posterior
postMn<-(N*tau*ybar + tau0*mu0)/(N*tau + tau0)
postPrec<-N*tau + tau0
postSd<-sqrt(1/postPrec)

## plot posterior and prior
x<-seq(4,9,0.01)
plot(x,dnorm(x=x,mean=postMn,sd=postSd),type='l')
lines(x,dnorm(x=x,mean=mu0,sd=sd0),lty=2)

## repeat with alternative priors, and add to the plot
tau0<-1e-3
sd0<-sqrt(1/tau0)

postMn<-(N*tau*ybar + tau0*mu0)/(N*tau + tau0)
postPrec<-N*tau + tau0
postSd<-sqrt(1/postPrec)

lines(x,dnorm(x=x,mean=postMn,sd=postSd),col="blue")
lines(x,dnorm(x=x,mean=mu0,sd=sd0),col="blue",lty=2)

tau0<-1
sd0<-sqrt(1/tau0)

postMn<-(N*tau*ybar + tau0*mu0)/(N*tau + tau0)
postPrec<-N*tau + tau0
postSd<-sqrt(1/postPrec)

lines(x,dnorm(x=x,mean=postMn,sd=postSd),col="green")
lines(x,dnorm(x=x,mean=mu0,sd=sd0),col="green",lty=2)

## now lets see what happens with very little data
N<-4
y<-rnorm(n=N,mean=mu,sd=sigma)
ybar<-mean(y)

mu0<-0
tau0<-1
sd0<-sqrt(1/tau0)

postMn<-(N*tau*ybar + tau0*mu0)/(N*tau + tau0)
postPrec<-N*tau + tau0
postSd<-sqrt(1/postPrec)

x<-seq(-5,15,0.01)
plot(x,dnorm(x=x,mean=postMn,sd=postSd),type='l')
lines(x,dnorm(x=x,mean=mu0,sd=sd0),lty=2)

## compare prior, sample, and posterior mean
mu0
ybar
postMn

## repeate N=400
N<-400
y<-rnorm(n=N,mean=mu,sd=sigma)
ybar<-mean(y)

mu0<-0
tau0<-1
sd0<-sqrt(1/tau0)

postMn<-(N*tau*ybar + tau0*mu0)/(N*tau + tau0)
postPrec<-N*tau + tau0
postSd<-sqrt(1/postPrec)

x<-seq(-5,15,0.01)
plot(x,dnorm(x=x,mean=postMn,sd=postSd),type='l')
lines(x,dnorm(x=x,mean=mu0,sd=sd0),lty=2)

## compare prior, sample, and posterior mean
mu0
ybar
postMn

## do it one more time with tau0 = 1e-6

N<-400
y<-rnorm(n=N,mean=mu,sd=sigma)
ybar<-mean(y)

mu0<-0
tau0<-1e-6
sd0<-sqrt(1/tau0)

postMn<-(N*tau*ybar + tau0*mu0)/(N*tau + tau0)
postPrec<-N*tau + tau0
postSd<-sqrt(1/postPrec)

x<-seq(-5,15,0.01)
plot(x,dnorm(x=x,mean=postMn,sd=postSd),type='l')
lines(x,dnorm(x=x,mean=mu0,sd=sd0),lty=2)

## compare prior, sample, and posterior mean
mu0
ybar
postMn
