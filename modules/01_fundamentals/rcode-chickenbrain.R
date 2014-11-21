## number of chickens in control and treatment
nc<-32
nt<-36

## sample means
xc<-1.013
xt<-1.173

## known precision
tau<-25

## prior on mu_c and mu_t
mu0<-0
tau0<-1e-6

## means of posterior distributions
muc<-(nc * tau * xc + tau0 * mu0)/(nc * tau + tau0)
mut<-(nt * tau * xt + tau0 * mu0)/(nt * tau + tau0)

## precisions of the posterior distributions
tauc<-nc * tau + tau0
taut<-nt * tau + tau0

## plot posterior distributions
mus<-seq(0.5,1.5,0.01)
plot(mus,dnorm(mus,mean=muc,sd=sqrt(1/tauc)),type='l')
plot(mus,dnorm(mus,mean=mut,sd=sqrt(1/taut)),type='l')

## compute summaries of posterior
qnorm(p=c(0.025,0.5,0.975),mean=muc,sd=sqrt(1/tauc))
qnorm(p=c(0.025,0.5,0.975),mean=mut,sd=sqrt(1/taut))

## samples from posteriors
samc<-rnorm(n=5000,mean=muc,sd=sqrt(1/tauc))
samt<-rnorm(n=5000,mean=mut,sd=sqrt(1/taut))

## note similarity to analytical quantiles
quantile(samc,probs=c(0.025,0.5,0.975))
quantile(samt,probs=c(0.025,0.5,0.975))

## samples make it easy to compute posteriors for derived parameters, such as the difference in means
dif<-samt - samc

## you can summarize that posterior too, isn't that cool!
quantile(dif,probs=c(0.025,0.5,0.975))
mean(dif > 0)
