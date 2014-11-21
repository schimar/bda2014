## create a sequence of possible parameter values
p<-seq(from=0,to=1,by=0.01)
## calculate and plot P(y | p, n) = c * L(p | y, n)
lik<-dbinom(x=24,size=50,prob=p)
plot(p,lik,type='l',xlab="p",ylab="likelihood")

## add the analytical MLE for p = 0.48
abline(v=0.48,col="orange")

## you can approximate the MLE by calculating the likelihood on a dense grid, and finding the maximum value
p[order(lik)]
p[lik == max(lik)]

## now plot the log-likelihood surface and compare it to the above result
plot(p,log(lik),type='l',xlab="p",ylab="lnLikelihood")

## and now plot the relative likelihood
rlik<-lik / dbinom(x=24,size=50,prob=0.48)
plot(p,rlik,type='l',xlab="p",ylab="Likelihood")

## we can plot the (approximate) 5% likelihood region,\\

abline(v=p[min(which(rlik >= 0.05))])
abline(v=p[max(which(rlik >= 0.05))])
