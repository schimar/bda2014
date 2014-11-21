## this code explores the Binomial probability distribution in R

## view the help page
help(dbinom)

## plot the PMF for n = 20 and p = 0.2
y<-0:20
plot(y,dbinom(x=y,size=20,prob=0.2),type='h',ylab="probability(Y = y)", xlab="y")

## add mean and 0.025, 0.5, and 0.975 quantiles
abline(v=(20*0.2),col="red",lty=2)
abline(v=qbinom(p=c(0.025,0.5,0.975),size=20,prob=0.2),col="green",lty=2)

## draw 1000 random variates from a Binomial
y<-rbinom(n=1000,size=20,prob=0.2)
hist(y,xlim=c(0,20))

## calculate mean and quantiles of sample
mean(y)
quantile(x=y,probs=c(0.025,0.5,0.975))
