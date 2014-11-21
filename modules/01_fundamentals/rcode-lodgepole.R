## read data into R
y<-read.table("data-lodgepole.txt",header=FALSE)

## total number of trees and sample size
sy<-sum(y)
n<-20

## set priors
a0<-0.001
b0<-0.001

## determine parameters of the posterior distribution for lambda
a<-sy + a0
b<-n + b0
    
## plot and summarize the posterior probability distribuiton for lambda
plot(seq(0,200,0.5),dgamma(seq(0,200,0.5),shape=a,rate=b),type='l')
mn<-a/b
qgamma(p=c(0.025,0.975),shape=a,rate=b)

## summaries of posterior for population size, a derived parameter
nTrees<- 40000 * a/b
40000 * qgamma(p=c(0.025,0.975),shape=a,rate=b)
