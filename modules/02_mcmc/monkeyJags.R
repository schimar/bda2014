library(rjags)

## read in data and create a list
monkey<-read.table("monkeys.txt",header=T)
data<-list(y=monkey$y,trt=monkey$T,n=dim(monkey)[1])

## define the model
modPois<-textConnection("model{

    for(i in 1:n){
        y[i] ~ dpois(lambda[trt[i]+1])
    }
    
    lambda[1] ~ dgamma(0.01,0.01)
    lambda[2] ~ dgamma(0.01,0.01)
    dif<-lambda[1]-lambda[2]


}")

## compile the model
mymodel<-jags.model(modPois,data=data,n.chains=2)
## run the model, burnin the samples
update(mymodel,n.iter=1000)
outc<-coda.samples(model=mymodel,variable.names=c("lambda","dif"),n.iter=10000,thin=3)

## examine the output
summary(outc)
plot(outc)
effectiveSize(outc)
autocorr.plot(outc,lag.max=100)

## compare results to analytical solutions, alpha = a0 + sum y_i, beta = b0 + n
aT1<-0.01+sum(monkey$y[monkey$T == 1])
bT1<-0.01+length(monkey$y[monkey$T == 1]) # 
qgamma(c(0.025,0.25,0.5,0.75,0.975),shape=aT1,rate=bT1)

aT0<-0.01+sum(monkey$y[monkey$T == 0])
bT0<-0.01+length(monkey$y[monkey$T == 0])
qgamma(c(0.025,0.25,0.5,0.75,0.975),shape=aT0,rate=bT0)

## assess convergence
gelman.diag(x=outc,transform=TRUE)
gelman.plot(x=outc,transform=TRUE)

## prob dif > 0
mean(c(outc[[1]][,1],outc[[2]][,1]) > 0)
