library(rjags)

## read in data and create a list
opH<-read.csv("op2_hnv.csv")
opV<-read.csv("op2_vic.csv")

dataH<-list(y=opH$M,n=apply(opH,1,sum),N=dim(opH)[1])
dataV<-list(y=opV$M,n=apply(opV,1,sum),N=dim(opV)[1])

## define the model
model<-textConnection("model{

    ## binomial likelihood
    for(i in 1:N){
        y[i] ~ dbinom(p[i],n[i])
    }        
        
    ## conditional prior on the ps
    for(i in 1:N){
        p[i] ~ dbeta(alpha, beta)
    }
    
    ## hyperpriors
    alpha ~ dunif(0.1,1000)
    beta ~ dunif(0.1,1000)
    
    ## mean and sd among butterflies
    mn<-alpha/(alpha+beta)
    sd<-sqrt((alpha * beta)/((alpha+beta)^2*(alpha+beta+1)))

}")



## compile the model
#modH<-jags.model(model,data=dataH,inits=list(alpha=runif(1,1,5),beta=runif(1,1,5)),n.chains=2)
modH<-jags.model(model,data=dataH,inits=list(list(alpha=0.5,beta=3),list(alpha=2,beta=5)),n.chains=2)
## run the model, burnin then samples
update(modH,n.iter=1000)
outH1<-coda.samples(model=modH,variable.names=c("p","mn","sd","alpha","beta"),n.iter=10000,thin=3)

par(ask=TRUE)
plot(outH1)
effectiveSize(outH1)


## define the model
model<-textConnection("model{

    ## binomial likelihood
    for(i in 1:N){
        y[i] ~ dbinom(p[i],n[i])
    }        
        
    ## conditional prior on the ps
    for(i in 1:N){
        p[i] ~ dbeta(pi * V,(1 - pi) *V)
    }
    
    ## hyperpriors
    pi ~ dbeta(1,1)
    V ~ dunif(0.1,1000)

}")



## compile the model
modH<-jags.model(model,data=dataH,n.chains=2)
## run the model, burnin then samples
update(modH,n.iter=1000)
outH2<-coda.samples(model=modH,variable.names=c("p","pi","V"),n.iter=10000,thin=3)

par(ask=TRUE)
plot(outH)

## define the model
model<-textConnection("model{

    ## binomial likelihood
    for(i in 1:N){
        y[i] ~ dbinom(p[i],n[i])
    }        
        
    ## conditional prior on the ps
    for(i in 1:N){
        p[i] ~ dbeta(pi * V,(1 - pi) *V)
    }
    
    ## hyperpriors
    pi ~ dbeta(1,1)
    V ~ dunif(1,1000)

}")



## compile the model
modV<-jags.model(model,data=dataV,n.chains=2)
## run the model, burnin then samples
update(modV,n.iter=1000)
outV<-coda.samples(model=modV,variable.names=c("p","pi","V"),n.iter=10000,thin=3)

## probability of natal host preference
mean(c(outH2[[1]][,27],outH2[[2]][,27]) < 0.5)
mean(c(outV[[1]][,26],outV[[2]][,26]) > 0.5)
## 95ETPI on diff in preference
dif<-c(outH2[[1]][,27],outH2[[2]][,27]) - c(outV[[1]][,26],outV[[2]][,26])
quantile(dif,probs=c(0.5,0.025,0.975))

