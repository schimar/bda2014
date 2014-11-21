library(rjags)

## read in data and create a list
papayas<-read.table("papayas.txt",header=F)
data<-list(y=papayas,Nfr=dim(papayas)[1],Ntr=dim(papayas)[2])

## define the model
modHnorm<-textConnection("model{

    ## normal likelihood
    for(i in 1:Nfr){
        for(j in 1:Ntr){
            y[i,j] ~ dnorm(mu[j],tau)  # changed tau[j] to tau
        }
    }
    
    ## independent priors on taus
    for(j in 1:Ntr){
        tau ~ dgamma(0.01,0.01)    ## changed tau[j] to tau
    }
    
    ## conditional priors on mus
    for(j in 1:Ntr){
        mu[j] ~ dnorm(nu,kappa)   ## nu is the "grand" mean and kappa is the precision
    }
    
    ## hyperpriors
    nu ~ dnorm(0,1e-6)
    kappa ~ dgamma(0.01,0.01)
    
    sdtrs<-sqrt(1/kappa)

}")

## compile the model
mymodel<-jags.model(modHnorm,data=data,n.chains=2)
## run the model, burnin then samples
update(mymodel,n.iter=1000)
out<-coda.samples(model=mymodel,variable.names=c("mu","tau","nu","kappa","sdtrs"),n.iter=10000,thin=3)

## simple summaries
summary(out)
par(ask=TRUE)
plot(out)
