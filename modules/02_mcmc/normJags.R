library(rjags)

## generate data and create a list
y<-rnorm(25,3,0.5)
data<-list(y=y,n=length(y))

## define the model
modNorm<-textConnection("model{

    for(i in 1:n){
        y[i] ~ dnorm(mu,tau)
    }
    
    mu ~ dnorm(0,1e-5)
    tau ~ dgamma(0.01,0.01)


}")

## compile the model
mymodel<-jags.model(modNorm,data=data,inits=list(mu=rnorm(1,3,2)),n.chains=2)
## run the model, burnin then samples
update(mymodel,n.iter=1000)
outc<-coda.samples(model=mymodel,variable.names=c("mu","tau"),n.iter=10000,thin=3)

## examine the output
summary(outc)
plot(outc)
effectiveSize(outc)
autocorr.plot(outc,lag.max=100)

## assess convergence
gelman.diag(x=outc) # if > 1.1 then unhappy; if < 1.1 then happy (in general, the closer to 1, the happier
gelman.plot(x=outc)





