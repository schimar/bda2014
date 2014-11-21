library(rjags)

## read and format data
fish<-read.table("fish.txt",header=T)
fishData<-list(y=fish$y,x=fish$x,N=length(fish$y))

## compile the model
fishModel<-jags.model(file="fishModel.txt",data=fishData,n.chains=3)
update(fishModel,n.iter=1000)
out<-coda.samples(model=fishModel,variable.names=c("beta","tau","mu"),n.iter=15000,thin=3)

## compare to OLS estiamtes
lmout<-lm(fish$y ~ fish$x)
