## read and format data
data(ToothGrowth)
toothData<-list(y=ToothGrowth$len,supp=as.numeric(ToothGrowth$supp),dose=ToothGrowth$dose,N=length(ToothGrowth$len))

## compile model
tooth<-jags.model(file="toothlmModel",data=toothData,n.chains=2)
update(tooth,n.iter=1000)
out<-coda.samples(model=tooth,n.iter=20000,variable.names=c("beta","tau"),thin=3)

## examine results, estimate treatment effect
plot(out)
summary(out)
quantile(c(out[[1]][,1],out[[2]][,1])-c(out[[1]][,2],out[[2]][,2]),probs=c(0.5,0.025,0.975))

## covariate centering
toothData<-list(y=ToothGrowth$len,supp=as.numeric(ToothGrowth$supp),dose=ToothGrowth$dose-mean(ToothGrowth$dose),N=length(ToothGrowth$len))

## compile model
tooth<-jags.model(file="toothlmModel",data=toothData,n.chains=2)
update(tooth,n.iter=1000)
out<-coda.samples(model=tooth,n.iter=20000,variable.names=c("beta","tau"),thin=3)

