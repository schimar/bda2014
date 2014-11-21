library(rjags)

## read and format data
data(ToothGrowth)
toothData<-list(y=ToothGrowth$len,supp=as.numeric(ToothGrowth$supp),dose=ToothGrowth$dose-mean(ToothGrowth$dose),N=length(ToothGrowth$len))

## compile model
tooth<-jags.model(file="toothlmModelV2",data=toothData,n.chains=2)
update(tooth,n.iter=1000)
out<-coda.samples(model=tooth,n.iter=20000,variable.names=c("beta","tau","ynew","resid","mu"),thin=3)

## examine results, estimate treatment effect
plot(out)
summary(out)

## extract interval estimatesd
estt<-summary(out)[[2]]
n<-60

## plot observed and predicted data
plot(1:n,toothData$y,ylim=c(1,38),xlab="observation number",ylab="data value")
segments(1:n,est[125:184,1],1:n,est[125:184,5],lty=3)

## plot expected vs. residual
plot(est[4:63,3],est[64:123,3],xlab="expected value",ylab="redsidual")
abline(h=0,lty=2)

## Bayesian analogue to r^2
## This can be interpreted as the proportional reduction of uncertainty concerning the response variable Y achieved by incorporating the explanatory variables X in the model
## 1 - sigma^2/Sy^2c
sigma2<-1/est[124,3]
r2<-1 - sigma2/var(toothData$y)
## you could calculate this over the posterior distribution

## cross-validation, leave out data for 20%
x<-sample(1:60,12,replace=F)
toothData$y[x]<-NA
## compile model
tooth<-jags.model(file="toothlmModelV2",data=toothData,n.chains=2)
update(tooth,n.iter=1000)
out<-coda.samples(model=tooth,n.iter=20000,variable.names=c("beta","tau","ynew","y"),thin=3)

## check with cross validation
estt<-summary(out)[[2]]
plot(x,ToothGrowth$len[x],ylim=c(1,38))
points(x,estt[x+4,3],col="darkgreen")
segments(x,estt[x+4,1],x,estt[x+4,5],col="darkgreen")
cor(ToothGrowth$len[x],est[x+4,3])
