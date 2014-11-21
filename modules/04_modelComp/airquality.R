## Combine data in list
data(airquality)
airData<-list(y=airquality$Ozone,temp=airquality$Temp-mean(airquality$Temp),wind=airquality$Wind-mean(airquality$Wind),N=length(airquality$Ozone))

## Read in model
airModel<-jags.model(file="airmodel",data=airData,n.chains=3)
## Run MCMC with JAGS
update(airModel,n.iter=1000)
airSamples<-coda.samples(model=airModel,n.iter=30000,variable.names=c("beta","mu", "tau","ynew", "resid"),thin=2)

## extract interval estimatesd
est <- summary(airSamples)[[2]]
n <- length(airquality$Ozone)

## plot observed and predicted data (y and ynew)
plot(1:n,airData$y,xlab="observation number",ylab="data value", ylim= c(-40, 170))
segments(1:n,est[312:464,1],1:n,est[312:464,5],lty=3)

## plot expected vs. residual (mu vs resid)
plot(est[5:157,3],est[158:310,3],xlab="expected value",ylab="residual")
abline(h=0,lty=2)


## Bayesian analogue to r^2
## This can be interpreted as the proportional reduction of uncertainty concerning the response variable Y achieved by incorporating the explanatory variables X in the model
## 1 - sigma^2/Sy^2c
sigma2<-1/est[311,3]
r2<-1 - sigma2/var(airData$y, na.rm=T)
## this could also be done in the model (see Zach's model)

## you could calculate this over the posterior distribution

## cross-validation, leave out data for 20%
x <- sample(1:153, 30, replace=F)
airData$y[x] <- NA

## compile model
airModel<-jags.model(file="airmodel",data=airData,n.chains=3)
## Run MCMC with JAGS
update(airModel,n.iter=1000)
airSamples<-coda.samples(model=airModel,n.iter=30000,variable.names=c("beta","mu", "tau","ynew", "resid"),thin=2)


## check with cross validation
est<-summary(airSamples)[[2]]
plot(x,airquality$Ozone[x],ylim=c(-20, 120))
points(x,est[x+4,3],col="darkgreen")
segment(x,est[x+4,1],x,est[x+4,5],col="darkgreen") # I need to track y, and then take those values[x] (didn't here, therefore the plot looks wacky)
cor.test(airquality$Ozone[x],est[x+4,3], na.rm=T)
