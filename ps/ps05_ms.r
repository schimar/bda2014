# ps 04

library(rjags)
###
# create data
bact <- read.csv("microbes.csv", header= T)


data <- list(y= bact$y, N= dim(bact)[1], temp= bact$temp, substrate= bact$sub)


model <- textConnection("model{
	for(i in 1:N){
		## normal likelihood on y
		y[i] ~ dnorm(mu[i], tau)
		## linear model with interaction and intercept
		mu[i] <- beta[1] + beta[2] * temp[i] + beta[3] * substrate[i] + beta[4] * temp[i] * substrate[i]
		## new data for model checking
		ynew[i] ~ dnorm(mu[i], tau)
		## residuals for model checking
		resid[i] <- y[i] - mu[i]
	}

	for(j in 1:4){
		## uninformative priors on reg coeff
		beta[j] ~ dnorm(0, 1e-6)
	}
	## prior on tau
	tau ~ dgamma(0.01, 0.01)

}")

####

bactModel <- jags.model(model, data= data, n.chains= 3)

update(bactModel, n.iter= 1000)

bactSamples <- coda.samples(model= bactModel, n.iter= 30000, variable.names= c("beta", "mu", "tau", "ynew", "resid"), thin= 3)



plot(bactSamples)

## 

est <- summary(bactSamples)[[2]]
n <- 75

## plot observed and predicted data (y and ynew)
plot(1:n,data$y,xlab="observation number",ylab="data value", ylim= c(-10, 150))
segments(1:n,est[156:230,1],1:n,est[156:230,5],lty=3)


## expected vs residual (mu vs resid)
plot(est[5:79,3], est[80:154,3], xlab= "expected value",ylab="residual")
abline(h=0, lty=2)

## r^2 (tau)

sigma <- 1/est[80,3]
r2 <- 1- sigma/var(bact$y) # 0.986 !!!


## cross-validation 

x <- sample(1:75, 15, replace=F)
data$y[x] <- NA

## compile model (remember to run the textConnection again!)
bactModel2 <- jags.model(model, data= data, n.chains=3)
update(bactModel2, n.iter= 1000)
bactSamples2 <- coda.samples(model= bactModel, n.iter= 500000, variable.names= c("beta", "mu", "tau","y", "ynew"), thin= 50)

## plot 
est <- summary(bactSamples2)[[2]]
ys <-est[81:155,]

plot(x,bact$y[x], ylim= c())
points(x, ys[x, 3], col= "darkgreen")
segments(x, ys[x, 1], x, ys[x, 5], col= "darkgreen")
