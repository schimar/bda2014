# ps 04

library(rjags)
###
# create data

mice <- read.csv("virusData.csv", header= T)


data <- list(y= mice$y, N= dim(mice)[1], strain= mice$strains, line= mice$lines, popSize = mice$popsize-mean(mice$popsize), Nstrains = 5, Nlines= 5)


model <- textConnection("model{
	for(i in 1:N){
		## Poisson likelihhod on y
		y[i] ~ dpois(lambda[i])
		## link function and linear model
		log(lambda[i]) <- beta1[strain[i]] + beta2[line[i]] + beta3 * popSize[strain[i]] + epsi[i]
	}

	for(j in 1:Nstrains){		
		beta1[j] ~ dnorm(alpha[1], nu[1])
	}
	
	for(k in 1:Nlines){
		## hierarch model for lines # with szc 
		beta2[k] ~ dnorm(alpha[2], nu[2])
	}
	## beta2[Nlines]<- 1-sum(beta2[2:Nlines])

	for(k in 1:N){
		## prior on epsilon
		epsi[k] ~ dnorm(0, tau)
	}
	
	for(d in 1:2){
		## Priors on mean and precision for line and strain
			alpha[d] ~ dnorm(0, 1e-6)
			nu[d] ~ dgamma(0.01, 0.01)
	}
	
	## Additional uninformative priors
	beta3 ~ dnorm(0,1e-6)
	tau ~ dgamma(0.01,0.01)
}")


## load and initialize model
virusModel<-jags.model(model,data=data,n.chains=2)
## run model
update(virusModel,n.iter=4000)

out<-coda.samples(model=virusModel,variable.names=c("beta1", "beta2", "beta3", "lambda", "tau", "alpha", "nu", "epsi"),n.iter=80000,thin=4)

######

plot(out)
##
est <- summary(out)[[2]]











