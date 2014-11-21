## Tue, Nov 18
####################
#####################
library(rjags)

bdata<-read.csv("beetledata.csv",header=T)
bcovar<-read.csv("beetlecovariates.csv",header=T)



data <- list(y= log(bdata$y/(1-bdata$y)), species = bdata$species, stand= bdata$stand, sd= bcovar$stand_density-mean(bcovar$stand_density), dsi= bcovar$drought_severity_index-mean(bcovar$drought_severity_index), N= length(bdata$y), Nstand= length(bcovar$stand_density), Nspecies= 5)


###################### model 
modBeetle <- textConnection("model{

	for(i in 1:N){
		## Normal likelihood on y
		y[i] ~ dnorm(mu[i], tau)
		## linear model 
		mu[i] <- beta1[species[i]] + beta2[stand[i]] + beta3 * sd[stand[i]] + beta4 * dsi[stand[i]]
	}
	for(j in 1:Nspecies){
		## hierarch model for species
		beta1[j] ~ dnorm(alpha[1], nu[1])
	}
	for(l in 1:Nstand){
		## hierarch model for stand intercepts
		beta2[l] ~ dnorm(alpha[2], nu[2])
	}
	for(k in 1:2){
		## priors on population mean and tau
		alpha[k] ~ dnorm(0, 1e-6)
		nu[k] ~ dgamma(0.01, 0.01)
	}
	## uninformative priors on the last 2 betas and tau
	beta3 ~ dnorm(0, 1e-6)
	beta4 ~ dnorm(0, 1e-6)
	tau ~ dgamma(0.01, 0.01)
}")



beetleModel <- jags.model(file= modBeetle, data= data, n.chains= 3) 

update(beetleModel, n.iter= 2000)


