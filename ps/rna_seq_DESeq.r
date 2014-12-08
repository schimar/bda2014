#library(DESeq2)
#library(GenomicRanges)
library(rjags)


##########

# create simulated data
#c19 <- read.delim("/home/schimar/Desktop/bda/ps/19c.counts", header= F)
#cnt19 <- c19[1:32678,]
#cnt19$V3 <- rnorm(n= 32678, 2000, 1600)

##########
# read express output of counts (in /home/mschilling/Desktop/retro/express_out) see bio.math.berkeley.edu/eXpress/manual.html for column descriptions
# retro ctrl
rc19 <- read.delim("19_results.xprs", header= T)
rc19 <- rc19[order(rc19$target_id),]
rc20 <- read.delim("20_results.xprs", header= T)
rc20 <- rc20[order(rc20$target_id),]
rc21 <- read.delim("21_results.xprs", header= T)
rc21 <- rc21[order(rc21$target_id),]

# m19_20 <- merge(x19, x20, by= "target_id")
#retro trtm
rt22 <- read.delim("22_results.xprs", header= T)
rt22 <- rt22[order(rt22$target_id),]
rt23 <- read.delim("23_results.xprs", header= T)
rt23 <- rt23[order(rt23$target_id),]
rt24 <- read.delim("24_results.xprs", header= T)
rt24 <- rt24[order(rt24$target_id),]


# preparing data for jags
# ctrl
y <- c(rc19$fpkm, rc20$fpkm, rc21$fpkm, rt22$fpkm, rt23$fpkm, rt24$fpkm)
rep <-  c(rep(1, 60709), rep(2, 60709), rep(3, 60709), rep(1, 60709), rep(2, 60709), rep(3, 60709))
trtm <- c(rep(0, 182127), rep(1, 182127))
transcript <- c(rc19$target_id, rc20$target_id, rc21$target_id, rt22$target_id, rt23$target_id, rt24$target_id)
#data <- list(y= y, trtm = trtm, N= length(trtm), Ntrtm = 2, Nrep= 3)
data <- list(y= y, trans= transcript, trtm= trtm, N= length(y), Ntrans= length(rc19$target_id))
##########

modDiffExp <- textConnection("model{

	## likelihood on y
	for(i in 1:N){
		y[i] ~ dpois(lambda[i])
		log(lambda[i]) <- beta1[trans[i]] + beta2[trans[i]] * trtm[i] 
	}
	
	## priors on betas 
	for(k in 1:Ntrans){
		beta1[k] ~ dnorm(mu, tau[1])
	}

	for(l in 1:Ntrans){
		beta2[l] ~ dnorm(0, tau[2])
	}
	beta2[Ntrans]<- sum(beta2[1:Ntrans-1)])	

	## hyperpriors
	for(m in 1:2){
		tau[m] ~ dgamma(0.01, 0.01)
	}
	mu ~ dnorm(0, 1e-6)
		
}")


# non-integer x = 8.672166Initializing model
# Deleting model

# Error in jags.model(modDiffExp, data = data, n.chains = 2) : 
#  Error in node y[1]
# Observed node inconsistent with unobserved parents at initialization.
# Try setting appropriate initial values.

##########

compMod <- jags.model(modDiffExp, data= data, n.chains= 2)

update(compMod, n.iter= 2000)

out <- coda.samples(model= compMod, variable.names= c("lambda"), n.iter= 50000, thin= 3)





