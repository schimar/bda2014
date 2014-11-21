# ps 04


library(rjags)

## read in data and create a list
geno <- read.csv("/home/schimar/Desktop/bda/ps/genotypes.csv", header= T)

bigAs <- (geno[,1]*2)+geno[,2]
data <- list(y= bigAs, n= 2*apply(geno,1,sum), N= dim(geno)[1])

## define the model
model<-textConnection("model{

	## binomial likelihood on pi
	for(i in 1:N){
		y[i] ~ dbinom(p[i], n[i])
	}

	## conditional on the ps
	for(i in 1:N){
		p[i] ~ dbeta(theta, theta)
	}

	## hyperprior
	theta ~ dunif(0.5, 1000)

	## population size
	size <- theta/(4*mu)
	mu <- 0.0000001

}")

## compile the model

modG <- jags.model(model, data= data, n.chains= 2)

## run the model, burnin then samples
update(modG,n.iter=1000)

out <- coda.samples(model= modG, variable.names= c("p", "size", "theta"), n.iter= 10000, thin= 3)

par(ask= T)
plot(out)
effectiveSize(out)

# size 
mean(c(out[[1]][,101], out[[2]][,101])) # 1834767

# theta
mean(c(out[[1]][,102], out[[2]][,102])) 

par(ask=T)
plot(out)
effectiveSize(out)

plot(p, dbeta(p, mean(c(out[[1]][,102], out[[2]][,102])), mean(c(out[[1]][,102], out[[2]][,102]))), type= 'l', xlab= "allele frequency", ylab= "prop. of loci")

