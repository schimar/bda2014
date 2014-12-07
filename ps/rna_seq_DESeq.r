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
y <- c(rc19$fpkm, rc20$fpkm, rc21$fpkm) #, rt22$fpkm, rt23$fpkm, rt24$fpkm)
rep <-  c(rep(1, 60709), rep(2, 60709), rep(3, 60709)) #, rep(1, 60709), rep(2, 60709), rep(3, 60709))
#trtm <- c(rep(0, 182127), rep(1, 182127))
#data <- list(y= y, trtm = trtm, N= length(trtm), Ntrtm = 2, Nrep= 3)
data <- list(y= y, N= length(y), Nrep= 3)
##########

modDiffExp <- textConnection("model{

	## Poisson likelihood on y
	for(i in 1:N){
		for(j in 1:Nrep){
			y[i,j] ~ dpois(lambda[i])
		}
	}
	
	## independent priors on lambda (replicates)
	for(k in 1:Nrep){
		lambda[k] ~ dgamma(0.01, 0.01) # set up hyperpriors for alpha & beta?
	}

	## priors on    (treatment)
	

	## hyperpriors
	# alpha ~ 
}")



##########

compMod <- jags.model(modDiffExp, data= data, n.chains= 2)

update(compMod, n.iter= 2000)

out <- coda.samples(model= compMod, variable.names= c("lambda"), n.iter= 50000, thin= 3)





