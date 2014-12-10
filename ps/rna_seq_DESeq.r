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
y <- c(rc19$tot_counts, rc20$tot_counts, rc21$tot_counts, rt22$tot_counts, rt23$tot_counts, rt24$tot_counts)
#rep <-  c(rep(1, 60709), rep(2, 60709), rep(3, 60709), rep(1, 60709), rep(2, 60709), rep(3, 60709))
trtm <- c(rep(0, 182127), rep(1, 182127))
transcript <- c(rc19$target_id, rc20$target_id, rc21$target_id, rt22$target_id, rt23$target_id, rt24$target_id)
#
data <- list(y= y, trans= transcript, trtm= trtm, N= length(y), Ntrans= length(rc19$target_id))

# subset for initial run
suby <- c(rc19$tot_counts[1:100], rc20$tot_counts[1:100], rc21$tot_counts[1:100], rt22$tot_counts[1:100], rt23$tot_counts[1:100], rt24$tot_counts[1:100])
subtrtm <- c(rep(0, 300), rep(1, 300))
subtrans <- c(rc19$target_id[1:100], rc20$target_id[1:100], rc21$target_id[1:100], rt22$target_id[1:100], rt23$target_id[1:100], rt24$target_id[1:100])
#
subdata <- list(y= suby, trans= subtrans, trtm= subtrtm, N= length(suby), Ntrans= 100)





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

	for(l in 1:Ntrans-1){
		beta2[l] ~ dnorm(0, tau[2])
	}
	beta2[Ntrans]<- -sum(beta2[1:(Ntrans-1)])	

	## hyperpriors
	for(m in 1:2){
		tau[m] ~ dgamma(0.01, 0.01)
	}
	mu ~ dnorm(0, 1e-6)
		
}")




##########

compMod <- jags.model(modDiffExp, data= data, n.chains= 2)

update(compMod, n.iter= 6000)
statK <- as.data.frame(est[[1]])
out <- coda.samples(model= compMod, variable.names= c("beta1", "beta2", "lambda", "tau", "mu"), n.iter= 10000, thin= 3)


# subset of 1000 random transcripts

sub19 <- droplevels(rc19[sample(nrow(rc19), 1000),])
sub20 <- droplevels(rc20[sub19$target_id,])
sub21 <- droplevels(rc21[sub19$target_id,])
sub22 <- droplevels(rt22[sub19$target_id,])
sub23 <- droplevels(rt23[sub19$target_id,])
sub24 <- droplevels(rt24[sub19$target_id,])

sub1000 <- c(sub19$tot_counts, sub20$tot_counts, sub21$tot_counts, sub22$tot_counts, sub23$tot_counts, sub24$tot_counts)
subtrtm <- c(rep(0, 3000), rep(1, 3000))
subtrans <- rep(sub19$target_id,6)
#
subdata1000 <- list(y= sub1000, trans= subtrans, trtm= subtrtm, N= length(sub1000), Ntrans= 1000)

## jags model w/ 1000

subMod1000 <- jags.model(modDiffExp, data= subdata1000, n.chains= 2)

update(subMod1000, n.iter= 6000)

out1000 <- coda.samples(model= subMod1000, variable.names= c("beta1", "beta2", "lambda", "tau", "mu"), n.iter= 60000, thin= 3)

est <- summary(out1000)

quaK <- as.data.frame(est[[2]])
statK <- as.data.frame(est[[1]])
colnames(quaK) <- c("q2.5", "q25", "q50", "q75", "q97.5")
beta1 <- quaK[1:1000,]
beta1$trans <- subdata1000$trans[1:1000]
beta1 <- beta1[order(beta1$trans),]
beta2 <- quaK[1001:2000,]
beta2$trans <- subdata1000$trans[1:1000]
beta2 <- beta2[order(beta2$trans),]   

beta2$exp <- rep(1, 1000)
beta2$exp[beta2$q2.5 >0] <- 2
beta2$exp[beta2$q97.5 < 0] <- 3

# proportions
# beta1 
table(beta1$q2.5 > 0)
table(beta1$q97.5 < 0)

# beta2 
table(beta2$q2.5 > 0)
table(beta2$q97.5 < 0)

overBeta2median <- beta2$q50[which(beta2$q2.5 >0)]
underBeta2median <- beta2$q50[which(beta2$q97.5 <0)]
###
overExpTrans <- sort(subdata1000$trans[which(beta2$q2.5 >0)])
underExpTrans <- sort(subdata1000$trans[which(beta2$q97.5 <0)])

write.csv(subdata1000$trans, "sub1000transcripts.csv")
write.csv(overExpTrans, "trans1000_overExp.csv")
write.csv(underExpTrans, "trans1000_underExp.csv")

# plot 

#plot(beta2$trans, beta2$q50, col= beta2$exp)
#abline(h=0)
#points(beta2$q50[beta2$q2.5>0], col= "green", pch= 20)
#points(beta2$q50[beta2$q97.5<0], col= "red", pch= 20)

# barplot
barx <- barplot(beta2$q50, col= as.vector(beta2$exp), ylim= c(-8, 15), border= F, ylab= "Estimated abundance", xlab= "Transcripts")
box()


# meh, these are way too big...
error.bar <- function(x, y, upper, lower, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    		stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
    }



