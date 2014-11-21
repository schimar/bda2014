# ps 01

#####  1.
# a)

a <- rbinom(20, 1, 0.5)

upper <- list()
lower <- list()

for (i in 1:1000){
	a <- rbinom(20, 1, 0.5)
	conf <- CI(a)
	upper[i] <- conf[1]
	lower[i] <- conf[3]
}

x <- 0

for (i in 1:1000){
	if (upper[i] < 0.5 | lower[i] > 0.5)
		x[i] <- "p not included"
	else 
		x[i] <- "p included"
}

table(x)	 # 955 p included, 45 not included

# How does this result compare with the expectation (which is based on the normal approximation)?


# b) 

b <- rbinom(20, 1, 0.02)

for (i in 1:1000){
	b <- rbinom(20, 1, 0.5)
	conf <- CI(b)
	upper[i] <- conf[1]
	lower[i] <- conf[3]
}


for (i in 1:1000){
	if (upper[i] <= 0.02 | lower[i] > 0.02)
		x[i] <- "p not included"
	else 
		x[i] <- "p included"
}

table(x)    # 6 p included, 994 not included

# Is the true coverage of the standard 95% confidence interval affected by the parameter value?
# absolutely. decreasing p moves the entire distribution towards 0.

# c) 

gbinom <- function(n, p, low=0, high=n,scale = F){
	sd <- sqrt(n * p * (1 - p))
	if(scale && (n > 10)) {
		low <- max(0, round(n * p - 4 * sd))
		high <- min(n, round(n * p + 4 * sd))
		}
	values <- low:high
	probs <- dbinom(values, n, p)
	plot(c(low,high), c(0,max(probs)), type = "n", xlab = "Possible Values", ylab = 
    "Probability", main = paste("Binomial Distribution \n", "n =", 
    n, ", p =", p))
	lines(values, probs, type = "h", col = 2)
	abline(h=0,col=3)
	return(invisible())
}

par(mfrow = c(2, 1))
gbinom(20, 0.5)
gbinom(20, 0.02)

dev.copy2pdf(file= "pmf_0.5_0.02")

# see tex file for ps01


##### 2. r
# Assume that the number of deaths in each hospital is independently and identically distributed

y <- c(6, 11, 9, 16, 12)
p <- seq(0, 40, 0.01)
# a)
# Write an equation that specifies the likelihood function for the death rate. Remember that the likelihood for a vector of identically and independently distributed random variables is the product of the likelihood for each random variable.

# b)

#lik6 <- dpois(x=y, lambda=6)
#lik11 <- dpois(x=y, lambda=11)
#lik16 <- dpois(x=y, lambda=16)
#lik12 <- dpois(x=y, lambda=12)
#lik9 <- dpois(x=y, lambda=9)

###
# par(mfrow=c(3,2))
#plot(p, lik6, type= 'l', ylab= 'Likelihood', xlab= 'y', main='H 1')
#abline(v=6, col='red')
#plot(p, lik9, type= 'l', ylab= 'Likelihood', xlab= 'y', main='H 3')
#abline(v=9, col='red')
#plot(p, lik16, type= 'l', ylab= 'Likelihood', xlab= 'y', main='H 4')
#abline(v=16, col='red')
#plot(p, lik12, type= 'l', ylab= 'Likelihood', xlab= 'y', main='H 5')
#abline(v=12, col='red')


lt<- dpois(x= y, lambda= p)

# b)
par(mfrow= c(1,2))
plot(p, lt, type= 'l', xlab= 'y', ylab= 'Likelihood')
p[1081] # order(lik_total) and then take the index for p
abline(v= 10.8, col= 'red')

# c)

#rell12 <- lt / dpois(x= y, lambda= 12)
#rlik_total  <- lik_total / dpois(x= y, lambda= )
#plot(p, rlik, type='l',xlab="p",ylab="Likelihood")
#dpois(x= y, lambda= p)


################################# Alberto

# create a vector for possible lambda values ##
lam <- seq(0.1,40,0.1)

# create a vector of the data 
dat <- c(6,11,9,16,12)

# compute poisson-likelihood surface 
r <- 1 
for(i in 1:5){
	r <- r*((lam^dat[i]*exp(1)^(-lam))/factorial(dat[i]))
	}

# 
plot(lam, r, type= 'l', xlab= "lambda", ylab= "Likelihood")

# print MLE 
lam[which.max(r)]

# relative likelihoods for 12 and 22
r[which(lam=='12')]/max(r)    # 0.7330568
r[which(lam=='22')]/max(r)    # 2.319802e-08 (pretty small)

# add lines to plot
abline(v= 10.8, col= "yellow")
abline(v= 12, col= 'blue')
abline(v= 22, col= 'orange')



r_ex <- exp(r)
which.max(r_ex) # 23
lam[23] # 2.2 ; is this MLE? 

####



# according to mle introduction (nyu):

mle <- 1/n * sum(y) # gives 10.8 




