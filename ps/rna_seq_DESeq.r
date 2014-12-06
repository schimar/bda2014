library(DESeq2)
library(GenomicRanges)


####

c19 <- read.delim("C268_19c.sam", header= F)
c20 <- read.delim("C268_20c.sam", header= F)
c21 <- read.delim("C268_20c.sam", header= F)
c22 <- read.delim("C268_20c.sam", header= F)
c23<- read.delim("C268_20c.sam", header= F)
c24 <- read.delim("C268_20c.sam", header= F)


##########
# create simulated data
c19 <- read.delim("/home/schimar/Desktop/bda/ps/19c.counts", header= F)

cnt19 <- c19[1:32678,]


cnt19$V3 <- rnorm(n= 32678, 2000, 1600)


# sample random subset of transcripts
s19 <- c19[sample(c19,100000),]

s20 <- c20[s19$V2]











