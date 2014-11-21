## Creating and indexing vectors
x<-seq(from=0,to=30,by=0.5)
length(x)
x[12]
x[12:15]
x[c(1,3,19:21)]
which(x == 33)
x[which(x == 33)]
x == 33
x > 5
x > 5 & x < 13
x[x > 5 & x < 13]
y<-seq(from=-30,to=0,by=0.5)
x[y > -5]
x[y == -3]

## Operations with vectors and matrixes
sam<-rbinom(n=100,size=50,prob=0.3)
phat<-sam/50
x<-1:100
y<-1:100
x * y
x<-1:5
x * y
z<-matrix(phat,nrow=10,ncol=10,byrow=TRUE)
z[1,3]
z[1:5,]
z > 0.2
z[z > 0.2]
z * 50
x<-1:5
z * x
z * z

## For loops
x<-vector(mode="numeric",length=100)
x<-numeric(100)
N<-100
for (i in 1:N){
	p<-rbeta(n=1,shape1=5,shape2=10)
	x[i]<-rbinom(n=1,size=20,prob=p)
}
x<-rbinom(n=100,size=20,prob=rbeta(n=1,shape1=5,shape2=10))

x[1]<-0
for(i in 2:N){
	x[i]<-x[i-1] + rpois(n=1,lambda=3)
}

z<-matrix(NA,nrow=10,ncol=10)
for(i in 1:10){
	for(j in 1:10){
		z[i,j]<-rbeta(n=1,shape1=i,shape2=j)
	}
}
