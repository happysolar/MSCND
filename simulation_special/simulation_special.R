source("functions.r")
source("generate_p.R")

### Exp 3, Sparse Case
n=10^6
beta=0.7
r=0.25
sigma=seq(0.2,2,by=0.2)
A=sqrt(2*r*log(n))
epsilon=n^(-beta)



### Exp 3, Dense Case
n=10^6
beta=0.2
r=0.4
sigma=seq(0.2,2,by=0.2)
A=n^(-r)
epsilon=n^(-beta)

error=array(NA,c(3,length(sigma)))
for(i in 1:length(sigma)) {
  si <- sigma[i]
  pval0 <- p.generation(0, n, epsilon, A, si, 100)
  pval1 <- p.generation(1, n, epsilon, A, si, 100)
  HC0=apply(pval0,2,higher.criticism)
  HC1=apply(pval1,2,higher.criticism)
  error[1,i]=min.error.2(rep(c(0,1),each=100),c(HC0,HC1),higher=TRUE)
  Max0=apply(pval0,2,Max)
  Max1=apply(pval1,2,Max)
  error[2,i]=min.error.2(rep(c(0,1),each=100),c(Max0,Max1),higher=FALSE)
  Ar0=apply(pval0,2,a.rOP.2)
  Ar1=apply(pval1,2,a.rOP.2)
  error[3,i]=min.error.2(rep(c(0,1),each=100),c(Ar0,Ar1),higher=TRUE)
  cat(i,'\n')
}
  
plot(error[1,]~sigma,col='red',type='l')
points(error[2,]~sigma,col='green',type='l')
points(error[3,]~sigma,col='blue',type='l')

### Exp 4, Sparse Case
n=10^6
beta=seq(0.55,1,by=0.05)
r=0.25
sigma=sqrt(0.5)
A=sqrt(2*r*log(n))
epsilon=n^(-beta)

### Exp 4, Dense Case
n=10^6
beta=seq(0.05,0.5,by=0.05)
r=0.3
sigma=1
A=n^(-r)
epsilon=n^(-beta)

for(i in 1:length(epsilon)) {
  ep <- epsilon[i]
  pval0 <- p.generation(0, n, epsilon, A, si, 100)
  pval1 <- p.generation(1, n, epsilon, A, si, 100)
  HC0=apply(pval0,2,higher.criticism)
  HC1=apply(pval1,2,higher.criticism)
  Max0=apply(pval0,2,Max)
  Max1=apply(pval1,2,Max)
  Ar0=apply(pval0,2,a.rOP.2)
  Ar1=apply(pval1,2,a.rOP.2)
}
