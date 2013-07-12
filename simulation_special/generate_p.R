
#p.generation<-function(hypothesis,n,epsilon,A,sigma,rep){
#  obs=array(rnorm(n*rep),dim=c(n,rep))
#  if (hypothesis==1){
#    n.signal=floor(n*epsilon)
#    obs[1:n.signal,]=obs[1:n.signal,]*sigma+A
#  }
#  p=1-pnorm(obs)
#  return(p)
#}

p.generation<-function(hypothesis,n,epsilon,A,sigma,rep){
  p=array(runif(n*rep),dim=c(n,rep))
  if (hypothesis==1){
    n.signal=floor(n*epsilon)
    obs=-abs(array(rnorm(n.signal*rep)*sigma+A,dim=c(n.signal,rep)))
    p[1:n.signal,]=2*pnorm(obs)
  }
  return(p)
}

