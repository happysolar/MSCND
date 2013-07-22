### Obtain the cutoff for U.WSum
U.WSum.cutoff <- function(pval, T, T0=1, T1, N, p0){
  ALPHA <- (1 - p0)/p0
  g <- function(u){
    u^2 * exp(u^2/2)/(ALPHA + exp(u^2/2))/2
  }
  g.moments <- computeMoments(g, 0)
  cutoff.scaled <- getCutoffMultisampleWeightedChisq(pval, T, T0/T, T1, N, ALPHA)
  cutoff.original <- cutoff.scaled * sqrt(N * g.moments$psidotdot) + g.moments$psidot * N
  return(cutoff.original)
}

getCutoffMultisampleWeightedChisq <- function(pval,m,delta,win,N,alpha){
    
    cat("Computing threshold for weighted chi-square...\n")
    THRES = 0.1*pval
    currb = 1
    prevsmallerb = currb
    prevlargerb = 200
    currpval=1
    
    while( abs(currpval-pval)>THRES ){
#        cat("pval =", pval, ", currpval = ",currpval,", THRES=", THRES,".\n",sep="")
        
        if( currpval>pval){
            # need to increase b.
            prevsmallerb = currb
            currb = currb+(prevlargerb-currb)/2
        } else {
            # need to decrease b.
            prevlargerb = currb
            currb = currb - (currb-prevsmallerb)/2
        }
    
#        cat("currb = ",currb,"\n")
        currpval = pvalueMultisampleWeightedChisq(currb,m,delta,win,alpha,N);
    }    
    currb
}

pvalueMultisampleWeightedChisq<-function(b,m,delta,win,ALPHA,N){
    delta1 = min(1, win/m);
#    if(msscan.debug.trace){
#        print(paste("m = ", m, "; b = ", b, "; delta = ", delta, "; delta1 = ", delta1))
#    }

    beta=computeBeta(ALPHA)

    integrand<-function(u){
        vu(sqrt(2)*b*sqrt(beta)/(sqrt(m)*sqrt(u*(1-u))))^2/(u^2*(1-u))
    }
    
    integral =  integrate(integrand,lower=delta,upper=delta1)
    pmarg = pmarg.sumweightedchisq(b,ALPHA,N)

    pval=b^3*beta^2*pmarg*integral$value
    pval
}


pmarg.sumweightedchisq<-function(b,ALPHA,N){
    g<-function(u){
         u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    g.moments=computeMoments(g,0);
    gnormed<-function(u){(g(u)-g.moments$psidot)/sqrt(g.moments$psidotdot)}
    theta0=sqrt(g.moments$psidotdot)/2   # need theta/sqrt(psidotdot) < 1/2 for psi(theta)<infty.
    b0=b/sqrt(N)
    THRESH=0.0001
    marg = computeTiltDirect(b0,gnormed,THRESH,theta0)
    thetabN = marg$theta*sqrt(N)
    pmarg = exp(-thetabN*b + N*marg$psi)/sqrt(2*pi*marg$psidotdot)
    pmarg
}

computeMoments<-function(g,theta){


    INTLIM.THRESH=0.001
    psidottop.int <-function(u){ g(u)*exp(theta*g(u)-(u)^2/2) }    
    for( INTLIM in 10:50 ){
        temp=psidottop.int(INTLIM)
        if (!is.finite(temp)){
            INTLIM=INTLIM-1
            break
        }
        if (temp<INTLIM.THRESH) break
    }
    psidottop =  integrate(psidottop.int,lower=-INTLIM,upper=INTLIM)
    psidottop = psidottop$value/sqrt(2*pi)

    psidotbot.int =  function(u){ exp(theta*g(u)-(u)^2/2)}  
    for( INTLIM in 10:50 ){
        temp=psidotbot.int(INTLIM)
        if( !is.finite(temp)){
            INTLIM=INTLIM-1
            break
        }
        if( temp<INTLIM.THRESH ) break
    }
    psidotbot = integrate(psidotbot.int, lower=-INTLIM, upper=INTLIM)
    psidotbot = psidotbot$value/sqrt(2*pi)
    psidot = psidottop/psidotbot

    EgUsq.int<-function(u){ g(u)^2*exp(theta*g(u))*exp(-(u)^2/2)/psidotbot }    
    
    # us = seq(-INTLIM,INTLIM,.1)
    # plot(us,EgUsq.int(us))
    
    for( INTLIM in 10:50 ){
        temp=EgUsq.int(INTLIM)
        if( !is.finite(temp)){
            INTLIM=INTLIM-1
            break
        }
        if( temp<INTLIM.THRESH) break
    }
    
    EgUsq =  integrate(EgUsq.int,lower=-INTLIM,upper=INTLIM)
    EgUsq = EgUsq$value/sqrt(2*pi)
    
    if(EgUsq<0) stop("Error, overshot theta: E[g^2 e^(th*g-psi(th))]-E[g e^(th*g)]^2 < 0")

    psidotdot = EgUsq - (psidottop/psidotbot)^2;
    psi = log(psidotbot)

    list(psi=psi,psidot=psidot,psidotdot=psidotdot)
}

computeBeta <- function(ALPHA){
    w<-function(u){
         exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    integrand<-function(u){
        u^2*w(u)^2*(2*u^4*w(u)*(1-w(u))+5*u^2*w(u)-3*u^2-2)*dchi(u,1)
    }
    
#    us=seq(0,10,0.1)
#    ys=integrand(us)
#    plot(us,ys)
    
    numerator=integrate(integrand,lower=0,upper=10)
    
    g<-function(u){
         u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    g.moments=computeMoments(g,0)
    numerator$value/(2*g.moments$psidotdot)
}


computeTiltDirect <-function(b,gr,THRESH,theta0, use.binomial=FALSE,binomial.winsize=NA) {
#    Computes, via Newton Raphson, the value of theta such that
#    E_{theta}[g(U)]= b, for U~N(0,1),
#    and given function handle g.
#    E_{theta} is defined as the expectation under the directly
#    tilted measure for g(U), and not for U.
#    
#    Used in importance sampling of:
#    P(max_t sum_1^N g(U_ti) / sqrt(N) > b'), where b'/sqrt(N) = b.

#    cat("Trying Newton Raphson with theta0=",theta0,"\n")
    
    theta = theta0 # if start theta at 0, then dh(theta)=0 at the first step for g(u)=u^2.
    prevtheta=Inf
    prevprevtheta=Inf
    thetarec = theta
    
    while( abs(theta-prevtheta)>THRESH && abs(theta-prevprevtheta)>THRESH ){
        
        if(use.binomial && binomial.winsize>0) {
            g.moments<-computeMomentsBinomial(gr,theta,binomial.winsize)
        } else { 
            g.moments<-computeMoments(gr,theta)
        }
        htheta=g.moments$psidot-b
        dhtheta = g.moments$psidotdot

        # cat("theta=", theta," htheta=",htheta,".\n",sep="")
        thetarec= c(thetarec, theta)
        prevprevtheta=prevtheta
        prevtheta=theta
        theta = prevtheta - htheta/dhtheta
    }

    theta=prevtheta
    psi = g.moments$psi

    list(theta=theta,psi=psi,psidot=g.moments$psidot,psidotdot=g.moments$psidotdot)
}
computeMomentsBinomial<-function(g,theta,m){
# m is the window size, assumes that u~binomial(m,1/2).

    u = (seq(0,m,1)-m/2)/sqrt(m*.5^2)
    p = (0.5^m)*choose(m,seq(0,m,1))
  
    psidottop =  sum( g(u)*exp(theta*g(u))*p)
    psidotbot = sum(exp(theta*g(u))*p)
    psidot = psidottop/psidotbot
    EgUsq = sum(g(u)^2*exp(theta*g(u))*p)
    psidotdot = (EgUsq*psidotbot - psidottop^2)/(psidotbot^2)
    psi = log(psidotbot)

    list(psi=psi,psidot=psidot,psidotdot=psidotdot)
}




getCutoffMultisample <- function(pval,m,m0,m1,N,g,gdot,p0){
    cat("\nComputing p-value cut off for m=",m,", m0=",m0,", m1=",m1,", N=",N,", p0=",p0,", pvalue=",pval,".\n",sep="")

    func<-function(z){ g(z,p0)*exp(-z^2/2)}
    Eg = integrate(func,lower=-10, upper=10)
    Eg = Eg$value/sqrt(2*pi)
    func<-function(z){ (g(z,p0)-Eg)^2*exp(-z^2/2)}
    Vg = integrate(func, lower=-10, upper=10)
    Vg = Vg$value/sqrt(2*pi)

    THRES = 0.1*pval
    currb = N*Eg+1*sqrt(Vg*N)
    prevsmallerb = currb
    prevlargerb = N*Eg+100*sqrt(Vg*N)
    currpval=1
    
    while( abs(currpval-pval)>THRES ){
        cat("     getCutoffMultisample: pval =", pval, ", currpval = ",currpval,", THRES=", THRES,".\n",sep="")
        
        if( currpval>pval){
            # need to increase b.
            prevsmallerb = currb
            currb = currb+(prevlargerb-currb)/2
        } else {
            # need to decrease b.
            prevlargerb = currb
            currb = currb - (currb-prevsmallerb)/2
        }
    
        cat("     currb = ",currb,"\n")
        currpval = pvalueMultisample(currb,m,m0,m1,N,g,gdot,p0)
    }    
    currb
}


pvalueMultisample<-function(b.3,m,m0,m1,N,g,gdot,p0){
# b: threshold
# m: total length
# m0: smallest window size
# m1: largest window size
# N: number of sequences
# g: the transformation function
# p0: weight factor.
# b/N must be larger than psidot(theta=0,b=0), which is the mean under the null.
# 
# Note: pvalueMultisample(b.3,..mllr,mllrdot,p0=1) should be equal to pvalueSumChisqNew(2*b.3,...)
# 
# Added 8/15: pvalueMultisample(b.3,...,g,gdot,p0) is the probability P(max_{s,t} \sum_{i=1}^N g(U_{ist}^2,p0) > b.3)
# 
    
    theta0=0.1
    gr <-function(u){ g(u,p0)}
    numtries=0
    MAX.TRIES=100
    while(numtries<MAX.TRIES){
        theta0 = sqrt(theta0)
        tilt=try(computeTiltDirect(b=b.3/N,g=gr,THRESH=0.00001,theta0=theta0, use.binomial=FALSE), silent=TRUE)
        if(!(class(tilt)=="try-error")) break
        numtries = numtries+1
    }    
    if(numtries==MAX.TRIES) stop("Newton Raphson could not determine correct tilt.")
    th = tilt$theta
    psi = tilt$psi
    psidot=tilt$psidot
    psidotdot = tilt$psidotdot
    
    mu = computeLocalIncrementMean(g=g,gdot=gdot,theta=tilt$theta,p0=p0)
    
    integrand<-function(t){
#        nuFunction(sqrt(2*N*mu)/sqrt(m*t*(1-t)))^2/(t^2*(1-t))
        (1-t)*nuFunction(sqrt(2*N*mu)/sqrt(m*t))^2/(t^2)
    }
    
    integral=integrate(integrand,lower=m0/m,upper=m1/m)
    
    cat("th=",th,", mu=",mu,", integral=", integral$value,"\n")
    
    # break the full expression below into two parts, to compare with pvalueSumChisqNew in the case p0=1:
    # exp(-N*(th*psidot-psi))*(2*pi*N*psidotdot)^(-1/2)*(th^(-1))*N^2*mu^2*integral$value
    part1=exp(-N*(th*psidot-psi))*(2*pi*N*psidotdot)^(-1/2)
    part2=(th^(-1))*N^2*mu^2*integral$value
    
    pval = part1*part2
    pval
}

nuFunction<-function(x){
    ((2/x)*(pnorm(x/2)-.5))/((x/2)*pnorm(x/2)+dnorm(x/2))
}
psiNormal<-function(theta,g,p0){

    INTLIM.THRESH=0.001
    psidotbot.int =  function(u){ exp(theta*g(u,p0)-(u)^2/2)}  
    for( INTLIM in 10:50 ){
        temp=psidotbot.int(INTLIM)
        if( !is.finite(temp)){
            INTLIM=INTLIM-1
            break
        }
        if( temp<INTLIM.THRESH ) break
    }
    psidotbot = integrate(psidotbot.int, lower=-INTLIM, upper=INTLIM)
    psidotbot = psidotbot$value/sqrt(2*pi)
   
    log(psidotbot)
}

computeLocalIncrementMean<-function(g,gdot,theta,p0, LOWER=-50, UPPER=50){
    
    psi = psiNormal(theta,g,p0)
    
    integrand<-function(z){
        gdot(z,p0)^2*exp(theta*g(z,p0)-psi-.5*z^2)
    }
    
#    LOWER=-50; UPPER=50
#    temp=seq(LOWER,UPPER,.1)
#    plot(temp, integrand(temp))
    
    
    integral = integrate(integrand, lower=LOWER,upper=UPPER)
    
    .5*theta^2*integral$value/sqrt(2*pi)
    
}

dchi<-function(y,N){
    fy <-(1-N/2)*log(2) + (N-1)*log(y) - y^2/2 - lgamma(N/2)
    exp(fy)
}


wchi<-function(u,p0){
    alpha=(1-p0)/p0
    eterm = exp(-u^2/2)*alpha
    u^2/(eterm+1)
}

wchidot<-function(u,p0){
    alpha=(1-p0)/p0
    eterm = exp(-u^2/2)*alpha
    (2*u*(eterm+1)+eterm*u^3/2)/(1+eterm)^2
}
vu<-function(x,maxn=1000,do.approx=(abs(x)<0.2)){
#    Evaluates the function vu(x) defined in (4.35) and (4.37) in 
#    Siegmund (1985) Sequential Analysis. 
#    
#    If x<0 then vu(x)= vu(-x).  x must be less than 0.
#    If approx = 1, uses the approximation:
#    v(x) = exp(-Px), where P=0.583.
#    maxn is the upper cap in evaluating the sum in the exponential.

    x=as.matrix(x)
    vux = matrix(0,nrow=nrow(x),ncol=ncol(x))

#    if(sum(is.na(do.approx)) > 0 && msscan.debug.trace){
#        #check for errors here
#        print(do.approx)
#        scan()
#        print(x)
#        scan()
#    }


    if(is.logical(do.approx)){
        if(sum(do.approx) > 0)   vux[do.approx] = exp(-x[do.approx]*0.583);
    }else if(length(do.approx) > 0){
        vux[do.approx] = exp(-x[do.approx]*0.583);
    }
  

    if (sum(do.approx)<length(x)){
        notdo.approx.ix = which(!do.approx)
        n=matrix(c(1:maxn),nrow=1,ncol=maxn)
        summands = pnorm(-0.5*matrix(x[notdo.approx.ix],nrow=length(notdo.approx.ix),ncol=1)%*%sqrt(n))/matrix(rep(n,length(notdo.approx.ix)),ncol=length(n), byrow=TRUE)
        expterm = -2*apply(summands,1,"sum");
        vux[notdo.approx.ix] = (2/x[notdo.approx.ix]^2)*exp(expterm);
    }

    vux
}


num.over <- function(x, threshold){
  sum(x > threshold, na.rm = TRUE)
}
