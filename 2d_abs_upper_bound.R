# This script calculate the upper bound for the selection of the 
# absolute maximum in 2D and generate figure 3.

setEPS(paper = "special", horizontal = T)
library(MASS)
library(lattice)
library(ggplot2)

## Global Constants
{
  alpha<-0.05; za<-qnorm(alpha/2)
}

f<-function(x,mu,sd=1){
  return(dnorm(x,mean=mu[1],sd=sd)*(pnorm(x,mu[2],sd)-pnorm(-x,mu[2],sd))*((-1)^(as.numeric(x<0))))
}
g<-function(x,mu,sd=1){
  return(dnorm(x,mean=mu[2],sd=sd)*(pnorm(x,mu[1],sd)-pnorm(-x,mu[1],sd))*((-1)^(as.numeric(x<0))))
}
calc_area<-function(mu,r=1){
  return(c(integrate(f,mu[1]+(r*za),mu[1]-(r*za),mu)$value,integrate(g,mu[2]+(r*za),mu[2]-(r*za),mu)$value))
}

optimal_r<-function(a,b,leng=300){
  r_tmp<-seq(1,1.15,length=301)  
  tmp<-sapply(X=r_tmp,FUN=calc_area,mu=c(a,b))
  tmp_sum<-apply(tmp,2,sum)
  ind<-which((tmp_sum-(1-alpha)>0)==T)[1]
  return(c(r_tmp[ind],tmp_sum[ind]))
}

mu<-c(2,1)
alpha<-0.05
za<-qnorm(alpha/2)
sum(calc_area(c(50,50)))
sum(calc_area(c(1,0.1)))

t<-seq(0,5,length=150)
tt<-sapply(t,optimal_r,0)
r<-tt[1,]
z_new<-(-r)*za

postscript("figures/upper_bound.pdf",width = 6, height = 6)
plot(z_new~t,type="l",xlab=expression(paste(mu[1])),ylab=expression(paste(c,"(",mu[1],",0)")),main=expression(paste(c,"(",mu[1],",0) as a function of ",mu[1])),lwd=2)
dev.off()


