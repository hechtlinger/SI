# This script generates a visualization of the acceptence regions for the 
# selection of the largest in absolute value (Figure 2). 

setEPS(paper = "special", horizontal = T)
trellis.device("postscript",color=TRUE)

mu=(c(0,0))
r=1
alpha<-0.05; za<-qnorm(alpha/2)

x<-seq(-8,8,length=201)
y<-seq(-8,8,length=201)

plot_2d<-function(mu,r=1){
  x<-seq(-8,8,length=201)
  y<-seq(-8,8,length=201)
  xy<-jitter(as.matrix(expand.grid(x,y)))
  
  max_ind<-max.col(abs(xy))
  max_value<-xy[cbind(1:length(max_ind),max_ind)]
  to_plot<-which((max_value>(mu[max_ind]+r*za))&(max_value<(mu[max_ind]-r*za)))
  tmp_col<-max_ind;tmp_col[tmp_col==1]<-1;tmp_col[tmp_col==2]<-2
  
  plot(y~x,type="l",lty=3,asp=1,main=bquote(mu=="("~.(mu[1])~","~.(mu[2])~")"),xlab=expression(x[1]),ylab=expression(x[2]))
  
  lines(x=(-x),y=y,lty=3,type="l")
  points(xy[to_plot,][tmp_col[to_plot]==1,], col=muted[1],cex=0.1)
  points(xy[to_plot,][tmp_col[to_plot]==2,], col=muted[3],cex=0.1)
  points(mu[1],mu[2],cex=1,pch=16)
}

muted=c("#4878CF", "#6ACC65", "#D65F5F","#B47CC7", "#C4AD66", "#77BEDB")

pdf("figures/AR_visualization.pdf",width = 9, height = 6)
par(mfrow=c(2,3),cex=0.6)
plot_2d(c(0,0))
plot_2d(c(1,0))
plot_2d(c(1.5,1))
plot_2d(c(3,1.25))
plot_2d(c(3.5,2.5))
plot_2d(c(5,5))
dev.off()