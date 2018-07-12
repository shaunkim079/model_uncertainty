bc_transform<-function(q,lambda){
  if(lambda[1]!=0){
    return(((((q)+lambda[2])^lambda[1])-1)/lambda[1])
  } else {
    return(log((q)+lambda[2]))
  }
}

bc_fitted_lm<-function(level,beta){
  beta[1]+beta[2]*(level)
}

sd_zero_mean<-function(x) sqrt(mean(x^2))

# gauge_rating_data_file<-"C:/Users/kim079/Documents/PhD/bates/from_justin/Bates gauging details.csv"
dat<-read.csv(gauge_rating_data_file,as.is=T)

# # change cumecs to m3/10min
dat$FLOW<-dat$FLOW*60*10

# change cumecs to m3/day
# dat$FLOW<-dat$FLOW*86400

layout(1:2)
plot(x=dat$M_GH,y=dat$FLOW)

# plot(x=log(dat$M_GH),y=log(dat$FLOW))
plot(x=dat$M_GH,y=dat$FLOW,log="xy")

lm<-lm(log(dat$FLOW)~log(dat$M_GH))

plot(x=log(dat$M_GH),y=log(dat$FLOW))
mod<-lm$coefficients[2]*log(dat$M_GH)+lm$coefficients[1]
lines(x=log(dat$M_GH),y=mod)


plot(x=dat$M_GH,y=dat$FLOW)
xx<-seq(9,11,length.out = 1000)
mod<-exp(lm$coefficients[1])*xx^lm$coefficients[2]
lines(x=xx,y=mod)

# box-cox fit
library(geoR)
bc2<-boxcoxfit(dat$FLOW,dat$M_GH,lambda2 = T)
bc2$lambda
cat(bc2$lambda,"\n")

# log_flow_with_hetero_bctrans<-(((log(flow_with_hetero)+bc2$lambda[2])^bc2$lambda[1])-1)/bc2$lambda[1]
flow_with_hetero_bctrans<-bc_transform(dat$FLOW,bc2$lambda)
layout(1:2)
plot(y=flow_with_hetero_bctrans,x=dat$M_GH)

# pred_flow<-bc2$beta.normal[1]+bc2$beta.normal[2]*log(all_level)
pred_flow<-bc_fitted_lm(dat$M_GH,bc2$beta.normal)
lines(y=range(pred_flow),range(dat$M_GH),col=2)
# bc2$optim.results
# stop()

# log_all_level<-log(all_level)
# lm<-lm(log_flow_with_hetero_bctrans ~ log_all_level)
# pred_flow<-predict(lm)
# lines(y=pred_flow,log_all_level,col=2)

residual_error<-pred_flow-flow_with_hetero_bctrans

hist(residual_error)

population_flow_sd<-sqrt(bc2$sigmasq.normal)

# check
sd_zero_mean(residual_error)


plot_for_paper<-function(){
  jpeg("C:/Users/kim079/Documents/PhD/state_uncertainty_paper1/streamflow_uncertainty.jpg",width=20000,height=10000,pointsize=20,res=1000)
  layout(t(1:2))
  par(mar=c(5.1, 5.5, 4.1, 2.1))
  plot(y=flow_with_hetero_bctrans,x=dat$M_GH,main="(a)",font.main=1,
       ylab="Box-Cox transformed flow",xlab="Stage (m)",cex.lab=1.5,cex.axis=1.5,cex.main=1.5) #expression("Box-Cox transformed flow"~~bgroup("(",m^{3}/s,")"))
  lines(y=range(pred_flow),range(dat$M_GH),col=2)
  legend("topleft",pch=c(1,NA),lty=c(NA,1),legend=c("Observed","Modelled"),col=c(1,2))
  
  par(mar=c(5.1, 5.1, 4.1, 2.1))
  # residual_error_cumecs<-pred_flow-flow_with_hetero_bctrans
  hist(residual_error,main="(b)",font.main=1,
       xlab="Residual",cex.lab=1.5,cex.axis=1.5,cex.main=1.5) # breaks=seq(min(residual_error),max(residual_error),length.out=9),mgp=c(3.5, 1, 0)
  abline(h=0)
  # hist(residual_error)
  # plot(density(residual_error))
  dev.off()
  
  
}












