dat<-read.csv("C:/Users/kim079/Documents/PhD/bates/request4_rainfall/107187/RainfallForSite.csv",as.is=T,skip=1)
coords<-read.csv("C:/Users/kim079/Documents/PhD/bates/request4_rainfall/107187/DR107187 Site Selections For the Request.csv",as.is=T)

image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}


head(dat)
sites<-unique(dat$Default.Site.Ref)

all_dates<-as.Date(unique(dat$Collect.Date),format="%d/%m/%Y")
range_all_dates<-range(all_dates)
all_dates<-seq(range_all_dates[1],range_all_dates[2],by=1)
all_sites_data<-NULL

# remove a bad one
sites<-sites[-which(sites==509400)]

for(i in 1:length(sites)){
  indices<-which(dat$Default.Site.Ref==sites[i])
  cur_site_dates<-as.Date(dat$Collect.Date[indices],format="%d/%m/%Y")
  cur_data<-rep(NA,length(all_dates))
  qual<-rep(1000,length(all_dates))
  matches<-match(cur_site_dates,all_dates)
  cur_data[matches]<-dat$Reading.Value[indices]
  qual[matches]<-dat$Quality[indices]
  qual_is_bad<-which(qual>1)
  cur_data[qual_is_bad]<-NA
  
  all_sites_data<-cbind(all_sites_data,cur_data)
}

all_sites_data<-cbind(as.character(all_dates),all_sites_data)

colnames(all_sites_data)<-c("date",sites)

all_sites_data<-all_sites_data[which(all_dates>="1988-01-01"),]

matplot(all_sites_data[,-1],type="l")
matplot(all_sites_data[,2:3],type="l")

# all(!is.na(all_sites_data[1,-1]))
row_all_good<-function(x){
  all(!is.na(x))
}
# row_all_good<-function(x){
#   ww<-which(!is.na(x))
#   return(length(ww)>=10)
# }

ok_rows<-apply(all_sites_data[,-1],1,row_all_good)
any(ok_rows)
length(which(ok_rows))
all_sites_data[ok_rows,1]
which(ok_rows)
which(ok_rows)[length(which(ok_rows))]-which(ok_rows)[1]
head(all_sites_data[ok_rows,])

ok_rows_indices<-which(ok_rows)
all_sites_data_good<-all_sites_data[ok_rows_indices[1]:ok_rows_indices[length(ok_rows_indices)],]

matplot(all_sites_data_good[,-1],type="l")

# which(is.na(all_sites_data_good[,13]))

ok_dates<-seq(as.Date(all_sites_data_good[1,1]),as.Date(all_sites_data_good[nrow(all_sites_data_good),1]),by=1)

length(ok_dates)
nrow(all_sites_data_good)

# View(all_sites_data_good)

# find the longest continuous period
counter<-0
period_indices<-c()
all_period_indices<-list()
all_counts<-c()
for(i in 1:nrow(all_sites_data_good)){
  ok<-all(!is.na(all_sites_data_good[i,]))
  if(ok){
    counter<-counter+1
    period_indices<-c(period_indices,i)
  } else {
    if(counter>0) all_counts<-c(all_counts,counter)
    counter<-0
    all_period_indices[[length(all_period_indices)+1]]<-period_indices
    period_indices<-c()
  }
}

longest_period_indices<-all_period_indices[which(all_counts==max(all_counts))][[1]]
all_sites_data_good_longest<-all_sites_data_good[longest_period_indices,]
matplot(all_sites_data_good_longest[,-1],type="l")


# calculate the variogram
library(geoR)

data_variog_precip<-apply(all_sites_data_good_longest[,-1],2,as.numeric)
indices<-match(colnames(data_variog_precip),coords$DEFAULT_SITE_REF)
locations_easting_northing<-cbind(coords$EASTING[indices],coords$NORTHING[indices])

layout(matrix(1:10,ncol=2))
for(i in 1:1){ #1:nrow(data_variog_precip)
  par(mfrow=c(2,2))
  v2<-variog(coords=locations_easting_northing,data=data_variog_precip[i,])
  plot(v2,type = "b")
  v2<-variog(coords=locations_easting_northing,data=data_variog_precip[i,],op="cloud")
  plot(v2, main="variogram cloud")
  v2<-variog(coords=locations_easting_northing,data=data_variog_precip[i,], bin.cloud=TRUE)
  plot(v2, bin.cloud=TRUE, main="clouds for binned variogram") 
  # v2<-variog(coords=locations_easting_northing,data=data_variog_precip[i,], op="sm", band=0.2)
  # plot(v2, main="smoothed variogram") 
  cat(data_variog_precip[i,],"\n")
  # browser()
}


#
# computing variograms:
#
# binned variogram
vario.b <- variog(s100, max.dist=1)
# variogram cloud
vario.c <- variog(s100, max.dist=1, op="cloud")
#binned variogram and stores the cloud
vario.bc <- variog(s100, max.dist=1, bin.cloud=TRUE)
# smoothed variogram
vario.s <- variog(s100, max.dist=1, op="sm", band=0.2)
#
#
# plotting the variograms:
par(mfrow=c(2,2))
plot(vario.b, main="binned variogram") 
plot(vario.c, main="variogram cloud")
plot(vario.bc, bin.cloud=TRUE, main="clouds for binned variogram")  
plot(vario.s, main="smoothed variogram") 

# computing a directional variogram
vario.0 <- variog(s100, max.dist=1, dir=0, tol=pi/8)
plot(vario.b, type="l", lty=2)
lines(vario.0)
legend("topleft", legend=c("omnidirectional", expression(0 * degree)), lty=c(2,1))




#############################################################
## ##
## Conditional simulation ## 
## ##
#############################################################
library(RandomFields)
# first we simulate some random values at a
# 100 random locations:
n <- 100
x <- runif(n=n, min=-1, max=1)
y <- runif(n=n, min=-1, max=1)
# model <- RMgauss(1)
model <- RMexp()
# plot(model)
# model<-RPtbm(RMexp())
RFoptions(seed=NA)
data <- RFsimulate(model = model, x=x, y=y, grid=FALSE,n=1)
# plot(data)
# x11(); plot(data)

# let simulate a field conditional on the above data
x.seq.cond <- y.seq.cond <- seq(-1.5, 1.5, length=n)
# model <- RMgauss()

cond <- RFsimulate(model, x=x.seq.cond, y=y.seq.cond, data=data)
x11(); plot(cond, data,zlim=range(cond$variable1))





# Examples
RFoptions(seed=0) ## *ANY* simulation will have the random seed 0; set
##                   RFoptions(seed=NA) to make them all random again


RFoptions(modus_operandi="sloppy")

#########################################################
## simulate some data first                            ## 
points <- 100

x <- runif(points, 0, 3)
y <- runif(points, 0, 3) ## random points in square [0, 3]^2
# x<-coords$EASTING[indices]
# y<-coords$NORTHING[indices]
model <- RMgencauchy(alpha=1, beta=2)
d <- RFsimulate(model, x=x, y=y, grid=FALSE, n=1) #1000
estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=2) +
  RMtrend(mean=NA)
eg.rff<-RFfit(estmodel, data=d)

# checking that the fitting looks ok
x11(); plot(eg.rff)

eg.variog<-variog(coords=cbind(x,y),data=d$variable1)
x11(); plot(eg.variog)
eg.variog<-variog(coords=cbind(x,y),data=d$variable1,op="cloud")
x11(); plot(eg.variog)
eg.variog<-variog(coords=cbind(x,y),data=d$variable1, bin.cloud=TRUE)
x11(); plot(eg.variog, bin.cloud=TRUE)

coords_RF<-cbind(coords$EASTING[indices],coords$NORTHING[indices])
colnames(coords_RF)<-c("coords.x1","coords.x2")
coords_with_id<-cbind(coords$DEFAULT_SITE_REF[indices],coords$EASTING[indices],coords$NORTHING[indices])

rotate <- function(x) t(apply(x, 2, rev))
data_variog_precip_rot<-rotate(rotate(rotate(data_variog_precip)))
data_variog_precip_rot<-data_variog_precip_rot[nrow(data_variog_precip_rot):1,]
# head(data_variog_precip_rot)
# data_variog_precip_rot[1,]
# data_variog_precip_rot[2,]
data_colnames<-paste0("variable1.n",1:ncol(data_variog_precip_rot))
rownames(data_variog_precip_rot)<-1:nrow(data_variog_precip_rot)
colnames(data_variog_precip_rot)<-data_colnames
data_df<-as.data.frame(data_variog_precip_rot)
rf_df<-RFspatialPointsDataFrame(coords=coords_RF,data=data_df,RFparams = list(n = ncol(data_df), vdim = 1))
# rf_df<-RFspatialPointsDataFrame(coords=coords_RF,data=data_df,RFparams = list(n = 1, vdim = ncol(data_df)))
# rf_df<-RFspatialPointsDataFrame(coords=coords_RF,data=data_df[,1],RFparams = list(n = 1, vdim = 1))

# as.data.frame(data_variog_precip_rot[,1:2])$variable1.n1

rf_df

#########################################################
## estimation; 'NA' means: "to be estimated"           ##
# estmodel <- RMexp(var=NA, scale=NA) +
#   RMtrend(mean=NA)
# rff<-RFfit(estmodel, data=rf_df)
# plot(rff)
# x11(); plot(rff)
# 
# estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=2) +
#   RMtrend(mean=NA)
# rff<-RFfit(estmodel, data=rf_df)
# plot(rff)
# x11(); plot(rff)
# 
# estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=NA) + 
#   RMtrend(NA)
# rff<-RFfit(estmodel, data=rf_df) ## just for information
# plot(rff)
# x11(); plot(rff)


estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=NA)
rff<-RFfit(estmodel, data=rf_df) ## just for information
plot(rff)


# loop through a few days to have a look
for(i in 1:1){
  if(length(which(data_df[,i]>0))>1){
    rf_df<-RFspatialPointsDataFrame(coords=coords_RF,data=data_df[,i],RFparams = list(n = 1, vdim = 1))
    estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=NA)
    rff<-RFfit(estmodel, data=rf_df)
    x11(); plot(rff)
    # v2<-variog(coords=locations_easting_northing,data=data_variog_precip[i,],op="cloud")
    v2<-variog(coords=locations_easting_northing,data=data_variog_precip[i,])
    v2<-variog(coords=locations_easting_northing,data=data_variog_precip[i,],uvec=seq(0,7000,length.out=20))
    x11(); plot(v2,type="b")
    
    # browser()
  }
}




rff_att<-attributes(rff)
gencauchy.s<-attributes(rff_att$ml)$param[1,1]
gencauchy.alpha<-attributes(rff_att$ml)$param[1,2]
gencauchy.beta<-attributes(rff_att$ml)$param[1,3]
globalvariance<-attributes(rff_att$ml)$globalvariance

model <- RMgencauchy(var=globalvariance, scale=gencauchy.s, alpha=gencauchy.alpha, beta=gencauchy.beta)
# model <- RMgencauchy(scale=gencauchy.s, alpha=gencauchy.alpha, beta=gencauchy.beta) # should give results same as above
# plot(model)

x.seq.cond <- seq(min(coords_RF[,1]),max(coords_RF[,1]),length.out = 90)
y.seq.cond <- seq(min(coords_RF[,2]),max(coords_RF[,2]),length.out = 90)

# first day
rf_df1<-RFspatialPointsDataFrame(coords=coords_RF,data=data_df[,1],RFparams = list(n = 1, vdim = 1))
cond <- RFsimulate(model, x=x.seq.cond, y=y.seq.cond, data=rf_df1, n=4)
x11(); plot(cond)
x11(); plot(rf_df1)


cond$variable1
length(cond$variable1)
attributes(cond)

RFoptions(seed=NA)
cond <- RFsimulate(model, x=x.seq.cond, y=y.seq.cond, data=rf_df1)
cond2 <- RFsimulate(model, x=x.seq.cond, y=y.seq.cond, data=rf_df1)
diff<-range(cond$variable1-cond2$variable1)

cond <- RFsimulate(model, x=x.seq.cond, y=y.seq.cond, data=rf_df1, n=4)
x11(); plot(cond)


# simulate ensembles for each day and calculate rainfall for just the catchment area
x.seq.cond.all <- seq(min(coords_RF[,1]),max(coords_RF[,1]),by = 400) #seq(min(coords_RF[,1]),max(coords_RF[,1]),length.out = 40)
y.seq.cond.all <- seq(min(coords_RF[,2]),max(coords_RF[,2]),by = 400) #seq(min(coords_RF[,2]),max(coords_RF[,2]),length.out = 40)

coords
library(sp)
library(rgdal)
# testing
xy <- data.frame(ID = 1, X = coords$LONGITUDE[1], Y = coords$LATITUDE[1])
coordinates(xy) <- c("X", "Y")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
rr<-spTransform(xy, CRS("+proj=utm +zone=50J +ellps=WGS84 +south"))
coords$EASTING[1]
coords$NORTHING[1]
spTransform(rr, CRS("+proj=longlat +ellps=WGS84 +south"))


library(sp)
library(rgdal)
xy <- data.frame(ID = 1:2, X = c(116.011,116.038), Y = c(-32.563,-32.583))
# xy <- data.frame(ID = 1:2, X = c(116.010238,116.049437), Y = c(-32.563613,-32.589544))
# xy <- data.frame(ID = 1:2, X = c(116.0273,116.0282), Y = c(-32.5734,-32.5742))
coordinates(xy) <- c("X", "Y")
proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
# proj4string(xy) <- CRS("+proj=longlat")  ## for example


res <- spTransform(xy, CRS("+proj=utm +zone=50J +ellps=WGS84 +south"))
spTransform(res, CRS("+proj=longlat +datum=WGS84"))

x.seq.cond <- seq(coordinates(res)[1,1],coordinates(res)[2,1],by = 100) #seq(407082.55,410788.24,length.out = 90)
y.seq.cond <- seq(coordinates(res)[2,2],coordinates(res)[1,2],by = 100) #seq(6393787.44,6396692.25,length.out = 90)

# x.seq.cond <- seq(407082.55,410788.24,by = 100) #seq(407082.55,410788.24,length.out = 90)
# y.seq.cond <- seq(6393787.44,6396692.25,by = 100) #seq(6393787.44,6396692.25,length.out = 90)

all_y.seq.cond<-sort(rep(y.seq.cond,length(x.seq.cond)))
all_x.seq.cond<-rep(x.seq.cond,length(y.seq.cond))

coords_points<-cbind(all_x.seq.cond,all_y.seq.cond)

# library(sp)
# library(rgdal)

xy <- data.frame(ID = 1:length(all_y.seq.cond), X = all_x.seq.cond, Y = all_y.seq.cond)
coordinates(xy) <- c("X", "Y")
proj4string(xy) <- CRS("+proj=utm +zone=50J +ellps=WGS84 +south")  ## for example

res <- spTransform(xy, CRS("+proj=longlat +datum=WGS84 +south"))
res

write.csv(coordinates(res),"C:/Users/kim079/Documents/PhD/bates/gis/coords_points.csv",quote=F,row.names=F)
# use the above to in arcgis and clip it is subcatchment boundary shapefile

library(foreign)
subcat_coords<-read.dbf("C:/Users/kim079/Documents/PhD/bates/gis/actual_bates_subcat/bates_cat_grid_rain_points_clipped_v4.dbf")
coordinates(subcat_coords) <- c("X", "Y")
proj4string(subcat_coords) <- CRS("+proj=longlat +datum=WGS84")  ## for example
subcat_coords_utm <- spTransform(subcat_coords, CRS("+proj=utm +zone=50J +ellps=WGS84 +south"))

coords_char<-paste(round(all_x.seq.cond,2),round(all_y.seq.cond,2),sep=",")
coords_char_mat<-matrix(coords_char,nrow=length(x.seq.cond),byrow=F)


subcat_coords_utm_char<-paste(round(coordinates(subcat_coords_utm)[,1],2),round(coordinates(subcat_coords_utm)[,2],2),sep=",")

index_subcat<-match(subcat_coords_utm_char,coords_char_mat)
coords_char_mat[index_subcat]
all(coords_char_mat[index_subcat]==subcat_coords_utm_char)

subcat_mask<-matrix(NA,nrow=nrow(coords_char_mat),ncol=ncol(coords_char_mat))
subcat_mask[index_subcat]<-1
x11(); image(subcat_mask,col = (heat.colors(12)))

# transforming the data perform a boxcox transform
# fit box-cox transform
data_df_agg<-unlist(data_df)
hist(data_df_agg)
library(geoR)
# bc2<-boxcoxfit(data_df_agg,data_df_agg,lambda2 = T)
bc2<-boxcoxfit(data_df_agg,lambda2 = T)
tmp<-rnorm(100,100,3)^3
hist(tmp)
bc2<-boxcoxfit(tmp,lambda2 = T)
# tmp_trans<-transform_data(tmp,bc2$lambda)
# hist(tmp_trans)

plot(bc2)
# layout(1:3)
# plot(bc2)
# lines(bc2)
bc2$lambda
cat(bc2$lambda,"\n")
transform_data<-function(q,lambda){
  if(lambda[1]!=0){
    return(((((q)+lambda[2])^lambda[1])-1)/lambda[1])
  } else {
    return(log((q)+lambda[2]))
  }
}
inverse_transform_data<-function(q_trans,lambda){
  if(lambda[1]!=0){
    return(((q_trans*lambda[1]+1)^(1/lambda[1]))-lambda[2])
  } else {
    return(exp(q_trans)-lambda[2])
  }
}
# inverse_transform_data(-4,c(0.3946199,0.0000114))
# inverse_transform_data(-4,c(-0.3946199,0.0000114))



data_bctrans<-transform_data(data_df_agg,bc2$lambda)
hist(data_bctrans)
layout(1:2)

plot(y=data_bctrans,x=data_df_agg)
pred<-rep(bc2$beta.normal,length(data_df_agg))
lines(y=pred,x=data_df_agg,col=2)
resid<-pred-data_bctrans

population_flow_sd<-sqrt(bc2$sigmasq.normal)
sd_zero_mean<-function(x) sqrt(mean(x^2))
population_flow_sd
sd_zero_mean(resid)

hist(data_bctrans)

graphics.off()
# check
# data_invbctrans<-inverse_transform_data(data_bctrans,bc2$lambda)

do_bctrans<-F
test_run<-F
if(test_run){
  nsim<-10
} else {
  nsim<-100
}

test_run2<-F
test_run3<-F
test_run4<-F

if(test_run4) nsim <- 50
paper_plot<-F
paper_plot2<-T
nsim<-6

RFoptions(seed=NA)

set.seed(as.integer(Sys.time()))

all_condvariable1_mat_invtrans_masked<-list()
all_v2<-list()
all_sim_semivar_all<-list()
all_dist_line<-list()
all_ensemble_rain<-NULL
all_simulated_dates<-c()
for(i in 6){ #1:ncol(data_df) #c(2,6,17,19)
  all_simulated_dates<-c(all_simulated_dates,all_sites_data_good_longest[i,1])
  daily_ensemble_rain<-c()
  all_sim_pars<-NULL
  all_opt_pars<-NULL
  all_sim_semivar<-NULL
  for(s in 1:nsim){
    cat("--------------------------------------------------------------- i =",i,"s =",s,"\n")
    if(length(which(data_df[,i]>0))>0){
      cat("data:",data_df[,i],"\n")
      cat("data range:",range(data_df[,i]),"\n")
      if(do_bctrans){
        # transform:
        bc2<-boxcoxfit(data_df[,i],lambda2 = T)
        data_df_trans<-transform_data(data_df[,i],bc2$lambda)
        inverse_transform_data(data_df_trans,bc2$lambda)
        # hist(data_df_trans)
      } else {
        # bc2<-list(lambda=NA)
        # transform_data<-function(x,y) x
        # inverse_transform_data<-function(x,y) x
        # diff_from_mean<-mean(data_df[,i])-data_df[,i]
        # dble_diff_from_mean<-diff_from_mean*2
        # # mean(data_df[,i])-diff_from_mean
        # # data_df_trans<-mean(data_df[,i])-dble_diff_from_mean
        # data_df_trans<- -dble_diff_from_mean
        # # data_df_trans<-data_df[,i]
        # 
        # rf_df<-RFspatialPointsDataFrame(coords=coords_RF,data=data_df[,i],RFparams = list(n = 1, vdim = 1))
        # estmodel <- RMspheric(var=NA, scale=NA)
        # rff_orig<-RFfit(estmodel, data=rf_df)
        # x11(); plot(rff_orig)
        
        bc2<-list(lambda=NA)
        # This now seems to cause an oversestimation of variance
        # try the built in boxcox transform but with limiting bounds - not done but not needed yet
        # transform_data<-function(x,y){
        #   sign<-sample(c(-1,1),length(x),replace = T)
        #   return(sqrt(x)*sign)
        # }
        # inverse_transform_data<-function(x,y) x^2
        transform_data<-function(x,y){
          return(x^(1/3))
          # return(base::log(x+1e-6,base=20))
          # return((x+101)^(1/1001))
          # return((x+1e-6)^(1/3))
          # return(log10(x+1e-6))
          # return(log(x+1e-6))
        }
        inverse_transform_data<-function(x,y){
          #(x^3)
          # browser()
          # if(any(x_inv_trans<0)){
          #   x_inv_trans<-x_inv_trans-min(x_inv_trans)
          # }
          return(x^3)
          # x_inv_trans<-(20^x)+1e-6
          # return(x_inv_trans)
          # return((x^1001)-101)
          # return((10^x)-1e-6)
          # return((x^3)-1e-6)
          # return(exp(x)-1e-6)
        }
        # transform_data<-function(x,y) log(x+1e-10)
        # inverse_transform_data<-function(x,y) exp(x)-1e-10
        data_to_use<-data_df[,i]
        # data_to_use[data_to_use==0]<-abs(rnorm(length(which(data_to_use==0)),0,1e-2))
        
        # Assume that the distribution is truncated at zero so adjust the zeros to extend into the negative domain
        # then get rid of the negatives later when resimulated
        
        # Assume below is truncated normal
        data_df_trans<-transform_data(data_to_use,bc2$lambda)
        
        if(test_run) hist(data_df_trans)
        
        replace_zeros<-function(data_df_trans){
          # get the parameters for the normal distribution
          optim_normal<-function(pars,data){
            # cat(pars,"\n")
            hh<-hist(data,plot=F)
            dens_actual<-hh$density[-1]
            dens_trial<-dnorm(hh$mids[-1],mean=pars[1],sd=pars[2])
            
            # dd<-density(data)
            # dens_actual<-dd$y[dd$x>hh$breaks[2]]
            # dens_trial<-dnorm(dd$x[dd$x>hh$breaks[2]],mean=pars[1],sd=pars[2])
            
            if(is.na(dens_trial[1])){
              ss<-Inf
            } else {
              ss<-sum((dens_actual-dens_trial)^2)
            }
            
            # ss<--sum(dnorm(data[data>0],mean=pars[1],sd=pars[2],log = T))
            # cat(ss,"\n")
            return(ss)
          }
          # library(hydromad)
          optim_normal(pars=c(mean(data_df_trans[data_df_trans>0]),0.01),data=data_df_trans)
          opt_normal<-hydromad::SCEoptim(par=c(mean(data_df_trans[data_df_trans>0]),0.01),FUN=optim_normal,data=data_df_trans,
                                         lower=c(-Inf,1e-30),upper=c(max(data_df_trans),diff(range(data_df_trans))),
                                         control=list(ncomplex=20))
          opt_normal<-hydromad::SCEoptim(par=opt_normal$par,FUN=optim_normal,data=data_df_trans,
                                         lower=c(-Inf,1e-30),upper=c(max(data_df_trans),diff(range(data_df_trans))),
                                         control=list(ncomplex=20))
          mean_normal<-opt_normal$par[1]
          sd_normal<-opt_normal$par[2]
          
          
          # library(DEoptim)
          # opt_normal<-DEoptim::DEoptim(fn=optim_normal,data=data_df_trans,
          #                     lower=c(0,1e-30),upper=c(max(data_df_trans),diff(range(data_df_trans))),
          #                     control=list(trace=F))
          # mean_normal<-opt_normal$optim$bestmem[1]
          # sd_normal<-opt_normal$optim$bestmem[2]
          

          if(test_run){
            tmp<-rnorm(10000,mean=mean_normal,sd=sd_normal)
            # tmp2<-hist(tmp,plot=F)
            x11()
            hist(data_df_trans,freq=F,xlim=c(-1,2),ylim=c(0,7))
            # lines(y=tmp2$density,x=tmp2$mids,col=2)
            # lines(density(tmp),col=2)
            # plot(y=tmp2$density,x=tmp2$mids,col=2,type="l")
            # plot(density(tmp),col=2)
            
            # plot(density(data_df_trans))
            # lines(density(tmp),col=2)
            # 
            # lines(density(data_df_trans[data_df_trans>0]))
            lines(density(tmp),col=2)
            # 
            if(test_run2) browser()
          }

          # sample and replace the zeros
          library(truncnorm)
          replacements<-rtruncnorm(length(which(data_df_trans<=0)),a=-Inf,b=0,mean=mean_normal,sd=sd_normal)
          data_df_trans[data_df_trans<=0]<-replacements
          
          if(test_run){
            x11()
            hist(data_df_trans,freq=F,xlim=c(-1,2),ylim=c(0,7))
            lines(density(tmp),col=2)
            # x11()
            # plot(density(data_df_trans))
            # lines(density(tmp),col=2)
            
            if(test_run2) browser()
          }

          
          return(data_df_trans)
        }
        
        data_df_trans<-replace_zeros(data_df_trans)
        
        # cat(data_df_trans,"\n")
        # hist(replace_zeros(data_df_trans))
        
        # replacements<-rtruncnorm(1,a=-Inf,b=0,mean=mean_data_df_trans,sd=sd_normal)
        # data_df_trans[data_df_trans<=0]<-rep(replacements,length(which(data_df_trans<=0)))
        
        # This is done for resimulation tests
        # TODO: something to think about: the replacement does not carry any information about the variance in terms of distance
        # TODO: Should be repeat a sample? might be ok because any negatives will be set to zero and covariance preserved for those
        # This is ok if we assume that the "real" distribution is the truncated normal
        
        if(test_run) hist(data_df_trans)
        
        
      }
      
      mean_trans_data<-mean(data_df_trans)
      rf_df<-RFspatialPointsDataFrame(coords=coords_RF,data=data_df_trans-mean_trans_data,RFparams = list(n = 1, vdim = 1))
      # rf_df<-RFspatialPointsDataFrame(coords=coords_RF,data=data_df_trans,RFparams = list(n = 1, vdim = 1))
      # estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=NA)
      # rff<-RFfit(estmodel, data=rf_df)
      # rff.tmp<-RFfit(estmodel, data=rf_df,upper=c(100000,20,100))
      
      RFoptions(bins=5)
      estmodel <- RMspheric(var=NA, scale=NA)
      rff<-RFfit(estmodel, data=rf_df)
      # plot(rff)
      
      # # checks
      if(test_run){
        x11(); plot(rff)
      }
      if(test_run2) browser()
      ## empirical variogram plots variog
      v2<-variog(coords=locations_easting_northing,data=data_df_trans-mean_trans_data,op="cloud")
      # v2<-variog(coords=locations_easting_northing,data=data_df_trans,op="cloud")
      if(test_run){
        x11(); plot(v2)
      }
      # v2<-variog(coords=locations_easting_northing,data=data_df_trans-mean_trans_data)
      # Fit without binning - binning is stupid
      v2$u
      v2$v
      
      spherical_semivariogram<-function(dist, range, scale){
        within_range<-dist<=range
        res<-rep(NA,length(dist))
        res[within_range]<-scale*(3/2*dist[within_range]/range - 1/2*(dist[within_range]/range)^3)
        res[!within_range]<-scale
        return(res)
      }
      
      # xx<-seq(0,3,length.out=100)
      # yy<-spherical_semivariogram(xx,1,4)
      # layout(1)
      # plot(x=xx,y=yy,type="l")
      # abline(v=1,col=2)
      
      # x11(); plot(rff,xlim=c(0,20000))
      rff_att<-attributes(rff)
      globalvariance<-attributes(rff_att$ml)$globalvariance
      scale<-attributes(rff_att$ml)$param[1,1]

      # xx<-seq(0,20000,length.out=20000)
      # yy<-spherical_semivariogram(xx,range=scale,scale=globalvariance)
      # layout(1)
      # plot(x=xx,y=yy,type="l")
      # abline(v=scale,col=2)
      
      spherical_semivariogram_optim_func<-function(par,dist,semivar_obs){
        semivar_sim<-spherical_semivariogram(dist, range=par[1], scale=par[2])
        ss<-sum((semivar_sim-semivar_obs)^2)
        return(ss)
      }
      
      spherical_semivariogram_optim_func(c(scale,globalvariance),v2$u,v2$v)
      # opt<-optim(par=c(scale,globalvariance),fn=spherical_semivariogram_optim_func,dist=v2$u,semivar_obs=v2$v)
      opt<-hydromad::SCEoptim(par=c(scale,globalvariance),FUN=spherical_semivariogram_optim_func,dist=v2$u,semivar_obs=v2$v,
                              lower=c(0,0))
      all_opt_pars<-rbind(all_opt_pars,opt$par)
      # opt<-DEoptim::DEoptim(fn=spherical_semivariogram_optim_func,dist=v2$u,semivar_obs=v2$v,
      #                         lower=c(0,0),upper=c(1e30,1e30))
      
      # layout(1)
      
      dist_line<-seq(0,max(v2$u),by=100)
      sim_semivar<-spherical_semivariogram(dist_line, opt$par[1], opt$par[2])
      all_sim_semivar<-rbind(all_sim_semivar,sim_semivar)
      
      if(test_run){
        x11()
        plot(x=v2$u,y=v2$v)
        lines(x=dist_line,y=sim_semivar,col=2)
        browser()

      }
      all_sim_semivar_all[[length(all_sim_semivar_all)+1]]<-sim_semivar
      all_dist_line[[length(all_dist_line)+1]]<-dist_line
      all_v2[[length(all_v2)+1]]<-v2
      # if(paper_plot){
      #   # binning for plotting
      #   x11()
      #   bin_intervals<-seq(0,max(v2$u),length.out=7)
      #   semivariance_bins<-list()
      #   for(bb in 2:length(bin_intervals)){
      #     bin_indices<-which(v2$u<bin_intervals[bb] & v2$u>bin_intervals[bb-1])
      #     v2$v[bin_indices]
      #     v2$u[bin_indices]
      #     # plot(x=v2$u[bin_indices],y=v2$v[bin_indices])
      #     semivariance_bins[[bb-1]]<-v2$v[bin_indices]
      #   }
      #   bin_intervals_mids<-bin_intervals[-1]
      #   bin_intervals_mids<-bin_intervals_mids-((bin_intervals[2]-bin_intervals[1])/2)
      #   semivariance_bins_mean<-sapply(semivariance_bins,mean)
      #   jpeg("C:/Users/kim079/Documents/PhD/state_uncertainty_paper1/variogram_fitting.jpg",width=20000,height=10000,pointsize=20,res=1000)
      #   layout(matrix(1:4,nrow=2,byrow=2))
      #   plot(x=bin_intervals_mids,y=semivariance_bins_mean,xlim=c(0,max(bin_intervals_mids)),ylim=c(0,max(semivariance_bins_mean)),type="b",
      #        ylab="Semivariance",xlab="Distance (m)",lwd=2)
      #   lines(x=dist_line,y=sim_semivar,col=2,lwd=2)
      #   dev.off()
      #   browser()
      # }

      
      if(test_run2) browser()

      # v2<-variog(coords=locations_easting_northing,data=data_df_trans-mean_trans_data,uvec=seq(0,18000,length.out=17))
      # x11(); plot(v2,type="b")
      
      ## empirical variogram plots RFempiricalvariogram
      # ev1 <- RFempiricalvariogram(data=rf_df, bin=seq(0,20000,length.out=10000))
      # ev1 <- RFempiricalvariogram(data=rf_df)
      # x11()
      # plot(ev1)
      
      
      # rff_att<-attributes(rff)
      
      # estmodel_att<-attributes(estmodel)
      
      # scale<-attributes(rff_att$ml)$param[1,1]
      # gencauchy.alpha<-attributes(rff_att$ml)$param[1,2]
      # gencauchy.beta<-attributes(rff_att$ml)$param[1,3]
      # globalvariance<-attributes(rff_att$ml)$globalvariance
      # scale<-attributes(rff_att$ml)$param[1,1]
      
      # model <- RMgencauchy(var=globalvariance, scale=scale, alpha=gencauchy.alpha, beta=gencauchy.beta)
      # model <- RMspheric(var=globalvariance, scale=scale)
      model <- RMspheric(var=opt$par[2], scale=opt$par[1])
      
      # modsim<-RFsimulate(model,x=seq(0,1000000,by=1))
      # plot(modsim)
      # plot(model,xlim=c(0,8000))
      condvariable1_mat_invtrans<-NaN
      counter<-0
      # create a while loop so that only stops when no NaNs are produced? good idea?
      while(any(is.nan(condvariable1_mat_invtrans) | any(condvariable1_mat_invtrans<0))){
        counter<-counter+1
        if(counter>1) {
          cat("seems to be taking too long to get a viable simulation\n")
          break
        }
        # RFoptions(seed=as.integer(Sys.time()))
        # set.seed(as.integer(Sys.time()))
        # cond <- RFsimulate(model, x=x.seq.cond, y=y.seq.cond, data=rf_df, n=1)
        cond <- RFsimulate(model, x=x.seq.cond, y=y.seq.cond, data=rf_df, n=1)
        # cond <- RFsimulate(model, x=x.seq.cond, y=y.seq.cond, data=rf_df, n=100)
        # x11(); plot(cond)
        
        # checks to see if observed gauge value is being preserved
        # cond$variable1.n1
        # condvariable1_mat<-matrix(cond$variable1.n1,nrow=length(y.seq.cond))
        cond$variable1<-cond$variable1+mean_trans_data
        # cond$variable1<-cond$variable1
        
        condvariable1_mat<-matrix(cond$variable1,nrow=length(x.seq.cond))
        condvariable1_mat<-t(apply(t(condvariable1_mat), 2, rev))
        # x11(); image(condvariable1_mat,col = (heat.colors(12)))
        # Remember: The matrices are rotated to the left when actually projected onto an image
        
        # inverse transform
        condvariable1_mat_tmp<-condvariable1_mat
        condvariable1_mat_tmp[condvariable1_mat_tmp<=0]<-0
        
        # condvariable1_mat_invtrans represents the final data (almost) i.e. the mean is added above
        condvariable1_mat_invtrans<-inverse_transform_data(condvariable1_mat_tmp,bc2$lambda)
        # inverse_transform_data(condvariable1_mat[20,20],bc2$lambda)
        # inverse_transform_data(condvariable1_mat,bc2$lambda)[20,20]
      }
      if(test_run){
        x11(); layout(1); plot(cond)
      }
      if(test_run){
        x11(); layout(1); image(condvariable1_mat_invtrans,col = (heat.colors(12)))
      }
      
      # this ensures there aren't any negatives
      index_nans<-which(is.nan(condvariable1_mat_invtrans),arr.ind = T)
      if(length(index_nans)>0) stop("there shouldn't be any NaNs")
      index_small<-which(condvariable1_mat[index_nans] < min(data_df_trans))
      condvariable1_mat_invtrans[index_nans][index_small] <- min(condvariable1_mat_invtrans,na.rm=T)
      index_big<-which(condvariable1_mat[index_nans] > max(data_df_trans))
      condvariable1_mat_invtrans[index_nans][index_big] <- max(condvariable1_mat_invtrans,na.rm=T)
      
      # condvariable1_mat_invtrans[is.nan(condvariable1_mat_invtrans)]<-0
      condvariable1_mat_invtrans[condvariable1_mat_invtrans<0]<-0
      
      if(test_run){
        x11(); image(condvariable1_mat_invtrans,col = (heat.colors(12)))
      }
      cat("sim range small:",range(condvariable1_mat_invtrans),"\n")
      if(test_run){
        x11(); hist(c(condvariable1_mat_invtrans))
      }
      if(test_run2) browser()
      
      index_GOI<-which(coords_with_id[,1]=="509579")
      abs_diff_x<-abs(coords_with_id[index_GOI,2]-x.seq.cond)
      index_min_abs_diff_x<-which(min(abs_diff_x)==abs_diff_x)
      
      abs_diff_y<-abs(coords_with_id[index_GOI,3]-y.seq.cond)
      index_min_abs_diff_y<-which(min(abs_diff_y)==abs_diff_y)
      
      # condvariable1_mat has mean of the transformed data
      sim_value_at_GOI<-condvariable1_mat[index_min_abs_diff_x,index_min_abs_diff_y]-mean_trans_data # adjusted for mean
      
      # rf_df$data has mean of zero
      actual_value_at_GOI<-rf_df$data[index_GOI]
      cat("sim_value_at_GOI =",sim_value_at_GOI,"| actual_value_at_GOI =",actual_value_at_GOI,"\n") 
      
      # cookie cut condvariable1_mat_invtrans
      # TODO: export some data for presentation
      if(test_run2) browser()
      
      condvariable1_mat_invtrans_masked<-matrix(NA,nrow=nrow(condvariable1_mat_invtrans),ncol=ncol(condvariable1_mat_invtrans))
      condvariable1_mat_invtrans_masked[index_subcat]<-condvariable1_mat_invtrans[index_subcat]
      if(test_run){
        x11(); image(condvariable1_mat_invtrans_masked,col = (heat.colors(12)))
      }
      
      # Plotting of masked rainfall field
      
      if(paper_plot2){
        all_condvariable1_mat_invtrans_masked[[length(all_condvariable1_mat_invtrans_masked)+1]]<-condvariable1_mat_invtrans_masked
      }
      if(test_run){
        png(paste0("C:/Users/kim079/Documents/PhD/bates/rainfall_generation/output/i",i,"_s",s,".png"),
            width = 570, height = 600)
        par(mar=c(0, 4.1, 4.1, 2.1))
        layout(1:2,heights=c(3,1))
        image(condvariable1_mat_invtrans_masked,col = (heat.colors(12)), axes=F)
        par(mar=c(5.1, 4.1, 4.1, 2.1))
        image.scale(condvariable1_mat_invtrans_masked,col = (heat.colors(12)), xlab="Rainfall (mm)", ylab="")
        dev.off()
      }
      
      mean_cat_rain<-mean(condvariable1_mat_invtrans_masked,na.rm=T)
      daily_ensemble_rain<-c(daily_ensemble_rain,mean_cat_rain)
      
      if(test_run | test_run4){
        # rf_df$data<-rf_df$data+mean(data_df[,i])
        RFoptions(seed=NA)
        cond.all <- RFsimulate(model, x=x.seq.cond.all, y=y.seq.cond.all, data=rf_df, n=4)
        cond.all$variable1.n1<-cond.all$variable1.n1+mean_trans_data
        if(test_run){
          x11(); plot(cond.all)
        }
        
        
        cond.all$variable1.n1
        # refit the semivariogram to see if we get the same answer back
        # Should do this for several realisations
        cond.all.variable1_mat<-matrix(cond.all$variable1.n1,nrow=length(x.seq.cond.all))
        
        cond.all.variable1_mat_adjusted<-t(apply(t(cond.all.variable1_mat), 2, rev))
        if(test_run){
          x11(); image(cond.all.variable1_mat_adjusted,col = (heat.colors(12)))
        }
        
        # get rid of negatives
        cond.all.variable1_mat_adjusted_tmp<-cond.all.variable1_mat_adjusted
        cond.all.variable1_mat_adjusted_tmp[cond.all.variable1_mat_adjusted_tmp<=0]<-0
        
        cond.all.variable1_mat_adjusted_invtrans<-inverse_transform_data(cond.all.variable1_mat_adjusted_tmp,bc2$lambda)
        if(test_run){
          x11(); image(cond.all.variable1_mat_adjusted_invtrans,col = (heat.colors(12)))
        }
        
        # TODO: check that below should be redundant
        # eliminate the negatives and NaNs
        # NaNs may be high numbers as well!!
        index_nans<-which(is.nan(cond.all.variable1_mat_adjusted_invtrans),arr.ind = T)
        if(length(index_nans)>0) stop("shouldn't be NaNs")
        
        index_small<-which(cond.all.variable1_mat_adjusted[index_nans] < min(data_df_trans))
        cond.all.variable1_mat_adjusted_invtrans[index_nans][index_small] <- min(cond.all.variable1_mat_adjusted_invtrans,na.rm=T)
        
        index_big<-which(cond.all.variable1_mat_adjusted[index_nans] > max(data_df_trans))
        cond.all.variable1_mat_adjusted_invtrans[index_nans][index_big] <- max(cond.all.variable1_mat_adjusted_invtrans,na.rm=T)
        
        # cond.all.variable1_mat_adjusted_invtrans[is.nan(cond.all.variable1_mat_adjusted_invtrans)]<-min(cond.all.variable1_mat_adjusted_invtrans,na.rm=T)
        if(length(which(cond.all.variable1_mat_adjusted_invtrans<0))>0) stop("shouldn't have negatives here")
        cond.all.variable1_mat_adjusted_invtrans[cond.all.variable1_mat_adjusted_invtrans<0]<-0
        if(test_run){
          x11(); image(cond.all.variable1_mat_adjusted_invtrans,col = (heat.colors(12)))
        }
        cat("all sim range:",range(cond.all.variable1_mat_adjusted_invtrans),"\n")
        
        # ind.bad<-which(cond.all.variable1_mat_adjusted_invtrans>20000,arr.ind=T)
        # cond.all.variable1_mat_adjusted_invtrans[ind.bad]
        # cond.all.variable1_mat_adjusted[ind.bad]
        if(test_run2) browser()
        
        # This re-transforms the data after eliminating the zeros to see if there is an effect on the re-fitted
        # variogram etc
        cond.all.variable1_mat_adjusted_invtrans_trans<-transform_data(cond.all.variable1_mat_adjusted_invtrans,bc2$lambda)
        
        # Remember: The matrices are rotated to the left when actually projected onto an image
        # matrices: x's top to bottom; y's left to right
        cond.all.variable1_v<-c(cond.all.variable1_mat_adjusted_invtrans_trans) # vectors are by column (all x's)
        
        # Perform the normal fitting and resampling of truncated normal to cond.all.variable1_v
        cond.all.variable1_v<-replace_zeros(cond.all.variable1_v)
        if(test_run){
          x11(); hist(cond.all.variable1_v)
        }
        
        cond.all.variable1_v_noelim<-c(cond.all.variable1_mat_adjusted) # vectors are by column (all x's)
        if(test_run){
          x11()
          hist(cond.all.variable1_v_noelim)
        }
        if(test_run2) browser()
        
        all_y.seq.cond.all<-sort(rep(y.seq.cond.all,length(x.seq.cond.all)))
        all_x.seq.cond.all<-rep(x.seq.cond.all,length(y.seq.cond.all))
        
        # # below is a different way of working out above but I can't remember how
        # cond.all.variable1_v<-as.vector(t(cond.all.variable1_mat))
        # all_y.seq.cond.all<-rep(y.seq.cond.all,length(x.seq.cond.all))
        # all_x.seq.cond.all<-sort(rep(x.seq.cond.all,length(y.seq.cond.all)))
        
        coords_all<-cbind(all_x.seq.cond.all,all_y.seq.cond.all)
        refit_rf_df<-RFspatialPointsDataFrame(coords=coords_all,data=cond.all.variable1_v-mean(cond.all.variable1_v),RFparams = list(n = 1, vdim = 1))
        refit_rf_df_noelim<-RFspatialPointsDataFrame(coords=coords_all,data=cond.all.variable1_v_noelim-mean(cond.all.variable1_v_noelim),RFparams = list(n = 1, vdim = 1))
        # estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=NA)
        # estmodel <- RMgencauchy(var=globalvariance, scale=NA, alpha=NA, beta=NA)
        # estmodel <- RMgencauchy(scale=NA, alpha=NA, beta=NA)
        
        try_error<-try({
          # RFoptions(bins=20)
          # RFoptions(bins=5)
          # refit_rff<-RFfit(estmodel, data=refit_rf_df) #,printlevel=4
          # refit_rff_att<-attributes(refit_rff)
          # refit_rff_att$Z
          # x11(); plot(refit_rff)
          emp.variog_rff<-variog(coords=coords_all,data=cond.all.variable1_v-mean(cond.all.variable1_v),op="cloud")
          # x11(); plot(emp.variog_rff)
          # opt_refit<-optim(par=c(scale,globalvariance),fn=spherical_semivariogram_optim_func,dist=emp.variog_rff$u,semivar_obs=emp.variog_rff$v)
          opt_refit<-hydromad::SCEoptim(par=c(scale,globalvariance),FUN=spherical_semivariogram_optim_func,dist=emp.variog_rff$u,semivar_obs=emp.variog_rff$v,
                                        lower=c(0,0))
          all_sim_pars<-rbind(all_sim_pars,opt_refit$par)
          
          # refit_rff.tmp<-RFfit(estmodel, data=refit_rf_df,upper=c(100000,20,100)) #,printlevel=4
          # refit_rff_noelim<-RFfit(estmodel, data=refit_rf_df_noelim)
          emp.variog_rff_noelim<-variog(coords=coords_all,data=cond.all.variable1_v_noelim-mean(cond.all.variable1_v_noelim),op="cloud")
          # x11(); plot(emp.variog_rff_noelim)
          # opt_refit_noelim<-optim(par=c(scale,globalvariance),fn=spherical_semivariogram_optim_func,dist=emp.variog_rff_noelim$u,semivar_obs=emp.variog_rff_noelim$v)
          opt_refit_noelim<-hydromad::SCEoptim(par=c(scale,globalvariance),FUN=spherical_semivariogram_optim_func,dist=emp.variog_rff_noelim$u,semivar_obs=emp.variog_rff_noelim$v,
                                               lower=c(0,0))
          
        }, silent = T)
        
        if(length(grep("Error",try_error))==0){ #class(try_error)[1]=="RFfit" 
          # if(test_run) x11(); plot(rff)
          # # x11(); plot(rff,ylim=c(0,0.1),xlim=c(0,20000))
          # if(test_run) x11(); plot(refit_rff_noelim)
          # # x11(); plot(refit_rff_noelim,ylim=c(0,0.2),xlim=c(0,20000))
          # # x11(); plot(rff_orig)
          # if(test_run) x11(); plot(refit_rff)
          # # x11(); plot(refit_rff,ylim=c(0,0.2),xlim=c(0,20000))
          
          # x11(); plot(v2)
          # lines(x=dist_line,y=sim_semivar,col=2)
          if(test_run){
            x11(); plot(variog(coords=locations_easting_northing,data=data_df_trans-mean_trans_data))
            lines(x=dist_line,y=sim_semivar,col=2)
          }

          # x11(); plot(emp.variog_rff_noelim)
          sim_semivar_rff_noelim<-spherical_semivariogram(dist_line, opt_refit_noelim$par[1], opt_refit_noelim$par[2])
          # lines(x=dist_line,y=sim_semivar_rff_noelim,col=2)
          if(test_run){
            x11(); plot(variog(coords=coords_all,data=cond.all.variable1_v_noelim-mean(cond.all.variable1_v_noelim)))
            lines(x=dist_line,y=sim_semivar_rff_noelim,col=2)
          }
                    
          # x11(); plot(emp.variog_rff)
          sim_semivar_rff<-spherical_semivariogram(dist_line, opt_refit$par[1], opt_refit$par[2])
          # lines(x=dist_line,y=sim_semivar_rff,col=2)
          if(test_run){
            x11(); plot(variog(coords=coords_all,data=cond.all.variable1_v-mean(cond.all.variable1_v)))
            lines(x=dist_line,y=sim_semivar_rff,col=2)
          }

          if(test_run2) browser()
            # attributes(refit_rff)
          # globalvariance_refit<-attributes(attributes(refit_rff)$ml)$globalvariance
          # scale_refit<-attributes(attributes(refit_rff)$ml)$param[1,1]
          globalvariance_refit<-opt_refit$par[2]
          scale_refit<-opt_refit$par[1]
          model_refit<-RMspheric(var=globalvariance_refit, scale=scale_refit)
          modsim_refit<-RFsimulate(model_refit,x=seq(0,1000000,by=1))
          
          # globalvariance_refit_noelim<-attributes(attributes(refit_rff_noelim)$ml)$globalvariance
          # scale_refit_noelim<-attributes(attributes(refit_rff_noelim)$ml)$param[1,1]
          globalvariance_refit_noelim<-opt_refit_noelim$par[2]
          scale_refit_noelim<-opt_refit_noelim$par[1]
          model_refit_noelim<-RMspheric(var=globalvariance_refit_noelim, scale=scale_refit_noelim)
          modsim_refit_noelim<-RFsimulate(model_refit_noelim,x=seq(0,1000000,by=1))
          if(test_run3){
            x11()
            layout(1:3)
            # original
            modsim_orig <- RFsimulate(model, x=seq(0,1000000,by=1))
            plot(attributes(modsim_orig)$data[,1],type="l")
            
            plot(attributes(modsim_refit_noelim)$data[,1],type="l")
            
            # refitted
            plot(attributes(modsim_refit)$data[,1],type="l")
            mean(attributes(modsim_refit)$data[,1])
          }

          
          # need to transform rainfall somehow so that only positives are produced (keep mean at zero?) 
          # then we can check if things are ok if the variograms can be reproduced (more or less) - done
          if(test_run2) browser()
          
        } else {
          cat("someting went wrong with refitting \n")
          stop()
        }
      }
    } else {
      # just take all zeros
      condvariable1_mat_invtrans<-matrix(0,nrow=length(x.seq.cond),ncol=length(y.seq.cond))
      condvariable1_mat_invtrans<-t(apply(t(condvariable1_mat_invtrans), 2, rev))
      if(test_run){
        x11(); image(condvariable1_mat_invtrans,col = (heat.colors(12)))
      }
      if(test_run2) browser()
      
      condvariable1_mat_invtrans_masked<-matrix(NA,nrow=nrow(condvariable1_mat_invtrans),ncol=ncol(condvariable1_mat_invtrans))
      condvariable1_mat_invtrans_masked[index_subcat]<-condvariable1_mat_invtrans[index_subcat]
      if(test_run){
        x11(); image(condvariable1_mat_invtrans_masked,col = (heat.colors(12)))
      }
      
      mean_cat_rain<-mean(condvariable1_mat_invtrans_masked,na.rm=T)
      daily_ensemble_rain<-c(daily_ensemble_rain,mean_cat_rain)
      
    }
    if(test_run) graphics.off()
  }
  
  all_ensemble_rain<-rbind(all_ensemble_rain,daily_ensemble_rain)
  
  if(length(which(data_df[,i]>0))>0){
    if(test_run4){
      
      #x11()
      png(paste0("C:/Users/kim079/Documents/PhD/bates/rainfall_generation/output/variog_i",i,".png"),
          width = 570, height = 600)
      plot(variog(coords=locations_easting_northing,data=data_df_trans-mean_trans_data,op="cloud"))
      all_sim_semivar_rff<-NULL
      for(pp in 1:nrow(all_sim_pars)){
        sim_semivar_rff<-spherical_semivariogram(dist_line, all_sim_pars[pp,1], all_sim_pars[pp,2])
        lines(x=dist_line,y=sim_semivar_rff,col=2)
        all_sim_semivar_rff<-rbind(all_sim_semivar_rff,sim_semivar_rff)
      }
      for(ass in 1:nrow(all_sim_semivar)){
        lines(x=dist_line,y=all_sim_semivar[ass,],col=1,lwd=3)
      }
      
      dev.off()
      if(test_run2) browser()
      # browser()
      all_sim_pars_export<-cbind(all_opt_pars,all_sim_pars)
      all_sim_pars_export<-cbind(1:nrow(all_sim_pars),all_sim_pars_export)
      colnames(all_sim_pars_export)<-c("ensemble","range.orig","scale.orig","range.sim","scale.sim")
      write.csv(all_sim_pars_export,paste0("C:/Users/kim079/Documents/PhD/bates/rainfall_generation/output/semivar_params_recalc_i",i,".csv"),
                quote = F,row.names = F)
      
      sim_semivar_rff_export<-cbind(all_sim_semivar,all_sim_semivar_rff)
      sim_semivar_rff_export<-cbind(1:nrow(sim_semivar_rff_export),sim_semivar_rff_export)
      colnames(sim_semivar_rff_export)<-c("ensemble",paste0("orig_dist_",dist_line),paste0("sim_dist_",dist_line))
      write.csv(sim_semivar_rff_export,paste0("C:/Users/kim079/Documents/PhD/bates/rainfall_generation/output/semivar_simulations_i",i,".csv"),
                quote = F,row.names = F)
      system(paste0("C:/Users/kim079/Documents/model_optimisation_framework/software/gzip/gzip.exe -f C:/Users/kim079/Documents/PhD/bates/rainfall_generation/output/semivar_simulations_i",i,".csv"))
    }

  }
}

if(paper_plot){
  jpeg("C:/Users/kim079/Documents/PhD/state_uncertainty_paper1/variogram_fitting.jpg",width=20000,height=10000,pointsize=20,res=1000)
  par(mar=c(5.1, 5.1, 3.1, 2.1))
  layout(rbind(cbind(c(5,5),matrix(1:4,nrow=2,byrow=2)),c(7,6,6)),
         widths = c(0.5,4,4),heights=c(3,3,1)) #layout(matrix(1:4,nrow=2,byrow=2))
  for(i in 1:length(all_v2)){
    v2<-all_v2[[i]]
    # binning for plotting
    bin_intervals<-seq(0,max(v2$u),length.out=7)
    semivariance_bins<-list()
    for(bb in 2:length(bin_intervals)){
      bin_indices<-which(v2$u<bin_intervals[bb] & v2$u>bin_intervals[bb-1])
      v2$v[bin_indices]
      v2$u[bin_indices]
      # plot(x=v2$u[bin_indices],y=v2$v[bin_indices])
      semivariance_bins[[bb-1]]<-v2$v[bin_indices]
    }
    bin_intervals_mids<-bin_intervals[-1]
    bin_intervals_mids<-bin_intervals_mids-((bin_intervals[2]-bin_intervals[1])/2)
    semivariance_bins_mean<-sapply(semivariance_bins,mean)
    if(i==1 | i==2){
      if(i==1){
        plot(x=bin_intervals_mids,y=semivariance_bins_mean,xlim=c(0,max(bin_intervals_mids)),ylim=c(0,max(semivariance_bins_mean)),type="b",
             ylab="",xlab="",lwd=2,main=all_simulated_dates[i],cex.main=1.7,font.main=1,cex.lab=1.5,cex.axis=1.5,pch=19)
        lines(x=all_dist_line[[i]],y=all_sim_semivar_all[[i]],col=2,lwd=2)
      } else {
        plot(x=bin_intervals_mids,y=semivariance_bins_mean,xlim=c(0,max(bin_intervals_mids)),ylim=c(0,max(semivariance_bins_mean)),type="b",
             ylab="",xlab="",lwd=2,main=all_simulated_dates[i],cex.main=1.7,font.main=1,cex.lab=1.5,cex.axis=1.5,pch=19)
        lines(x=all_dist_line[[i]],y=all_sim_semivar_all[[i]],col=2,lwd=2)
      }

    } else {
      if(i==3){
        plot(x=bin_intervals_mids,y=semivariance_bins_mean,xlim=c(0,max(bin_intervals_mids)),ylim=c(0,max(semivariance_bins_mean)),type="b",
             ylab="",xlab="",lwd=2,main=all_simulated_dates[i],cex.main=1.7,font.main=1,cex.lab=1.5,cex.axis=1.5,pch=19)
        lines(x=all_dist_line[[i]],y=all_sim_semivar_all[[i]],col=2,lwd=2)
      } else {
        plot(x=bin_intervals_mids,y=semivariance_bins_mean,xlim=c(0,max(bin_intervals_mids)),ylim=c(0,max(semivariance_bins_mean)),type="b",
             ylab="",xlab="",lwd=2,main=all_simulated_dates[i],cex.main=1.7,font.main=1,cex.lab=1.5,cex.axis=1.5,pch=19)
        lines(x=all_dist_line[[i]],y=all_sim_semivar_all[[i]],col=2,lwd=2)
      }

    }

  }
  par(mar=c(0,0,0,0))
  plot(0,0,type="n",axes=F,ann=F,xlim=c(-1,1))
  text(y=0,x=0,labels="Semivariance",cex=4,pos=1,srt=90)
  
  plot(0,0,type="n",axes=F,ann=F,xlim=c(-1,1))
  text(y=0,x=0,labels="Distance (m)",cex=4,pos=3)
  legend(x="topright",legend=c("Fitted","Empirical"),lty=c(1,1),col=c(2,1),pch=c(NA,19),inset=c(0.03,0),cex=1.5,lwd=2)
  
  dev.off()
  
}


if(paper_plot2){
  library(RColorBrewer)
  # col.ramp<-colorRampPalette(brewer.pal(9,"Blues")[2:9])(100)
  # col.ramp<-colorRampPalette(brewer.pal(11,"RdYlBu")[2:9])(100)
  # col.ramp<-colorRampPalette(brewer.pal(11,"RdYlBu"))(100)
  # col.ramp<-colorRampPalette(brewer.pal(11,"RdYlBu")[c(3,11)])(100)
  # col.ramp<-colorRampPalette(brewer.pal(11,"RdYlBu"))(11)
  # col.ramp<-brewer.pal(11,"RdYlBu")
  # col.ramp<-colorRampPalette(brewer.pal(11,"RdYlBu"))(20)
  # col.ramp<-colorRampPalette(rainbow(100))
  # col.ramp<-rainbow(100)
  # col.ramp<-colorRampPalette(brewer.pal(11,"RdYlBu"))(100)
  
  col.ramp<-rev(colorRamps::matlab.like(200))[18:170]
  png_width<-20000
  png_height<-png_width * 0.666666 #png_width/ncol(all_condvariable1_mat_invtrans_masked[[1]])*nrow(all_condvariable1_mat_invtrans_masked[[1]])*0.666666
  png(paste0("C:/Users/kim079/Documents/PhD/bates/rainfall_generation/output/rainfall_simulation_example_",all_simulated_dates,".png"),
      width=png_width,height=png_height,pointsize=20,res=1000)
  par(mar=c(0, 0, 0, 0)) #par(mar=c(0, 4.1, 4.1, 2.1))
  layout(rbind(matrix(1:6,nrow=2,byrow = T),c(7,7,7)),height=c(1,1,0.3))
  zlim<-range(unlist(all_condvariable1_mat_invtrans_masked),na.rm = T)
  for(i in 1:length(all_condvariable1_mat_invtrans_masked)){
    image(all_condvariable1_mat_invtrans_masked[[i]],col = col.ramp, axes=F, zlim=zlim) #(heat.colors(12))
    text(x=0.5,y=0.99,labels=paste("Member",i),cex=2)
  }
  par(mar=c(5, 2.2, 0, 2.2))
  # par(mar=c(0,0,0,0))
  image.scale(zlim,col = col.ramp, xlab="Rainfall (mm)", ylab="",cex.lab=1.7,cex.axis=1.5)
  dev.off()
}


stop()


all_ensemble_rain_mean<-apply(all_ensemble_rain,1,mean)
apply(all_ensemble_rain,2,mean)

plot(all_ensemble_rain_mean,type="l")
hist(all_ensemble_rain[1,])
hist(all_ensemble_rain[2,])
# hist(all_ensemble_rain[19,])

out_filename<-"C:/Users/kim079/Documents/PhD/bates/rainfall_generation/output/ensemble_rainfall.csv"
all_ensemble_rain<-cbind(all_sites_data_good_longest[1:nrow(all_ensemble_rain),1],all_ensemble_rain)

write.csv(all_ensemble_rain,out_filename,quote = F,row.names = F)


all_ensemble_rain<-read.csv(out_filename,as.is=T)

rownames(all_ensemble_rain)<-NULL
sds_each_ts<-apply(all_ensemble_rain[,-1],1,function(x) sd(as.numeric(x)))
plot(sds_each_ts,type="l")
plot(sds_each_ts,type="l",log="y")
hist(sds_each_ts)

all_ensemble_rain_numeric<-t(apply(all_ensemble_rain[,-1],1,function(x) as.numeric(x)))
hist(c(all_ensemble_rain_numeric))
hist(c(all_ensemble_rain_numeric)[which(c(all_ensemble_rain_numeric)>0)])

hist(c(all_ensemble_rain_numeric)^(1/3))
hist(log(c(all_ensemble_rain_numeric)))

# boxcox transform
bc_transform_data<-function(q,lambda){
  if(lambda[1]!=0){
    return(((((q)+lambda[2])^lambda[1])-1)/lambda[1])
  } else {
    return(log((q)+lambda[2]))
  }
}
inverse_bc_transform_data<-function(q_trans,lambda){
  if(lambda[1]!=0){
    return(((q_trans*lambda[1]+1)^(1/lambda[1]))-lambda[2])
  } else {
    return(exp(q_trans)-lambda[2])
  }
}

# including zeros
bc2<-boxcoxfit(c(all_ensemble_rain_numeric),lambda2 = T)
all_ensemble_rain_numeric_trans<-bc_transform_data(c(all_ensemble_rain_numeric),bc2$lambda)
inverse_bc_transform_data(all_ensemble_rain_numeric_trans,bc2$lambda)

hist(all_ensemble_rain_numeric_trans)
# check
range(c(all_ensemble_rain_numeric)-inverse_bc_transform_data(all_ensemble_rain_numeric_trans,bc2$lambda))
inverse_bc_transform_data(3,bc2$lambda)

# excluding zeros
all_ensemble_rain_numeric_v_exczero<-c(all_ensemble_rain_numeric)[c(all_ensemble_rain_numeric)>0]
bc2<-boxcoxfit(all_ensemble_rain_numeric_v_exczero,lambda2 = T)
all_ensemble_rain_numeric_trans<-bc_transform_data(all_ensemble_rain_numeric_v_exczero,bc2$lambda)
inverse_bc_transform_data(all_ensemble_rain_numeric_trans,bc2$lambda)

hist(all_ensemble_rain_numeric_trans)
# check
range(all_ensemble_rain_numeric_v_exczero-inverse_bc_transform_data(all_ensemble_rain_numeric_trans,bc2$lambda))
inverse_bc_transform_data(-2,bc2$lambda)

# get stats for each day
daily_ensemble_trans_mean<-c()
daily_ensemble_trans_sd<-c()
daily_ensemble_trans_bcparam1<-c()
daily_ensemble_trans_bcparam2<-c()
for(i in 1:nrow(all_ensemble_rain_numeric)){
  cat(i,"/",nrow(all_ensemble_rain_numeric),"\n")
  if(all(all_ensemble_rain_numeric[i,]<1e-9)){
    daily_ensemble_trans_mean<-c(daily_ensemble_trans_mean,0)
    daily_ensemble_trans_sd<-c(daily_ensemble_trans_sd,0)
    daily_ensemble_trans_bcparam1<-c(daily_ensemble_trans_bcparam1,1)
    daily_ensemble_trans_bcparam2<-c(daily_ensemble_trans_bcparam2,1)
  } else {
    # hist(all_ensemble_rain_numeric[i,])
    try<-try({bc2<-list(lambda=c(boxcoxfit(all_ensemble_rain_numeric[i,])$lambda,0))},silent=T)
    if(length(grep("Error",try))>0){
      bc2<-boxcoxfit(all_ensemble_rain_numeric[i,],lambda2 = T)
      
    }
    ensemble_rain_numeric_trans<-bc_transform_data(all_ensemble_rain_numeric[i,],bc2$lambda)
    # hist(ensemble_rain_numeric_trans)
    daily_ensemble_trans_mean<-c(daily_ensemble_trans_mean,mean(ensemble_rain_numeric_trans))
    daily_ensemble_trans_sd<-c(daily_ensemble_trans_sd,sd(ensemble_rain_numeric_trans))
    daily_ensemble_trans_bcparam1<-c(daily_ensemble_trans_bcparam1,bc2$lambda[1])
    daily_ensemble_trans_bcparam2<-c(daily_ensemble_trans_bcparam2,bc2$lambda[2])
  }

}

all_stats<-cbind(all_ensemble_rain[,1],daily_ensemble_trans_mean,daily_ensemble_trans_sd,daily_ensemble_trans_bcparam1,daily_ensemble_trans_bcparam2)
colnames(all_stats)<-c("date","mean","sd","bcparam1","bcparam2")
write.csv(all_stats,"C:/Users/kim079/Documents/PhD/bates/rainfall_generation/output/ensemble_rainfall_stats.csv",quote=F,row.names = F)



# some check plots
rain_0.05<-apply(all_ensemble_rain_numeric,1,quantile,probs=0.05)
rain_0.95<-apply(all_ensemble_rain_numeric,1,quantile,probs=0.95)
plot(rain_0.05,type="l")
plot(rain_0.95,type="l")
summary_ts<-cbind(rain_0.05,as.numeric(data_df[index_GOI,]),rain_0.95)

rain_min<-apply(all_ensemble_rain_numeric,1,min)
rain_max<-apply(all_ensemble_rain_numeric,1,max)
summary_ts<-cbind(rain_min,as.numeric(data_df[index_GOI,]),rain_max)
poly_y<-c(rain_max,rev(rain_min))
poly_x<-c(1:length(rain_max),length(rain_max):1)

plot_index<-1:nrow(summary_ts)
# plot_index<-1:2
to_plot<-summary_ts[plot_index,]
matplot(to_plot,type="l",col="white")
poly_y<-c(rain_max[plot_index],rev(rain_min[plot_index]))
poly_x<-c(1:length(rain_max[plot_index]),length(rain_max[plot_index]):1)
polygon(y=poly_y,x=poly_x,col="grey",border="grey")

# lines(as.numeric(data_df[index_GOI,plot_index]),col=2,lty=2)
lines(apply(all_ensemble_rain_numeric,1,median)[plot_index],col=2,lty=2)


stop()

model <- RMgencauchy(alpha=1, beta=2)
d <- RFsimulate(model, x=x, y=y, grid=FALSE, n=100) #1000
d$variable1.n1

d <- RFsimulate(model, x=x, y=y, grid=FALSE, n=2) #1000
d$variable1


x11(); plot(d)

#########################################################
## estimation; 'NA' means: "to be estimated"           ##
estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=2) +
  RMtrend(mean=NA)
RFfit(estmodel, data=d)


#########################################################
## coupling alpha and beta                             ##
estmodel <- RMgencauchy(var=NA, scale=NA, alpha=NA, beta=NA) + 
  RMtrend(NA)
RFfit(estmodel, data=d, transform = NA) ## just for information
trafo <- function(a) c(a[1], rep(a[2], 2))
fit <- RFfit(estmodel, data=d,
             transform = list(c(TRUE, TRUE, FALSE), trafo))
print(fit)
print(fit, full=TRUE)












