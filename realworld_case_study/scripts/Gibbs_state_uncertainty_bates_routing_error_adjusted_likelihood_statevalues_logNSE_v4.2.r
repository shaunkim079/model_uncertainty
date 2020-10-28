
## MCMC AM Algorithm Shell: READ ME --------------------------------------#
#This script provides a basic shell of the Adaptive Metropolis
#algorithm (Haario et al, 2001). There are several locations where
#user-specific entries must be made (e.g., global variables, input
#data, parameter initialization, priors, simulation model function,
#likelihood function, etc.) to replace current place-holder code. Note
#that any parameter constraints (e.g., PAR1 > 0) must be built into the
#simulation model directly.


# perform_restart<-F
# if(perform_restart){
#   # RRparam_file<-"output/state_uncertainty/AM/bates/RRparameters_bates_30_1515467127.7433.csv.gz"
#   init_cov_file<-"output/state_uncertainty/AM/bates/CovPar_bates_30_1515467127.7433.csv.gz"
#   # SD2<-0.00857755438888274
#   SD2_file<-"output/state_uncertainty/AM/bates/SD2_bates_30_1515467127.7433.csv.gz"
#   theta_file<-"output/state_uncertainty/AM/bates/theta_bates_30_1515467127.7433.csv.gz"
#   restart_number<-500000 # this is the number of previous iterations to use
# } else {
#   rm(RRparam_file)
#   rm(init_cov_file)
#   rm(SD2_file)
#   rm(theta_file)
#   rm(restart_number)
#   
# }


if(length(commandArgs(TRUE))!=0){
  arguments <- commandArgs(TRUE)
  cat("arguments:",arguments,"\n")
  
  start_ts<-as.numeric(arguments[[1]])
  length_ts<-as.numeric(arguments[[2]])
  if(length(arguments)>2) initial_seed<-as.numeric(arguments[[3]])
  if(length(arguments)>3) output_dir<-as.character(arguments[[4]])
  if(length(arguments)>4) initial_state_errors_initial_file<-as.character(arguments[[5]])
  if(length(arguments)>5) state_errors_initial_file<-as.character(arguments[[6]])
  
  if(length(arguments)>6) init_cov_file<-as.character(arguments[[7]])
  if(length(arguments)>7) SD2_file<-as.character(arguments[[8]])
  if(length(arguments)>8) theta_file<-as.character(arguments[[9]])
  if(length(arguments)>9) restart_number<-as.numeric(arguments[[10]])
  if(length(arguments)>10) seed<-as.numeric(arguments[[11]])
  if(length(arguments)>11) stop_update_covariance<-as.numeric(arguments[[12]])
  if(length(arguments)>12) ITER<-as.numeric(arguments[[13]])
}

if(exists("initial_state_errors_initial_file")) rm("initial_state_errors_initial_file")
if(exists("state_errors_initial_file")) rm("state_errors_initial_file")

if(exists("init_cov_file")){
  if(is.na(init_cov_file) | init_cov_file=="NA") rm("init_cov_file")
}

if(exists("SD2_file")){
  if(is.na(SD2_file) | SD2_file=="NA") rm("SD2_file")
}

if(exists("restart_number")){
  if(is.na(restart_number) | restart_number=="NA") rm("restart_number")
}

if(Sys.info()["sysname"]=="Linux"){
  wd<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/model_optimisation_framework_v2"
  ensemble_rainfall_file<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/PhD/bates/rainfall_generation/output/ensemble_rainfall.csv"
  ensemble_rainfall_stats_file<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/PhD/bates/rainfall_generation/output/ensemble_rainfall_stats.csv"
  PET_file<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/PhD/bates/silo/116053255.csv"
  stage_data_file<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/PhD/bates/from_justin/Bates Stage data/Bates Stage data.csv"
  gauge_rating_data_file<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/PhD/bates/from_justin/Bates gauging details.csv"
} else {
  # wd<-"C:/Users/kim079/Documents/model_optimisation_framework"
  # ensemble_rainfall_file<-"C:/Users/kim079/Documents/PhD/bates/rainfall_generation/output/ensemble_rainfall.csv"
  # ensemble_rainfall_stats_file<-"C:/Users/kim079/Documents/PhD/bates/rainfall_generation/output/ensemble_rainfall_stats.csv"
  # PET_file<-"C:/Users/kim079/Documents/PhD/bates/silo/116053255.csv"
  # stage_data_file<-"C:/Users/kim079/Documents/PhD/bates/from_justin/Bates Stage data/Bates Stage data.csv"
  # gauge_rating_data_file<-"C:/Users/kim079/Documents/PhD/bates/from_justin/Bates gauging details.csv"
  wd<-"Y:/"
  ensemble_rainfall_file<-"C:/Users/kim079/Documents/PhD/bates/rainfall_generation/output/ensemble_rainfall.csv"
  ensemble_rainfall_stats_file<-"C:/Users/kim079/Documents/PhD/bates/rainfall_generation/output/ensemble_rainfall_stats.csv"
  PET_file<-"C:/Users/kim079/Documents/PhD/bates/silo/116053255.csv"
  stage_data_file<-"C:/Users/kim079/Documents/PhD/bates/from_justin/Bates Stage data/Bates Stage data.csv"
  gauge_rating_data_file<-"C:/Users/kim079/Documents/PhD/bates/from_justin/Bates gauging details.csv"
}
# initial_state_errors_initial_file<-"output/state_uncertainty/AM/bates/initial_state_errors_initial.csv"
# state_errors_initial_file<-"output/state_uncertainty/AM/bates/state_errors_initial.csv"
# if(!exists("output_dir")) output_dir<-"output/state_uncertainty/AM/bates/adjusted_likelihood_v9.0"
if(!exists("output_dir")) output_dir<-"Y:/output/state_uncertainty/AM/bates/adjusted_likelihood_statevalues_logNSE_v4.1"
RRparam_file<-"output/state_uncertainty/AM/bates/gr4j_params_logNSE.csv"


if(!exists("length_ts")) length_ts<-30 #365
if(!exists("start_ts")) start_ts<-1
if(!exists("initial_seed")) initial_seed<-1
initial_state_normaliser<-1 #0.001
initial_state_normaliser_R<-0.01
# seed<-1496902861.98095 
if(!exists("ITER")) ITER<-1000000
stop_update_covariance<-1 # stop_update_covariance<-0

export_diagnostic<-F
if(Sys.info()["sysname"]=="Linux"){
  export_diagnostic<-T
}

setwd(wd)
dir.create(output_dir,showWarnings = F)

if(exists("SD2_file")) SD2<-as.numeric(readLines(SD2_file))

# if(perform_restart | exists("theta_file")){
#   rm(initial_state_errors_initial_file)
#   rm(state_errors_initial_file)
# }

# source("scripts/gr4j_sma.r")
# source("packages/gr4j/R/gr4j_sma.r")
# if(!is.loaded("sma_gr4j_sk")){
#   if(Sys.info()["sysname"]=="Linux"){
#     dyn.load("packages/gr4j/src/gr4j.so")
#   } else {
#     dyn.load("packages/gr4j/src/gr4j.dll")
#   }
# }
source("packages/gr4j_with_routing_error/R/gr4j_sma_routing.r")
if(!is.loaded("routing_gr4j_sk")){
  if(Sys.info()["sysname"]=="Linux"){
    dyn.load("packages/gr4j_with_routing_error/src/gr4j_with_routing_error.so")
  } else {
    dyn.load("packages/gr4j_with_routing_error/src/gr4j_with_routing_error.dll")
  }
}
# source("scripts/state_uncertainty_trial_gr4j_sma_likelihood_prior_AM.r")
# source("scripts/generate.rating.curve.error.r")
library(mvtnorm)
library(MASS)
library(lattice)
#library(mvnfast)
#library(parallel)

library(condMVNorm)


sd_zero_mean<-function(x) sqrt(mean(x^2))
normalise<-function(x,factor=1){
  x*factor
  #log10(x+1000)
  #(x^factor)-1)/factor
  #trans<-log10(x+)
  #return(trans)
}
inv.normalise<-function(x,factor=1){
  x/factor
  #(10^x)-1000
  #(x*factor+1)^(1/factor)
  #(10^x)
}

# # normalisation working
# normalise(c(actual_initial_state,rnorm(length(input_trial)-1,mean=0,sd=1e-30),input_error))
# aa<-1e-30
# aa<-log10(1e-30)
# 10^aa
# bb<-1e-30
# aa<-normalise(bb)
# inv.normalise(aa)
# bb<-normalise(aa,0)
# inv.normalise(bb,1)


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


# from PyMC3 (metropolis.py)
tune<-function(scale, acc_rate){
  if(acc_rate < 0.001){
    # reduce by 90%
    scale<-scale*0.1
  } else if(acc_rate < 0.05){
    # reduce by 50%
    scale<-scale*0.5
  } else if(acc_rate < 0.2){
    # reduce by 10%
    scale<-scale*0.9
  } else if(acc_rate > 0.95){
    # increase by factor of 10
    scale<-scale*10
  } else if(acc_rate > 0.75){
    # increase by double
    scale<-scale*2
  } else if(acc_rate > 0.5){
    # increase by 10%
    scale<-scale*1.1
  }
  return(scale)
}

set.seed(12321)

#length_ts<-7
ensemble_rainfall<-read.csv(ensemble_rainfall_file,as.is=T)
ensemble_rainfall_stats<-read.csv(ensemble_rainfall_stats_file,as.is=T)


# stage data at 10 min intervals
# TODO: need to convert this to daily
stage_data<-read.csv(stage_data_file,as.is=T,skip=3)
stage_data_dates<-as.POSIXct(stage_data$X,format="%d/%m/%Y %H:%M")

PET<-read.csv(PET_file,as.is=T)



dates_of_interest<-ensemble_rainfall_stats$date[start_ts:(start_ts+length_ts-1)]
PET$Date[match(dates_of_interest,PET$Date)]

# combined<-read.csv(ts_file,as.is=T,comment="#")
# true_input_trial<-combined$rainfall.river[start_ts:(start_ts+length_ts-1)]
E_input<-PET$APET_mm[match(dates_of_interest,PET$Date)]
#true_input_trial<-runif(length_ts,0,20) #runif(length_ts,100,112)
#E_input<-runif(length_ts,0,20) #runif(length_ts,100,105)
#input_error_sd<-1e-1 #1e-1



first_index<-which(as.character(as.Date(dates_of_interest[1]))==substr(as.character(stage_data_dates),1,10))[1]
stage_data[first_index,]

last_index<-which(as.character(as.Date(dates_of_interest[length_ts]))==substr(as.character(stage_data_dates),1,10))
last_index<-last_index[length(last_index)]

stage_data_of_interest<-stage_data[first_index:last_index,]


set.seed(12321)
# input_error<-rnorm(length(true_input_trial),mean=0,sd=input_error_sd) #sd=0.5
# input_trial<-true_input_trial+input_error
# input_trial<-pmax(input_trial,0)

input_trial<-c()
for(i in 1:nrow(ensemble_rainfall_stats)){
  input_trial<-c(input_trial,inverse_bc_transform_data(ensemble_rainfall_stats$mean[i],c(ensemble_rainfall_stats$bcparam1[i],ensemble_rainfall_stats$bcparam2[i])))
}


median_ensemble_rain<-apply(ensemble_rainfall[,-1],1,median)
layout(1:2)
plot(median_ensemble_rain,type="l")
plot(input_trial,type="l")

plot(x=median_ensemble_rain,y=input_trial)
lines(x=c(0,1000),y=c(0,1000),col=2)

plot(x=median_ensemble_rain,y=input_trial,ylim=c(0,1),xlim=c(0,1))
lines(x=c(0,1000),y=c(0,1000),col=2)

# cut down period
input_trial<-input_trial[start_ts:(start_ts+length_ts-1)]
ensemble_rainfall_stats<-ensemble_rainfall_stats[start_ts:(start_ts+length_ts-1),]

index_sampled_inputs<-which(input_trial!=0)
number_of_sampled_inputs<-length(index_sampled_inputs)


set.seed(12321)
# actual_state_error<-rnorm(length(input_trial)-1,mean=0,sd=state_error_sd) #1e-6 # minimum should be max_state_error_bound/max_normaliser_value

#discharge_error_sd<-0.01 #0.01

source("scripts/determine_rating_curve_uncertainty_bates.r")
discharge_error_sd<-population_flow_sd # from determine_rating_curve_uncertainty_bates.r
discharge_beta_normal<-bc2$beta.normal
discharge_lambda<-bc2$lambda

actual_obs_discharge_trans<-bc_fitted_lm(stage_data_of_interest$Inst,bc2$beta.normal)

actual_obs_discharge_untrans<-inverse_bc_transform_data(actual_obs_discharge_trans,discharge_lambda)


# actual_discharge_error<-rnorm(length(input_trial),mean=0,sd=discharge_error_sd) #sd=1e-6
# library(hydrodiy)

# true_input<-data.frame(P=true_input_trial,E=E_input)

# need to get parameters for the model
library(hydroGOF)
# optim_gr4j.sma<-function(x,input,obs,area){
#   sim<-gr4j.sma(x[1],x[2],rep(0,nrow(input)-1),input)/1000*area/86400
#   if(!is.na(sim[1])){
#     # calculate for each 10 min
#     sim<-sim*60*10
#     sim_10min<-c()
#     for(i in 1:length(sim)){
#       sim_10min<-c(sim_10min,rep(sim[i],24*6))
#     }
# 
#     out<-NSE(sim=sim_10min,obs=obs)
#   } else {
#     out<--1e30
#   }
# 
#   return(out)
# }
# 
rrm_input<-data.frame(P=input_trial,E=E_input)
# 
# # actual area for Bates catchment
area_m2<-2254067
# area_m2<-3051394 # This is wrong!!
# optim_gr4j.sma(c(1000,500),rrm_input,actual_obs_discharge_untrans,area=area_m2)
# if(!exists("rrm_opt")){
#   rrm_opt<-hydromad::SCEoptim(FUN=optim_gr4j.sma,par=c(1000,500),input=rrm_input,obs=actual_obs_discharge_untrans,area=area_m2,
#                               control=list(fnscale=-1,trace=1))
#                               # lower=c(hydromad::hydromad.options()$gr4j$x1[1],hydromad::hydromad.options()$gr4j$x1[1]),
#                               # upper=c(hydromad::hydromad.options()$gr4j$x1[2],hydromad::hydromad.options()$gr4j$x1[2]))
# }
# rrm_opt$value
# rrm_opt$par
# optim_gr4j.sma(rrm_opt$par,rrm_input,actual_obs_discharge_untrans,area=area_m2)

optim_gr4j.run<-function(x,input,obs,area,transformed=F){
  sim<-gr4j.run(x[1:4],x[5],x[6],rep(0,nrow(input)-1),input,run_compiled=T,state_error_R=rep(0,nrow(input)-1),
                transformed = transformed)/1000*area/86400
  if(!is.na(sim[1])){
    # # calculate for each 10 min
    # sim<-sim*60*10
    # sim_10min<-c()
    # for(i in 1:length(sim)){
    #   sim_10min<-c(sim_10min,rep(sim[i],24*6))
    # }
    # # calculate for each day
    # sim<-sim*86400
    # indices<-sort(rep(1:length(sim),24*6))
    # obs_agg<-aggregate(obs,by=list(indices),sum)[,2]
    # # resid_trans<-log(obs_agg+1e-30)-log(sim+1e-30)
    # resid_trans<-bc_transform_data(obs_agg,discharge_lambda)-bc_transform_data(sim,discharge_lambda)
    # out<--sum(abs(resid_trans))
    # # out<-NSE(sim=sim,obs=obs_agg)
    
    # calculate for each 10 min
    sim_flow<-sim*60*10
    sim_flow_trans<-bc_transform_data(sim_flow,discharge_lambda)
    if(all(!is.na(sim_flow_trans))){
      # agregate obs flow - assume that instantaneous flow is equivalent to 10 min volume!
      sim_flow_per_day<-sim_flow_trans*24*6
      indices<-sort(rep(1:length(sim_flow),24*6))
      obs_trans_agg<-aggregate(obs,by=list(indices),sum)[,2]
      resid_trans<-sim_flow_per_day-obs_trans_agg
      out<--sum(abs(resid_trans))
    } else {
      out<--1e30
    }
    
    
  } else {
    out<--1e30
  }
  
  return(out)
}
plot_gr4j.run<-function(x,input,obs,area){
  sim<-gr4j.run(x[1:4],x[5],x[6],rep(0,nrow(input)-1),input,run_compiled=T,state_error_R=rep(0,nrow(input)-1))/1000*area/86400
  if(!is.na(sim[1])){
    # # calculate for each 10 min
    # sim<-sim*60*10
    # sim_10min<-c()
    # for(i in 1:length(sim)){
    #   sim_10min<-c(sim_10min,rep(sim[i],24*6))
    # }
    # calculate for each day
    sim<-sim*86400
    indices<-sort(rep(1:length(sim),24*6))
    obs_agg<-aggregate(obs,by=list(indices),sum)[,2]
    out<-NSE(sim=sim,obs=obs_agg)
    matplot(cbind(obs_agg,sim),type="l")
  }
  
}


optim_gr4j.run(c(1000,1,100,2,500,100),rrm_input,actual_obs_discharge_untrans,area=area_m2)
optim_gr4j.run(c(log(1000),asinh(1),log(100),log(2-0.5),log(500),log(100)),rrm_input,actual_obs_discharge_untrans,area=area_m2,
               transformed = T)

# lower_gr4j<-c(
#   hydromad::hydromad.options()$gr4j$x1[1],
#   hydromad::hydromad.options()$gr4jrouting$x2[1],
#   hydromad::hydromad.options()$gr4jrouting$x3[1],
#   hydromad::hydromad.options()$gr4jrouting$x4[1],
#   0,
#   0
# )
# upper_gr4j<-c(
#   hydromad::hydromad.options()$gr4j$x1[2],
#   hydromad::hydromad.options()$gr4jrouting$x2[2],
#   hydromad::hydromad.options()$gr4jrouting$x3[2],
#   hydromad::hydromad.options()$gr4jrouting$x4[2],
#   hydromad::hydromad.options()$gr4j$x1[2],
#   hydromad::hydromad.options()$gr4jrouting$x3[2]
# )
lower_gr4j<-c(
  10,
  -25,
  1,
  0.1,
  10,
  1
)
upper_gr4j<-c(
  4000,
  5,
  1000,
  5,
  4000,
  500
)

if(exists("RRparam_file")){
  rrparams<-read.csv(RRparam_file,as.is=T)
  rrm_opt<-list(par=unlist(rrparams[1,]))
}

# rrm_opt<-list(par=c(1000,1,100,2,500,10))

if(!exists("rrm_opt")){
  lower_trans<-c(log(lower_gr4j[1]),asinh(lower_gr4j[2]),log(lower_gr4j[3]),log(0.5001-0.5),log(lower_gr4j[5]),log(lower_gr4j[6]))
  upper_trans<-c(log(upper_gr4j[1]),asinh(upper_gr4j[2]),log(upper_gr4j[3]),log(upper_gr4j[4]-0.5),log(upper_gr4j[5]),log(upper_gr4j[6]))

  rrm_opt<-hydromad::SCEoptim(FUN=optim_gr4j.run,par=c(1000,1,100,2,500,100),input=rrm_input,obs=actual_obs_discharge_trans,area=area_m2,
                              control=list(fnscale=-1,trace=1,ncomplex=10),
                              lower=lower_gr4j,
                              upper=upper_gr4j)
  rrm_opt2<-hydromad::SCEoptim(FUN=optim_gr4j.run,par=c(1000,1,100,2,500,100),input=rrm_input,obs=actual_obs_discharge_trans,area=area_m2,
                              control=list(fnscale=-1,trace=1,ncomplex=10),
                              lower=lower_gr4j,
                              upper=upper_gr4j)
  if(rrm_opt2$value<rrm_opt$value) rrm_opt<-rrm_opt2
  rrm_opt3<-hydromad::SCEoptim(FUN=optim_gr4j.run,par=c(log(1000),asinh(1),log(100),log(2-0.5),log(500),log(100)),input=rrm_input,obs=actual_obs_discharge_trans,area=area_m2,
                               control=list(fnscale=-1,trace=1,ncomplex=10),
                               lower=lower_trans,
                               upper=upper_trans,
                               transformed=T)
  if(rrm_opt3$value<rrm_opt$value) rrm_opt<-rrm_opt3
}

cat(rrm_opt$par,"\n")
optim_gr4j.run(rrm_opt$par,rrm_input,actual_obs_discharge_untrans,area=area_m2)
plot_gr4j.run(rrm_opt$par,rrm_input,actual_obs_discharge_untrans,area=area_m2)


actual_model_param<-rrm_opt$par[1]
actual_initial_state<-rrm_opt$par[5]
actual_initial_state_R<-rrm_opt$par[6]
all_model_params<-rrm_opt$par[1:4]

if(exists("initial_state_errors_initial_file")){
  initial_state_errors_initial<-read.csv(initial_state_errors_initial_file,as.is=T)
  
  
  actual_model_param<-rrm_opt$par[1]
  actual_initial_state<-rrm_opt$par[5]+initial_state_errors_initial[,1]
  actual_initial_state_R<-rrm_opt$par[6]+initial_state_errors_initial[,2]
  all_model_params<-rrm_opt$par[1:4]
}

if(exists("state_errors_initial_file")){
  state_errors_initial<-read.csv(state_errors_initial_file,as.is=T)
}


# tt<-gr4j.sma(actual_model_param,actual_initial_state,actual_state_error,true_input)


# rating_error_uncertainty<-generate_rating_error(length(tt),
#                                                 input_error_sd=input_error_sd,
#                                                 state_error_sd=state_error_sd)


# true_flow<-rating_error_uncertainty$flow_data_without_error
# actual_discharge_error<--(rating_error_uncertainty$flow_with_hetero-true_flow)

# actual_obs_discharge<-true_flow-actual_discharge_error

# actual_obs_discharge+actual_discharge_error


population_flow_sd_trans<-population_flow_sd
bc_lambdas<-discharge_lambda
bc_betas<-discharge_beta_normal

true_level_data<-stage_data_of_interest$Inst
rating_flow_trans<-bc_fitted_lm(true_level_data,bc_betas)


# TODO: use rating_flow_trans instead of actual_obs_flow_trans

# actual_obs_flow_trans<-log_bc_transform(actual_obs_discharge,bc_lambdas)
# plot(actual_obs_trans)

# true_residual_trans<-rating_flow_trans-actual_obs_flow_trans
# hist(true_residual_trans)
# plot(y=actual_obs_trans,x=true_obs_flow_log_bc)

cat("population_flow_sd_trans = ",population_flow_sd_trans,"\n")
# sd_zero_mean(true_residual_trans)


# cat("actual_state_error_sd=",sd_zero_mean(actual_state_error),"\n")
# cat("input_error_sd=",sqrt(mean(input_error^2)),"\n")
# cat("actual_discharge_error_sd=",sqrt(mean(actual_discharge_error^2)),"\n")

# input_error
# normalise(actual_state_error,0.01)
# normalise(actual_initial_state,0.0001)
# normalisers<-c(0.0001,rep(0.01,length_ts-1),rep(1,length_ts))


if(!exists("seed")){
  seed<-as.numeric(Sys.time())
}
# seed<-1122
# seed<-as.numeric(Sys.time())
set.seed(seed)

#initial_state_normaliser<-10^(runif(1,1e-3,1))-1
#initial_state_normaliser<-50
# ff<-c()
# for(uu in 1:1000){
#   initial_state_normaliser<-10^(runif(1,1e-3,1.65))-1
#   ff<-c(ff,initial_state_normaliser)
# }
# hist(ff)

# normalisers<-c(0.0001,rep(initial_state_normaliser,length_ts-1),rep(1,length_ts))
# normalisers<-c(0.001,rep(initial_state_normaliser,length_ts-1),rep(1,length_ts))

# TODO: change the input normaliser to the inverse sd transformed values. Do i need to abs this?
# normalisers<-c(0.0001,rep(initial_state_normaliser,length_ts-1),rep(1,length_ts),0.0001)
# normalisers<-c(0.0001,rep(initial_state_normaliser,length_ts-1),1/ensemble_rainfall_stats$sd[index_sampled_inputs],0.0001)

# The scaling assumes that the samples are in the range of -1 and +1 so the magnitude should be around 0.1
# working of the input scaling
for(i in 1:length(ensemble_rainfall_stats$sd[index_sampled_inputs])){
  trial_input<-seq(-1,1,by=0.1)/(1/ensemble_rainfall_stats$sd[index_sampled_inputs][i]/3)
  cat(trial_input,"\n")
  range_error<-ensemble_rainfall_stats$sd[index_sampled_inputs][i]*3
  cat(range_error,"\n")
  
}

# normalisers #######################################
# Super important!!
# These scaling factors are chosen so that the actual sampled range is around -1 to +1
# If the correlation matrix is showing very high correlations then this is likely the problem!!
# initial production store, production store errors, input forcing errors, initial routing store, routing store errors
normalisers<-c(1/all_model_params[1],
               rep(initial_state_normaliser,length_ts-1),
               1/ensemble_rainfall_stats$sd[index_sampled_inputs]/3,
               1/all_model_params[3])
#normalisers<-c(0.0001,rep(1e+29,length_ts-1),rep(1,length_ts))

# error_discharge_variance<-sum(actual_discharge_error^2)/(length(actual_discharge_error)-1)
# error_input_variance<-sum(input_error^2)/(length(input_error)-1)

# error_discharge_variance<-discharge_error_sd^2
population_flow_variance_trans<-population_flow_sd_trans^2
# error_input_variance<-input_error_sd^2
error_input_variance<-ensemble_rainfall_stats$sd^2

logprior_fun<-function(params,min,max,data,likelihood=NULL){
  
  return(0)
  
}

loglikelihood_fun<-function(data,params){
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  
  length_ts<-length(data$input)
  # state_normaliser<-data$normalisers[2] #10^(params[length_ts+data$number_of_sampled_inputs+2]*10)
  
  # params_unnorm<-inv.normalise(params[1:(length_ts+data$number_of_sampled_inputs+1)],
  #                              c(data$normalisers[1],
  #                                rep(state_normaliser,length_ts-1),
  #                                data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)],
  #                                data$normalisers[length_ts+data$number_of_sampled_inputs+1]))
  params_unnorm<-params[1:(length_ts+data$number_of_sampled_inputs+1)]
  initial_state<-params_unnorm[1]
  state_values<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[length_ts+data$number_of_sampled_inputs+1]
  
  
  population_flow_variance_trans<-data$population_flow_variance_trans
  # error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  # browser()
  # TODO: REALLY IMPORTANT!!!! have to have box-cox untransformed rainfall data here before the simulation!!!
  
  # transform rainfall data, introduce the error, then untransform again to simulate
  input_ts<-data$input
  input_trans<-rep(NA,data$number_of_sampled_inputs)
  input_trans_with_error<-rep(NA,data$number_of_sampled_inputs)
  input_untrans_with_error<-rep(NA,data$number_of_sampled_inputs)
  log_likelihoods_input<-rep(NA,data$number_of_sampled_inputs)
  input_not_zero<-which(input_ts!=0)
  for(ii in 1:data$number_of_sampled_inputs){
    # TODO: move trans below to data input to improve efficiency
    input_trans[ii]<-bc_transform_data(data$input[input_not_zero[ii]],
                                       c(data$ensemble_rainfall_stats$bcparam1[input_not_zero[ii]],data$ensemble_rainfall_stats$bcparam2[input_not_zero[ii]]))
    input_trans_with_error[ii]<-input_trans[ii]-input_error[ii]
    
    input_untrans_with_error[ii]<-inverse_bc_transform_data(input_trans_with_error[ii],
                                                            c(data$ensemble_rainfall_stats$bcparam1[input_not_zero[ii]],
                                                              data$ensemble_rainfall_stats$bcparam2[input_not_zero[ii]]))
    
    # # WTF is below??!!
    # log_likelihoods_input[ii]<-dnorm(data$error_input_variance[input_not_zero[ii]],
    #                                  mean=input_trans[ii],
    #                                  sd=sqrt(data$error_input_variance[input_not_zero[ii]]),log = T)
    log_likelihoods_input[ii]<-dnorm(input_trans_with_error[ii],
                                     mean=input_trans[ii],
                                     sd=sqrt(data$error_input_variance[input_not_zero[ii]]),log = T)
  }
  # log_likelihoods_input<-log_likelihoods_input[is.finite(log_likelihoods_input)]
  
  input_ts[which(input_ts!=0)]<-input_untrans_with_error
  
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.run(param=data$all_model_params, initial_state_S=initial_state, initial_state_R=init_state_R, 
                      state_error=rep(0,length_ts-1), state_S_ts = state_values, input=input, state_error_R=rep(0,length_ts-1))/1000*data$area_m2/86400
  if(!is.na(model_run[1])){
    
    # calculate for each 10 min
    sim_flow<-model_run*60*10
    # sim_10min<-c()
    # for(i in 1:length(sim_flow)){
    #   sim_10min<-c(sim_10min,rep(sim_flow[i],24*6))
    # }
    
    
    sim_flow_trans<-bc_transform_data(sim_flow,data$rating_error_uncertainty$lambda)
    # sim_flow_trans<-log_bc_transform(model_run,data$rating_error_uncertainty$bc_lambdas)
    if(all(!is.na(sim_flow_trans))){
      # agregate obs flow - assume that instantaneous flow is equivalent to 10 min volume!
      sim_flow_per_day<-sim_flow_trans*24*6
      indices<-sort(rep(1:length(sim_flow),24*6))
      obs_trans_agg<-aggregate(data$rating_flow_trans,by=list(indices),sum)[,2]
      
      residual_trans<-sim_flow_per_day-obs_trans_agg
      
      # assume that variances are additive:
      population_flow_variance_trans_agg<-population_flow_variance_trans*24*6
      
      # sd_zero_mean(residual_trans)
      # error_discharge<-model_run-data$obs_discharge
      #var_zero_mean(error_discharge)
      error_discharge_variance_calc<-sum(residual_trans^2)/(length(residual_trans))

      log_likelihood_error_discharge<-sum(dnorm(residual_trans,0,sqrt(population_flow_variance_trans_agg),log=T))
      
      # input error mean likelihood assuming true mean is zero
      combined_log_likehood<-log_likelihood_error_discharge+sum(log_likelihoods_input)
      
      return(list(combined_log_likehood=combined_log_likehood,error_discharge=residual_trans,error_discharge_variance_calc=error_discharge_variance_calc))
      
    } else {
      return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
    }
    
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
  }
  
}

plot_function<-function(){
  
  if(export_diagnostic){
    prefix<-paste("bates",start_ts,length_ts,initial_seed,seed,sep="_")
    # write parameters
    output_name<-paste0(output_dir,"/diagnostic_",prefix,".png")
    png(output_name,height=1000)
  }
  layout(matrix(1:10,nrow=5,byrow=T))
  # state_normaliser<-10^(theta[1:i,length_ts+data$number_of_sampled_inputs+2]*10)
  state_normaliser<-rep(data$normalisers[2],i)
  #all_state_normaliser<-cbind(state_normaliser,state_normaliser,state_normaliser,state_normaliser,state_normaliser,state_normaliser)
  all_state_normaliser<-matrix(rep(state_normaliser,length_ts-1),ncol=length_ts-1)
  all_normaliser<-cbind(rep(data$normalisers[1],length(state_normaliser)),
                        all_state_normaliser,
                        matrix(rep(data$normalisers[(length_ts+1):(length_ts+number_of_sampled_inputs)],length(state_normaliser)),
                               nrow=length(state_normaliser),ncol=number_of_sampled_inputs,byrow=T),
                        rep(data$normalisers[length_ts+number_of_sampled_inputs+1],length(state_normaliser)))
  
  # theta_unnorm<-cbind(inv.normalise(theta[1:i,1:(length_ts+number_of_sampled_inputs+1)],all_normaliser),state_normaliser)
  
  theta_unnorm<-theta[1:i,1:(length_ts+number_of_sampled_inputs+1)]
  
  #theta_unnorm<-t(apply(theta[!is.na(theta[,1]),],1,inv.normalise,factor=normalisers))
  #t(theta_unnorm)[1,]
  #inv.normalise(theta[1,],normalisers)
  plot(theta_unnorm[,1],type="l")
  plot(theta_unnorm[,2],type="l")
  plot(theta_unnorm[,length_ts+1],type="l")
  plot(theta_unnorm[,length_ts+number_of_sampled_inputs+1],type="l")
  # plot(theta_unnorm[,length_ts*2+1],type="l")
  #plot(theta_unnorm[,15],type="l")
  #plot(theta[!is.na(theta[,31]),31],type="l")
  #plot(theta[!is.na(theta[,61]),61],type="l")
  # plot(state_normaliser,type="l",log="y",main=state_normaliser[length(state_normaliser)])
  logposterior<-L+Pr
  if(!all(is.infinite(logposterior[!is.na(logposterior)]))){
    #plot(logposterior[max(1,i-20000):i],type="l")
    plot(logposterior[1:i],type="l")
    state_sd_calc<-apply(theta_unnorm[1:i,2:length_ts],1,sd_zero_mean)
    plot(x=state_sd_calc,y=logposterior[1:i],log="x")
    points(x=state_sd_calc[(length(state_sd_calc)-1999):length(state_sd_calc)],y=logposterior[(length(state_sd_calc)-1999):length(state_sd_calc)],col=3)
    points(x=state_sd_calc[1],y=logposterior[1],col=2,cex=4)
    
    
  }
  
  # input error sd plot
  # input_error_sd_calc<-apply(theta_unnorm[1:i,(length_ts+1):(length_ts+number_of_sampled_inputs)],1,sd_zero_mean)
  # plot(input_error_sd_calc,type="l")
  
  # transform rainfall data, introduce the error, then untransform again to simulate
  input_error_added<-data$input
  input_trans<-rep(NA,data$number_of_sampled_inputs)
  input_trans_with_error<-rep(NA,data$number_of_sampled_inputs)
  input_untrans_with_error<-rep(NA,data$number_of_sampled_inputs)
  input_not_zero<-which(input_error_added!=0)
  for(ii in 1:data$number_of_sampled_inputs){
    # TODO: move trans below to data input to improve efficiency
    input_trans[ii]<-bc_transform_data(data$input[input_not_zero[ii]],
                                       c(data$ensemble_rainfall_stats$bcparam1[input_not_zero[ii]],data$ensemble_rainfall_stats$bcparam2[input_not_zero[ii]]))
    input_trans_with_error[ii]<-input_trans[ii]-theta_unnorm[i,(length_ts+1):(length_ts+number_of_sampled_inputs)][ii]
    
    input_untrans_with_error[ii]<-inverse_bc_transform_data(input_trans_with_error[ii],
                                                            c(data$ensemble_rainfall_stats$bcparam1[input_not_zero[ii]],
                                                              data$ensemble_rainfall_stats$bcparam2[input_not_zero[ii]]))
    
  }
  
  input_error_added[which(input_error_added!=0)]<-input_untrans_with_error
  matplot(cbind(data$input,input_error_added),type="l")
  
  model_run<-gr4j.run(param=data$all_model_params, initial_state_S=theta_unnorm[i,1],
                      initial_state_R=theta_unnorm[i,length_ts+number_of_sampled_inputs+1], 
                      state_error=rep(0,length_ts-1),state_S_ts = theta_unnorm[i,2:length_ts], input=data.frame(P=input_error_added,E=E_input),
                      state_error_R=rep(0,length_ts-1))/1000*data$area_m2/86400
  
  if(!is.na(model_run[1])){
    # calculate for each day
    sim_flow<-model_run*86400
    indices<-sort(rep(1:length(model_run),24*6))
    obs_agg<-aggregate(actual_obs_discharge_untrans,by=list(indices),sum)[,2]
    
    matplot(cbind(obs_agg,sim_flow),type="l",main=paste("disch trans sd =",
                                                        round(sqrt(logL$error_discharge_variance_calc),2),
                                                        "( expected = ",round(sqrt((population_flow_sd_trans^2)*24*6),2),")"))
  }
  
  
  
  # if(i%%cor_plot_interval==0){
  #   # correlations
  #   cor_theta<-cor(theta[1:i,])
  #   lp<-levelplot(cor_theta,ylim=c(nrow(cor_theta)+0.5,0.5),at=seq(-1,1,length.out=51))
  #   print(lp)
  #   index_all_cor_theta<-i/cor_plot_interval
  #   all_cor_theta[index_all_cor_theta,]<-c(cor_theta)
  #   if(index_all_cor_theta>1 & i%%(cor_plot_interval*2)==0){
  #     layout(1)
  #     matplot(all_cor_theta[1:index_all_cor_theta,],type="l")
  #   }
  # }
  
  
  acceptance_rate<-length(which(Jump[!is.na(Jump)]==1))/length(Jump[!is.na(Jump)])
  cat("total acceptance_rate =",acceptance_rate,"\n")
  
  cat("discharge trans sd =",sqrt(logL$error_discharge_variance_calc),"( expected = ",sqrt((population_flow_sd_trans^2)*24*6),") \n")
  if(export_diagnostic){
    dev.off()
  }
}

#min_par<-c(0,0,1e-6)
#max_par<-c(100,100,100)

#initial_params<-c(500,rep(0,length(actual_state_error)),rep(1e-6,length(input_error)))
#initial_params<-c(600,rep(0,length(actual_state_error)),rnorm(length(input_error),0,0.1))
#initial_params<-c(actual_initial_state,actual_state_error,input_error)
#initial_params<-c(normalise(c(actual_initial_state,actual_state_error,input_error),normalisers),log10(initial_state_normaliser)) #log10(0.01)
# initial_params<-c(normalise(
#   c(rnorm(length(actual_initial_state),0,1000),
#     rnorm(length(actual_state_error),0,sd_zero_mean(actual_state_error)),
#     rnorm(length(input_error),0,input_error_sd)),
#   normalisers),log10(initial_state_normaliser)) #log10(0.01)
#initial_params<-c(600,rnorm(length(actual_state_error),0,3),rnorm(length(input_error),0,0.1))
#params_min<-c(normalise(c(-1000,rep(-1000,length(actual_state_error)),rep(-1000,length(input_trial))),normalisers),log10(1e-12))
#params_max<-c(normalise(c(1000,rep(1000,length(actual_state_error)),rep(1000,length(input_trial))),normalisers),log10(1e12))
# params_min<-c(-1000,rep(-1,length(actual_state_error)),rep(-1000,length(input_trial)),log10(1e-12)/10,rep(-1000,length(actual_state_error)))
# params_max<-c(1000,rep(1,length(actual_state_error)),rep(1000,length(input_trial)),log10(1e12)/10,rep(1000,length(actual_state_error)))

# params_min<-c(-1000,rep(-1,length(actual_state_error)),rep(-0.3,length(input_trial)),log10(1e-12)/10,rep(-1000,length(actual_state_error)))
# params_max<-c(1000,rep(1,length(actual_state_error)),rep(0.3,length(input_trial)),log10(1e12)/10,rep(1000,length(actual_state_error)))
# TODO: double check these, update initial production state bounds and input error bounds
# Scaled values are within the sampling range nominally about -1 to 1.
# Below bounds are as follows:
# actual production store initial state (unscaled),
# production store state samples before unscaling,
# input error samples after scaling but before boxcox untransform,
# actual routing initial state (unscaled),
# scaled scaling factor for production store (before itself is unscaled),
# production store state samples after unscaling
# routing state samples after unscaling
params_min<-c(0,
              rep(0,length(input_trial)-1),
              -ensemble_rainfall_stats$sd[index_sampled_inputs]*5,
              0,
              log10(1e-12)/10,
              rep(-4000,length(input_trial)-1) # used for unnormalised state S bounds
              ) 
params_max<-c(actual_model_param,
              rep(actual_model_param,length(input_trial)-1),
              ensemble_rainfall_stats$sd[index_sampled_inputs]*5,
              all_model_params[3],
              log10(1e12)/10,
              rep(4000,length(input_trial)-1) # used for unnormalised state S bounds
              )

data<-list(input=input_trial,E=E_input,model_param=actual_model_param,
           obs_discharge=actual_obs_discharge_untrans,
           # error_discharge_variance=error_discharge_variance,
           population_flow_variance_trans=population_flow_variance_trans,error_input_variance=error_input_variance,
           error_input_variance_min=1e-20,error_input_variance_max=100,
           error_discharge_variance_min=1e-20,error_discharge_variance_max=100,
           state_sds_min=1e-60,state_sds_max=20,
           normalisers=normalisers,
           initial_state_R=actual_initial_state_R,
           all_model_params=all_model_params,
           rating_error_uncertainty=bc2,
           rating_flow_trans=rating_flow_trans,
           area_m2=area_m2,
           ensemble_rainfall_stats=ensemble_rainfall_stats,
           number_of_sampled_inputs=number_of_sampled_inputs)

# data<-list(input=input_trial,E=E_input,model_param=actual_model_param,
#            obs_discharge=true_flow-actual_discharge_error,
#            # error_discharge_variance=error_discharge_variance,
#            population_flow_variance_trans=population_flow_variance_trans,error_input_variance=error_input_variance,
#            error_input_variance_min=1e-20,error_input_variance_max=100,
#            error_discharge_variance_min=1e-20,error_discharge_variance_max=100,
#            state_sds_min=1e-60,state_sds_max=20,
#            normalisers=normalisers,
#            rating_error_uncertainty=rating_error_uncertainty,
#            rating_flow_trans=rating_flow_trans,
#            actual_obs_flow_trans=actual_obs_flow_trans)

# browser()
# check initialise params and test
# Initialisation ############################################

if(exists("initial_seed")) set.seed(initial_seed)

# init_counter<-0
# init_fail<-T
# while(init_fail){
#   init_counter<-init_counter+1
#   if((init_counter%%1000)==0) cat("init_counter",init_counter,"\n")
#   # initial_params<-c(normalise(
#   #   c(rnorm(length(actual_initial_state),0,1000),
#   #     rnorm(length(actual_state_error),0,sd_zero_mean(actual_state_error)),
#   #     rnorm(length(input_error),0,input_error_sd)),
#   #   normalisers),log10(initial_state_normaliser)) #log10(0.01)
#   # logprior_fun<-trial_log_prior4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_adjusted_v5.8
#   # loglikelihood_fun<-log_likelihood_trial4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error
#   
#   # initial_params<-c(normalise(c(actual_initial_state,
#   #                                  rep(0,length(input_trial)-1),
#   #                                  rep(0,length(input_trial)),
#   #                                  actual_initial_state_R),normalisers),log10(initial_state_normaliser)/10)
#   if(exists("state_errors_initial_file")){
#     initial_params<-c(normalise(c(actual_initial_state,
#                                   state_errors_initial[,1],
#                                   rep(0,number_of_sampled_inputs),
#                                   actual_initial_state_R),normalisers))
#     
#     
#   } else {
#     # initial_params<-c(normalise(c(actual_initial_state,
#     #                               rep(0,length(input_trial)-1),
#     #                               rep(0,number_of_sampled_inputs),
#     #                               actual_initial_state_R),normalisers),
#     #                   log10(initial_state_normaliser)/10)
#     input_params_min<-params_min[(length_ts+1):(length_ts+number_of_sampled_inputs)]
#     input_params_max<-params_max[(length_ts+1):(length_ts+number_of_sampled_inputs)]
#     all_init_input<-c()
#     for(iinput in 1:number_of_sampled_inputs){
#       all_init_input<-c(all_init_input,runif(1,input_params_min[iinput],input_params_max[iinput]))
#     }
#     
#     # initial_params<-c(normalise(c(runif(1,params_min[1],params_max[1]),
#     #                               runif(length(input_trial)-1,params_min[2],params_max[2]),
#     #                               all_init_input,
#     #                               runif(1,params_min[length_ts+number_of_sampled_inputs+1],params_max[length_ts+number_of_sampled_inputs+1])),normalisers)
#     #                   )
#     initial_params<-c(runif(1,params_min[1],params_max[1]),
#                                   runif(length(input_trial)-1,params_min[2],params_max[2]),
#                                   all_init_input,
#                                   runif(1,params_min[length_ts+number_of_sampled_inputs+1],params_max[length_ts+number_of_sampled_inputs+1]))
#   }
#   logprior_init<-logprior_fun(initial_params,params_min,params_max,data)
#   loglike_init<-loglikelihood_fun(data,initial_params)[[1]]
#   if(!is.infinite(logprior_init) & !is.infinite(loglike_init)) init_fail<-F
#   if(init_counter>1000000) stop("Initialisation failed!!")
# }
if(exists("theta_file")){
  initial_params<-c(rrm_opt$par[5],rep(rrm_opt$par[5],length(input_trial)-1),
                                          rep(0,number_of_sampled_inputs),
                                          rrm_opt$par[6])
} else {
  init_opt_fun<-function(par){
    logprior_init<-logprior_fun(par,params_min,params_max,data)
    loglike_init<-loglikelihood_fun(data,par)[[1]]
    return(-(logprior_init+loglike_init))
  }
  
  library(DEoptim)
  init_opt<-DEoptim(init_opt_fun,
                    lower = params_min[1:(length_ts+number_of_sampled_inputs+1)],
                    upper = params_max[1:(length_ts+number_of_sampled_inputs+1)],
                    control = list(trace=1))
  bestval<-round(init_opt$optim$bestval,2)
  prev_bestval<-bestval+1
  while(bestval!=prev_bestval){
    prev_bestval<-bestval
    init_opt<-DEoptim(init_opt_fun,
                      lower = params_min[1:(length_ts+number_of_sampled_inputs+1)],
                      upper = params_max[1:(length_ts+number_of_sampled_inputs+1)],
                      control = list(trace=1,itermax=2000,initialpop=init_opt$member$pop))
    bestval<-round(init_opt$optim$bestval,2)
  }
  
  initial_params<-init_opt$optim$bestmem
}


# # check initialise params and test
# init_counter<-0
# init_fail<-T
# while(init_fail){
#   init_counter<-init_counter+1
#   # initial_params<-c(normalise(
#   #   c(rnorm(length(actual_initial_state),0,1000),
#   #     rnorm(length(actual_state_error),0,sd_zero_mean(actual_state_error)),
#   #     rnorm(length(input_error),0,input_error_sd)),
#   #   normalisers),log10(initial_state_normaliser)/10) #log10(0.01)
#   initial_params<-c(normalise(c(actual_initial_state,actual_state_error,input_error),normalisers),log10(initial_state_normaliser)/10)
#   logprior_init<-logprior_fun(initial_params,params_min,params_max,data)
#   loglike_init<-loglikelihood_fun(data,initial_params)[[1]]
#   if(!is.infinite(logprior_init) & !is.infinite(loglike_init)) init_fail<-F
#   if(init_counter>1000) stop("Initialisation failed!!")
# }


num_timesteps_cov<-Inf #100000
cor_plot_interval<-20000
tune_scale<-T # tunes the scaling factor (SD2) according to the acceptance rate over the last tune_interval. From pymc3 (metropolis.py)
tune_interval<-2000 # this should be divisible by update_covariance_interval
update_covariance_interval<-100 #20 #100 # this should be divide into tune_interval
start_collecting_CovPar<-Inf #1000001
start_sampling_CovPar<-Inf #1000001
resample_thetas<-F
start_gibbs_iter<-1

# if(!exists("ITER")) ITER = 2000000  #number of iterations to perform
i0   = 0.001 #0.001  #percentage of initial iterations before adaptation
# rm("SD2")
if(exists("SD2")){
  SD1  = SD2  #initial covariance matrix scaling factor (i0)
} else {
  SD1  = 0.50  #initial covariance matrix scaling factor (i0)
  SD2  = (2.4^2)/length(initial_params) #from Haario #0.30 #0.009 #0.15  #adaptive covariance matrix scaling factor (1-i0) (lower scaling is higher acceptance) )
}
if(!exists("stop_update_covariance")) stop_update_covariance<-0.5

if(export_diagnostic){
  diagnostic_plot_interval<-ITER
} else {
  diagnostic_plot_interval<-2000
}

#ncores<-detectCores(logical=F)

### 2 - Simulation model parameters
#PAR1 = 0.5	#place holder
#PAR2 = 15.00	#place holder
### 3 - Likelihood function parameters
#VARP = 0.185	#variance parameter
### 4 - Define parameter matrix
INIT_PAR<-initial_params
if(exists("theta_file")){
  theta<-as.matrix(read.csv(theta_file,as.is=T))
  prev_theta<-theta[!is.na(theta[,1]),]
  if(is.null(nrow(prev_theta))){
    INIT_PAR<-prev_theta
    prev_theta<-t(prev_theta)
  } else if(nrow(prev_theta)==1) {
    INIT_PAR<-prev_theta[nrow(prev_theta),]
  } else {
    INIT_PAR<-prev_theta[nrow(prev_theta)-1,]
  }
}

theta<-matrix(NA,nrow=ITER,ncol=length(INIT_PAR))
if(exists("restart_number") & exists("prev_theta")){
  indices_to_replace<-(nrow(prev_theta)-restart_number):(nrow(prev_theta)-1)
  theta[1:length(indices_to_replace),]<-prev_theta[indices_to_replace,]
  start_iter<-length(indices_to_replace)+1
} else {
  start_iter<-2
  theta[1,] = INIT_PAR #first row using initial values
}

nPAR = ncol(theta) #determines number of parameters
### 5 - Define parameter names for output
#ParName = c('Slope','Intercept','Variance') #redefine as appropriate

## Define parameter priors -----------------------------------------------#
# Define (log) prior to be used for each calibrated parameter.
# Pr1 = dbeta(PAR1,0.5,1.0,log=T)
# Pr2 = dgamma(PAR2,0.5,scale=1.0,log=T)
# Pr3 = dexp(VARP,0.5,log=T)
# Pr<-rep(NA,ITER)
# Pr[1] = Pr1+Pr2+Pr3
Pr<-rep(NA,ITER)
Pr[start_iter-1] = logprior_fun(theta[start_iter-1,],params_min,params_max,data)
if(is.infinite(Pr[start_iter-1])) stop("Initial prior is infinite!!")

## Initialize covariance matrix ------------------------------------------#
### 1 - Define best guess variances for each calibrated parameter
VarPar = rep(1E-6,nPAR)
# A <- matrix(VarPar,ncol=nPAR) # allocates some memory for rmvn function

### 2 - Populate initial covariance matrix
if(!exists("init_cov_file")){
  CovPar<-matrix(0,nrow=length(VarPar),ncol=length(VarPar))
  for(i in 1:nPAR){
    CovPar[i,i] = SD1*VarPar[i]
  }
} else {
  CovPar<-as.matrix(read.csv(init_cov_file,as.is=T,header=F))
  if(!all(dim(CovPar)==c(length(VarPar),length(VarPar)))){
    stop("dimensions of initial covariance matrix is wrong")
    CovPar<-matrix(0,nrow=length(VarPar),ncol=length(VarPar))
  }
}

### 3 - Define other covariance terms needed later
epsilon = 1e-20 #1e-20 #small number preventing CovPar from becoming singular during adaptation
Id = diag(nPAR)  #identity matrix used in CovPar adaptation

recalculate_CovPar<-T
if(recalculate_CovPar & exists("prev_theta")){
  if(nrow(prev_theta)>1){
    CovPar_unscaled<-cov(theta[max(1,start_iter-num_timesteps_cov):(start_iter-1),])
    CovPar = SD2*CovPar_unscaled+SD2*epsilon*Id
  }
}


## Initial run of simulation model ---------------------------------------#
# Replace the following line with your function as appropriate, where
# y is the output, SimModel is the function name, PAR1 & PAR2 are
# simulation model parameters
# SimModel<-function(PAR1,PAR2){
#   1:10*PAR1+PAR2
# }
# y = SimModel(theta[1,],data)

# Replace the following line with the appropriate log likelihood
# function, where obs is the observed data, pred is the predicted data,
# and varp is the variance parameter.

# normal_loglikelihood<-function(obs,pred,varp){
#   sum(dnorm(obs-pred,mean=0,sd=sqrt(varp),log=T))
# }
# obs<-15:24
# varp<-VARP

logL = loglikelihood_fun(data,theta[start_iter-1,])
L<-rep(NA,ITER)

L[start_iter-1] = logL[[1]] #logL$combined_log_likehood #initial value for the log likelihood





## Run AM algorithm ####################################################
# theta_pro<-matrix(NA,nrow=ITER,ncol=nPAR)
L_pro<-rep(NA,ITER)
Pr_pro<-rep(NA,ITER)
psi<-rep(NA,ITER)
Jump<-rep(NA,ITER)
Jumps = 0 #jump counter
beginning_time<-Sys.time()
all_cor_theta<-matrix(NA,nrow=ITER/cor_plot_interval,ncol=nPAR*nPAR)
all_SD2<-rep(NA,ITER/tune_interval)
all_interval_acceptance<-rep(NA,ITER/tune_interval)
all_discharge_error_sd<-c()
first_covpar_calculated<-F
iii<-0
param_indices<-1:nPAR

for(i in start_iter:ITER){
  if(i%%1000==0) cat(i,"/",ITER,"|",start_ts,"|",length_ts,"|",initial_seed,"\n")
  # if(i%%1000==0) browser()
  iii<-iii+1
  if(iii>nPAR) iii<-1
  
  ### 1 - Adapt covariance matrix
  if (i > i0*ITER & i%%update_covariance_interval==0 & i<stop_update_covariance*ITER){  #adaptive covariance matrix routine
    thin_for_covar<-F
    if(thin_for_covar){
      theta_tmp<-theta[max(1,i-num_timesteps_cov):(i-1),]
      theta_tmp<-theta_tmp[seq(1,nrow(theta_tmp),by=100),]
      CovPar_unscaled<-cov(theta_tmp)
    } else {
      use_updateCovariance<-T
      if(use_updateCovariance & exists("CovPar_unscaled") & first_covpar_calculated & is.infinite(num_timesteps_cov)){
        library(onlinePCA)
        CovPar_unscaled_new<-updateCovariance(CovPar_unscaled,
                                              theta[max(1,i-update_covariance_interval):(i-1),],
                                              i-update_covariance_interval-1,
                                              colMeans(theta[1:max(1,i-update_covariance_interval-1),]))
        
        if(i%%(update_covariance_interval*10000)==0){
          actual_CovPar_unscaled<-cov(theta[1:(i-1),])
          if(!all.equal(CovPar_unscaled_new,actual_CovPar_unscaled)){
            cat("covariance update problem \n")
            browser()
            stop("covariance update problem")
            CovPar_unscaled_new<-actual_CovPar_unscaled
          }
        }
        CovPar_unscaled<-CovPar_unscaled_new
        
      } else {
        CovPar_unscaled<-cov(theta[max(1,i-num_timesteps_cov):(i-1),])
        first_covpar_calculated<-T
      }
      
      
    }
    
    CovPar = SD2*CovPar_unscaled+SD2*epsilon*Id #updates the covariance
    if(i>=start_collecting_CovPar){
      all_CovPar[[length(all_CovPar)+1]]<-CovPar
      all_CovPar_means<-rbind(all_CovPar_means,theta[i-1,])
    }
    if(i>=start_sampling_CovPar & i%%(update_covariance_interval*100)==0){
      sample_CovPar<-sample(length(all_CovPar),1)
      CovPar<-all_CovPar[[sample_CovPar]]
      if(resample_thetas){
        theta[i-1,]<-all_CovPar_means[sample_CovPar,]
      }
      
    }
  }
  # Delayed Rejection ###########
  number_delayed_rejections<-4
  attempt<-0
  DR_scaling_factor<-1
  while(attempt<=number_delayed_rejections){
    attempt<-attempt+1
    ### 2 - Generate parameter proposals
    if(i>=start_gibbs_iter){
      
      theta[i,]<-theta[i-1,]
      
      # check positive definiteness
      # if(!is.positive.definite(CovPar)){
      #   stop("CovPar not positive definite")
      # }
      try_rcmvnorm<-try({
        theta[i,iii]<-rcmvnorm(n=1,
                               mean=theta[i-1,],
                               sigma=CovPar*DR_scaling_factor,
                               dependent.ind = iii,
                               given.ind = param_indices[-iii],
                               X.given = theta[i-1,-iii],
                               method="chol",check.sigma = F)
      },silent = T)
      if(length(grep("Error",try_rcmvnorm))>0){
        theta[i,]<-mvrnorm(1,theta[i-1,],CovPar*DR_scaling_factor) #proposed values - more robust faster
      }

    } else {
      theta[i,]<-mvrnorm(1,theta[i-1,],CovPar*DR_scaling_factor) #proposed values - more robust faster
    }
    #theta[i,] = rmvnorm(1,theta[i-1,],CovPar) #proposed values - less robust slower
    
    DR_scaling_factor<-DR_scaling_factor*0.5 # for next delayed rejection stage
    
    # # adjust covpar
    # input_error_indices<-(length_ts+1):(length_ts*2)
    # for(ii in input_error_indices){
    #   # CovPar[ii,ii] = SD2*(input_error_sd^2)+SD2*epsilon
    #   CovPar[ii,ii] = (input_error_sd^2)
    # }
    # prop_mean<-theta[i-1,]
    # prop_mean[input_error_indices]<-0
    # theta[i,]<-mvrnorm(1,prop_mean,CovPar) #proposed values - more robust faster
    # # plot(mvrnorm(1,theta[i-1,],CovPar))
    
    
    #system.time({for(jj in 1:1000000) theta[i,]<-mvrnorm(1,theta[i-1,],CovPar)})
    #system.time({for(jj in 1:1000000) theta[i,]<-rmvn(1,theta[i-1,],CovPar,ncores=ncores)})
    #system.time({for(jj in 1:1000000) rmvn(1,theta[i-1,],CovPar,ncores=ncores,A=A); theta[i,]<-A})
    #rmvn(1,theta[i-1,],CovPar,ncores=ncores,isChol=T,A=A); theta[i,]<-A #proposed values - faster
    #rmvn(1,theta[i-1,],CovPar,ncores=ncores,isChol=F,A=A); theta[i,]<-A #proposed values - faster
    
    # # replace the input error proposal with sample from actual distribution
    # input_error_indices<-(length_ts+1):(length_ts*2)
    # input_error_sample<-rnorm(length(input_error_indices),mean=0,sd=input_error_sd)
    # # sd_zero_mean(input_error_sample)
    # # TODO: scale this
    # theta[i,input_error_indices]<-input_error_sample
    
    # Ensure within range
    #theta[i,theta[i,]<params_min]<-params_min[theta[i,]<params_min]
    #theta[i,theta[i,]>params_max]<-params_max[theta[i,]>params_max]
    #theta_restrict<-theta[i,]
    #theta_restrict[theta[i,]<params_min]<-params_min[theta[i,]<params_min]
    #theta_restrict[theta[i,]>params_max]<-params_max[theta[i,]>params_max]
    
    # theta_pro[i,] <- theta[i,] #matrix of proposed theta's (rejected & accepted)
    
    ### 3 - Run simulation model
    #y = SimModel(theta[i,],data)
    
    ### 5 - Compute log prior
    #   Pr1 = dbeta(theta[i,1],0.5,1.0,log=T)
    #   Pr2 = dgamma(theta[i,2],0.5,scale=1.0,log=T)
    #   Pr3 = dexp(theta[i,3],0.5,log=T)
    #   Pr[i] = Pr1+Pr2+Pr3
    Pr[i] = logprior_fun(theta[i,],params_min,params_max,data)
    #Pr[i] = logprior_fun(theta_restrict,params_min,params_max,data)
    Pr_pro[i] = Pr[i] #proposed prior
    
    ### 4 - Compute log likelihood
    if(is.infinite(Pr[i])){
      logL = -Inf
    } else {
      logL = loglikelihood_fun(data,theta[i,])
    }
    # logL = loglikelihood_fun(data,theta[i,])
    #logL = loglikelihood_fun(data,theta_restrict)
    
    L[i] = logL[[1]] #logL$combined_log_likehood
    L_pro[i] = L[i] #proposed likelihood
    
    ### 6 - Compute Metropolis ratio
    psi[i] = exp((L[i]-L[i-1])+(Pr[i]-Pr[i-1])) #Metropolis ratio
    if(is.nan(psi[i])) psi[i]<-0
    ### 7 - Determine accept/reject of proposal
    z = runif(1,0,1)
    if(z <= psi[i]){
      #theta[i,] = theta[i,] #jump to next theta
      Jumps = Jumps+1
      Jump[i] = 1
      break
    } else{
      theta[i,] = theta[i-1,] #remain at current theta
      L[i] = L[i-1]
      Pr[i] = Pr[i-1]
      Jump[i] = 0
    }
  }
  
  ## 8 - Iterate
  if(i%%diagnostic_plot_interval==0){
    plot_function()
  }
  
  if(tune_scale & (i%%tune_interval)==0 & i > i0*ITER){
    Jump_interval<-Jump[(i-tune_interval+1):i]
    interval_acc_rate<-length(which(Jump_interval==1))/length(Jump_interval)
    
    SD2<-tune(SD2,interval_acc_rate)
    cat("interval acceptance_rate =",interval_acc_rate,"SD2 =",SD2,"\n")
    all_SD2[ITER/tune_interval]<-SD2
    all_interval_acceptance[ITER/tune_interval]<-interval_acc_rate
    CovPar = SD2*CovPar_unscaled+SD2*epsilon*Id
  }
  
  
  
}





elapsed_time<-difftime(Sys.time(),beginning_time,units="hours")
cat("Elapsed time: ",elapsed_time," ",units(elapsed_time),"\n")

# write the output
prefix<-paste("bates",start_ts,length_ts,initial_seed,seed,sep="_")

if(Sys.info()[1]=="Windows"){
  gzip<-"software/gzip/gzip.exe -f"
} else {
  gzip<-"gzip -f"
}

# write parameters
if(!exists("output_dir")){
  output_dir<-"output/state_uncertainty/AM"
}
output_name<-paste0(output_dir,"/theta_",prefix,".csv")
write.csv(theta,output_name,row.names=F,quote=F)
system(paste(gzip,output_name))

# write SD2
output_name<-paste0(output_dir,"/SD2_",prefix,".csv")
write.table(SD2,output_name,row.names=F,quote=F,col.names=F,sep=",")
system(paste(gzip,output_name))


# output_name<-paste0(output_dir,"/scaling_S_",prefix,".png")
# png(output_name)
# layout(1)
# plot(theta[1:i,2*length_ts+number_of_sampled_inputs+1],type="l")
# dev.off()

# output_name<-paste0(output_dir,"/scaling_R_",prefix,".png")
# png(output_name)
# layout(1)
# plot(theta[1:i,2*length_ts+number_of_sampled_inputs+2],type="l")
# dev.off()


# write Covariance
output_name<-paste0(output_dir,"/CovPar_",prefix,".csv")
write.table(CovPar,output_name,row.names=F,quote=F,col.names=F,sep=",")
system(paste(gzip,output_name))

# write rainfall-runoff model parameters
output_name<-paste0(output_dir,"/RRparameters_",prefix,".csv")
write.csv(t(rrm_opt$par),output_name,row.names=F,quote=F)
system(paste(gzip,output_name))

# if(Sys.info()[1]=="Linux"){
#   subject<-prefix
#   message<-prefix
#   email_command<-sprintf("echo \"%s\" | mail -s \"%s\" shaunsanghokim@gmail.com",message,subject)
#   system(email_command)
# }

# # write log likelihood
# output_name<-paste0(output_dir,"/loglikelihood_",prefix,".csv")
# write.table(L,output_name,row.names=F,quote=F,col.names=F,sep=",")
# system(paste(gzip,output_name))

# # write log prior
# output_name<-paste0(output_dir,"/logprior_",prefix,".csv")
# write.table(Pr,output_name,row.names=F,quote=F,col.names=F,sep=",")
# system(paste(gzip,output_name))

#write log posterior
output_name<-paste0(output_dir,"/logposterior_",prefix,".csv")
write.table(L+Pr,output_name,row.names=F,quote=F,col.names=F,sep=",")
system(paste(gzip,output_name))

# # write jump
# output_name<-paste0(output_dir,"/Jump_",prefix,".csv")
# write.table(Jump,output_name,row.names=F,quote=F,col.names=F,sep=",")
# system(paste(gzip,output_name))


# # write correlation timeseries
# output_name<-paste0(output_dir,"/cor_timeseries_",prefix,".csv")
# write.table(all_cor_theta,output_name,row.names=F,quote=F,col.names=F,sep=",")
# system(paste(gzip,output_name))


# fraction of accepted jumps
acceptance_rate<-length(which(Jump[!is.na(Jump)]==1))/length(Jump[!is.na(Jump)])
cat("acceptance_rate =",acceptance_rate,"\n")
cat("Finished!! \n")

# 
# # stop()
# # write.table(CovPar,"scripts/CovPar.csv",row.names=F,quote=F,col.names=F,sep=",")
# # 
# # write.csv(theta,"output/state_uncertainty/AM/theta.csv",row.names=F,quote=F)
# # 
# # # fraction of accepted jumps
# # acceptance_rate<-length(which(Jump[!is.na(Jump)]==1))/length(Jump[!is.na(Jump)])
# # cat("acceptance_rate =",acceptance_rate,"\n")
# # 
# # 
# # correlations
# library(lattice)
# cor_theta<-cor(theta)
# layout(1)
# lp<-levelplot(cor_theta,ylim=c(nrow(cor_theta)+0.5,0.5),at=seq(-1,1,length.out=51))
# print(lp)
# # 
# # # library(car)
# # # scatterplotMatrix(theta[1:100,1:3])
# # # pairs(theta[,1:3])
# # 
# # library(PerformanceAnalytics)
# # chart.Correlation(theta,method="pearson",histogram=T,pch=16)
# # 
# # # get rid of nas
# # theta_nona<-theta[!is.na(theta[,1]),]
# # 
# # indices of different data types
# groups_indices<-list(1,2:(length(actual_state_error)+1),
#                      (length(actual_state_error)+2):(length(actual_state_error)+length(input_trial)+1))
# # 
# # # plot intitial state, state error 1, input error 1
# # layout(1:3)
# # plot(theta_nona[,groups_indices[[1]]],type="l")
# # plot(theta_nona[,groups_indices[[2]][1]],type="l",main="2nd time step state error")
# # plot(theta_nona[,groups_indices[[3]][1]],type="l")
# # 
# # 
# # # Likelihood removed nas
# # layout(1:3)
# L_no_na<-L[!is.na(L)]
# # plot(L_no_na,type="l",main="log likelihood")
# # 
# # plot(Pr,type="l",main="log prior")
# # 
# # Calculate posterior
# Pr_no_na<-Pr[!is.na(Pr)]
# logposterior<-L_no_na+Pr_no_na
# # plot(logposterior,type="l",main="log posterior")
# # write.csv(logposterior,"output/state_uncertainty/AM/posterior.csv",row.names=F,quote=F)
# # 
# # unnormalise theta
# #theta_unnorm<-t(apply(theta[!is.na(theta[,1]),],1,inv.normalise,factor=normalisers))
# state_normaliser<-10^theta[1:i,length_ts+length_ts+1]
# all_state_normaliser<-matrix(rep(state_normaliser,length_ts-1),ncol=length_ts-1)
# all_normaliser<-cbind(rep(data$normalisers[1],length(state_normaliser)),
#                       all_state_normaliser,
#                       matrix(1,nrow=length(state_normaliser),ncol=length_ts))
# 
# theta_unnorm<-cbind(inv.normalise(theta[1:i,-ncol(theta)],all_normaliser),state_normaliser)
# 
# # layout(1:3)
# # calculate state error sd
# state_sd_calc<-apply(theta_unnorm[,groups_indices[[2]]],1,sd_zero_mean)
# plot(state_sd_calc,type="l",main="state error sd",log="y")
# 
# burn.in<-100000
# thin<-100
# theta_p<-theta_unnorm[seq(from=min(nrow(theta_unnorm),burn.in),to=nrow(theta_unnorm),by=thin),]
# #acf(theta_p[,15])
# state_sd_calc_p<-apply(theta_p[,groups_indices[[2]]],1,sd_zero_mean)
# layout(1:2)
# boxplot(state_sd_calc_p,main="state error sd burn-in thinned")
# boxplot(state_sd_calc_p,main="state error sd burn-in thinned (log y)",log="y")
# 
# 
# 
# # 
# # # calculate input error sd
# # input_sd_calc<-apply(theta_unnorm[,groups_indices[[3]]],1,sd_zero_mean)
# # plot(input_sd_calc,type="l",main="input error sd")
# # 
# # # calculate discharge error sd
# # discharge_sd_calc<-c()
# # for(ii in 1:nrow(theta_unnorm)){
# #   if(ii%%1000==0) cat(ii,"/",nrow(theta_unnorm),"\n")
# #   input_run<-data.frame(P=input_trial-as.numeric(theta_unnorm[ii,groups_indices[[3]]]),E=E_input)
# #   model_run<-gr4j.sma(data$model_param,as.numeric(theta_unnorm[ii,1]),as.numeric(theta_unnorm[ii,groups_indices[[2]]]),input_run)
# #   error_discharge<-model_run-data$obs_discharge
# #   discharge_sd_calc<-c(discharge_sd_calc,sd_zero_mean(error_discharge))
# # }
# # plot(discharge_sd_calc,type="l",main="discharge error sd")
# # 
# # layout(1)
# # boxplot(state_sd_calc,main="state sd calc",log="y")
# # abline(h=sd_zero_mean(actual_state_error),col=2,lty=2)
# # 
# # cat("actual_state_error_sd=",sd_zero_mean(actual_state_error),"\n")
# # cat("input_error_sd=",sqrt(mean(input_error^2)),"\n")
# # cat("actual_discharge_error_sd=",sqrt(mean(actual_discharge_error^2)),"\n")
# 
# # colormap
# library(RColorBrewer)
# #library(colorRamps)
# #col.ramp<-colorRampPalette(c("white","blue","red"))(100)
# col.ramp<-rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))
# layout(1)
# smoothScatter(x=state_sd_calc,y=logposterior[1:i],nrpoints=0,colramp=colorRampPalette(col.ramp))
# 
# 
# #theta_unnorm<-t(apply(prev_theta[!is.na(prev_theta[,1]),],1,inv.normalise,factor=normalisers))
# #ppp<-prev_theta
