
## MCMC AM Algorithm Shell: READ ME --------------------------------------#
#This script provides a basic shell of the Adaptive Metropolis
#algorithm (Haario et al, 2001). There are several locations where
#user-specific entries must be made (e.g., global variables, input
#data, parameter initialization, priors, simulation model function,
#likelihood function, etc.) to replace current place-holder code. Note
#that any parameter constraints (e.g., PAR1 > 0) must be built into the
#simulation model directly.

export_diagnostic<-F
diagnostic_plot_interval<-20000

if(length(commandArgs(TRUE))!=0){
  arguments <- commandArgs(TRUE)
  cat("arguments:",arguments,"\n")
  wd<-as.character(arguments[[1]])
  input_error_sd<-as.numeric(arguments[[2]])
  state_error_sd<-as.numeric(arguments[[3]])
  discharge_error_sd<-as.numeric(arguments[[4]])
  length_ts<-as.numeric(arguments[[5]])
  if(length(arguments)>5) initial_state_normaliser<-as.numeric(arguments[[6]])
  if(length(arguments)>6) output_dir<-as.character(arguments[[7]])
  if(length(arguments)>7) seed<-as.numeric(arguments[[8]])
  if(length(arguments)>8) ITER<-as.numeric(arguments[[9]])
  if(length(arguments)>9) stop_update_covariance<-as.numeric(arguments[[10]])
  
  if(length(arguments)>10) restart_number<-as.numeric(arguments[[11]])
  if(length(arguments)>11) init_cov_file<-as.character(arguments[[12]])
  if(length(arguments)>12) SD2_file<-as.character(arguments[[13]])
  if(length(arguments)>13) theta_file<-as.character(arguments[[14]])
  
  export_diagnostic<-T
  diagnostic_plot_interval<-1e6

}

if(!exists("wd")) wd<-"C:/Users/kim079/Documents/model_optimisation_framework"
if(!exists("input_error_sd")) input_error_sd<-0.1
if(!exists("state_error_sd")) state_error_sd<-1
if(!exists("discharge_error_sd")) discharge_error_sd<-0.01
if(!exists("length_ts")) length_ts<-60
# initial_state_normaliser<-0.01 
# init_cov_file<-"//pearceydata.csiro.au/data/kim079/model_optimisation_framework/output/state_uncertainty/AM/restarts/ts60_gr4jsmarouting_smallobserr/CovPar_0.1_1_0.01_60_1496644561.44581.csv.gz"
# SD2<-0.0160691998016972 
# theta_file<-"//pearceydata.csiro.au/data/kim079/model_optimisation_framework/output/state_uncertainty/AM/restarts/ts60_gr4jsmarouting_smallobserr/theta_0.1_1_0.01_60_1496644561.44581.csv.gz"
# seed<-1496902861.98095 
# output_dir<-"//pearceydata.csiro.au/data/kim079/model_optimisation_framework/output/state_uncertainty/AM/5.6_ts60_gr4jsmarouting_smallobserr"
# ITER<-1e+06 
# stop_update_covariance<-1 
# restart_number<-5e+05 


setwd(wd)
#source("scripts/gr4j_sma.r")
source("packages/gr4j_with_routing/R/gr4j_sma_routing.r")
if(!is.loaded("routing_gr4j_sk")){
  if(Sys.info()["sysname"]=="Linux"){
    dyn.load("packages/gr4j_with_routing/src/gr4j_with_routing.so")
  } else {
    dyn.load("packages/gr4j_with_routing/src/gr4j_with_routing.dll")
  }
}
source("scripts/state_uncertainty_trial_gr4j_sma_likelihood_prior_AM.r")
library(mvtnorm)
library(MASS)
library(lattice)
#library(mvnfast)
#library(parallel)
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
ts_file<-"data/MDB-Murray/combined/combined_401012_sitev5.001.csv.gz"
start_ts<-60

combined<-read.csv(ts_file,as.is=T,comment="#")
true_input_trial<-combined$rainfall.river[start_ts:(start_ts+length_ts-1)]
E_input<-combined$evap.river[start_ts:(start_ts+length_ts-1)]
#true_input_trial<-runif(length_ts,0,20) #runif(length_ts,100,112)
#E_input<-runif(length_ts,0,20) #runif(length_ts,100,105)
#input_error_sd<-1e-1 #1e-1
input_error<-rnorm(length(true_input_trial),mean=0,sd=input_error_sd) #sd=0.5
input_trial<-true_input_trial+input_error
input_trial<-pmax(input_trial,0)
actual_state_error<-rnorm(length(input_trial)-1,mean=0,sd=state_error_sd) #1e-6 # minimum should be max_state_error_bound/max_normaliser_value
#discharge_error_sd<-0.01 #0.01
actual_discharge_error<-rnorm(length(input_trial),mean=0,sd=discharge_error_sd) #sd=1e-6
true_input<-data.frame(P=true_input_trial,E=E_input)
actual_model_param<-1000
actual_initial_state<-500
actual_initial_state_R<-100
all_model_params<-c(actual_model_param,1,200,2)

initial_condition_UH1<-rep(0,ceiling(all_model_params[4])-1)
initial_condition_UH2<-rep(0,floor(all_model_params[4])+length(ceiling(all_model_params[4]):ceiling(all_model_params[4]*2))-1)


tt<-gr4j.run(param=all_model_params, initial_state_S=actual_initial_state, initial_state_R=actual_initial_state_R, 
              state_error=actual_state_error, input=true_input,initial_condition_UH1=initial_condition_UH1,initial_condition_UH2=initial_condition_UH2)

#tt<-gr4j.sma(actual_model_param,actual_initial_state,actual_state_error,true_input)

cat("actual_state_error_sd=",sd_zero_mean(actual_state_error),"\n")
cat("input_error_sd=",sqrt(mean(input_error^2)),"\n")
cat("actual_discharge_error_sd=",sqrt(mean(actual_discharge_error^2)),"\n")

# input_error
# normalise(actual_state_error,0.01)
# normalise(actual_initial_state,0.0001)
# normalisers<-c(0.0001,rep(0.01,length_ts-1),rep(1,length_ts))

if(exists("SD2_file")) SD2<-as.numeric(readLines(SD2_file))

if(!exists("seed")){
  seed<-as.numeric(Sys.time())
}
#seed<-1122
set.seed(seed)

#initial_state_normaliser<-10^(runif(1,1e-3,1))-1
#initial_state_normaliser<-50
# ff<-c()
# for(uu in 1:1000){
#   initial_state_normaliser<-10^(runif(1,1e-3,1.65))-1
#   ff<-c(ff,initial_state_normaliser)
# }
# hist(ff)
initial_state_normaliser<-1/max(abs(actual_state_error))
normalisers<-c(0.001,rep(initial_state_normaliser,length_ts-1),rep(1,length_ts),0.001)
#normalisers<-c(0.0001,rep(1e+29,length_ts-1),rep(1,length_ts))

# error_discharge_variance<-sum(actual_discharge_error^2)/(length(actual_discharge_error)-1)
# error_input_variance<-sum(input_error^2)/(length(input_error)-1)
error_discharge_variance<-discharge_error_sd^2
error_input_variance<-input_error_sd^2



logprior_fun<-trial_log_prior4_gr4jwithrouting_allinitstates4
loglikelihood_fun<-log_likelihood_trial4_gr4jwithrouting_allinitstates4

#min_par<-c(0,0,1e-6)
#max_par<-c(100,100,100)

#initial_params4.7<-c(500,rep(0,length(actual_state_error)),rep(1e-6,length(input_error)))
#initial_params4.7<-c(600,rep(0,length(actual_state_error)),rnorm(length(input_error),0,0.1))
#initial_params4.7<-c(actual_initial_state,actual_state_error,input_error)
#initial_params4.7<-c(normalise(c(actual_initial_state,actual_state_error,input_error),normalisers),log10(initial_state_normaliser)) #log10(0.01)
# initial_params4.7<-c(normalise(
#   c(rnorm(length(actual_initial_state),0,1000),
#     rnorm(length(actual_state_error),0,sd_zero_mean(actual_state_error)),
#     rnorm(length(input_error),0,input_error_sd)),
#   normalisers),log10(initial_state_normaliser)) #log10(0.01)
#initial_params4.7<-c(600,rnorm(length(actual_state_error),0,3),rnorm(length(input_error),0,0.1))
#params_min<-c(normalise(c(-1000,rep(-1000,length(actual_state_error)),rep(-1000,length(input_trial))),normalisers),log10(1e-12))
#params_max<-c(normalise(c(1000,rep(1000,length(actual_state_error)),rep(1000,length(input_trial))),normalisers),log10(1e12))
params_min<-c(-1000,rep(-1,length(actual_state_error)),rep(-1000,length(input_trial)),0,log10(1e-12)/10,rep(-1000,length(actual_state_error)))
params_max<-c(1000,rep(1,length(actual_state_error)),rep(1000,length(input_trial)),all_model_params[3],log10(1e12)/10,rep(1000,length(actual_state_error)))

data<-list(input=input_trial,E=E_input,model_param=actual_model_param,
           obs_discharge=tt-actual_discharge_error,
           error_discharge_variance=error_discharge_variance,error_input_variance=error_input_variance,
           error_input_variance_min=1e-20,error_input_variance_max=100,
           error_discharge_variance_min=1e-20,error_discharge_variance_max=100,
           state_sds_min=1e-60,state_sds_max=20,
           normalisers=normalisers,
           initial_state_R=actual_initial_state_R,
           all_model_params=all_model_params)

# check initialise params and test
init_counter<-0
init_fail<-T
while(init_fail){
  init_counter<-init_counter+1
  # initial_params4.7<-c(normalise(
  #   c(rnorm(length(actual_initial_state),0,1000),
  #     rnorm(length(actual_state_error),0,sd_zero_mean(actual_state_error)),
  #     rnorm(length(input_error),0,input_error_sd)),
  #   normalisers),log10(initial_state_normaliser)) #log10(0.01)
  initial_params4.7<-c(normalise(c(actual_initial_state,actual_state_error,input_error,actual_initial_state_R),normalisers),log10(initial_state_normaliser)/10)
  logprior_init<-logprior_fun(initial_params4.7,params_min,params_max,data)
  loglike_init<-loglikelihood_fun(data,initial_params4.7)[[1]]
  if(!is.infinite(logprior_init) & !is.infinite(loglike_init)) init_fail<-F
  if(init_counter>1000) stop("Initialisation failed!!")
}

cor_plot_interval<-20000
tune_scale<-T # tunes the scaling factor (SD2) according to the acceptance rate over the last tune_interval. From pymc3 (metropolis.py)
tune_interval<-2000
update_covariance_interval<-50
if(!exists("ITER")) ITER = 2000000  #number of iterations to perform
i0   = 0.001 #0.10  #percentage of initial iterations before adaptation
if(exists("SD2")){
  SD1  = SD2  #initial covariance matrix scaling factor (i0)
} else {
  SD1  = 0.50  #initial covariance matrix scaling factor (i0)
  SD2  = (2.4^2)/(length_ts*2+1) #from Haario #0.30 #0.009 #0.15  #adaptive covariance matrix scaling factor (1-i0) (lower scaling is higher acceptance) )
}
if(!exists("stop_update_covariance")) stop_update_covariance<-0.5
#ncores<-detectCores(logical=F)

### 2 - Simulation model parameters
#PAR1 = 0.5	#place holder
#PAR2 = 15.00	#place holder
### 3 - Likelihood function parameters
#VARP = 0.185	#variance parameter
### 4 - Define parameter matrix
INIT_PAR<-initial_params4.7
if(exists("theta_file")){
  theta<-as.matrix(read.csv(theta_file,as.is=T))
  prev_theta<-theta[!is.na(theta[,1]),]
  if(is.null(nrow(prev_theta))){
    INIT_PAR<-prev_theta
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
  CovPar_unscaled<-cov(prev_theta)
  CovPar = SD2*CovPar_unscaled+SD2*epsilon*Id
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


## Run AM algorithm ------------------------------------------------------#
theta_pro<-matrix(NA,nrow=ITER,ncol=nPAR)
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
prefix<-paste(input_error_sd,state_error_sd,discharge_error_sd,length_ts,seed,sep="_")
first_covpar_calculated<-F
for(i in start_iter:ITER){
  if(i%%1000==0) cat(i,"/",ITER,"|",prefix,"\n")
  ### 1 - Adapt covariance matrix
  if (i > i0*ITER & i%%update_covariance_interval==0 & i<stop_update_covariance*ITER){  #adaptive covariance matrix routine
    # CovPar_unscaled<-cov(theta[1:(i-1),])
    
    use_updateCovariance<-T
    if(use_updateCovariance & exists("CovPar_unscaled") & first_covpar_calculated){
      library(onlinePCA)
      CovPar_unscaled_new<-updateCovariance(CovPar_unscaled,
                                            theta[max(1,i-update_covariance_interval):(i-1),],
                                            i-update_covariance_interval-1,
                                            colMeans(theta[1:max(1,i-update_covariance_interval-1),]))
      
      if(i%%update_covariance_interval*10000==0){
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
      CovPar_unscaled<-cov(theta[1:(i-1),])
      first_covpar_calculated<-T
    }
    
    CovPar = SD2*CovPar_unscaled+SD2*epsilon*Id #updates the covariance
    
  }
  ### 2 - Generate parameter proposals
  #theta[i,] = rmvnorm(1,theta[i-1,],CovPar) #proposed values - less robust slower
  theta[i,]<-mvrnorm(1,theta[i-1,],CovPar) #proposed values - more robust faster
  #system.time({for(jj in 1:1000000) theta[i,]<-mvrnorm(1,theta[i-1,],CovPar)})
  #system.time({for(jj in 1:1000000) theta[i,]<-rmvn(1,theta[i-1,],CovPar,ncores=ncores)})
  #system.time({for(jj in 1:1000000) rmvn(1,theta[i-1,],CovPar,ncores=ncores,A=A); theta[i,]<-A})
  #rmvn(1,theta[i-1,],CovPar,ncores=ncores,isChol=T,A=A); theta[i,]<-A #proposed values - faster
  #rmvn(1,theta[i-1,],CovPar,ncores=ncores,isChol=F,A=A); theta[i,]<-A #proposed values - faster
  
  # Ensure within range
  #theta[i,theta[i,]<params_min]<-params_min[theta[i,]<params_min]
  #theta[i,theta[i,]>params_max]<-params_max[theta[i,]>params_max]
  #theta_restrict<-theta[i,]
  #theta_restrict[theta[i,]<params_min]<-params_min[theta[i,]<params_min]
  #theta_restrict[theta[i,]>params_max]<-params_max[theta[i,]>params_max]
  
  theta_pro[i,] = theta[i,] #matrix of proposed theta's (rejected & accepted)
  
  ### 3 - Run simulation model
  #y = SimModel(theta[i,],data)
  
  ### 4 - Compute log likelihood
  logL = loglikelihood_fun(data,theta[i,])
  #logL = loglikelihood_fun(data,theta_restrict)
  
  L[i] = logL[[1]] #logL$combined_log_likehood
  L_pro[i] = L[i] #proposed likelihood
  ### 5 - Compute log prior
  #   Pr1 = dbeta(theta[i,1],0.5,1.0,log=T)
  #   Pr2 = dgamma(theta[i,2],0.5,scale=1.0,log=T)
  #   Pr3 = dexp(theta[i,3],0.5,log=T)
  #   Pr[i] = Pr1+Pr2+Pr3
  Pr[i] = logprior_fun(theta[i,],params_min,params_max,data)
  #Pr[i] = logprior_fun(theta_restrict,params_min,params_max,data)
  Pr_pro[i] = Pr[i] #proposed prior
  ### 6 - Compute Metropolis ratio
  psi[i] = exp((L[i]-L[i-1])+(Pr[i]-Pr[i-1])) #Metropolis ratio
  if(is.nan(psi[i])) psi[i]<-0
  ### 7 - Determine accept/reject of proposal
  z = runif(1,0,1)
  if(z <= psi[i]){
    #theta[i,] = theta[i,] #jump to next theta
    Jumps = Jumps+1
    Jump[i] = 1
  } else{
    theta[i,] = theta[i-1,] #remain at current theta
    L[i] = L[i-1]
    Pr[i] = Pr[i-1]
    Jump[i] = 0
  }
  
  ## 8 - Iterate
  if(i%%2000==0){
    if(z > psi[i]){
      cur_discharge_error_sd<-sd_zero_mean(logL[[2]])
      all_discharge_error_sd<-c(all_discharge_error_sd,cur_discharge_error_sd)
    }
  }
  if(i%%diagnostic_plot_interval==0){
    if(export_diagnostic){
      # write parameters
      if(!exists("output_dir")){
        output_dir<-"output/state_uncertainty/AM"
      }
      output_name<-paste0(output_dir,"/diagnostic_",prefix,".png")
      png(output_name,height=1000)
    }
    cat("loglikelihood:",L[i],"\n")
    layout(matrix(1:10,nrow=5,byrow=T))
    state_normaliser<-10^(theta[1:i,length_ts+length_ts+2]*10)
    #all_state_normaliser<-cbind(state_normaliser,state_normaliser,state_normaliser,state_normaliser,state_normaliser,state_normaliser)
    all_state_normaliser<-matrix(rep(state_normaliser,length_ts-1),ncol=length_ts-1)
    all_normaliser<-cbind(rep(data$normalisers[1],length(state_normaliser)),
                          all_state_normaliser,
                          matrix(1,nrow=length(state_normaliser),ncol=length_ts),
                          rep(data$normalisers[length_ts*2+1],length(state_normaliser)))

    theta_unnorm<-cbind(inv.normalise(theta[1:i,-ncol(theta)],all_normaliser),state_normaliser)

    #theta_unnorm<-t(apply(theta[!is.na(theta[,1]),],1,inv.normalise,factor=normalisers))
    #t(theta_unnorm)[1,]
    #inv.normalise(theta[1,],normalisers)
    plot(theta_unnorm[,1],type="l")
    plot(theta_unnorm[,2],type="l")
    plot(theta_unnorm[,length_ts+1],type="l")
    plot(theta_unnorm[,length_ts*2+1],type="l")
    #plot(theta_unnorm[,15],type="l")
    #plot(theta[!is.na(theta[,31]),31],type="l")
    #plot(theta[!is.na(theta[,61]),61],type="l")
    plot(state_normaliser,type="l",log="y",main=state_normaliser[length(state_normaliser)])
    logposterior<-L+Pr
    if(!all(is.infinite(logposterior[!is.na(logposterior)]))){
      #plot(logposterior[max(1,i-20000):i],type="l")
      plot(logposterior[1:i],type="l")
      state_sd_calc<-apply(theta_unnorm[1:i,2:length_ts],1,sd_zero_mean)
      plot(x=state_sd_calc,y=logposterior[1:i],log="x")
      points(x=state_sd_calc[(length(state_sd_calc)-1999):length(state_sd_calc)],y=logposterior[(length(state_sd_calc)-1999):length(state_sd_calc)],col=3)
      points(x=state_sd_calc[1],y=logposterior[1],col=2,cex=4)
    }

    # if(i%%10000==0){
    boxplot(all_discharge_error_sd,main=paste("disch err sd:",cur_discharge_error_sd,"\nexpected:",discharge_error_sd))
    abline(h=discharge_error_sd,col=2,lty=2)
    boxplot(apply(theta_unnorm[1:i,(length_ts+1):(length_ts*2)],1,sd_zero_mean),main=paste("input err:",sd_zero_mean(theta_unnorm[i,(length_ts+1):(length_ts*2)]),"\nexpected:",input_error_sd))
    abline(h=input_error_sd,col=2,lty=2)
    #   # correlations
    #   library(lattice)
    #   cor_theta<-cor(theta[1:i,])
    #   lp<-levelplot(cor_theta,ylim=c(nrow(cor_theta)+0.5,0.5),at=seq(-1,1,length.out=51))
    #   print(lp)
    # }
    if(export_diagnostic){
      dev.off()
    }


    acceptance_rate<-length(which(Jump[!is.na(Jump)]==1))/length(Jump[!is.na(Jump)])
    cat("total acceptance_rate =",acceptance_rate,"\n")
  }
  # if(i%%cor_plot_interval==0){
  #   # correlations
  #   cor_theta<-cor(theta[1:i,])
  #   index_all_cor_theta<-i/cor_plot_interval
  #   all_cor_theta[index_all_cor_theta,]<-c(cor_theta)
  # }
  
  if(tune_scale & i%%tune_interval==0 & exists("CovPar_unscaled")){
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

output_name<-paste0(output_dir,"/scaling_",prefix,".png")
png(output_name)
layout(1)
plot(theta[1:i,length_ts+length_ts+2],type="l")
dev.off()

# write Covariance
output_name<-paste0(output_dir,"/CovPar_",prefix,".csv")
write.table(CovPar,output_name,row.names=F,quote=F,col.names=F,sep=",")
system(paste(gzip,output_name))

# write log likelihood
output_name<-paste0(output_dir,"/loglikelihood_",prefix,".csv")
write.table(L,output_name,row.names=F,quote=F,col.names=F,sep=",")
system(paste(gzip,output_name))

# write log prior
output_name<-paste0(output_dir,"/logprior_",prefix,".csv")
write.table(Pr,output_name,row.names=F,quote=F,col.names=F,sep=",")
system(paste(gzip,output_name))

#write log posterior
output_name<-paste0(output_dir,"/logposterior_",prefix,".csv")
write.table(L+Pr,output_name,row.names=F,quote=F,col.names=F,sep=",")
system(paste(gzip,output_name))

# write jump
output_name<-paste0(output_dir,"/Jump_",prefix,".csv")
write.table(Jump,output_name,row.names=F,quote=F,col.names=F,sep=",")
system(paste(gzip,output_name))


# # write correlation timeseries
# output_name<-paste0(output_dir,"/cor_timeseries_",prefix,".csv")
# write.table(all_cor_theta,output_name,row.names=F,quote=F,col.names=F,sep=",")
# system(paste(gzip,output_name))


# fraction of accepted jumps
acceptance_rate<-length(which(Jump[!is.na(Jump)]==1))/length(Jump[!is.na(Jump)])
cat("acceptance_rate =",acceptance_rate,"\n")
cat("Finished!! \n")

if(Sys.info()[1]=="Linux"){
  subject<-prefix
  message<-prefix
  email_command<-sprintf("echo \"%s\" | mail -s \"%s\" shaunsanghokim@gmail.com",message,subject)
  system(email_command)
}

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
