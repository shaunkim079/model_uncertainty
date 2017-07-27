# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial<-function(data,params){
  library(Brobdingnag)
  browser()
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  
  length_ts<-length(data$input)
  params_unnorm<-inv.normalise(params,data$normalisers)
  
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  var_zero_mean(input_error)
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_error,input=input)
  
  error_discharge<-model_run-data$obs_discharge
  var_zero_mean(error_discharge)
  
  error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
  log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
  #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
  if(is.infinite(log_likelihood_error_discharge_variance)){
    browser()
  }
  
  error_input_variance_calc<-sum(input_error^2)/(length(input_error)-1)
  #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
  #if(error_input_variance_calc<error_input_variance) browser()
  log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
  #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
  # for testing
  #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
  if(is.infinite(log_likelihood_error_input_variance)){
    browser()
  }
  
  combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance
  #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
  
  if(is.infinite(combined_log_likehood)) browser()
  
  return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  
  
}


trial_log_prior<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  params_unnorm<-inv.normalise(params,data$normalisers)
  
  initial_state<-params_unnorm[1]
  initial_state_min<-min[1]
  initial_state_max<-max[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-min[2:length_ts]
  state_errors_max<-max[2:length_ts]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #   log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-min[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-max[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_state_errors)+sum(log_prior_input_errors)+log_prior_initial_state #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}


# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial2<-function(data,params){
  library(Brobdingnag)
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  
  length_ts<-length(data$input)
  params_unnorm<-inv.normalise(params,data$normalisers)
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  var_zero_mean(input_error)
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_error,input=input)
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    if(is.infinite(log_likelihood_error_discharge_variance)){
      browser()
    }
    
    error_input_variance_calc<-sum(input_error^2)/(length(input_error)-1)
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    if(is.infinite(log_likelihood_error_input_variance)){
      browser()
    }
    
    state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance #+log_likelihood_state_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior2<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  params_unnorm<-inv.normalise(params,data$normalisers)
  params_min_unnorm<-inv.normalise(min,data$normalisers)
  params_max_unnorm<-inv.normalise(max,data$normalisers)
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[2:length_ts]
  state_errors_max<-params_max_unnorm[2:length_ts]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  log_prior_state_sds<-log(1/state_errors_sd_calc)
  if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_state_errors)+sum(log_prior_input_errors)+log_prior_initial_state #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}


# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial2_state_error_only<-function(data,params){
  library(Brobdingnag)
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  params_unnorm<-inv.normalise(params,data$normalisers)
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-data$input_error #params_unnorm[(length_ts+1):(length_ts+length_ts)]
  var_zero_mean(input_error)
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_error,input=input)
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)
    #log_likelihood_error_discharge_variance<-dnorm(error_discharge_variance_calc,mean=error_discharge_variance,sd=1e-6,log=T)
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    if(is.infinite(log_likelihood_error_discharge_variance)){
      browser()
    }
    
    error_input_variance_calc<-sum(input_error^2)/(length(input_error)-1)
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    if(is.infinite(log_likelihood_error_input_variance)){
      browser()
    }
    
    state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    combined_log_likehood<-log_likelihood_error_discharge_variance #+log_likelihood_state_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior2_state_error_only<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  params_unnorm<-inv.normalise(params,data$normalisers)
  params_min_unnorm<-inv.normalise(min,data$normalisers)
  params_max_unnorm<-inv.normalise(max,data$normalisers)
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[2:length_ts]
  state_errors_max<-params_max_unnorm[2:length_ts]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  log_prior_state_sds<-log(1/state_errors_sd_calc)
  if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_state_errors)+log_prior_initial_state #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}

# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial3<-function(data,params){
  library(Brobdingnag)
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  
  length_ts<-length(data$input)
  params_unnorm<-inv.normalise(params,data$normalisers)
  initial_state<-params_unnorm[1]
  state_with_error<-params_unnorm[2:length_ts]
  input_minus_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_error<-data$input-input_minus_error
  var_zero_mean(input_error)
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-input_minus_error #data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.sma2(param=data$model_param,initial_state=initial_state,state_with_error=state_with_error,input=input)
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    if(is.infinite(log_likelihood_error_discharge_variance)){
      browser()
    }
    
    error_input_variance_calc<-sum(input_error^2)/(length(input_error)-1)
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    if(is.infinite(log_likelihood_error_input_variance)){
      browser()
    }
    
    #state_with_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    #log_likelihood_state_with_error_sd<-sum(dnorm(state_with_error,mean=0,sd=state_with_error_sd,log=T))
    #if(is.nan(log_likelihood_state_with_error_sd)) log_likelihood_state_with_error_sd<- -Inf
    
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance #+log_likelihood_state_with_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior3<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  params_unnorm<-inv.normalise(params,data$normalisers)
  params_min_unnorm<-inv.normalise(min,data$normalisers)
  params_max_unnorm<-inv.normalise(max,data$normalisers)
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_with_errors<-params_unnorm[2:length_ts]
  state_with_errors_min<-params_min_unnorm[2:length_ts]
  state_with_errors_max<-params_max_unnorm[2:length_ts]
  log_prior_state_with_errors<-dunif(state_with_errors,min=state_with_errors_min,max=state_with_errors_max,log=T)
  #log_prior_state_with_errors<-dnorm(state_with_errors,mean=0,sd=0.001,log=T)
  
  state_with_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_with_errors_sd_calc<-sqrt(mean(state_with_errors^2))
  #   state_with_errors_sd_calc<-max(state_with_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_with_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  log_prior_state_sds<-log(1/state_with_errors_sd_calc)
  if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_minus_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_minus_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_minus_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_minus_errors<-dunif(input_minus_errors,min=input_minus_errors_min,max=input_minus_errors_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma2(param=data$model_param,initial_state=initial_state,state_error=state_with_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_state_with_errors)+sum(log_prior_input_minus_errors)+log_prior_initial_state #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}



# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4<-function(data,params){
  #library(Brobdingnag,lib.loc="packages")
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    # TODO: brobdingnag
    #frac<-as.brob(variance)/as.brob(expected_variance)
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^params[length_ts+length_ts+1]
  params_unnorm<-inv.normalise(params[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  #var_zero_mean(input_error)
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_error,input=input)
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    #var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    error_input_variance_calc<-sum(input_error^2)/(length(input_error)-1)
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    #     if(is.infinite(log_likelihood_error_input_variance)){
    #       browser()
    #     }
    #     state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    #     log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    #     if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance #+log_likelihood_state_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior4<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^params[length_ts+length_ts+1]
  state_normaliser_min<-10^min[length_ts+length_ts+1]
  state_normaliser_max<-10^max[length_ts+length_ts+1]
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[2:length_ts]
  state_errors_max<-params_max_unnorm[2:length_ts]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on normalised state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_state_errors)+sum(log_prior_input_errors)+log_prior_initial_state+state_normaliser_prior+log_prior_state_errors_normalised #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}



# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_without_normaliser<-function(data,params){
  #library(Brobdingnag,lib.loc="packages")
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    # TODO: brobdingnag
    #frac<-as.brob(variance)/as.brob(expected_variance)
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  #   state_normaliser<-10^params[length_ts+length_ts+1]
  #   params_unnorm<-inv.normalise(params[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  initial_state<-params[1]
  state_error<-params[2:length_ts]
  input_error<-params[(length_ts+1):(length_ts+length_ts)]
  var_zero_mean(input_error)
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_error,input=input)
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    error_input_variance_calc<-sum(input_error^2)/(length(input_error)-1)
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    #     if(is.infinite(log_likelihood_error_input_variance)){
    #       browser()
    #     }
    #     state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    #     log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    #     if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance #+log_likelihood_state_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior4_without_normaliser<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  #   state_normaliser<-10^params[length_ts+length_ts+1]
  #   state_normaliser_min<-10^min[length_ts+length_ts+1]
  #   state_normaliser_max<-10^max[length_ts+length_ts+1]
  #   state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  #   
  #   params_unnorm<-inv.normalise(params[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  #   params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  #   params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params[1]
  initial_state_min<-min[1]
  initial_state_max<-max[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params[2:length_ts]
  state_errors_min<-min[2:length_ts]
  state_errors_max<-max[2:length_ts]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on normalised state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-min[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-max[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_state_errors)+sum(log_prior_input_errors)+log_prior_initial_state #+state_normaliser_prior+log_prior_state_errors_normalised #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}



# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithrouting<-function(data,params){
  #library(Brobdingnag,lib.loc="packages")
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    # TODO: brobdingnag
    #frac<-as.brob(variance)/as.brob(expected_variance)
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^params[length_ts+length_ts+1]
  params_unnorm<-inv.normalise(params[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  #var_zero_mean(input_error)
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.run(param=data$all_model_params, initial_state_S=initial_state, initial_state_R=data$initial_state_R, 
           state_error=state_error, input=input)
  
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    #var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    error_input_variance_calc<-sum(input_error^2)/(length(input_error)-1)
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    #     if(is.infinite(log_likelihood_error_input_variance)){
    #       browser()
    #     }
    #     state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    #     log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    #     if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance #+log_likelihood_state_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior4_gr4jwithrouting<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^params[length_ts+length_ts+1]
  state_normaliser_min<-10^min[length_ts+length_ts+1]
  state_normaliser_max<-10^max[length_ts+length_ts+1]
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[2:length_ts]
  state_errors_max<-params_max_unnorm[2:length_ts]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on normalised state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_state_errors)+sum(log_prior_input_errors)+log_prior_initial_state+state_normaliser_prior+log_prior_state_errors_normalised #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}


# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithrouting_allinitstates<-function(data,params){
  #library(Brobdingnag,lib.loc="packages")
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    # TODO: brobdingnag
    #frac<-as.brob(variance)/as.brob(expected_variance)
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^params[length_ts+length_ts+2]
  params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),data$normalisers[(length_ts+1):(length_ts*2)],data$normalisers[length_ts*2+1]))
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[length_ts*2+1]
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.run(param=data$all_model_params, initial_state_S=initial_state, initial_state_R=init_state_R, 
                      state_error=state_error, input=input)
  
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    #var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    error_input_variance_calc<-sum(input_error^2)/(length(input_error)-1)
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    #     if(is.infinite(log_likelihood_error_input_variance)){
    #       browser()
    #     }
    #     state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    #     log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    #     if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance #+log_likelihood_state_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior4_gr4jwithrouting_allinitstates<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^params[length_ts+length_ts+2]
  state_normaliser_min<-10^min[length_ts+length_ts+2]
  state_normaliser_max<-10^max[length_ts+length_ts+2]
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),data$normalisers[(length_ts+1):(length_ts*2)],data$normalisers[length_ts*2+1]))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[2:length_ts]
  state_errors_max<-params_max_unnorm[2:length_ts]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on normalised state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  # prior on initial R store
  initial_state_R<-params_unnorm[length_ts*2+1]
  initial_state_R_min<-params_min_unnorm[length_ts*2+1]
  initial_state_R_max<-params_max_unnorm[length_ts*2+1]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_state_errors)+sum(log_prior_input_errors)+log_prior_initial_state+state_normaliser_prior+log_prior_state_errors_normalised+log_prior_initial_state_R #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}


# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithrouting_allinitstates2<-function(data,params){
  #library(Brobdingnag,lib.loc="packages")
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    # TODO: brobdingnag
    #frac<-as.brob(variance)/as.brob(expected_variance)
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts+length_ts+2]*10)
  params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),data$normalisers[(length_ts+1):(length_ts*2)],data$normalisers[length_ts*2+1]))
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[length_ts*2+1]
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.run(param=data$all_model_params, initial_state_S=initial_state, initial_state_R=init_state_R, 
                      state_error=state_error, input=input)
  
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    #var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    error_input_variance_calc<-sum(input_error^2)/(length(input_error)-1)
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    #     if(is.infinite(log_likelihood_error_input_variance)){
    #       browser()
    #     }
    #     state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    #     log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    #     if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance #+log_likelihood_state_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior4_gr4jwithrouting_allinitstates2<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts+length_ts+2]*10)
  state_normaliser_min<-10^(min[length_ts+length_ts+2]*10)
  state_normaliser_max<-10^(max[length_ts+length_ts+2]*10)
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),data$normalisers[(length_ts+1):(length_ts*2)],data$normalisers[length_ts*2+1]))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[(length_ts*2+3):(length_ts*3+1)]
  state_errors_max<-params_max_unnorm[(length_ts*2+3):(length_ts*3+1)]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)

  # prior on normalised state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  # prior on initial R store
  initial_state_R<-params_unnorm[length_ts*2+1]
  initial_state_R_min<-params_min_unnorm[length_ts*2+1]
  initial_state_R_max<-params_max_unnorm[length_ts*2+1]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_input_errors)+log_prior_initial_state+state_normaliser_prior+log_prior_state_errors_normalised+log_prior_initial_state_R+sum(log_prior_state_errors) #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}

# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithrouting_allinitstates3<-function(data,params){
  #library(Brobdingnag,lib.loc="packages")
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    # TODO: brobdingnag
    #frac<-as.brob(variance)/as.brob(expected_variance)
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts+length_ts+2]*10)
  params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),data$normalisers[(length_ts+1):(length_ts*2)],data$normalisers[length_ts*2+1]))
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[length_ts*2+1]
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.run(param=data$all_model_params, initial_state_S=initial_state, initial_state_R=init_state_R, 
                      state_error=state_error, input=input)
  
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    #var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge))
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    error_input_variance_calc<-sum(input_error^2)/(length(input_error))
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    #     if(is.infinite(log_likelihood_error_input_variance)){
    #       browser()
    #     }
    #     state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    #     log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    #     if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    # discharge error mean likelihood assuming true mean is zero
    error_discharge_mean_calc<-mean(error_discharge)
    log_likelihood_error_discharge_mean_calc<-dnorm(error_discharge_mean_calc,mean=0,sd=sqrt(error_discharge_variance)/sqrt(length(error_discharge)),log=T)
    
    # input error mean likelihood assuming true mean is zero
    error_input_mean_calc<-mean(input_error)
    log_likelihood_error_input_mean_calc<-dnorm(error_input_mean_calc,mean=0,sd=sqrt(error_input_variance)/sqrt(length(input_error)),log=T)
    
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance+log_likelihood_error_discharge_mean_calc+log_likelihood_error_input_mean_calc #+log_likelihood_state_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior4_gr4jwithrouting_allinitstates3<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts+length_ts+2]*10)
  state_normaliser_min<-10^(min[length_ts+length_ts+2]*10)
  state_normaliser_max<-10^(max[length_ts+length_ts+2]*10)
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),data$normalisers[(length_ts+1):(length_ts*2)],data$normalisers[length_ts*2+1]))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[(length_ts*2+3):(length_ts*3+1)]
  state_errors_max<-params_max_unnorm[(length_ts*2+3):(length_ts*3+1)]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on normalised state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  # prior on initial R store
  initial_state_R<-params_unnorm[length_ts*2+1]
  initial_state_R_min<-params_min_unnorm[length_ts*2+1]
  initial_state_R_max<-params_max_unnorm[length_ts*2+1]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_input_errors)+log_prior_initial_state+state_normaliser_prior+log_prior_state_errors_normalised+log_prior_initial_state_R+sum(log_prior_state_errors) #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}


# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithoutrouting_allinitstates2<-function(data,params){
  #library(Brobdingnag,lib.loc="packages")
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    # TODO: brobdingnag
    #frac<-as.brob(variance)/as.brob(expected_variance)
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts+length_ts+1]*10)
  params_unnorm<-inv.normalise(params[1:(length_ts*2)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),data$normalisers[(length_ts+1):(length_ts*2)]))
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  #var_zero_mean(input_error)
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_error,input=input)
  
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    #var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    error_input_variance_calc<-sum(input_error^2)/(length(input_error)-1)
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    #     if(is.infinite(log_likelihood_error_input_variance)){
    #       browser()
    #     }
    #     state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    #     log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    #     if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance #+log_likelihood_state_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior4_gr4jwithoutrouting_allinitstates2<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts+length_ts+1]*10)
  state_normaliser_min<-10^(min[length_ts+length_ts+1]*10)
  state_normaliser_max<-10^(max[length_ts+length_ts+1]*10)
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts*2)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),data$normalisers[(length_ts+1):(length_ts*2)]))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[(length_ts*2+2):(length_ts*3)]
  state_errors_max<-params_max_unnorm[(length_ts*2+2):(length_ts*3)]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on normalised state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_input_errors)+log_prior_initial_state+state_normaliser_prior+log_prior_state_errors_normalised+sum(log_prior_state_errors) #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}


# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithoutrouting_allinitstates3<-function(data,params){
  #library(Brobdingnag,lib.loc="packages")
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    # TODO: brobdingnag
    #frac<-as.brob(variance)/as.brob(expected_variance)
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts+length_ts+1]*10)
  params_unnorm<-inv.normalise(params[1:(length_ts*2)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),data$normalisers[(length_ts+1):(length_ts*2)]))
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  #var_zero_mean(input_error)
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_error,input=input)
  
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    #var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge))
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    error_input_variance_calc<-sum(input_error^2)/(length(input_error))
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    #     if(is.infinite(log_likelihood_error_input_variance)){
    #       browser()
    #     }
    #     state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    #     log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    #     if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    # discharge error mean likelihood assuming true mean is zero
    error_discharge_mean_calc<-mean(error_discharge)
    log_likelihood_error_discharge_mean_calc<-dnorm(error_discharge_mean_calc,mean=0,sd=sqrt(error_discharge_variance)/sqrt(length(error_discharge)),log=T)
    
    # input error mean likelihood assuming true mean is zero
    error_input_mean_calc<-mean(input_error)
    log_likelihood_error_input_mean_calc<-dnorm(error_input_mean_calc,mean=0,sd=sqrt(error_input_variance)/sqrt(length(input_error)),log=T)
    
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance+log_likelihood_error_discharge_mean_calc+log_likelihood_error_input_mean_calc #+log_likelihood_state_error_sd
    
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior4_gr4jwithoutrouting_allinitstates3<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts+length_ts+1]*10)
  state_normaliser_min<-10^(min[length_ts+length_ts+1]*10)
  state_normaliser_max<-10^(max[length_ts+length_ts+1]*10)
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts*2)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),data$normalisers[(length_ts+1):(length_ts*2)]))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[(length_ts*2+2):(length_ts*3)]
  state_errors_max<-params_max_unnorm[(length_ts*2+2):(length_ts*3)]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on normalised state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_input_errors)+log_prior_initial_state+state_normaliser_prior+log_prior_state_errors_normalised+sum(log_prior_state_errors) #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}

# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithrouting_allinitstates_samplestatelogspace<-function(data,params){
  #library(Brobdingnag,lib.loc="packages")
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    # TODO: brobdingnag
    #frac<-as.brob(variance)/as.brob(expected_variance)
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],data$normalisers[2:length_ts],data$normalisers[(length_ts+1):(length_ts*2)],data$normalisers[length_ts*2+1]))
  params_unnorm[2:length_ts]<-sapply(params_unnorm[2:length_ts],untranform_state_error)

  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[length_ts*2+1]
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.run(param=data$all_model_params, initial_state_S=initial_state, initial_state_R=init_state_R, 
                      state_error=state_error, input=input)
  
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    #var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    error_input_variance_calc<-sum(input_error^2)/(length(input_error)-1)
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    #     if(is.infinite(log_likelihood_error_input_variance)){
    #       browser()
    #     }
    #     state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    #     log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    #     if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance #+log_likelihood_state_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior4_gr4jwithrouting_allinitstates_samplestatelogspace<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  
  params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],data$normalisers[2:length_ts],data$normalisers[(length_ts+1):(length_ts*2)],data$normalisers[length_ts*2+1]))
  params_unnorm[2:length_ts]<-sapply(params_unnorm[2:length_ts],untranform_state_error)
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[2:length_ts]
  state_errors_max<-params_max_unnorm[2:length_ts]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on normalised state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  # prior on initial R store
  initial_state_R<-params_unnorm[length_ts*2+1]
  initial_state_R_min<-params_min_unnorm[length_ts*2+1]
  initial_state_R_max<-params_max_unnorm[length_ts*2+1]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_state_errors)+sum(log_prior_input_errors)+log_prior_initial_state+log_prior_state_errors_normalised+log_prior_initial_state_R #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}



# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithrouting_allinitstates_samplealllogspace<-function(data,params){
  #library(Brobdingnag,lib.loc="packages")
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    # TODO: brobdingnag
    #frac<-as.brob(variance)/as.brob(expected_variance)
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  params_unnorm<-sapply(params,untranform_state_error)
  
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[length_ts*2+1]
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.run(param=data$all_model_params, initial_state_S=initial_state, initial_state_R=init_state_R, 
                      state_error=state_error, input=input)
  
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    #var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    error_input_variance_calc<-sum(input_error^2)/(length(input_error)-1)
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    #     if(is.infinite(log_likelihood_error_input_variance)){
    #       browser()
    #     }
    #     state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    #     log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    #     if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance #+log_likelihood_state_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior4_gr4jwithrouting_allinitstates_samplealllogspace<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  params_unnorm<-sapply(params,untranform_state_error)
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[2:length_ts]
  state_errors_max<-params_max_unnorm[2:length_ts]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on normalised state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  # prior on initial R store
  initial_state_R<-params_unnorm[length_ts*2+1]
  initial_state_R_min<-params_min_unnorm[length_ts*2+1]
  initial_state_R_max<-params_max_unnorm[length_ts*2+1]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_state_errors)+sum(log_prior_input_errors)+log_prior_initial_state+log_prior_state_errors_normalised+log_prior_initial_state_R #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}


# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithoutrouting_allinitstates3_play<-function(data,params){
  #library(Brobdingnag,lib.loc="packages")
  # Cochran's theorem shows that s^2 follows a scaled chi-squared distribution
  logd_sample_variance<-function(variance,N,expected_variance){
    frac<-variance/expected_variance
    # TODO: brobdingnag
    #frac<-as.brob(variance)/as.brob(expected_variance)
    return(dens=dchisq(frac*(N-1),N-1,log=T)+log(N-1))
  }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  var_zero_mean<-function(x) sum(x^2)/(length(x)-1)
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts+length_ts+1]*10)
  params_unnorm<-inv.normalise(params[1:(length_ts*2)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),data$normalisers[(length_ts+1):(length_ts*2)]))
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  #var_zero_mean(input_error)
  
  error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_error,input=input)
  
  if(!is.na(model_run[1])){
    error_discharge<-model_run-data$obs_discharge
    #var_zero_mean(error_discharge)
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge))
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    error_input_variance_calc<-sum(input_error^2)/(length(input_error))
    #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
    #if(error_input_variance_calc<error_input_variance) browser()
    log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
    #log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=100000,expected_variance=error_input_variance)
    # for testing
    #log_likelihood_error_input_variance<-dnorm(error_input_variance_calc/error_input_variance,1,sd=0.14,log=T)
    #     if(is.infinite(log_likelihood_error_input_variance)){
    #       browser()
    #     }
    #     state_error_sd<-params_unnorm[(length_ts+length_ts)+1]
    #     log_likelihood_state_error_sd<-sum(dnorm(state_error,mean=0,sd=state_error_sd,log=T))
    #     if(is.nan(log_likelihood_state_error_sd)) log_likelihood_state_error_sd<- -Inf
    
    # discharge error mean likelihood assuming true mean is zero
    error_discharge_mean_calc<-mean(error_discharge)
    log_likelihood_error_discharge_mean_calc<-dnorm(error_discharge_mean_calc,mean=0,sd=sqrt(error_discharge_variance)/sqrt(length(error_discharge)),log=T)
    
    # input error mean likelihood assuming true mean is zero
    error_input_mean_calc<-mean(input_error)
    log_likelihood_error_input_mean_calc<-dnorm(error_input_mean_calc,mean=0,sd=sqrt(error_input_variance)/sqrt(length(input_error)),log=T)
    combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance+log_likelihood_error_discharge_mean_calc+log_likelihood_error_input_mean_calc #+log_likelihood_state_error_sd
    
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}



trial_log_prior4_gr4jwithoutrouting_allinitstates3_play<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts+length_ts+1]*10)
  state_normaliser_min<-10^(min[length_ts+length_ts+1]*10)
  state_normaliser_max<-10^(max[length_ts+length_ts+1]*10)
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts*2)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),data$normalisers[(length_ts+1):(length_ts*2)]))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[(length_ts*2+2):(length_ts*3)]
  state_errors_max<-params_max_unnorm[(length_ts*2+2):(length_ts*3)]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on normalised state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+length_ts)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+length_ts)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  #   error_input_variance_calc<-sum(input_errors^2)/(length(input_errors)-1)
  #   #if(error_input_variance_calc<data$error_input_variance_min | error_input_variance_calc>data$error_input_variance_max) browser()
  #   log_prior_error_input_variance<-log(1/(error_input_variance_calc*(log(data$error_input_variance_max)-log(data$error_input_variance_min))))
  
  # discharge error sd calc prior?
  #   input_ts<-data$input-input_errors
  #   E_ts<-data$E
  #   if(is.null(likelihood)){
  #     input<-data.frame(P=input_ts,E=E_ts)
  #     model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_errors,input=input)
  #     error_discharge<-model_run-data$obs_discharge
  #   } else {
  #     error_discharge<-likelihood$error_discharge
  #   }
  #   error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge)-1)
  #   #if(error_discharge_variance_calc<data$error_discharge_variance_min | error_discharge_variance_calc>data$error_discharge_variance_max) browser()
  #   log_prior_error_discharge_variance<-log(1/(error_discharge_variance_calc*(log(data$error_discharge_variance_max)-log(data$error_discharge_variance_min))))
  
  
  combined_log_prior<-sum(log_prior_input_errors)+log_prior_initial_state+state_normaliser_prior+log_prior_state_errors_normalised+sum(log_prior_state_errors) #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}
