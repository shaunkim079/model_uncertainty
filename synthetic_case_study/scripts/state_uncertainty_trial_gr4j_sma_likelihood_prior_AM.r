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
log_likelihood_trial4_gr4jwithrouting_allinitstates4<-function(data,params){
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
    log_likelihood_error_discharge<-sum(dnorm(error_discharge,0,sqrt(error_discharge_variance),log=T))
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge))
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    
    log_likelihood_error_input<-sum(dnorm(input_error,0,sqrt(error_input_variance),log=T))
    
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
    
    combined_log_likehood<-log_likelihood_error_input_variance+
      log_likelihood_error_discharge_variance+
      log_likelihood_error_discharge_mean_calc+
      log_likelihood_error_input_mean_calc+
      log_likelihood_error_input+
      log_likelihood_error_discharge #+log_likelihood_state_error_sd
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior4_gr4jwithrouting_allinitstates4<-function(params,min,max,data,likelihood=NULL){
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
log_likelihood_trial4_gr4jwithoutrouting_allinitstates4<-function(data,params){
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
    log_likelihood_error_discharge<-sum(dnorm(error_discharge,0,sqrt(error_discharge_variance),log=T))
    
    error_discharge_variance_calc<-sum(error_discharge^2)/(length(error_discharge))
    #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
    log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(error_discharge),expected_variance=error_discharge_variance)  
    #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
    #     if(is.infinite(log_likelihood_error_discharge_variance)){
    #       browser()
    #     }
    
    
    log_likelihood_error_input<-sum(dnorm(input_error,0,sqrt(error_input_variance),log=T))
    
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
    
    combined_log_likehood<-log_likelihood_error_input_variance+
      log_likelihood_error_discharge_variance+
      log_likelihood_error_discharge_mean_calc+
      log_likelihood_error_input_mean_calc+
      log_likelihood_error_input+
      log_likelihood_error_discharge #+log_likelihood_state_error_sd
    
    #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
    
    #if(is.infinite(combined_log_likehood)) browser()
    # TODO: check that small probabilities are not given -Inf likelihoods
    
    return(list(combined_log_likehood=combined_log_likehood,error_discharge=error_discharge))
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA))
  }
  
}


trial_log_prior4_gr4jwithoutrouting_allinitstates4<-function(params,min,max,data,likelihood=NULL){
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
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  
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
  
  # input_errors_sample_sd<-sd_zero_mean(input_errors)
  # 1-dnorm(input_errors_sample_sd,0.17,sd=0.01)
  
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
log_likelihood_trial4_gr4jwithoutrouting_allinitstates3_hs_play<-function(data,params){
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
  
  population_flow_variance_trans<-data$population_flow_variance_trans
  # error_discharge_variance<-data$error_discharge_variance
  error_input_variance<-data$error_input_variance
  
  input_ts<-data$input-input_error
  E_ts<-data$E
  input<-data.frame(P=input_ts,E=E_ts)
  model_run<-gr4j.sma(param=data$model_param,initial_state=initial_state,state_error=state_error,input=input)
  if(!is.na(model_run[1])){
    sim_flow_trans<-log_bc_transform(model_run,data$rating_error_uncertainty$bc_lambdas)
    if(all(!is.na(sim_flow_trans))){
      # residual_trans<-data$rating_flow_trans-sim_flow_trans # ??
      # sd_zero_mean(residual_trans)
      # residual_trans<-sim_flow_trans-data$rating_flow_trans
      # sd_zero_mean(residual_trans)
      
      # residual_trans<-sim_flow_trans-log_bc_transform(data$obs_discharge,data$rating_error_uncertainty$bc_lambdas)
      residual_trans<-sim_flow_trans-data$actual_obs_flow_trans
      # sd_zero_mean(residual_trans)

      # error_discharge<-model_run-data$obs_discharge
      #var_zero_mean(error_discharge)
      error_discharge_variance_calc<-sum(residual_trans^2)/(length(residual_trans))
      # sqrt(error_discharge_variance_calc)
      #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
      log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(residual_trans),expected_variance=population_flow_variance_trans)  
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
      error_discharge_mean_calc<-mean(residual_trans)
      log_likelihood_error_discharge_mean_calc<-dnorm(error_discharge_mean_calc,mean=0,sd=sqrt(population_flow_variance_trans)/sqrt(length(residual_trans)),log=T)
      
      # input error mean likelihood assuming true mean is zero
      error_input_mean_calc<-mean(input_error)
      log_likelihood_error_input_mean_calc<-dnorm(error_input_mean_calc,mean=0,sd=sqrt(error_input_variance)/sqrt(length(input_error)),log=T)
      combined_log_likehood<-log_likelihood_error_input_variance+log_likelihood_error_discharge_variance+log_likelihood_error_discharge_mean_calc+log_likelihood_error_input_mean_calc #+log_likelihood_state_error_sd
      
      #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
      
      #if(is.infinite(combined_log_likehood)) browser()
      # TODO: check that small probabilities are not given -Inf likelihoods
      
      return(list(combined_log_likehood=combined_log_likehood,error_discharge=residual_trans,error_discharge_variance_calc=error_discharge_variance_calc))
      
    } else {
      return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
    }
    
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
  }
  
}



trial_log_prior4_gr4jwithoutrouting_allinitstates3_hs_play<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  
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
  
  # input_errors_sample_sd<-sd_zero_mean(input_errors)
  # 1-dnorm(input_errors_sample_sd,0.17,sd=0.01)
  
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
log_likelihood_trial4_gr4jwithrouting_allinitstates3_hs_play_bates<-function(data,params){
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
  state_normaliser<-10^(params[length_ts+data$number_of_sampled_inputs+2]*10)
  # params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],
  #                                                          rep(state_normaliser,length_ts-1),
  #                                                          data$normalisers[(length_ts+1):(length_ts*2)],
  #                                                          data$normalisers[length_ts*2+1]))
  params_unnorm<-inv.normalise(params[1:(length_ts+data$number_of_sampled_inputs+1)],c(data$normalisers[1],
                                                                                       rep(state_normaliser,length_ts-1),
                                                                                       data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)],
                                                                                       data$normalisers[length_ts+data$number_of_sampled_inputs+1]))
  
  
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
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
                      state_error=state_error, input=input)/1000*data$area_m2/86400
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
      # residual_trans<-data$rating_flow_trans-sim_flow_trans # ??
      # sd_zero_mean(residual_trans)
      # residual_trans<-sim_flow_trans-data$rating_flow_trans
      # sd_zero_mean(residual_trans)
      
      # residual_trans<-sim_flow_trans-log_bc_transform(data$obs_discharge,data$rating_error_uncertainty$bc_lambdas)
      
      # agregate obs flow - assume that instantaneous flow is equivalent to 10 min volume!
      sim_flow_per_day<-sim_flow_trans*24*6
      indices<-sort(rep(1:length(sim_flow),24*6))
      obs_trans_agg<-aggregate(data$rating_flow_trans,by=list(indices),sum)[,2]
      
      residual_trans<-sim_flow_per_day-obs_trans_agg
      
      # plot(sim_flow_per_day,type="l")
      # lines(obs_trans_agg,col="2",lty=2)
      # plot(residual_trans,type="l")
      # 
      # plot(sim_flow_trans,type="l")
      # lines(data$rating_flow_trans,col="2",lty=2)
      # 
      # plot(inverse_bc_transform_data(sim_flow_trans,data$rating_error_uncertainty$lambda),type="l")
      # lines(inverse_bc_transform_data(data$rating_flow_trans,data$rating_error_uncertainty$lambda),col=2)
      
      # assume that variances are additive:
      population_flow_variance_trans_agg<-population_flow_variance_trans*24*6
      
      # sd_zero_mean(residual_trans)
      
      # error_discharge<-model_run-data$obs_discharge
      #var_zero_mean(error_discharge)
      error_discharge_variance_calc<-sum(residual_trans^2)/(length(residual_trans))
      # sqrt(error_discharge_variance_calc)
      #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
      log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(residual_trans),expected_variance=population_flow_variance_trans_agg)  
      #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
      #     if(is.infinite(log_likelihood_error_discharge_variance)){
      #       browser()
      #     }
      # log_likelihood_error_discharge<-sum(dnorm(residual_trans,0,sqrt(population_flow_variance_trans_agg),log=T))
      
      # input_trans<-rep(NA,data$number_of_sampled_inputs)
      # input_trans_with_error<-rep(NA,data$number_of_sampled_inputs)
      # log_likelihoods_input<-rep(NA,data$number_of_sampled_inputs)
      # input_not_zero<-which(input_ts!=0)
      # for(ii in 1:data$number_of_sampled_inputs){
      #   # TODO: move trans below to data input to improve efficiency
      #   input_trans[ii]<-bc_transform_data(data$input[input_not_zero[ii]],
      #                                      c(data$ensemble_rainfall_stats$bcparam1[input_not_zero[ii]],data$ensemble_rainfall_stats$bcparam2[input_not_zero[ii]]))
      #   input_trans_with_error[ii]<-input_trans[ii]-input_error[ii]
      #   log_likelihoods_input[ii]<-dnorm(data$error_input_variance[input_not_zero[ii]],
      #                                    mean=input_trans[ii],
      #                                    sd=sqrt(data$error_input_variance[input_not_zero[ii]]),log = T)
      # }
      # log_likelihoods_input<-log_likelihoods_input[is.finite(log_likelihoods_input)]
      
      # error_input_variance_calc<-sum(input_error^2)/(length(input_error))
      #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
      #if(error_input_variance_calc<error_input_variance) browser()
      # log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
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
      error_discharge_mean_calc<-mean(residual_trans)
      log_likelihood_error_discharge_mean_calc<-dnorm(error_discharge_mean_calc,mean=0,sd=sqrt(population_flow_variance_trans_agg)/sqrt(length(residual_trans)),log=T)
      
      # input error mean likelihood assuming true mean is zero
      # error_input_mean_calc<-mean(input_error)
      # log_likelihood_error_input_mean_calc<-dnorm(error_input_mean_calc,mean=0,sd=sqrt(error_input_variance)/sqrt(length(input_error)),log=T)
      combined_log_likehood<-log_likelihood_error_discharge_variance+log_likelihood_error_discharge_mean_calc+sum(log_likelihoods_input) #+log_likelihood_error_discharge #+log_likelihood_error_input_mean_calc #+log_likelihood_state_error_sd #+log_likelihood_error_input_variance
      
      #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
      
      #if(is.infinite(combined_log_likehood)) browser()
      # TODO: check that small probabilities are not given -Inf likelihoods
      
      return(list(combined_log_likehood=combined_log_likehood,error_discharge=residual_trans,error_discharge_variance_calc=error_discharge_variance_calc))
      
    } else {
      return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
    }
    
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
  }
  
}



trial_log_prior4_gr4jwithrouting_allinitstates3_hs_play_bates<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts+data$number_of_sampled_inputs+2]*10)
  state_normaliser_min<-10^(min[length_ts+data$number_of_sampled_inputs+2]*10)
  state_normaliser_max<-10^(max[length_ts+data$number_of_sampled_inputs+2]*10)
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts+data$number_of_sampled_inputs+1)],c(data$normalisers[1],
                                                                                       rep(state_normaliser,length_ts-1),
                                                                                       data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)],
                                                                                       data$normalisers[length_ts+data$number_of_sampled_inputs+1]))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[(length_ts+data$number_of_sampled_inputs+3):(length_ts*2+data$number_of_sampled_inputs+1)]
  state_errors_max<-params_max_unnorm[(length_ts+data$number_of_sampled_inputs+3):(length_ts*2+data$number_of_sampled_inputs+1)]
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
  # TODO: maybe I should add this to Box-Cox transformed data, then untransform the whole thing to assess the prior probability?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  
  # prior on initial R store
  initial_state_R<-params_unnorm[length_ts+data$number_of_sampled_inputs+1]
  initial_state_R_min<-params_min_unnorm[length_ts+data$number_of_sampled_inputs+1]
  initial_state_R_max<-params_max_unnorm[length_ts+data$number_of_sampled_inputs+1]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  # input_errors_sample_sd<-sd_zero_mean(input_errors)
  # 1-dnorm(input_errors_sample_sd,0.17,sd=0.01)
  
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
log_likelihood_trial4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error<-function(data,params){
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
  state_normaliser<-10^(params[2*length_ts+data$number_of_sampled_inputs+1]*10)
  
  state_normaliser_R<-10^(params[2*length_ts+data$number_of_sampled_inputs+2]*10)
  # params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],
  #                                                          rep(state_normaliser,length_ts-1),
  #                                                          data$normalisers[(length_ts+1):(length_ts*2)],
  #                                                          data$normalisers[length_ts*2+1]))

  params_unnorm<-inv.normalise(params[1:(length_ts*2+data$number_of_sampled_inputs)],
                               c(data$normalisers[1],
                                 rep(state_normaliser,length_ts-1),
                                 data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)],
                                 data$normalisers[length_ts+data$number_of_sampled_inputs+1],
                                 rep(state_normaliser_R,length_ts-1)))
  
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[length_ts+data$number_of_sampled_inputs+1]
  state_error_R<-params_unnorm[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  
  
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
                      state_error=state_error, input=input, state_error_R=state_error_R)/1000*data$area_m2/86400
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
      # residual_trans<-data$rating_flow_trans-sim_flow_trans # ??
      # sd_zero_mean(residual_trans)
      # residual_trans<-sim_flow_trans-data$rating_flow_trans
      # sd_zero_mean(residual_trans)
      
      # residual_trans<-sim_flow_trans-log_bc_transform(data$obs_discharge,data$rating_error_uncertainty$bc_lambdas)
      
      # agregate obs flow - assume that instantaneous flow is equivalent to 10 min volume!
      sim_flow_per_day<-sim_flow_trans*24*6
      indices<-sort(rep(1:length(sim_flow),24*6))
      obs_trans_agg<-aggregate(data$rating_flow_trans,by=list(indices),sum)[,2]
      
      residual_trans<-sim_flow_per_day-obs_trans_agg
      # plot(sim_flow_per_day,type="l")
      # lines(obs_trans_agg,col="2",lty=2)
      # plot(residual_trans,type="l")
      # 
      # plot(sim_flow_trans,type="l")
      # lines(data$rating_flow_trans,col="2",lty=2)
      # 
      # plot(inverse_bc_transform_data(sim_flow_trans,data$rating_error_uncertainty$lambda),type="l")
      # lines(inverse_bc_transform_data(data$rating_flow_trans,data$rating_error_uncertainty$lambda),col=2)
      
      # assume that variances are additive:
      population_flow_variance_trans_agg<-population_flow_variance_trans*24*6
      
      # sd_zero_mean(residual_trans)
      # error_discharge<-model_run-data$obs_discharge
      #var_zero_mean(error_discharge)
      error_discharge_variance_calc<-sum(residual_trans^2)/(length(residual_trans))
      # sqrt(error_discharge_variance_calc)
      #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
      log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(residual_trans),expected_variance=population_flow_variance_trans_agg)  
      #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
      #     if(is.infinite(log_likelihood_error_discharge_variance)){
      #       browser()
      #     }
      
      # input_trans<-rep(NA,data$number_of_sampled_inputs)
      # input_trans_with_error<-rep(NA,data$number_of_sampled_inputs)
      # log_likelihoods_input<-rep(NA,data$number_of_sampled_inputs)
      # input_not_zero<-which(input_ts!=0)
      # for(ii in 1:data$number_of_sampled_inputs){
      #   # TODO: move trans below to data input to improve efficiency
      #   input_trans[ii]<-bc_transform_data(data$input[input_not_zero[ii]],
      #                                      c(data$ensemble_rainfall_stats$bcparam1[input_not_zero[ii]],data$ensemble_rainfall_stats$bcparam2[input_not_zero[ii]]))
      #   input_trans_with_error[ii]<-input_trans[ii]-input_error[ii]
      #   log_likelihoods_input[ii]<-dnorm(data$error_input_variance[input_not_zero[ii]],
      #                                    mean=input_trans[ii],
      #                                    sd=sqrt(data$error_input_variance[input_not_zero[ii]]),log = T)
      # }
      # log_likelihoods_input<-log_likelihoods_input[is.finite(log_likelihoods_input)]
      
      # error_input_variance_calc<-sum(input_error^2)/(length(input_error))
      #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
      #if(error_input_variance_calc<error_input_variance) browser()
      # log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
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
      error_discharge_mean_calc<-mean(residual_trans)
      log_likelihood_error_discharge_mean_calc<-dnorm(error_discharge_mean_calc,mean=0,sd=sqrt(population_flow_variance_trans_agg)/sqrt(length(residual_trans)),log=T)
      
      # input error mean likelihood assuming true mean is zero
      # error_input_mean_calc<-mean(input_error)
      # log_likelihood_error_input_mean_calc<-dnorm(error_input_mean_calc,mean=0,sd=sqrt(error_input_variance)/sqrt(length(input_error)),log=T)
      combined_log_likehood<-log_likelihood_error_discharge_variance+
        log_likelihood_error_discharge_mean_calc+
        sum(log_likelihoods_input) #+log_likelihood_error_input_mean_calc #+log_likelihood_state_error_sd #+log_likelihood_error_input_variance
      
      #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
      
      #if(is.infinite(combined_log_likehood)) browser()
      # TODO: check that small probabilities are not given -Inf likelihoods
      
      return(list(combined_log_likehood=combined_log_likehood,error_discharge=residual_trans,error_discharge_variance_calc=error_discharge_variance_calc))
      
    } else {
      return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
    }
    
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
  }
  
}



trial_log_prior4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts*2+data$number_of_sampled_inputs+1]*10)
  state_normaliser_min<-10^(min[length_ts*2+data$number_of_sampled_inputs+1]*10)
  state_normaliser_max<-10^(max[length_ts*2+data$number_of_sampled_inputs+1]*10)
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  state_normaliser_R<-10^(params[length_ts*2+data$number_of_sampled_inputs+2]*10)
  state_normaliser_R_min<-10^(min[length_ts*2+data$number_of_sampled_inputs+2]*10)
  state_normaliser_R_max<-10^(max[length_ts*2+data$number_of_sampled_inputs+2]*10)
  state_normaliser_R_prior<-dunif(state_normaliser_R,min=state_normaliser_R_min,max=state_normaliser_R_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts*2+data$number_of_sampled_inputs)],
                               c(data$normalisers[1],
                                 rep(state_normaliser,length_ts-1),
                                 data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)],
                                 data$normalisers[length_ts+data$number_of_sampled_inputs+1],
                                 rep(state_normaliser_R,length_ts-1)))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  # prior on unnormalised S state values
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  state_errors_max<-params_max_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on unnormalised R state values
  state_errors_R<-params_unnorm[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_min<-params_min_unnorm[(3*length_ts+data$number_of_sampled_inputs+2):(4*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_max<-params_max_unnorm[(3*length_ts+data$number_of_sampled_inputs+2):(4*length_ts+data$number_of_sampled_inputs)]
  log_prior_state_errors_R<-dunif(state_errors_R,min=state_errors_R_min,max=state_errors_R_max,log=T)

  # prior on normalised S state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  # prior on normalised R state values
  state_errors_R_normalised<-params[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_min_normalised<-min[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_max_normalised<-max[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  log_prior_state_errors_R_normalised<-sum(dunif(state_errors_R_normalised,min=state_errors_R_min_normalised,max=state_errors_R_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  # TODO: maybe I should add this to Box-Cox transformed data, then untransform the whole thing to assess the prior probability?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  
  # prior on initial R store
  initial_state_R<-params_unnorm[length_ts+data$number_of_sampled_inputs+1]
  initial_state_R_min<-params_min_unnorm[length_ts+data$number_of_sampled_inputs+1]
  initial_state_R_max<-params_max_unnorm[length_ts+data$number_of_sampled_inputs+1]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  # input_errors_sample_sd<-sd_zero_mean(input_errors)
  # 1-dnorm(input_errors_sample_sd,0.17,sd=0.01)
  
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
  
  combined_log_prior<-sum(log_prior_input_errors)+
    log_prior_initial_state+
    state_normaliser_prior+
    log_prior_state_errors_normalised+
    log_prior_initial_state_R+
    sum(log_prior_state_errors)+
    sum(log_prior_state_errors_R)+
    state_normaliser_R_prior+
    log_prior_state_errors_R_normalised
    #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}



# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_adjusted<-function(data,params){
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
  state_normaliser<-10^(params[2*length_ts+data$number_of_sampled_inputs+1]*10)
  
  state_normaliser_R<-10^(params[2*length_ts+data$number_of_sampled_inputs+2]*10)
  # params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],
  #                                                          rep(state_normaliser,length_ts-1),
  #                                                          data$normalisers[(length_ts+1):(length_ts*2)],
  #                                                          data$normalisers[length_ts*2+1]))
  
  params_unnorm<-inv.normalise(params[1:(length_ts*2+data$number_of_sampled_inputs)],
                               c(data$normalisers[1],
                                 rep(state_normaliser,length_ts-1),
                                 data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)],
                                 data$normalisers[length_ts+data$number_of_sampled_inputs+1],
                                 rep(state_normaliser_R,length_ts-1)))
  
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[length_ts+data$number_of_sampled_inputs+1]
  state_error_R<-params_unnorm[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  
  
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
                      state_error=state_error, input=input, state_error_R=state_error_R)/1000*data$area_m2/86400
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
      # residual_trans<-data$rating_flow_trans-sim_flow_trans # ??
      # sd_zero_mean(residual_trans)
      # residual_trans<-sim_flow_trans-data$rating_flow_trans
      # sd_zero_mean(residual_trans)
      
      # residual_trans<-sim_flow_trans-log_bc_transform(data$obs_discharge,data$rating_error_uncertainty$bc_lambdas)
      
      # agregate obs flow - assume that instantaneous flow is equivalent to 10 min volume!
      sim_flow_per_day<-sim_flow_trans*24*6
      indices<-sort(rep(1:length(sim_flow),24*6))
      obs_trans_agg<-aggregate(data$rating_flow_trans,by=list(indices),sum)[,2]
      
      residual_trans<-sim_flow_per_day-obs_trans_agg
      # plot(sim_flow_per_day,type="l")
      # lines(obs_trans_agg,col="2",lty=2)
      # plot(residual_trans,type="l")
      # 
      # plot(sim_flow_trans,type="l")
      # lines(data$rating_flow_trans,col="2",lty=2)
      # 
      # plot(inverse_bc_transform_data(sim_flow_trans,data$rating_error_uncertainty$lambda),type="l")
      # lines(inverse_bc_transform_data(data$rating_flow_trans,data$rating_error_uncertainty$lambda),col=2)
      
      # assume that variances are additive:
      population_flow_variance_trans_agg<-population_flow_variance_trans*24*6
      
      # sd_zero_mean(residual_trans)
      # error_discharge<-model_run-data$obs_discharge
      #var_zero_mean(error_discharge)
      error_discharge_variance_calc<-sum(residual_trans^2)/(length(residual_trans))
      # sqrt(error_discharge_variance_calc)
      #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
      log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(residual_trans),expected_variance=population_flow_variance_trans_agg)
      #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
      #     if(is.infinite(log_likelihood_error_discharge_variance)){
      #       browser()
      #     }
      log_likelihood_error_discharge<-sum(dnorm(residual_trans,0,sqrt(population_flow_variance_trans_agg),log=T))
      
      # input_trans<-rep(NA,data$number_of_sampled_inputs)
      # input_trans_with_error<-rep(NA,data$number_of_sampled_inputs)
      # log_likelihoods_input<-rep(NA,data$number_of_sampled_inputs)
      # input_not_zero<-which(input_ts!=0)
      # for(ii in 1:data$number_of_sampled_inputs){
      #   # TODO: move trans below to data input to improve efficiency
      #   input_trans[ii]<-bc_transform_data(data$input[input_not_zero[ii]],
      #                                      c(data$ensemble_rainfall_stats$bcparam1[input_not_zero[ii]],data$ensemble_rainfall_stats$bcparam2[input_not_zero[ii]]))
      #   input_trans_with_error[ii]<-input_trans[ii]-input_error[ii]
      #   log_likelihoods_input[ii]<-dnorm(data$error_input_variance[input_not_zero[ii]],
      #                                    mean=input_trans[ii],
      #                                    sd=sqrt(data$error_input_variance[input_not_zero[ii]]),log = T)
      # }
      # log_likelihoods_input<-log_likelihoods_input[is.finite(log_likelihoods_input)]
      
      # error_input_variance_calc<-sum(input_error^2)/(length(input_error))
      #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
      #if(error_input_variance_calc<error_input_variance) browser()
      # log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
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
      error_discharge_mean_calc<-mean(residual_trans)
      log_likelihood_error_discharge_mean_calc<-dnorm(error_discharge_mean_calc,mean=0,sd=sqrt(population_flow_variance_trans_agg)/sqrt(length(residual_trans)),log=T)
      
      # input error mean likelihood assuming true mean is zero
      # error_input_mean_calc<-mean(input_error)
      # log_likelihood_error_input_mean_calc<-dnorm(error_input_mean_calc,mean=0,sd=sqrt(error_input_variance)/sqrt(length(input_error)),log=T)
      combined_log_likehood<-log_likelihood_error_discharge+
        log_likelihood_error_discharge_variance+
        log_likelihood_error_discharge_mean_calc+
        sum(log_likelihoods_input) #+log_likelihood_error_input_mean_calc #+log_likelihood_state_error_sd #+log_likelihood_error_input_variance
      
      #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
      
      #if(is.infinite(combined_log_likehood)) browser()
      # TODO: check that small probabilities are not given -Inf likelihoods
      
      return(list(combined_log_likehood=combined_log_likehood,error_discharge=residual_trans,error_discharge_variance_calc=error_discharge_variance_calc))
      
    } else {
      return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
    }
    
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
  }
  
}



trial_log_prior4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_adjusted<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts*2+data$number_of_sampled_inputs+1]*10)
  state_normaliser_min<-10^(min[length_ts*2+data$number_of_sampled_inputs+1]*10)
  state_normaliser_max<-10^(max[length_ts*2+data$number_of_sampled_inputs+1]*10)
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  state_normaliser_R<-10^(params[length_ts*2+data$number_of_sampled_inputs+2]*10)
  state_normaliser_R_min<-10^(min[length_ts*2+data$number_of_sampled_inputs+2]*10)
  state_normaliser_R_max<-10^(max[length_ts*2+data$number_of_sampled_inputs+2]*10)
  state_normaliser_R_prior<-dunif(state_normaliser_R,min=state_normaliser_R_min,max=state_normaliser_R_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts*2+data$number_of_sampled_inputs)],
                               c(data$normalisers[1],
                                 rep(state_normaliser,length_ts-1),
                                 data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)],
                                 data$normalisers[length_ts+data$number_of_sampled_inputs+1],
                                 rep(state_normaliser_R,length_ts-1)))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  # prior on unnormalised S state values
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  state_errors_max<-params_max_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on unnormalised R state values
  state_errors_R<-params_unnorm[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_min<-params_min_unnorm[(3*length_ts+data$number_of_sampled_inputs+2):(4*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_max<-params_max_unnorm[(3*length_ts+data$number_of_sampled_inputs+2):(4*length_ts+data$number_of_sampled_inputs)]
  log_prior_state_errors_R<-dunif(state_errors_R,min=state_errors_R_min,max=state_errors_R_max,log=T)
  
  # prior on normalised S state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  # prior on normalised R state values
  state_errors_R_normalised<-params[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_min_normalised<-min[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_max_normalised<-max[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  log_prior_state_errors_R_normalised<-sum(dunif(state_errors_R_normalised,min=state_errors_R_min_normalised,max=state_errors_R_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  # TODO: maybe I should add this to Box-Cox transformed data, then untransform the whole thing to assess the prior probability?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  
  # prior on initial R store
  initial_state_R<-params_unnorm[length_ts+data$number_of_sampled_inputs+1]
  initial_state_R_min<-params_min_unnorm[length_ts+data$number_of_sampled_inputs+1]
  initial_state_R_max<-params_max_unnorm[length_ts+data$number_of_sampled_inputs+1]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  # input_errors_sample_sd<-sd_zero_mean(input_errors)
  # 1-dnorm(input_errors_sample_sd,0.17,sd=0.01)
  
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
  
  combined_log_prior<-sum(log_prior_input_errors)+
    log_prior_initial_state+
    state_normaliser_prior+
    log_prior_state_errors_normalised+
    log_prior_initial_state_R+
    sum(log_prior_state_errors)+
    sum(log_prior_state_errors_R)+
    state_normaliser_R_prior+
    log_prior_state_errors_R_normalised
  #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}



# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_only_adjusted<-function(data,params){
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
  # state_normaliser<-10^(params[2*length_ts+data$number_of_sampled_inputs+1]*10)
  
  state_normaliser_R<-10^(params[length_ts+data$number_of_sampled_inputs+2]*10)
  # params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],
  #                                                          rep(state_normaliser,length_ts-1),
  #                                                          data$normalisers[(length_ts+1):(length_ts*2)],
  #                                                          data$normalisers[length_ts*2+1]))
  
  params_unnorm<-inv.normalise(params[1:(length_ts+data$number_of_sampled_inputs+1)],
                               c(data$normalisers[1],
                                 rep(state_normaliser_R,length_ts-1),
                                 data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)],
                                 data$normalisers[length_ts+data$number_of_sampled_inputs+1]))
  
  initial_state<-params_unnorm[length_ts+data$number_of_sampled_inputs+1]
  # state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[1]
  state_error_R<-params_unnorm[2:length_ts]
  
  
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
                      state_error=rep(0,length_ts), input=input, state_error_R=state_error_R)/1000*data$area_m2/86400
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
      # residual_trans<-data$rating_flow_trans-sim_flow_trans # ??
      # sd_zero_mean(residual_trans)
      # residual_trans<-sim_flow_trans-data$rating_flow_trans
      # sd_zero_mean(residual_trans)
      
      # residual_trans<-sim_flow_trans-log_bc_transform(data$obs_discharge,data$rating_error_uncertainty$bc_lambdas)
      
      # agregate obs flow - assume that instantaneous flow is equivalent to 10 min volume!
      sim_flow_per_day<-sim_flow_trans*24*6
      indices<-sort(rep(1:length(sim_flow),24*6))
      obs_trans_agg<-aggregate(data$rating_flow_trans,by=list(indices),sum)[,2]
      
      residual_trans<-sim_flow_per_day-obs_trans_agg
      # plot(sim_flow_per_day,type="l")
      # lines(obs_trans_agg,col="2",lty=2)
      # plot(residual_trans,type="l")
      # 
      # plot(sim_flow_trans,type="l")
      # lines(data$rating_flow_trans,col="2",lty=2)
      # 
      # plot(inverse_bc_transform_data(sim_flow_trans,data$rating_error_uncertainty$lambda),type="l")
      # lines(inverse_bc_transform_data(data$rating_flow_trans,data$rating_error_uncertainty$lambda),col=2)
      
      # assume that variances are additive:
      population_flow_variance_trans_agg<-population_flow_variance_trans*24*6
      
      # sd_zero_mean(residual_trans)
      # error_discharge<-model_run-data$obs_discharge
      #var_zero_mean(error_discharge)
      error_discharge_variance_calc<-sum(residual_trans^2)/(length(residual_trans))
      # sqrt(error_discharge_variance_calc)
      #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
      log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(residual_trans),expected_variance=population_flow_variance_trans_agg)
      #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
      #     if(is.infinite(log_likelihood_error_discharge_variance)){
      #       browser()
      #     }
      log_likelihood_error_discharge<-sum(dnorm(residual_trans,0,sqrt(population_flow_variance_trans_agg),log=T))
      
      # input_trans<-rep(NA,data$number_of_sampled_inputs)
      # input_trans_with_error<-rep(NA,data$number_of_sampled_inputs)
      # log_likelihoods_input<-rep(NA,data$number_of_sampled_inputs)
      # input_not_zero<-which(input_ts!=0)
      # for(ii in 1:data$number_of_sampled_inputs){
      #   # TODO: move trans below to data input to improve efficiency
      #   input_trans[ii]<-bc_transform_data(data$input[input_not_zero[ii]],
      #                                      c(data$ensemble_rainfall_stats$bcparam1[input_not_zero[ii]],data$ensemble_rainfall_stats$bcparam2[input_not_zero[ii]]))
      #   input_trans_with_error[ii]<-input_trans[ii]-input_error[ii]
      #   log_likelihoods_input[ii]<-dnorm(data$error_input_variance[input_not_zero[ii]],
      #                                    mean=input_trans[ii],
      #                                    sd=sqrt(data$error_input_variance[input_not_zero[ii]]),log = T)
      # }
      # log_likelihoods_input<-log_likelihoods_input[is.finite(log_likelihoods_input)]
      
      # error_input_variance_calc<-sum(input_error^2)/(length(input_error))
      #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
      #if(error_input_variance_calc<error_input_variance) browser()
      # log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
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
      error_discharge_mean_calc<-mean(residual_trans)
      log_likelihood_error_discharge_mean_calc<-dnorm(error_discharge_mean_calc,mean=0,sd=sqrt(population_flow_variance_trans_agg)/sqrt(length(residual_trans)),log=T)
      
      # input error mean likelihood assuming true mean is zero
      # error_input_mean_calc<-mean(input_error)
      # log_likelihood_error_input_mean_calc<-dnorm(error_input_mean_calc,mean=0,sd=sqrt(error_input_variance)/sqrt(length(input_error)),log=T)
      combined_log_likehood<-log_likelihood_error_discharge+
        log_likelihood_error_discharge_variance+
        log_likelihood_error_discharge_mean_calc+
        sum(log_likelihoods_input) #+log_likelihood_error_input_mean_calc #+log_likelihood_state_error_sd #+log_likelihood_error_input_variance
      
      #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
      
      #if(is.infinite(combined_log_likehood)) browser()
      # TODO: check that small probabilities are not given -Inf likelihoods
      
      return(list(combined_log_likehood=combined_log_likehood,error_discharge=residual_trans,error_discharge_variance_calc=error_discharge_variance_calc))
      
    } else {
      return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
    }
    
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
  }
  
}



trial_log_prior4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_only_adjusted<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  length_ts<-length(data$input)
  
  state_normaliser_R<-10^(params[length_ts+data$number_of_sampled_inputs+2]*10)
  state_normaliser_R_min<-10^(min[length_ts+data$number_of_sampled_inputs+2]*10)
  state_normaliser_R_max<-10^(max[length_ts+data$number_of_sampled_inputs+2]*10)
  state_normaliser_R_prior<-dunif(state_normaliser_R,min=state_normaliser_R_min,max=state_normaliser_R_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts+data$number_of_sampled_inputs+1)],
                               c(data$normalisers[1],
                                 rep(state_normaliser_R,length_ts-1),
                                 data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)],
                                 data$normalisers[length_ts+data$number_of_sampled_inputs+1]))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[length_ts+data$number_of_sampled_inputs+1]
  initial_state_min<-params_min_unnorm[length_ts+data$number_of_sampled_inputs+1]
  initial_state_max<-params_max_unnorm[length_ts+data$number_of_sampled_inputs+1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  # prior on unnormalised S state values
  # state_errors<-params_unnorm[2:length_ts]
  # state_errors_min<-params_min_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  # state_errors_max<-params_max_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  # log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on unnormalised R state values
  state_errors_R<-params_unnorm[2:length_ts]
  state_errors_R_min<-params_min_unnorm[(length_ts+data$number_of_sampled_inputs+3):(2*length_ts+data$number_of_sampled_inputs+1)]
  state_errors_R_max<-params_max_unnorm[(length_ts+data$number_of_sampled_inputs+3):(2*length_ts+data$number_of_sampled_inputs+1)]
  log_prior_state_errors_R<-dunif(state_errors_R,min=state_errors_R_min,max=state_errors_R_max,log=T)
  
  # prior on normalised S state values
  # state_errors_normalised<-params[2:length_ts]
  # state_errors_min_normalised<-min[2:length_ts]
  # state_errors_max_normalised<-max[2:length_ts]
  # log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  # prior on normalised R state values
  state_errors_R_normalised<-params[2:length_ts]
  state_errors_R_min_normalised<-min[2:length_ts]
  state_errors_R_max_normalised<-max[2:length_ts]
  log_prior_state_errors_R_normalised<-sum(dunif(state_errors_R_normalised,min=state_errors_R_min_normalised,max=state_errors_R_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  # TODO: maybe I should add this to Box-Cox transformed data, then untransform the whole thing to assess the prior probability?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  
  # prior on initial R store
  initial_state_R<-params_unnorm[1]
  initial_state_R_min<-params_min_unnorm[1]
  initial_state_R_max<-params_max_unnorm[1]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  # input_errors_sample_sd<-sd_zero_mean(input_errors)
  # 1-dnorm(input_errors_sample_sd,0.17,sd=0.01)
  
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
  
  combined_log_prior<-sum(log_prior_input_errors)+
    log_prior_initial_state+
    # state_normaliser_prior+
    # log_prior_state_errors_normalised+
    log_prior_initial_state_R+
    # sum(log_prior_state_errors)+
    sum(log_prior_state_errors_R)+
    state_normaliser_R_prior+
    log_prior_state_errors_R_normalised
  #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}



# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_AR_laplace<-function(data,params){
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
  
  state_normaliser<-10^(params[11+data$number_of_sampled_input]*10)
  
  state_normaliser_R<-10^(params[12+data$number_of_sampled_input]*10)
  
  params_unnorm<-inv.normalise(params[1:(10+data$number_of_sampled_inputs)],
                               c(data$normalisers[1],
                                 c(1,1,state_normaliser,state_normaliser), # TODO: need more scaling factors
                                 data$normalisers[6:(5+data$number_of_sampled_inputs)],
                                 data$normalisers[6+data$number_of_sampled_inputs],
                                 c(1,1,state_normaliser_R,state_normaliser_R)))
  
  initial_state<-params_unnorm[1]
  
  library(rmutil)
  random_state_error<-rlaplace(n=length_ts-1,m=params_unnorm[4],s=params_unnorm[5])
  initial_AR_state_error<-params_unnorm[3]
  state_error<-rep(NA,length_ts-1)
  all_initial_AR_state_error<-c()
  

  input_error<-params_unnorm[6:(5+data$number_of_sampled_inputs)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[6+data$number_of_sampled_inputs]
  random_state_error_R<-rlaplace(n=length_ts-1,
                                 m=params_unnorm[9+data$number_of_sampled_inputs],
                                 s=params_unnorm[10+data$number_of_sampled_inputs])
  initial_AR_state_error_R<-params_unnorm[8+data$number_of_sampled_inputs]
  state_error_R<-rep(NA,length_ts-1)
  all_initial_AR_state_error_R<-c()
  
  for(ii in 1:length(state_error)){
    all_initial_AR_state_error<-c(all_initial_AR_state_error,initial_AR_state_error)
    state_error[ii]<-initial_AR_state_error+random_state_error[ii]
    # for next time step
    initial_AR_state_error<-state_error[ii]*params_unnorm[2]+params_unnorm[3]
    
    all_initial_AR_state_error_R<-c(all_initial_AR_state_error_R,initial_AR_state_error_R)
    state_error_R[ii]<-initial_AR_state_error_R+random_state_error_R[ii]
    # for next time step
    initial_AR_state_error_R<-state_error_R[ii]*
      params_unnorm[7+data$number_of_sampled_inputs]+
      params_unnorm[8+data$number_of_sampled_inputs]
  }
  
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
                      state_error=state_error, input=input, state_error_R=state_error_R)/1000*data$area_m2/86400
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
      # residual_trans<-data$rating_flow_trans-sim_flow_trans # ??
      # sd_zero_mean(residual_trans)
      # residual_trans<-sim_flow_trans-data$rating_flow_trans
      # sd_zero_mean(residual_trans)
      
      # residual_trans<-sim_flow_trans-log_bc_transform(data$obs_discharge,data$rating_error_uncertainty$bc_lambdas)
      
      # agregate obs flow - assume that instantaneous flow is equivalent to 10 min volume!
      sim_flow_per_day<-sim_flow_trans*24*6
      indices<-sort(rep(1:length(sim_flow),24*6))
      obs_trans_agg<-aggregate(data$rating_flow_trans,by=list(indices),sum)[,2]
      
      residual_trans<-sim_flow_per_day-obs_trans_agg
      # plot(sim_flow_per_day,type="l")
      # lines(obs_trans_agg,col="2",lty=2)
      # plot(residual_trans,type="l")
      # 
      # plot(sim_flow_trans,type="l")
      # lines(data$rating_flow_trans,col="2",lty=2)
      # 
      # plot(inverse_bc_transform_data(sim_flow_trans,data$rating_error_uncertainty$lambda),type="l")
      # lines(inverse_bc_transform_data(data$rating_flow_trans,data$rating_error_uncertainty$lambda),col=2)
      
      # assume that variances are additive:
      population_flow_variance_trans_agg<-population_flow_variance_trans*24*6
      
      # sd_zero_mean(residual_trans)
      # error_discharge<-model_run-data$obs_discharge
      #var_zero_mean(error_discharge)
      error_discharge_variance_calc<-sum(residual_trans^2)/(length(residual_trans))
      # sqrt(error_discharge_variance_calc)
      #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
      log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(residual_trans),expected_variance=population_flow_variance_trans_agg)  
      #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
      #     if(is.infinite(log_likelihood_error_discharge_variance)){
      #       browser()
      #     }
      
      log_likelihood_error_discharge<-sum(dnorm(residual_trans,0,sqrt(population_flow_variance_trans_agg),log=T))
      
      # input_trans<-rep(NA,data$number_of_sampled_inputs)
      # input_trans_with_error<-rep(NA,data$number_of_sampled_inputs)
      # log_likelihoods_input<-rep(NA,data$number_of_sampled_inputs)
      # input_not_zero<-which(input_ts!=0)
      # for(ii in 1:data$number_of_sampled_inputs){
      #   # TODO: move trans below to data input to improve efficiency
      #   input_trans[ii]<-bc_transform_data(data$input[input_not_zero[ii]],
      #                                      c(data$ensemble_rainfall_stats$bcparam1[input_not_zero[ii]],data$ensemble_rainfall_stats$bcparam2[input_not_zero[ii]]))
      #   input_trans_with_error[ii]<-input_trans[ii]-input_error[ii]
      #   log_likelihoods_input[ii]<-dnorm(data$error_input_variance[input_not_zero[ii]],
      #                                    mean=input_trans[ii],
      #                                    sd=sqrt(data$error_input_variance[input_not_zero[ii]]),log = T)
      # }
      # log_likelihoods_input<-log_likelihoods_input[is.finite(log_likelihoods_input)]
      
      # error_input_variance_calc<-sum(input_error^2)/(length(input_error))
      #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
      #if(error_input_variance_calc<error_input_variance) browser()
      # log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
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
      error_discharge_mean_calc<-mean(residual_trans)
      log_likelihood_error_discharge_mean_calc<-dnorm(error_discharge_mean_calc,mean=0,sd=sqrt(population_flow_variance_trans_agg)/sqrt(length(residual_trans)),log=T)
      
      # input error mean likelihood assuming true mean is zero
      # error_input_mean_calc<-mean(input_error)
      # log_likelihood_error_input_mean_calc<-dnorm(error_input_mean_calc,mean=0,sd=sqrt(error_input_variance)/sqrt(length(input_error)),log=T)
      combined_log_likehood<-log_likelihood_error_discharge+
        log_likelihood_error_discharge_variance+
        log_likelihood_error_discharge_mean_calc+
        sum(log_likelihoods_input) #+log_likelihood_error_input_mean_calc #+log_likelihood_state_error_sd #+log_likelihood_error_input_variance
      
      #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
      
      #if(is.infinite(combined_log_likehood)) browser()
      # TODO: check that small probabilities are not given -Inf likelihoods
      
      return(list(combined_log_likehood=combined_log_likehood,error_discharge=residual_trans,error_discharge_variance_calc=error_discharge_variance_calc))
      
    } else {
      return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
    }
    
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
  }
  
}



trial_log_prior4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_AR_laplace<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  
  # TODO: look at below for prior
  # https://stats.stackexchange.com/questions/84340/does-there-exist-a-conjugate-prior-for-the-laplace-distribution
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  length_ts<-length(data$input)
  state_normaliser<-10^(params[11+data$number_of_sampled_input]*10)
  state_normaliser_min<-10^(min[11+data$number_of_sampled_input]*10)
  state_normaliser_max<-10^(max[11+data$number_of_sampled_input]*10)
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  state_normaliser_R<-10^(params[12+data$number_of_sampled_input]*10)
  state_normaliser_R_min<-10^(min[12+data$number_of_sampled_input]*10)
  state_normaliser_R_max<-10^(max[12+data$number_of_sampled_input]*10)
  state_normaliser_R_prior<-dunif(state_normaliser_R,min=state_normaliser_R_min,max=state_normaliser_R_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(10+data$number_of_sampled_inputs)],
                               c(data$normalisers[1],
                                 c(1,1,state_normaliser,state_normaliser), # TODO: need more scaling factors
                                 data$normalisers[6:(5+data$number_of_sampled_inputs)],
                                 data$normalisers[6+data$number_of_sampled_inputs],
                                 c(1,1,state_normaliser_R,state_normaliser_R)))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  # # prior on unnormalised S state values
  # state_errors<-params_unnorm[2:length_ts]
  # state_errors_min<-params_min_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  # state_errors_max<-params_max_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  # log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  # #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # # prior on unnormalised R state values
  # state_errors_R<-params_unnorm[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  # state_errors_R_min<-params_min_unnorm[(3*length_ts+data$number_of_sampled_inputs+2):(4*length_ts+data$number_of_sampled_inputs)]
  # state_errors_R_max<-params_max_unnorm[(3*length_ts+data$number_of_sampled_inputs+2):(4*length_ts+data$number_of_sampled_inputs)]
  # log_prior_state_errors_R<-dunif(state_errors_R,min=state_errors_R_min,max=state_errors_R_max,log=T)
  
  # # prior on normalised S state values
  # state_errors_normalised<-params[2:length_ts]
  # state_errors_min_normalised<-min[2:length_ts]
  # state_errors_max_normalised<-max[2:length_ts]
  # log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  # # prior on normalised R state values
  # state_errors_R_normalised<-params[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  # state_errors_R_min_normalised<-min[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  # state_errors_R_max_normalised<-max[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  # log_prior_state_errors_R_normalised<-sum(dunif(state_errors_R_normalised,min=state_errors_R_min_normalised,max=state_errors_R_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  # TODO: maybe I should add this to Box-Cox transformed data, then untransform the whole thing to assess the prior probability?
  input_errors<-params_unnorm[6:(5+data$number_of_sampled_inputs)]
  input_errors_min<-params_min_unnorm[6:(5+data$number_of_sampled_inputs)]
  input_errors_max<-params_max_unnorm[6:(5+data$number_of_sampled_inputs)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  
  # prior on initial R store
  initial_state_R<-params_unnorm[data$number_of_sampled_inputs+6]
  initial_state_R_min<-params_min_unnorm[data$number_of_sampled_inputs+6]
  initial_state_R_max<-params_max_unnorm[data$number_of_sampled_inputs+6]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  # input_errors_sample_sd<-sd_zero_mean(input_errors)
  # 1-dnorm(input_errors_sample_sd,0.17,sd=0.01)
  
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
  
  combined_log_prior<-sum(log_prior_input_errors)+
    log_prior_initial_state+
    state_normaliser_prior+
    # log_prior_state_errors_normalised+
    log_prior_initial_state_R+
    # sum(log_prior_state_errors)+
    # sum(log_prior_state_errors_R)+
    state_normaliser_R_prior
    # log_prior_state_errors_R_normalised
  #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}



# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_AR_laplace2<-function(data,params){
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
  
  state_normaliser<-10^(params[2*length_ts+data$number_of_sampled_inputs+1]*10)
  
  state_normaliser_R<-10^(params[2*length_ts+data$number_of_sampled_inputs+2]*10)
  # params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],
  #                                                          rep(state_normaliser,length_ts-1),
  #                                                          data$normalisers[(length_ts+1):(length_ts*2)],
  #                                                          data$normalisers[length_ts*2+1]))
  
  params_unnorm<-inv.normalise(params[1:(length_ts*2+data$number_of_sampled_inputs)],
                               c(data$normalisers[1],
                                 rep(state_normaliser,length_ts-1),
                                 data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)],
                                 data$normalisers[length_ts+data$number_of_sampled_inputs+1],
                                 rep(state_normaliser_R,length_ts-1)))
  
  initial_state<-params_unnorm[1]
  state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[length_ts+data$number_of_sampled_inputs+1]
  state_error_R<-params_unnorm[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  
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
                      state_error=state_error, input=input, state_error_R=state_error_R)/1000*data$area_m2/86400
  if(!is.na(model_run[1])){
    # AR and laplace for state errors
    AR_p1_S<-params[2*length_ts+data$number_of_sampled_inputs+3]
    AR_int_scale_S<-10^(params[2*length_ts+data$number_of_sampled_inputs+5]*10)
    AR_p2_S<-params[2*length_ts+data$number_of_sampled_inputs+4]/AR_int_scale_S
    laplace_p2_S<-params[2*length_ts+data$number_of_sampled_inputs+6]/state_normaliser
    
    AR_p1_R<-params[2*length_ts+data$number_of_sampled_inputs+7]
    AR_int_scale_R<-10^(params[2*length_ts+data$number_of_sampled_inputs+9]*10)
    AR_p2_R<-params[2*length_ts+data$number_of_sampled_inputs+8]/AR_int_scale_R
    laplace_p2_R<-params[2*length_ts+data$number_of_sampled_inputs+10]/state_normaliser_R
    
    AR_component_S<-state_error[-length(state_error)]*AR_p1_S+AR_p2_S
    state_error_AR_removed_S<-state_error[-1]-AR_component_S
    AR_component_R<-state_error_R[-length(state_error_R)]*AR_p1_R+AR_p2_R
    state_error_AR_removed_R<-state_error_R[-1]-AR_component_R
    library(rmutil)
    loglike_S<-sum(dlaplace(state_error_AR_removed_S,m=0,s=laplace_p2_S,log=T))
    loglike_R<-sum(dlaplace(state_error_AR_removed_R,m=0,s=laplace_p2_R,log=T))
    
    
    # calculate for each 10 min
    sim_flow<-model_run*60*10
    # sim_10min<-c()
    # for(i in 1:length(sim_flow)){
    #   sim_10min<-c(sim_10min,rep(sim_flow[i],24*6))
    # }
    
    
    sim_flow_trans<-bc_transform_data(sim_flow,data$rating_error_uncertainty$lambda)
    # sim_flow_trans<-log_bc_transform(model_run,data$rating_error_uncertainty$bc_lambdas)
    if(all(!is.na(sim_flow_trans))){
      # residual_trans<-data$rating_flow_trans-sim_flow_trans # ??
      # sd_zero_mean(residual_trans)
      # residual_trans<-sim_flow_trans-data$rating_flow_trans
      # sd_zero_mean(residual_trans)
      
      # residual_trans<-sim_flow_trans-log_bc_transform(data$obs_discharge,data$rating_error_uncertainty$bc_lambdas)
      
      # agregate obs flow - assume that instantaneous flow is equivalent to 10 min volume!
      sim_flow_per_day<-sim_flow_trans*24*6
      indices<-sort(rep(1:length(sim_flow),24*6))
      obs_trans_agg<-aggregate(data$rating_flow_trans,by=list(indices),sum)[,2]
      
      residual_trans<-sim_flow_per_day-obs_trans_agg
      # plot(sim_flow_per_day,type="l")
      # lines(obs_trans_agg,col="2",lty=2)
      # plot(residual_trans,type="l")
      # 
      # plot(sim_flow_trans,type="l")
      # lines(data$rating_flow_trans,col="2",lty=2)
      # 
      # plot(inverse_bc_transform_data(sim_flow_trans,data$rating_error_uncertainty$lambda),type="l")
      # lines(inverse_bc_transform_data(data$rating_flow_trans,data$rating_error_uncertainty$lambda),col=2)
      
      # assume that variances are additive:
      population_flow_variance_trans_agg<-population_flow_variance_trans*24*6
      
      # sd_zero_mean(residual_trans)
      # error_discharge<-model_run-data$obs_discharge
      #var_zero_mean(error_discharge)
      error_discharge_variance_calc<-sum(residual_trans^2)/(length(residual_trans))
      # sqrt(error_discharge_variance_calc)
      #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
      log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(residual_trans),expected_variance=population_flow_variance_trans_agg)  
      #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
      #     if(is.infinite(log_likelihood_error_discharge_variance)){
      #       browser()
      #     }
      
      log_likelihood_error_discharge<-sum(dnorm(residual_trans,0,sqrt(population_flow_variance_trans_agg),log=T))
      
      # input_trans<-rep(NA,data$number_of_sampled_inputs)
      # input_trans_with_error<-rep(NA,data$number_of_sampled_inputs)
      # log_likelihoods_input<-rep(NA,data$number_of_sampled_inputs)
      # input_not_zero<-which(input_ts!=0)
      # for(ii in 1:data$number_of_sampled_inputs){
      #   # TODO: move trans below to data input to improve efficiency
      #   input_trans[ii]<-bc_transform_data(data$input[input_not_zero[ii]],
      #                                      c(data$ensemble_rainfall_stats$bcparam1[input_not_zero[ii]],data$ensemble_rainfall_stats$bcparam2[input_not_zero[ii]]))
      #   input_trans_with_error[ii]<-input_trans[ii]-input_error[ii]
      #   log_likelihoods_input[ii]<-dnorm(data$error_input_variance[input_not_zero[ii]],
      #                                    mean=input_trans[ii],
      #                                    sd=sqrt(data$error_input_variance[input_not_zero[ii]]),log = T)
      # }
      # log_likelihoods_input<-log_likelihoods_input[is.finite(log_likelihoods_input)]
      
      # error_input_variance_calc<-sum(input_error^2)/(length(input_error))
      #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
      #if(error_input_variance_calc<error_input_variance) browser()
      # log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
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
      error_discharge_mean_calc<-mean(residual_trans)
      log_likelihood_error_discharge_mean_calc<-dnorm(error_discharge_mean_calc,mean=0,sd=sqrt(population_flow_variance_trans_agg)/sqrt(length(residual_trans)),log=T)
      
      # input error mean likelihood assuming true mean is zero
      # error_input_mean_calc<-mean(input_error)
      # log_likelihood_error_input_mean_calc<-dnorm(error_input_mean_calc,mean=0,sd=sqrt(error_input_variance)/sqrt(length(input_error)),log=T)
      combined_log_likehood<-log_likelihood_error_discharge_variance+
        log_likelihood_error_discharge_mean_calc+
        sum(log_likelihoods_input)+
        loglike_S+
        loglike_R+
        log_likelihood_error_discharge #+log_likelihood_error_input_mean_calc #+log_likelihood_state_error_sd #+log_likelihood_error_input_variance
      
      #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
      
      #if(is.infinite(combined_log_likehood)) browser()
      # TODO: check that small probabilities are not given -Inf likelihoods
      
      return(list(combined_log_likehood=combined_log_likehood,error_discharge=residual_trans,error_discharge_variance_calc=error_discharge_variance_calc))
      
    } else {
      return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
    }
    
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
  }
  
}



trial_log_prior4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_AR_laplace2<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  length_ts<-length(data$input)
  state_normaliser<-10^(params[length_ts*2+data$number_of_sampled_inputs+1]*10)
  state_normaliser_min<-10^(min[length_ts*2+data$number_of_sampled_inputs+1]*10)
  state_normaliser_max<-10^(max[length_ts*2+data$number_of_sampled_inputs+1]*10)
  state_normaliser_prior<-dunif(state_normaliser,min=state_normaliser_min,max=state_normaliser_max,log=T)
  
  state_normaliser_R<-10^(params[length_ts*2+data$number_of_sampled_inputs+2]*10)
  state_normaliser_R_min<-10^(min[length_ts*2+data$number_of_sampled_inputs+2]*10)
  state_normaliser_R_max<-10^(max[length_ts*2+data$number_of_sampled_inputs+2]*10)
  state_normaliser_R_prior<-dunif(state_normaliser_R,min=state_normaliser_R_min,max=state_normaliser_R_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts*2+data$number_of_sampled_inputs)],
                               c(data$normalisers[1],
                                 rep(state_normaliser,length_ts-1),
                                 data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)],
                                 data$normalisers[length_ts+data$number_of_sampled_inputs+1],
                                 rep(state_normaliser_R,length_ts-1)))
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-params_unnorm[1]
  initial_state_min<-params_min_unnorm[1]
  initial_state_max<-params_max_unnorm[1]
  log_prior_initial_state<-dunif(initial_state,min=initial_state_min,max=initial_state_max,log=T)
  
  # prior on unnormalised S state values
  state_errors<-params_unnorm[2:length_ts]
  state_errors_min<-params_min_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  state_errors_max<-params_max_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on unnormalised R state values
  state_errors_R<-params_unnorm[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_min<-params_min_unnorm[(3*length_ts+data$number_of_sampled_inputs+2):(4*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_max<-params_max_unnorm[(3*length_ts+data$number_of_sampled_inputs+2):(4*length_ts+data$number_of_sampled_inputs)]
  log_prior_state_errors_R<-dunif(state_errors_R,min=state_errors_R_min,max=state_errors_R_max,log=T)
  
  # prior on normalised S state values
  state_errors_normalised<-params[2:length_ts]
  state_errors_min_normalised<-min[2:length_ts]
  state_errors_max_normalised<-max[2:length_ts]
  log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  # prior on normalised R state values
  state_errors_R_normalised<-params[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_min_normalised<-min[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_max_normalised<-max[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  log_prior_state_errors_R_normalised<-sum(dunif(state_errors_R_normalised,min=state_errors_R_min_normalised,max=state_errors_R_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  # TODO: maybe I should add this to Box-Cox transformed data, then untransform the whole thing to assess the prior probability?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  
  # prior on initial R store
  initial_state_R<-params_unnorm[length_ts+data$number_of_sampled_inputs+1]
  initial_state_R_min<-params_min_unnorm[length_ts+data$number_of_sampled_inputs+1]
  initial_state_R_max<-params_max_unnorm[length_ts+data$number_of_sampled_inputs+1]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  # input_errors_sample_sd<-sd_zero_mean(input_errors)
  # 1-dnorm(input_errors_sample_sd,0.17,sd=0.01)
  
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
  
  # AR_p1_S<-params[2*length_ts+data$number_of_sampled_inputs+3]
  # log_prior_AR_p1_S<-dunif(AR_p1_S,min=-2,max=2,log=T)
  # 
  # AR_p2_S<-params[2*length_ts+data$number_of_sampled_inputs+4]
  # log_prior_AR_p2_S<-dunif(AR_p1_S,min=-10,max=10,log=T)
  # 
  # laplace_p1_S<-params[2*length_ts+data$number_of_sampled_inputs+5]
  # laplace_p2_S<-params[2*length_ts+data$number_of_sampled_inputs+6]/state_normaliser
  # 
  # AR_p1_R<-params[2*length_ts+data$number_of_sampled_inputs+7]
  # log_prior_AR_p1_R<-dunif(AR_p1_R,min=-2,max=2,log=T)
  # 
  # AR_p2_R<-params[2*length_ts+data$number_of_sampled_inputs+8]
  # laplace_p1_R<-params[2*length_ts+data$number_of_sampled_inputs+9]
  # laplace_p2_R<-params[2*length_ts+data$number_of_sampled_inputs+10]/state_normaliser_R
  
  AR_p1_S<-params[2*length_ts+data$number_of_sampled_inputs+3]
  log_prior_AR_p1_S<-dunif(AR_p1_S,min=-2,max=2,log=T)
  
  AR_int_scale_S<-10^(params[2*length_ts+data$number_of_sampled_inputs+5]*10)
  log_prior_AR_int_scale_S<-dunif(params[2*length_ts+data$number_of_sampled_inputs+5],min=-1,max=1,log=T)
  
  AR_p2_S<-params[2*length_ts+data$number_of_sampled_inputs+4]/AR_int_scale_S
  log_prior_AR_p2_S<-dunif(params[2*length_ts+data$number_of_sampled_inputs+4],min=-1,max=1,log=T)
  
  laplace_p2_S<-params[2*length_ts+data$number_of_sampled_inputs+6]/state_normaliser
  log_prior_laplace_p2_S<-dunif(params[2*length_ts+data$number_of_sampled_inputs+6],min=-1,max=1,log=T)
  
  AR_p1_R<-params[2*length_ts+data$number_of_sampled_inputs+7]
  log_prior_AR_p1_R<-dunif(AR_p1_R,min=-2,max=2,log=T)
  
  AR_int_scale_R<-10^(params[2*length_ts+data$number_of_sampled_inputs+9]*10)
  log_prior_AR_int_scale_R<-dunif(params[2*length_ts+data$number_of_sampled_inputs+9],min=-1,max=1,log=T)
  
  AR_p2_R<-params[2*length_ts+data$number_of_sampled_inputs+8]/AR_int_scale_R
  log_prior_AR_p2_R<-dunif(params[2*length_ts+data$number_of_sampled_inputs+8],min=-1,max=1,log=T)
  
  laplace_p2_R<-params[2*length_ts+data$number_of_sampled_inputs+10]/state_normaliser_R
  log_prior_laplace_p2_R<-dunif(params[2*length_ts+data$number_of_sampled_inputs+10],min=-1,max=1,log=T)
  
  combined_log_prior<-sum(log_prior_input_errors)+
    log_prior_initial_state+
    state_normaliser_prior+
    log_prior_state_errors_normalised+
    log_prior_initial_state_R+
    sum(log_prior_state_errors)+
    sum(log_prior_state_errors_R)+
    state_normaliser_R_prior+
    log_prior_state_errors_R_normalised+
    log_prior_AR_p1_S+
    log_prior_AR_int_scale_S+
    log_prior_AR_p2_S+
    log_prior_laplace_p2_S+
    log_prior_AR_p1_R+
    log_prior_AR_int_scale_R+
    log_prior_AR_p2_R+
    log_prior_laplace_p2_R
    
  #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}


# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_only_adjusted_no_init_soil<-function(data,params){
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
  # state_normaliser<-10^(params[2*length_ts+data$number_of_sampled_inputs+1]*10)
  
  state_normaliser_R<-10^(params[length_ts+data$number_of_sampled_inputs+1]*10)
  # params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],
  #                                                          rep(state_normaliser,length_ts-1),
  #                                                          data$normalisers[(length_ts+1):(length_ts*2)],
  #                                                          data$normalisers[length_ts*2+1]))
  
  params_unnorm<-inv.normalise(params[1:(length_ts+data$number_of_sampled_inputs)],
                               c(data$normalisers[1],
                                 rep(state_normaliser_R,length_ts-1),
                                 data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]))
  
  initial_state<-data$initial_state_S
  # state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[1]
  state_error_R<-params_unnorm[2:length_ts]
  
  
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
                      state_error=rep(0,length_ts), input=input, state_error_R=state_error_R)/1000*data$area_m2/86400
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
      # residual_trans<-data$rating_flow_trans-sim_flow_trans # ??
      # sd_zero_mean(residual_trans)
      # residual_trans<-sim_flow_trans-data$rating_flow_trans
      # sd_zero_mean(residual_trans)
      
      # residual_trans<-sim_flow_trans-log_bc_transform(data$obs_discharge,data$rating_error_uncertainty$bc_lambdas)
      
      # agregate obs flow - assume that instantaneous flow is equivalent to 10 min volume!
      sim_flow_per_day<-sim_flow_trans*24*6
      indices<-sort(rep(1:length(sim_flow),24*6))
      obs_trans_agg<-aggregate(data$rating_flow_trans,by=list(indices),sum)[,2]
      
      residual_trans<-sim_flow_per_day-obs_trans_agg
      # plot(sim_flow_per_day,type="l")
      # lines(obs_trans_agg,col="2",lty=2)
      # plot(residual_trans,type="l")
      # 
      # plot(sim_flow_trans,type="l")
      # lines(data$rating_flow_trans,col="2",lty=2)
      # 
      # plot(inverse_bc_transform_data(sim_flow_trans,data$rating_error_uncertainty$lambda),type="l")
      # lines(inverse_bc_transform_data(data$rating_flow_trans,data$rating_error_uncertainty$lambda),col=2)
      
      # assume that variances are additive:
      population_flow_variance_trans_agg<-population_flow_variance_trans*24*6
      
      # sd_zero_mean(residual_trans)
      # error_discharge<-model_run-data$obs_discharge
      #var_zero_mean(error_discharge)
      error_discharge_variance_calc<-sum(residual_trans^2)/(length(residual_trans))
      # sqrt(error_discharge_variance_calc)
      #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
      log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(residual_trans),expected_variance=population_flow_variance_trans_agg)
      #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
      #     if(is.infinite(log_likelihood_error_discharge_variance)){
      #       browser()
      #     }
      log_likelihood_error_discharge<-sum(dnorm(residual_trans,0,sqrt(population_flow_variance_trans_agg),log=T))
      
      # input_trans<-rep(NA,data$number_of_sampled_inputs)
      # input_trans_with_error<-rep(NA,data$number_of_sampled_inputs)
      # log_likelihoods_input<-rep(NA,data$number_of_sampled_inputs)
      # input_not_zero<-which(input_ts!=0)
      # for(ii in 1:data$number_of_sampled_inputs){
      #   # TODO: move trans below to data input to improve efficiency
      #   input_trans[ii]<-bc_transform_data(data$input[input_not_zero[ii]],
      #                                      c(data$ensemble_rainfall_stats$bcparam1[input_not_zero[ii]],data$ensemble_rainfall_stats$bcparam2[input_not_zero[ii]]))
      #   input_trans_with_error[ii]<-input_trans[ii]-input_error[ii]
      #   log_likelihoods_input[ii]<-dnorm(data$error_input_variance[input_not_zero[ii]],
      #                                    mean=input_trans[ii],
      #                                    sd=sqrt(data$error_input_variance[input_not_zero[ii]]),log = T)
      # }
      # log_likelihoods_input<-log_likelihoods_input[is.finite(log_likelihoods_input)]
      
      # error_input_variance_calc<-sum(input_error^2)/(length(input_error))
      #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
      #if(error_input_variance_calc<error_input_variance) browser()
      # log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
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
      error_discharge_mean_calc<-mean(residual_trans)
      log_likelihood_error_discharge_mean_calc<-dnorm(error_discharge_mean_calc,mean=0,sd=sqrt(population_flow_variance_trans_agg)/sqrt(length(residual_trans)),log=T)
      
      # input error mean likelihood assuming true mean is zero
      # error_input_mean_calc<-mean(input_error)
      # log_likelihood_error_input_mean_calc<-dnorm(error_input_mean_calc,mean=0,sd=sqrt(error_input_variance)/sqrt(length(input_error)),log=T)
      combined_log_likehood<-log_likelihood_error_discharge+
        log_likelihood_error_discharge_variance+
        log_likelihood_error_discharge_mean_calc+
        sum(log_likelihoods_input) #+log_likelihood_error_input_mean_calc #+log_likelihood_state_error_sd #+log_likelihood_error_input_variance
      
      #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
      
      #if(is.infinite(combined_log_likehood)) browser()
      # TODO: check that small probabilities are not given -Inf likelihoods
      
      return(list(combined_log_likehood=combined_log_likehood,error_discharge=residual_trans,error_discharge_variance_calc=error_discharge_variance_calc))
      
    } else {
      return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
    }
    
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
  }
  
}



trial_log_prior4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_only_adjusted_no_init_soil<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  length_ts<-length(data$input)
  
  state_normaliser_R<-10^(params[length_ts+data$number_of_sampled_inputs+1]*10)
  state_normaliser_R_min<-10^(min[length_ts+data$number_of_sampled_inputs+1]*10)
  state_normaliser_R_max<-10^(max[length_ts+data$number_of_sampled_inputs+1]*10)
  state_normaliser_R_prior<-dunif(state_normaliser_R,min=state_normaliser_R_min,max=state_normaliser_R_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts+data$number_of_sampled_inputs)],
                               c(data$normalisers[1],
                                 rep(state_normaliser_R,length_ts-1),
                                 data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]))
  
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-data$initial_state_S
  
  # prior on unnormalised S state values
  # state_errors<-params_unnorm[2:length_ts]
  # state_errors_min<-params_min_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  # state_errors_max<-params_max_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  # log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on unnormalised R state values
  state_errors_R<-params_unnorm[2:length_ts]
  state_errors_R_min<-params_min_unnorm[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_max<-params_max_unnorm[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  log_prior_state_errors_R<-dunif(state_errors_R,min=state_errors_R_min,max=state_errors_R_max,log=T)
  
  # prior on normalised S state values
  # state_errors_normalised<-params[2:length_ts]
  # state_errors_min_normalised<-min[2:length_ts]
  # state_errors_max_normalised<-max[2:length_ts]
  # log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  # prior on normalised R state values
  state_errors_R_normalised<-params[2:length_ts]
  state_errors_R_min_normalised<-min[2:length_ts]
  state_errors_R_max_normalised<-max[2:length_ts]
  log_prior_state_errors_R_normalised<-sum(dunif(state_errors_R_normalised,min=state_errors_R_min_normalised,max=state_errors_R_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  # TODO: maybe I should add this to Box-Cox transformed data, then untransform the whole thing to assess the prior probability?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  
  # prior on initial R store
  initial_state_R<-params_unnorm[1]
  initial_state_R_min<-params_min_unnorm[1]
  initial_state_R_max<-params_max_unnorm[1]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  # input_errors_sample_sd<-sd_zero_mean(input_errors)
  # 1-dnorm(input_errors_sample_sd,0.17,sd=0.01)
  
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
  
  combined_log_prior<-sum(log_prior_input_errors)+
    # log_prior_initial_state+
    # state_normaliser_prior+
    # log_prior_state_errors_normalised+
    log_prior_initial_state_R+
    # sum(log_prior_state_errors)+
    sum(log_prior_state_errors_R)+
    state_normaliser_R_prior+
    log_prior_state_errors_R_normalised
  #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}


# only sample state error and input error using discharge error sd sd and input error sd sd, using distribution of sample variance #############################
# https://en.wikipedia.org/wiki/Variance#Distribution_of_the_sample_variance
log_likelihood_trial4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_only_adjusted_no_init_soil_set_normaliser<-function(data,params){
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
  # state_normaliser<-10^(params[2*length_ts+data$number_of_sampled_inputs+1]*10)
  state_normaliser_R<-data$state_normaliser_R
  # state_normaliser_R<-10^(params[length_ts+data$number_of_sampled_inputs+1]*10)
  # params_unnorm<-inv.normalise(params[1:(length_ts*2+1)],c(data$normalisers[1],
  #                                                          rep(state_normaliser,length_ts-1),
  #                                                          data$normalisers[(length_ts+1):(length_ts*2)],
  #                                                          data$normalisers[length_ts*2+1]))
  
  params_unnorm<-inv.normalise(params[1:(length_ts+data$number_of_sampled_inputs)],
                               c(data$normalisers[1],
                                 rep(state_normaliser_R,length_ts-1),
                                 data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]))
  
  initial_state<-data$initial_state_S
  # state_error<-params_unnorm[2:length_ts]
  input_error<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  #var_zero_mean(input_error)
  init_state_R<-params_unnorm[1]
  state_error_R<-params_unnorm[2:length_ts]
  
  
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
                      state_error=rep(0,length_ts), input=input, state_error_R=state_error_R)/1000*data$area_m2/86400
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
      # residual_trans<-data$rating_flow_trans-sim_flow_trans # ??
      # sd_zero_mean(residual_trans)
      # residual_trans<-sim_flow_trans-data$rating_flow_trans
      # sd_zero_mean(residual_trans)
      
      # residual_trans<-sim_flow_trans-log_bc_transform(data$obs_discharge,data$rating_error_uncertainty$bc_lambdas)
      
      # agregate obs flow - assume that instantaneous flow is equivalent to 10 min volume!
      sim_flow_per_day<-sim_flow_trans*24*6
      indices<-sort(rep(1:length(sim_flow),24*6))
      obs_trans_agg<-aggregate(data$rating_flow_trans,by=list(indices),sum)[,2]
      
      residual_trans<-sim_flow_per_day-obs_trans_agg
      # plot(sim_flow_per_day,type="l")
      # lines(obs_trans_agg,col="2",lty=2)
      # plot(residual_trans,type="l")
      # 
      # plot(sim_flow_trans,type="l")
      # lines(data$rating_flow_trans,col="2",lty=2)
      # 
      # plot(inverse_bc_transform_data(sim_flow_trans,data$rating_error_uncertainty$lambda),type="l")
      # lines(inverse_bc_transform_data(data$rating_flow_trans,data$rating_error_uncertainty$lambda),col=2)
      
      # assume that variances are additive:
      population_flow_variance_trans_agg<-population_flow_variance_trans*24*6
      
      # sd_zero_mean(residual_trans)
      # error_discharge<-model_run-data$obs_discharge
      #var_zero_mean(error_discharge)
      error_discharge_variance_calc<-sum(residual_trans^2)/(length(residual_trans))
      # sqrt(error_discharge_variance_calc)
      #cat("error_discharge_variance_calc =",error_discharge_variance_calc,"\n")
      log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=length(residual_trans),expected_variance=population_flow_variance_trans_agg)
      #log_likelihood_error_discharge_variance<-logd_sample_variance(variance=error_discharge_variance_calc,N=100000,expected_variance=error_discharge_variance)  
      #     if(is.infinite(log_likelihood_error_discharge_variance)){
      #       browser()
      #     }
      log_likelihood_error_discharge<-sum(dnorm(residual_trans,0,sqrt(population_flow_variance_trans_agg),log=T))
      
      # input_trans<-rep(NA,data$number_of_sampled_inputs)
      # input_trans_with_error<-rep(NA,data$number_of_sampled_inputs)
      # log_likelihoods_input<-rep(NA,data$number_of_sampled_inputs)
      # input_not_zero<-which(input_ts!=0)
      # for(ii in 1:data$number_of_sampled_inputs){
      #   # TODO: move trans below to data input to improve efficiency
      #   input_trans[ii]<-bc_transform_data(data$input[input_not_zero[ii]],
      #                                      c(data$ensemble_rainfall_stats$bcparam1[input_not_zero[ii]],data$ensemble_rainfall_stats$bcparam2[input_not_zero[ii]]))
      #   input_trans_with_error[ii]<-input_trans[ii]-input_error[ii]
      #   log_likelihoods_input[ii]<-dnorm(data$error_input_variance[input_not_zero[ii]],
      #                                    mean=input_trans[ii],
      #                                    sd=sqrt(data$error_input_variance[input_not_zero[ii]]),log = T)
      # }
      # log_likelihoods_input<-log_likelihoods_input[is.finite(log_likelihoods_input)]
      
      # error_input_variance_calc<-sum(input_error^2)/(length(input_error))
      #cat("error_input_variance_calc =",error_input_variance_calc,"\n")
      #if(error_input_variance_calc<error_input_variance) browser()
      # log_likelihood_error_input_variance<-logd_sample_variance(variance=error_input_variance_calc,N=length(input_error),expected_variance=error_input_variance)
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
      error_discharge_mean_calc<-mean(residual_trans)
      log_likelihood_error_discharge_mean_calc<-dnorm(error_discharge_mean_calc,mean=0,sd=sqrt(population_flow_variance_trans_agg)/sqrt(length(residual_trans)),log=T)
      
      # input error mean likelihood assuming true mean is zero
      # error_input_mean_calc<-mean(input_error)
      # log_likelihood_error_input_mean_calc<-dnorm(error_input_mean_calc,mean=0,sd=sqrt(error_input_variance)/sqrt(length(input_error)),log=T)
      combined_log_likehood<-log_likelihood_error_discharge+
        log_likelihood_error_discharge_variance+
        log_likelihood_error_discharge_mean_calc+
        sum(log_likelihoods_input) #+log_likelihood_error_input_mean_calc #+log_likelihood_state_error_sd #+log_likelihood_error_input_variance
      
      #cat("log_likelihood_error_input_variance =",log_likelihood_error_input_variance,"\n")
      
      #if(is.infinite(combined_log_likehood)) browser()
      # TODO: check that small probabilities are not given -Inf likelihoods
      
      return(list(combined_log_likehood=combined_log_likehood,error_discharge=residual_trans,error_discharge_variance_calc=error_discharge_variance_calc))
      
    } else {
      return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
    }
    
  } else {
    return(list(combined_log_likehood=-Inf,error_discharge=NA,error_discharge_variance_calc=NA))
  }
  
}



trial_log_prior4_gr4jwithrouting_allinitstates3_hs_play_bates_routing_error_only_adjusted_no_init_soil_set_normaliser<-function(params,min,max,data,likelihood=NULL){
  #   inv.normalise<-function(x,factor){
  #     x/factor
  #   }
  sd_zero_mean<-function(x) sqrt(mean(x^2))
  length_ts<-length(data$input)
  
  state_normaliser_R<-data$state_normaliser_R
  # state_normaliser_R<-10^(params[length_ts+data$number_of_sampled_inputs+1]*10)
  # state_normaliser_R_min<-10^(min[length_ts+data$number_of_sampled_inputs+1]*10)
  # state_normaliser_R_max<-10^(max[length_ts+data$number_of_sampled_inputs+1]*10)
  # state_normaliser_R_prior<-dunif(state_normaliser_R,min=state_normaliser_R_min,max=state_normaliser_R_max,log=T)
  
  params_unnorm<-inv.normalise(params[1:(length_ts+data$number_of_sampled_inputs)],
                               c(data$normalisers[1],
                                 rep(state_normaliser_R,length_ts-1),
                                 data$normalisers[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]))
  
  params_min_unnorm<-min #inv.normalise(min[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  params_max_unnorm<-max #inv.normalise(max[1:(length_ts+length_ts)],c(data$normalisers[1],rep(state_normaliser,length_ts-1),rep(1,length_ts)))
  
  initial_state<-data$initial_state_S
  
  # prior on unnormalised S state values
  # state_errors<-params_unnorm[2:length_ts]
  # state_errors_min<-params_min_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  # state_errors_max<-params_max_unnorm[(2*length_ts+data$number_of_sampled_inputs+3):(3*length_ts+data$number_of_sampled_inputs+1)]
  # log_prior_state_errors<-dunif(state_errors,min=state_errors_min,max=state_errors_max,log=T)
  #log_prior_state_errors<-dnorm(state_errors,mean=0,sd=0.001,log=T)
  
  # prior on unnormalised R state values
  state_errors_R<-params_unnorm[2:length_ts]
  state_errors_R_min<-params_min_unnorm[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  state_errors_R_max<-params_max_unnorm[(length_ts+data$number_of_sampled_inputs+2):(2*length_ts+data$number_of_sampled_inputs)]
  log_prior_state_errors_R<-dunif(state_errors_R,min=state_errors_R_min,max=state_errors_R_max,log=T)
  
  # prior on normalised S state values
  # state_errors_normalised<-params[2:length_ts]
  # state_errors_min_normalised<-min[2:length_ts]
  # state_errors_max_normalised<-max[2:length_ts]
  # log_prior_state_errors_normalised<-sum(dunif(state_errors_normalised,min=state_errors_min_normalised,max=state_errors_max_normalised,log=T))
  
  # prior on normalised R state values
  state_errors_R_normalised<-params[2:length_ts]
  state_errors_R_min_normalised<-min[2:length_ts]
  state_errors_R_max_normalised<-max[2:length_ts]
  log_prior_state_errors_R_normalised<-sum(dunif(state_errors_R_normalised,min=state_errors_R_min_normalised,max=state_errors_R_max_normalised,log=T))
  
  #   state_errors_sd_calc<-params_unnorm[length_ts+length_ts+1]
  #   state_errors_sd_calc<-sqrt(mean(state_errors^2))
  #   state_errors_sd_calc<-max(state_errors_sd_calc,1e-60)
  #log_prior_state_sds<-log(1/(state_errors_sd_calc*(log(data$state_sds_max)-log(data$state_sds_min))))
  #   log_prior_state_sds<-log(1/state_errors_sd_calc)
  #   if(is.nan(log_prior_state_sds)) log_prior_state_sds<- -Inf
  
  # input error sd calc prior?
  # TODO: maybe I should add this to Box-Cox transformed data, then untransform the whole thing to assess the prior probability?
  input_errors<-params_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_min<-params_min_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  input_errors_max<-params_max_unnorm[(length_ts+1):(length_ts+data$number_of_sampled_inputs)]
  log_prior_input_errors<-dunif(input_errors,min=input_errors_min,max=input_errors_max,log=T)
  
  
  # prior on initial R store
  initial_state_R<-params_unnorm[1]
  initial_state_R_min<-params_min_unnorm[1]
  initial_state_R_max<-params_max_unnorm[1]
  log_prior_initial_state_R<-dunif(initial_state_R,min=initial_state_R_min,max=initial_state_R_max,log=T)
  
  # input_errors_sample_sd<-sd_zero_mean(input_errors)
  # 1-dnorm(input_errors_sample_sd,0.17,sd=0.01)
  
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
  
  combined_log_prior<-sum(log_prior_input_errors)+
    # log_prior_initial_state+
    # state_normaliser_prior+
    # log_prior_state_errors_normalised+
    log_prior_initial_state_R+
    # sum(log_prior_state_errors)+
    sum(log_prior_state_errors_R)+
    # state_normaliser_R_prior+
    log_prior_state_errors_R_normalised
  #+log_prior_state_sds #+log_prior_error_input_variance+log_prior_error_discharge_variance
  
  #if(is.infinite(combined_log_prior)) browser()
  
  return(combined_log_prior)
  
}



