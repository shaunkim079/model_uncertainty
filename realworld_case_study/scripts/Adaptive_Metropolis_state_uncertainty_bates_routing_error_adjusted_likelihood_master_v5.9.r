

if(Sys.info()["sysname"]=="Linux"){
  wd<-"/data/kim079/model_optimisation_framework_v2"
  rscript_runner<-"/data/kim079/model_optimisation_framework_v2/scripts/Adaptive_Metropolis_state_uncertainty_bates_routing_error_adjusted_likelihood_runner_v5.9.r"
} else {
  wd<-"C:/Users/kim079/Documents/model_optimisation_framework"
  rscript_runner<-"C:/Users/kim079/Documents/model_optimisation_framework/scripts/Adaptive_Metropolis_state_uncertainty_bates_routing_error_adjusted_likelihood_runner_v5.9.r"
}

all_start_ts<-c(1,61,121,181,241,301)


for(start_ts_index in 1:length(all_start_ts)){
  system2("Rscript",c(rscript_runner,all_start_ts[start_ts_index]),wait=F,invisible = F)
}
