

if(length(commandArgs(TRUE))!=0){
  args <- commandArgs(TRUE)
  cat("arguments:",args,"\n")
  start_ts<-as.numeric(args[[1]])
}


if(Sys.info()["sysname"]=="Linux"){
  wd<-"/data/kim079/model_optimisation_framework_v2"
  rscript<-"/data/kim079/model_optimisation_framework_v2/scripts/Adaptive_Metropolis_state_uncertainty_bates_routing_error_adjusted_likelihood_v5.8.r"
} else {
  wd<-"C:/Users/kim079/Documents/model_optimisation_framework"
  rscript<-"C:/Users/kim079/Documents/model_optimisation_framework/scripts/Adaptive_Metropolis_state_uncertainty_bates_routing_error_adjusted_likelihood_v5.8.r"
}

setwd(wd)

# output_dir<-"//pearceydata.csiro.au/data/kim079/model_optimisation_framework_v2/output/state_uncertainty/AM/bates"

output_dir<-"//lw-osm-03-cdc.it.csiro.au/OSM_CBR_LW_TVD_MDBA_work/8_Working/7_Shaun/state_uncertainty/AM/bates/adjusted_likelihood_v5.8"
# start_ts<-1
if(!exists("start_ts")) start_ts<-1 #1,61,121,181,241,301
length_ts<-60
restart_number<-500000

replace_this<-"//lw-osm-03-cdc.it.csiro.au/OSM_CBR_LW_TVD_MDBA_work"
replace_with<-"/OSM/CBR/LW_TVD_MDBA/work"

if(Sys.info()["sysname"]=="Linux"){
  output_dir<-gsub(replace_this,replace_with,output_dir)
}

while(T){
  new_seed<-as.character(as.numeric(Sys.time()))
  initial_state_errors_initial_file<-list.files(output_dir,pattern=paste0("^initial_state_errors_initial_",start_ts,"_",length_ts,".csv"),full.names = T)
  if(length(initial_state_errors_initial_file)==0){
    cat("Error: there is no initial_state_errors_initial_file file \n")
    break
  }
  state_errors_initial_file<-list.files(output_dir,pattern=paste0("^state_errors_initial_",start_ts,"_",length_ts,".csv"),full.names = T)
  if(length(state_errors_initial_file)==0){
    cat("Error: there is no state_errors_initial_file file \n")
    break
  }
  # find most recent init cov file
  all_init_cov_file<-list.files(output_dir,pattern=paste0("CovPar_bates_",start_ts,"_",length_ts,"_.*.csv.gz"))
  if(length(all_init_cov_file)==0){

    arguments<-c(start_ts,length_ts,output_dir,initial_state_errors_initial_file,state_errors_initial_file)
  } else {
    old_seeds<-gsub(paste0("CovPar_bates_",start_ts,"_",length_ts,"_|.csv.gz"),"",all_init_cov_file)
    old_seed<-as.character(max(as.numeric(old_seeds)))
    
    init_cov_file<-paste0(output_dir,"/CovPar_bates_",start_ts,"_",length_ts,"_",old_seed,".csv.gz")
    SD2_file<-paste0(output_dir,"/SD2_bates_",start_ts,"_",length_ts,"_",old_seed,".csv.gz")
    theta_file<-paste0(output_dir,"/theta_bates_",start_ts,"_",length_ts,"_",old_seed,".csv.gz")
    if(!file.exists(init_cov_file) | !file.exists(SD2_file) | !file.exists(theta_file)) stop("can't find file")
    
    arguments<-c(start_ts,length_ts,output_dir,initial_state_errors_initial_file,state_errors_initial_file,
                 init_cov_file,SD2_file,theta_file,restart_number,new_seed)
    seed<-new_seed
  }

  if(Sys.info()["sysname"]=="Linux"){
    arguments<-gsub(replace_this,replace_with,arguments)
  }
  
  # arguments<-paste(start_ts,length_ts,init_cov_file,SD2_file,theta_file,restart_number,new_seed)
  # command<-paste("Rscript,rscript,arguments)
  # system(command)
  # stop()
  system2("Rscript",c(rscript,arguments))
}








