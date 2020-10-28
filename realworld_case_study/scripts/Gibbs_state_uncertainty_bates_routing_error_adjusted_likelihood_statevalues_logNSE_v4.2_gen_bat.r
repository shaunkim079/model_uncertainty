batch.file.initial.lines<-c(
  "#!/bin/bash",
  "#SBATCH --job-name=batsvlN4.2",
  "#SBATCH --time=72:00:00",
  "#SBATCH --nodes=1",
  #"#SBATCH --ntasks-per-node=20",
  "#SBATCH --cpus-per-task=1",
  "#SBATCH --mem=5000",
  "",
  "module load R/3.4.0" #"module load R/3.4.0"
)
# #SBATCH --time= 96:00:00


if(Sys.info()["sysname"]=="Linux"){
  wd<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/model_optimisation_framework_v2"
  rscript_runner<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/model_optimisation_framework_v2/scripts/Gibbs_state_uncertainty_bates_routing_error_adjusted_likelihood_runner_statevalues_logNSE_v4.2.r"
} else {
  # wd<-"C:/Users/kim079/Documents/model_optimisation_framework"
  wd<-"Y:/"
  # rscript_runner<-"C:/Users/kim079/Documents/model_optimisation_framework/scripts/Gibbs_state_uncertainty_bates_routing_error_adjusted_likelihood_runner_v10.r"
  rscript_runner<-"Y:/scripts/Gibbs_state_uncertainty_bates_routing_error_adjusted_likelihood_runner_statevalues_logNSE_v4.2.r"

}

remove_for_linux<-"//pearceydata.csiro.au|//pearceyflush1.csiro.au|//pearceyflush2.csiro.au"

# replace_this<-"//gpfs2-cbr.san.csiro.au/lw_tvd_mdba_work" #"//lw-osm-03-cdc.it.csiro.au/OSM_CBR_LW_TVD_MDBA_work"
# replace_with<-"/datasets/work/LW_TVD_MDBA_WORK" #"/OSM/CBR/LW_TVD_MDBA/work"
replace_this<-"Y:"
replace_with<-"/datasets/work/LW_TVD_MDBA_WORK/8_Working/7_Shaun/data_backup/kim079/model_optimisation_framework_v2"

output_dir<-"output/state_uncertainty/AM/bates/adjusted_likelihood_statevalues_logNSE_v4.2"
batch.write.dir<-"scripts/bates_statevalues_logNSE_v4.2"


setwd(wd)

# dir.create(output_dir, showWarnings = F)
dir.create(batch.write.dir, showWarnings = F)

source("scripts/Gibbs_state_uncertainty_bates_routing_error_adjusted_likelihood_statevalues_gen_init_theta_v4.2.r")
use_init_DEoptim_parameters(output_dir,batch.write.dir)

all_start_ts<-c(1,31,61,91,121,151,181,211,241,271,301,331)
# all_start_ts<-c(1,31)
all_seeds<-1:4

all_batch_names<-c()
for(start_ts_index in 1:length(all_start_ts)){
  for(iseeds in all_seeds){
    
    
    command<-paste("Rscript",rscript_runner,all_start_ts[start_ts_index],iseeds)
    batch.file.lines<-c(batch.file.initial.lines,command)
    batch_fn<-paste0(batch.write.dir,"/",all_start_ts[start_ts_index],"_",iseeds)
    all_batch_names<-c(all_batch_names,basename(batch_fn))
    batch.file.lines<-gsub(remove_for_linux,"",batch.file.lines)
    batch.file.lines<-gsub(replace_this,replace_with,batch.file.lines)
    writeLines(batch.file.lines,batch_fn)
    
  }
}





batch_runner_fn<-paste0(batch.write.dir,"/batch_runner.sh")
all_cluster_lines<-c(paste("dos2unix",all_batch_names),paste("sbatch",all_batch_names))
writeLines(all_cluster_lines,batch_runner_fn)

