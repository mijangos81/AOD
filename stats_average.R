  folders <- paste0(getwd(),"/",dir(path=getwd(),pattern="^SIM"))
    # folders <- paste0(getwd(),"/",dir(path=getwd(),pattern="^AOD_fly"))
    # folders <- paste0(getwd(),"/",dir(path=getwd(),pattern="^first_sim"))

  variables_number <- 36

for(i in 1:length(folders)){
path.folder_stats_average <- folders[i]    

files_ave <-paste0(path.folder_stats_average,"/",dir(path.folder_stats_average,pattern = "generations"))
# files_ave <-paste0(path.folder_stats_average,"/",dir(path.folder_stats_average,pattern = "\\.csv$"))
pop_history <- read.table(files_ave[1],header = F,sep=",",nrows = variables_number, colClasses = "character")
simulation_type <- pop_history[1,2]

if(simulation_type=="general"){
 variables <-  paste0( 
     "GRAL",
    "_B", pop_history[6,2], # population_size_dispersal,
     "_C", pop_history[14,2], # c_gral,
     "_D", pop_history[26,2], # h_gral,
      "_F", pop_history[27,2] #  s_gral,
    # "_G", pop_history[11,2], #  msats,
     # "_H", pop_history[3,2], #  pre_adaptation,
    # "_I", pop_history[22,2], #  gamma_scale,
    # "_J", pop_history[23,2], #  gamma_shape,
    # "_K", pop_history[21,2], #  rate,
    # "_L", pop_history[20,2], #  intercept,
    # "_M", pop_history[24,2], #  mutation_rate,
    # "_N", pop_history[27,2], # s_gral
    # "_O", pop_history[26,2], # h_gral,
    # "_P", pop_history[28,2], #q_gral,
    # "_Q", pop_history[14,2], #c_gral
    # "_R",pop_history[15,2] # windows_gral,
    #  "_S",pop_history[31,2], # experiment_freq,
    #  "_T",pop_history[29,2], # fitness_for_repulsion,
    # "_rep_", iteration,
    # ".csv"
    )
}
if(simulation_type=="fly"){
    variables <-  paste0(
      "FLY",
       "_B_", pop_history[7,2], # rate,
      "_F_", pop_history[10,2] #  mutation_rate,
      #  "_B_", pop_history[21,2], # rate,
      # "_C_", pop_history[23,2], # gamma_scale,
      # "_D_", pop_history[24,2], # gamma_shape,
      # "_F_", pop_history[18,2] #  mutation_rate,
      # "_G", fdist_pop_history[11,2], #  msats,
      # "_H", fdist_pop_history[3,2], #  pre_adaptation,
      # "_I", fdist_pop_history[22,2], #  gamma_scale,
      # "_J", fdist_pop_history[23,2], #  gamma_shape,
      # "_K", fdist_pop_history[21,2], #  rate,
      # # "_L", fdist_pop_history[20,2], #  intercept,
      # "_M", fdist_pop_history[24,2], #  mutation_rate,
      # # "_N", fdist_pop_history[27,2], # s_gral
      # "_O", fdist_pop_history[26,2], # h_gral,
      # "_P", fdist_pop_history[28,2], #q_gral,
      # "_Q", fdist_pop_history[14,2] #c_gral
      # "_R",fdist_pop_history[15,2], # windows_gral,
      # "_S",fdist_pop_history[31,2] # experiment_freq
      # "_T",fdist_pop_history[29,2] # fitness_for_repulsion,
      # "_rep_", iteration,
      # ".csv"
    )
  }


if (length(files_ave)>1) {
  myfiles <- lapply(files_ave, read.table,header=F,sep=",",skip= variables_number,fill=T)
  mean_df <- aaply(laply(myfiles, as.matrix), c(2, 3), mean,na.rm=T)
  # sd_df <- aaply(laply(myfiles, as.matrix), c(2, 3), sd,na.rm=T)
}else{
  myfiles <- read.table(files_ave,header=F,sep=",",skip=variables_number,fill=T)
  mean_df <- myfiles
  sd_df <- myfiles
}

# if (neutral_simulations==T){
#   # This is the final dataframe 
# write.table(pop_history,file = paste0(getwd(),"/","neutral_simulations","/","final_stats_average_neutral_",variables,".csv"),sep = ",",row.names = F,col.names = F)
# suppressWarnings(write.table(mean_df,file = paste0(getwd(),"/","neutral_simulations","/","final_stats_average_neutral_",variables,".csv"),sep = ",",row.names = F,append = T,col.names=F))
# 
#  write.table(pop_history,file = paste0(getwd(),"/","neutral_simulations","/","final_stats_sd_neutral_",variables,".csv"),sep = ",",row.names = F,col.names = F)
#  suppressWarnings(write.table(sd_df,file = paste0(getwd(),"/","neutral_simulations","/","final_stats_sd_neutral_",variables,".csv"),sep = ",",row.names = F,append = T,col.names=F))
#  
# }
# else{
# This is the final dataframe 
write.table(pop_history,file = paste0(getwd(),"/","final_stats_average_",variables,".csv"),sep = ",",row.names = F,col.names = F)
suppressWarnings(write.table(mean_df,file = paste0(getwd(),"/","final_stats_average_",variables,".csv"),sep = ",",row.names = F,append = T,col.names=F))

 # write.table(pop_history,file = paste0(getwd(),"/","final_stats_sd_",variables,".csv"),sep = ",",row.names = F,col.names = F)
 # suppressWarnings(write.table(sd_df,file = paste0(getwd(),"/","final_stats_sd_",variables,".csv"),sep = ",",row.names = F,append = T,col.names=F))
}
# }

