folders <- dir(path=getwd(),pattern="^SIM_Neutral_FLY")
for (i in 1:length(folders)){
  path.folder_sim <- folders[i]   
  
  files_pop_history <- dir(path.folder_sim,pattern = "^fdist_neutral")
  pop_history <- read.table(paste0(path.folder_sim,"/",files_pop_history[1]),header = F,sep=",",nrows=variables_number, colClasses = "character")
  simulation_type <- pop_history[1,2]
  
  if(simulation_type=="general"){
    variables <-  paste0( 
      "GRAL",
      "_B", pop_history[6,2], # population_size_dispersal,
      "_C", pop_history[14,2] # c_gral,
     # "_D", pop_history[26,2], # h_gral,
      #"_F", pop_history[27,2] #  s_gral,
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
       "_B", pop_history[24,2] # mutation_rate,
      # "_C", fdist_pop_history[19,2], # natural_selection_model,
      # "_D", fdist_pop_history[25,2], # loci_number,
      # # "_F", fdist_pop_history[10,2], #  same_line,
      # "_G", fdist_pop_history[11,2], #  msats,
      # "_H", fdist_pop_history[3,2], #  pre_adaptation,
      # "_I", fdist_pop_history[22,2], #  gamma_scale,
      # "_J", fdist_pop_history[23,2], #  gamma_shape,
      # "_K", fdist_pop_history[21,2], #  rate,
      # # "_L", fdist_pop_history[20,2], #  intercept,
      # "_M", fdist_pop_history[24,2], #  mutation_rate,
      # "_N", fdist_pop_history[27,2], # s_gral
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

if(B_N_test==T){
    files_B_N_neutral <- dir(path.folder_sim,pattern = paste0("^fdist_neutral_"))
    B_N_files_neutral <- lapply(paste0(path.folder_sim,"/",files_B_N_neutral),read.table,header=T,sep=",",skip=nrow(variables_file))
    B_N_files_neutral_b <- do.call(rbind, B_N_files_neutral)
    write.table(B_N_files_neutral_b,file = paste0(getwd(),"/","neutral_simulations","/","fdist_",variables,".csv"),sep = ",",row.names = F)
  }
  if(W_L_test==T){ 
      files_W_L_neutral <- dir(path.folder_sim,pattern = paste0("^outflank_neutral_"))
    W_L_files_neutral <- lapply(paste0(path.folder_sim,"/",files_W_L_neutral),read.table,header=T,sep=",",skip=nrow(variables_file))
   W_L_files_neutral_a <- do.call(rbind, W_L_files_neutral)
    outflank_pop_history_neutral <- read.table(paste0(path.folder_sim,"/",files_W_L_neutral[1]),header = F,sep=",",nrows = nrow(variables_file),colClasses = "character")
    outflank_msats <- as.numeric(outflank_pop_history_neutral[11,2])
    if(outflank_msats>2){
      He_min <- 0.2
    }else{
      He_min <- 0.1 
    }
   W_L_files_neutral_a$FSTNoCorr[W_L_files_neutral_a$FSTNoCorr==1] <- NA
  W_L_files_neutral_a <- W_L_files_neutral_a[complete.cases(W_L_files_neutral_a$FSTNoCorr),]
  W_L_files_neutral_b <- OutFLANK(FstDataFrame=W_L_files_neutral_a, LeftTrimFraction = 0.05, RightTrimFraction = 0.05, Hmin = He_min, NumberOfSamples=2, qthreshold = 0.05)

  Fstbar_neutral <- W_L_files_neutral_b$FSTbar
  df_Inferred <- W_L_files_neutral_b$dfInferred
    write.table(c(Fstbar_neutral,df_Inferred),file = paste0(getwd(),"/","neutral_simulations","/","outflank_",variables,".csv"),sep = ",",row.names = F)
  }

}