

source('functions.R')
source('hierfstat.R')
source('Fdist_functions.R')
source('Bayescan_functions.R')
source('Outflank_functions.R')
folders <- dir(path=getwd(),pattern="^SIM")
for (i in 1:length(folders)){
  path.folder_sim <- folders[i]
  # VARIABLES
  B_N_test <- T # Beaumont and Nichols' test (FDIST2)
  B_B_test <- T # Balding and Beaumont's test (BAYESCAN)
  W_L_test <- T # Whitlock and Lotterhos' test (OUTFLANK)
  # CONSTANT VARIABLES
  variables_number <- 36
  bayescan.path <- "/usr/local/bin/bayescan" # this is the location of the binary file of the program BAYESCAN
  number_cores <- 7 # this is the number of cores that the program BAYESCAN uses
  # Bayescan settings. these setting have been tested to ensure that convergence of the
  # Bayescan's RJ-MCMC algorithm has been reached. For this we used the function gelman.diag
  # from the R package coda.
  n_iter <- 3000 # number of outputted iterations
  thinning <- 10 # thinning interval size
  n_pilot <- 10 # number of pilot runs
  l_pilot <- 5000 # length of pilot runs
  burnin <- 5000 # Burn-in length
  # the width of the bins (number of observations) used to estimate the density of
  # the Johnson distribution.
  bins_Johnson_distribution <- 400
  # DATAFRAMES
  selection_tests <- as.data.frame(matrix(nrow = length(folders),ncol = 16))
  colnames(selection_tests) <- c("FDIST2_FST_balancing","FDIST2_Shua_balancing","FDIST2_FST_Sha_diff_balancing","FDIST2_Sha_pop_balancing","BAYESCAN_A_balancing","BAYESCAN_B_balancing","OUTFLANK_A_balancing","OUTFLANK_B_balancing","FDIST2_FST_directional","FDIST2_Shua_directional","FDIST2_FST_Sha_diff_directional","FDIST2_Sha_pop_directional","BAYESCAN_A_directional","BAYESCAN_B_directional","OUTFLANK_A_directional","OUTFLANK_B_directional")

  files_B_N_emp <- dir(path.folder_sim,pattern = "^fdist_selection")
  fdist_pop_history <- read.table(paste0(path.folder_sim,"/",files_B_N_emp[1]),header = F,sep=",",nrows=variables_number, colClasses = "character")
  simulation_type <- fdist_pop_history[1,2]
  pop_history<- fdist_pop_history
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
       "_B", pop_history[21,2], # rate,
      "_C", pop_history[22,2], # gamma_scale,
      "_D", pop_history[23,2], # gamma_shape,
      "_F", pop_history[24,2] #  mutation_rate,
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
  ##### Beaumont and Nichol's method (B_N_method or FDIST2) #####
  if(B_N_test==TRUE){
    message("FDIST")
    files_B_N_emp <- dir(path.folder_sim,pattern = "^fdist_selection")
    B_N_files_emp <- lapply(paste0(path.folder_sim,"/",files_B_N_emp), read.table,header=T,sep=",",skip=variables_number)
   
    
    #  # nrow(variables_file))
    # # merging replicates in groups of 10
    # if (length(files_B_N_emp) %% 10 != 0){
    #   files_to_exclude <- length(files_B_N_emp) %% 10
    #   files_B_N_emp <- files_B_N_emp[1:(length(files_B_N_emp)-files_to_exclude)]
    # }
    # seq_start <- seq(1,length(files_B_N_emp),10)
    # seq_end <- seq(10,length(files_B_N_emp),10)
    # 
    # B_N_files_emp_2 <- list()
    # for(i in 1:length(seq_start) ){
    #   B_N_files_emp_temp <- rbindlist(B_N_files_emp[seq_start[i]:seq_end[i]])
    #   B_N_files_emp_2[[i]] <- B_N_files_emp_temp
    # }
    
    
B_N_files_emp_2 <- B_N_files_emp
    # fdist_pop_history <- read.table(paste0(path.folder_sim,"/",files_B_N_emp[1]),header = F,sep=",",nrows=variables_number, colClasses = "character")
    # nrow(variables_file))
    simulation_type <- fdist_pop_history[1,2]
    if(simulation_type=="fly"){
      fdist_pop_history_b <- fdist_pop_history[24,2]
      fdist_pop_history_c <- paste0("GRAL_","B",fdist_pop_history_b[1],collapse = "")
      files_B_N_neutral <- dir(paste0(getwd(),"/","neutral_simulations"),pattern = paste0("fdist_",fdist_pop_history_c,"\\b"))
      }
    if(simulation_type=="general"){
      fdist_pop_history_b <- fdist_pop_history[c(6,14),2]
      fdist_pop_history_c <- paste0("GRAL_","B",fdist_pop_history_b[1],"_C",fdist_pop_history_b[2],collapse = "")
      files_B_N_neutral <- dir(paste0(getwd(),"/","neutral_simulations"),pattern = paste0("fdist_",fdist_pop_history_c,"\\b"))
    }
    B_N_files_neutral <- read.table(paste0(getwd(),"/","neutral_simulations","/",files_B_N_neutral),header=T,sep=",")
    # removing the cases when one allele is fixed
    B_N_files_neutral$het[B_N_files_neutral$het==1] <- NA
    B_N_files_neutral <- B_N_files_neutral[complete.cases(B_N_files_neutral),]
    B_N_files_neutral <- B_N_files_neutral[sample(x=1:nrow(B_N_files_neutral) ,size = 1500),]


    FST.simDF_res <- B_N_files_neutral[,c("het","FST")]
    colnames(FST.simDF_res) <- c("alpha","beta")
    # ranking by heterozygosity as in Beaumont & Nichols 1996
    FST.simDF_res <- FST.simDF_res[order(FST.simDF_res$alpha),]
    FST.simDF_res <- FST.simDF_res[complete.cases(FST.simDF_res),]

    # Shannon.simDF_res <- B_N_files_neutral[,1:7]
    # Shannon.simDF_res[which(Shannon.simDF_res$sha_diff==1),] <- NA
    # Shannon.simDF_res[which(Shannon.simDF_res$sha_diff==0),] <- NA
    # Shannon.simDF_res <- Shannon.simDF_res[complete.cases(Shannon.simDF_res),]
    # 
    # Shannon_mean_sim <- cbind(Shannon.simDF_res$sha_pop1,Shannon.simDF_res$sha_pop2)
    # Shannon_mean_sim <- apply(Shannon_mean_sim,1,mean)
    # Shua.simDF_res <- as.data.frame(cbind(Shannon_mean_sim,Shannon.simDF_res$shua))
    # colnames(Shua.simDF_res) <- c("alpha","beta")
    # Shua.simDF_res <- Shua.simDF_res[order(Shua.simDF_res$alpha),]
    # 
    # D_Shannon_mean_sim <- cbind(Shannon.simDF_res$one_D_alpha_pop1,Shannon.simDF_res$one_D_alpha_pop2)
    # D_Shannon_mean_sim <- apply(D_Shannon_mean_sim,1,mean)
    # Sha_diff.simDF_res <- as.data.frame(cbind(D_Shannon_mean_sim,Shannon.simDF_res$sha_diff))
    # colnames(Sha_diff.simDF_res) <- c("alpha","beta")
    # Sha_diff.simDF_res <- Sha_diff.simDF_res[order(Sha_diff.simDF_res$alpha),]
    # 
    # Sha_pops.simDF_res <- as.data.frame(cbind(Shannon.simDF_res$sha_pops,Shannon.simDF_res$shua))
    # colnames(Sha_pops.simDF_res) <- c("alpha","beta")
    # Sha_pops.simDF_res <- Sha_pops.simDF_res[order(Sha_pops.simDF_res$alpha),]

    for(i in 1:length(B_N_files_emp_2)){
      FST.empDF_res <- B_N_files_emp_2[[i]][,c("het","FST")]
      FST.empDF_res$selection <- T
      FST.empDF_res[as.numeric(neutral_loci_location),"selection"] <- F
      colnames(FST.empDF_res) <- c("alpha","beta","selection")
      FST.empDF_res <- FST.empDF_res[complete.cases(FST.empDF_res),]
      FST_fdist_res <- get.pval.dataset(beta.empDF=FST.empDF_res, beta.simDF=FST.simDF_res, beta_max=1.01)
      FST_fdist_res_b <- correct.pval.dataframe(dataframe=FST_fdist_res, p.colName="p.val.cum")
      FST_fdist_res_balancing <- FST_fdist_res_b[FST_fdist_res_b$tail=="L" & FST_fdist_res_b$FDR.05==TRUE, ]
      selection_tests[i,1] <- nrow(FST_fdist_res_balancing)
      FST_fdist_res_directional <- FST_fdist_res_b[FST_fdist_res_b$tail=="R" & FST_fdist_res_b$FDR.05==TRUE, ]
      selection_tests[i,9] <- nrow(FST_fdist_res_directional)

 
    }
  }
  ##### Balding and Beaumont's method (B_B_method or BAYESCAN)#####
  if(B_B_test==TRUE){
    message("BAYESCAN")

    files_B_B <- dir(path.folder_sim,pattern = "^bayescan")
    files_B_B_2 <- lapply(paste0(path.folder_sim,"/",files_B_B),read.table,header=T,sep=",")
    files_B_B_2 <- lapply(files_B_B_2,"[",,-1)
   
    
    #  # merging replicates in groups of 10
    #    if (length(files_B_B_2) %% 10 != 0){
    #   files_to_exclude <- length(files_B_B_2) %% 10
    #   files_B_B_2 <- files_B_B_2[1:(length(files_B_B_2)-files_to_exclude)]
    # }
    # seq_start <- seq(1,length(files_B_B_2),10)
    # seq_end <- seq(10,length(files_B_B_2),10)
    # 
    # files_B_B_3 <- list()
    # for(i in 1:length(seq_start) ){
    #   files_B_B_3_temp <- list.cbind(files_B_B_2[seq_start[i]:seq_end[i]])
    #   files_B_B_3[[i]] <- files_B_B_3_temp
    # }
    
    
    files_B_B_3 <- files_B_B_2
    files_B_B_3 <- lapply(files_B_B_3,function(x){colnames(x) <- c(1:ncol(files_B_B_3[[1]]));x})
    pops <- c(rep(1,nrow(files_B_B_3[[1]])/2),rep(2,nrow(files_B_B_3[[1]])/2))
    for(i in 1:length(files_B_B_3)){
      files_B_B_3[[i]] <- cbind(pops,files_B_B_3[[i]])
      write.bayescan(files_B_B_3[[i]],file_name = paste0(path.folder_sim,"/","group_bayescan",variables,"_rep_", i,".txt"))
    }
    files_B_B_group_temp <- dir(path.folder_sim,pattern = "^group_bayescan")
    # running bayescan in just 10% of the groups
    files_B_B_group <- sample(files_B_B_group_temp,size = ceiling(length(files_B_B_group_temp)*0.1))
    # Files should be in the working directory for bayescan to work.
    working_dir <- getwd()
    setwd(path.folder_sim)
    for(i in 1:length(files_B_B_group)){
      B_B_res_1 <- run_bayescan(data = files_B_B_group[i],n=n_iter,thin=thinning,nbp=n_pilot,pilot=l_pilot,burn=burnin,pr_odds=1,parallel.core=number_cores,bayescan.path = bayescan.path,time_report = T)

      #B_B_res_1 <- run_bayescan(data = files_B_B_group[i],n=n_iter,thin=thinning,nbp=n_pilot,pilot=l_pilot,burn=burnin,pr_odds=1,parallel.core=number_cores,bayescan.path = bayescan.path,time_report = T)
      ### BayeScan outputs a file with q-values, but not p-values.
      # ### this function returns corrected p values for a Bayescan outfile
      # # Function taken from Whitlock and Lotterhos 2014
      B_B_res_1_b <- ReturnCorrectedPVal.BS(bayescan_ouput=B_B_res_1$bayescan)
      B_B_res_1_balancing <- B_B_res_1_b[B_B_res_1_b$SELECTION == "balancing",]
      selection_tests[i,5] <- sum(B_B_res_1_balancing$FDR.05)
      B_B_res_1_directional <- B_B_res_1_b[B_B_res_1_b$SELECTION == "diversifying",]
      selection_tests[i,13] <- sum(B_B_res_1_directional$FDR.05)
      B_B_res_2 <- run_bayescan(data= files_B_B_group[i],n=n_iter,thin=thinning,nbp=n_pilot,pilot=l_pilot,burn=burnin,pr_odds=10,parallel.core=number_cores,bayescan.path = bayescan.path,time_report = T)
      # ### BayeScan outputs a file with q-values, but not p-values.
      # ### this function returns corrected p values for a Bayescan outfile
      # # Function taken from Whitlock and Lotterhos 2014
      B_B_res_2_b <- ReturnCorrectedPVal.BS(bayescan_ouput=B_B_res_2$bayescan)
      B_B_res_2_balancing <- B_B_res_2_b[B_B_res_2_b$SELECTION == "balancing",]
      selection_tests[i,6] <- sum(B_B_res_2_balancing$FDR.05)
      B_B_res_2_directional <- B_B_res_2_b[B_B_res_2_b$SELECTION == "diversifying",]
      selection_tests[i,14] <- sum(B_B_res_2_directional$FDR.05)
      #DELETING BAYESCAN FOLDER
      unlink("radiator_bayescan_*",recursive = T)
    }
    # After all the files have been analysed the working directory is set back
    # to the original directory
    setwd(working_dir)
  }
  ##### Whitlock and Lotterhos' method (W_L_method or OUTFLANK)#####
  if (W_L_test==TRUE){
    message("OUTFLANK")
    files_W_L <- dir(path.folder_sim,pattern = "^outflankGRAL")
    W_L_files <- lapply(paste0(path.folder_sim,"/",files_W_L), read.table,header=T,sep=",",skip=variables_number )
    # nrow(variables_file))
    outflank_pop_history <- read.table(paste0(path.folder_sim,"/",files_W_L[1]),header = F,sep=",",nrows = variables_number,colClasses = "character")
    # nrow(variables_file))
    if(simulation_type=="fly"){
           outflank_pop_history_b <- outflank_pop_history[24,2]
      outflank_pop_history_c <-paste0("GRAL_","B",outflank_pop_history_b[1],collapse = "")
      files_W_L_neutral <- dir(paste0(getwd(),"/","neutral_simulations"),pattern = paste0("outflank_",outflank_pop_history_c,"\\b"))
      }
    if(simulation_type=="general"){
      outflank_pop_history_b <- outflank_pop_history[c(6,14),2]
      outflank_pop_history_c <-paste0("GRAL_","B",outflank_pop_history_b[1],"_C",outflank_pop_history_b[2],collapse = "")
      files_W_L_neutral <- dir(paste0(getwd(),"/","neutral_simulations"),pattern = paste0("outflank_",outflank_pop_history_c,"\\b"))
    }
    outflank_neutral <- read.table(paste0(getwd(),"/","neutral_simulations","/",files_W_L_neutral),header=T,sep=",")
    Fstbar_neut <- outflank_neutral[1,1]
    df_infer <- outflank_neutral[2,1]
    outflank_msats <- as.numeric(outflank_pop_history[11,2])
    if(outflank_msats>2){
      He_min <- 0.2
    }else{
      He_min <- 0.1
    }
   
    
    #  # merging replicates in groups of 10
    #       if (length(W_L_files) %% 10 != 0){
    #   files_to_exclude <- length(W_L_files) %% 10
    #   W_L_files <- W_L_files[1:(length(W_L_files)-files_to_exclude)]
    # }
    # seq_start <- seq(1,length(W_L_files),10)
    # seq_end <- seq(10,length(W_L_files),10)
    # 
    # W_L_files_2 <- list()
    # for(i in 1:length(seq_start) ){
    #   W_L_files_2_temp <- rbindlist(W_L_files[seq_start[i]:seq_end[i]],fill=TRUE)
    #   W_L_files_2[[i]] <- W_L_files_2_temp
    # }

    
    W_L_files_2 <- W_L_files
    for(i in 1:length(W_L_files_2)){
      outflank_input <- W_L_files_2[[i]]
      outflank_input$FSTNoCorr[outflank_input$FSTNoCorr==1] <- NA
      outflank_input <- outflank_input[complete.cases(outflank_input$FSTNoCorr),]
      outflank_input_a <- outflank_input[complete.cases(outflank_input$FST),]
      outflank_input_b <- cbind(outflank_input_a,qvalues=NA, OutlierFlag=NA)
      W_L_res_a <- OutFLANK(outflank_input_a,LeftTrimFraction = 0.05, RightTrimFraction = 0.05,Hmin = He_min, NumberOfSamples=2, qthreshold = 0.05)
      W_L_res_b <- pOutlierFinderChiSqNoCorr(DataList=outflank_input_b, Fstbar=Fstbar_neut, dfInferred=df_infer, qthreshold = 0.05,Hmin = 0.1)

      if(is.list(W_L_res_a)){selection_tests[i,7] <- length(which(W_L_res_a$OutlierFlag_left==T))
      }else{selection_tests[i,7] <- NA}
      if(is.list(W_L_res_b)){selection_tests[i,8] <- length(which(W_L_res_b$OutlierFlag_left==T))
      }else{selection_tests[i,8] <- NA}

      if(is.list(W_L_res_a)){selection_tests[i,15] <- length(which(W_L_res_a$OutlierFlag==T))
      }else{selection_tests[i,15] <- NA}
      if(is.list(W_L_res_b)){selection_tests[i,16] <- length(which(W_L_res_b$OutlierFlag==T))
      }else{selection_tests[i,16] <- NA}

    }
  }
  # This is the final dataframe
  # final_selection_tests <- colMeans(selection_tests,na.rm = T)
  # write.table(fdist_pop_history,file = paste0(getwd(),"/","final_selection_tests_",variables,".csv"),sep = ",",row.names = F,col.names = F)
  # suppressWarnings(write.table(final_selection_tests,file = paste0(getwd(),"/","final_selection_tests_",variables,".csv"),sep = ",",row.names = T,append = T))
}


