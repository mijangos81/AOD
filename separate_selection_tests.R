
# folders <- dir(path=getwd(),pattern="^final_selection")
# path.folder_sim <- folders[1]    
variables_number <-36

 path.folder_sim <- paste0(getwd(),"/","gral_sim_res_copy/")
files_ave <- paste0(path.folder_sim,"/",dir(path.folder_sim,pattern = "^final_selection_tests.*csv$"))

final_res <- as.data.frame(matrix(nrow = length(files_ave) , ncol = 16 +1))

for(i in 1:length(files_ave)){
 
pop_history <- read.table(files_ave[i],header = F,sep=",",nrows = variables_number, colClasses = "character")
vars <- pop_history[c(6,14,26,27,18),2]
vars_2 <- paste0("Ne_",vars[1],"_c_",vars[2],"_h_",vars[3],"_s_",vars[4],"_sel_",vars[5])
#vars <- pop_history[c(14,26,27,28),2]
#vars_2 <- paste0("c_",vars[1],"_h_",vars[2],"_s_",vars[3],"_q_",vars[4],"_d_",vars[5])
selection_df <- read.table(files_ave[i],header=F,sep=",",skip= variables_number+1,fill=T)
tests <- as.character(selection_df$V1)
selection_df <- as.numeric(t(selection_df)[2,])
final_res[i,1] <- vars_2
final_res[i,2:17] <- selection_df
}

colnames(final_res) <- c( 
  "vars",
    tests
#   "FDIST2_FST",
# "FDIST2_Shua",
# "FDIST2_FST_Sha_diff",
# "FDIST2_Sha_pop",
# "BAYESCAN_A",
# "BAYESCAN_B",
# "OUTFLANK_A",
# "OUTFLANK_B"
  )

write.csv(final_res,file = "final_res_selection.csv")
