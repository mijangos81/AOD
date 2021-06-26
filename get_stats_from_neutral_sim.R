


folders <- paste0(getwd(),"/",dir(path=getwd(),pattern="^SIM"))
variables_number <- 36
number_loci <- 231
number_of_stats_calculated <- 25

final_res_2 <- as.data.frame(matrix(ncol = 701))

for(fol in 1:length(folders)){
path.folder_stats_average <- folders[fol]    

files_ave <-paste0(path.folder_stats_average,"/",dir(path.folder_stats_average,pattern = "generations"))
pop_history <- read.table(files_ave[1],header = F,sep=",",nrows = variables_number, colClasses = "character")
simulation_type <- pop_history[9,2]

if(simulation_type=="1"){
  number_of_generations <- 12
  dis <- "high"
  }
if(simulation_type=="2"){
  number_of_generations <- 26
  dis <- "mod"
  }
if(simulation_type=="8"){
  number_of_generations <- 34
  dis <- "low"
  }

number_of_lines <- number_of_stats_calculated * number_of_generations
initial_breaks <- seq(1,number_of_lines,number_of_generations)
end_breaks <- c((initial_breaks-1)[-1],number_of_lines)

final_res <- as.data.frame(matrix(nrow = 6))
final_res[,1] <- dis 
rownames(final_res) <- c("exp_het_pop1","exp_het_pop2","het_pop1","het_pop2","Fst","Fstp")

for(i in 1:length(files_ave)){
 
# pop_history <- read.table(files_ave[i],header = F,sep=",",nrows = variables_number, colClasses = "character")
# vars <- pop_history[c(9,10),2]
# vars_2 <- paste0("size_",vars[1],"_c_",vars[2])

# vars <- pop_history[c(3,19,31,24),2]
# vars_2 <- paste0("adap_",vars[1],"_sel_",vars[2],"_exp_",vars[3],"_mut_",vars[4])

mean_df <- read.table(files_ave[i],header=F,sep=",",skip= variables_number,fill=T)

# mean_AFD <- mean_df[initial_breaks[1]:end_breaks[1],1:number_loci]
# mean_sha_diff <- mean_df[initial_breaks[2]:end_breaks[2],1:number_loci]
# mean_one_D_alpha_pop1  <-  mean_df[initial_breaks[3]:end_breaks[3],1:number_loci]
# mean_one_D_alpha_pop2 <-  mean_df[initial_breaks[4]:end_breaks[4],1:number_loci]
# mean_two_D_alpha_pop1 <-  mean_df[initial_breaks[5]:end_breaks[5],1:number_loci]
# mean_two_D_alpha_pop2 <-  mean_df[initial_breaks[6]:end_breaks[6],1:number_loci]
# mean_one_D_beta <-  mean_df[initial_breaks[7]:end_breaks[7],1:number_loci]
# mean_two_D_beta <-  mean_df[initial_breaks[8]:end_breaks[8],1:number_loci]
# mean_sha_pop1 <-  mean_df[initial_breaks[9]:end_breaks[9],1:number_loci]
# mean_sha_pop2 <-  mean_df[initial_breaks[10]:end_breaks[10],1:number_loci]
# mean_pops <- mean_df[initial_breaks[11]:end_breaks[11],1:number_loci]
# mean_shua <-  mean_df[initial_breaks[12]:end_breaks[12],1:number_loci]
mean_expected_het_pop1 <-  mean_df[end_breaks[13],c(30,50,84,126,151,180,208)]
mean_expected_het_pop2 <-  mean_df[end_breaks[14],c(30,50,84,126,151,180,208)]
mean_het_pop1 <-  mean_df[end_breaks[15],c(30,50,84,126,151,180,208)]
mean_het_pop2 <-  mean_df[end_breaks[16],c(30,50,84,126,151,180,208)]
# mean_Ht <-  mean_df[initial_breaks[17]:end_breaks[17],1:number_loci]
mean_Fst <-  mean_df[end_breaks[18],c(30,50,84,126,151,180,208)]
mean_Fstp <-  mean_df[end_breaks[19],c(30,50,84,126,151,180,208)]
# mean_Dest <-  mean_df[initial_breaks[20]:end_breaks[20],1:number_loci]
# mean_deleterious_eliminated <-  mean_df[initial_breaks[21]:end_breaks[21],1]
# mean_mean_s_final <-  mean_df[initial_breaks[22]:end_breaks[22],1]
# mean_mean_h_final <-  mean_df[initial_breaks[23]:end_breaks[23],1]
# mean_g_load_a_final <-  mean_df[initial_breaks[24]:end_breaks[24],1]
# mean_g_load_m_final <-  mean_df[initial_breaks[25]:end_breaks[25],1]

res_temp <- rbind(mean_expected_het_pop1,mean_expected_het_pop2,mean_het_pop1,mean_het_pop2,mean_Fst,mean_Fstp)

# final_res[i,1] <- vars_2
# final_res[i,2] <- rowMeans(mean_AFD, na.rm = T)[number_of_generations]
# final_res[i,3] <- rowMeans(mean_sha_diff, na.rm = T)[number_of_generations]
# final_res[i,4] <- rowMeans(mean_one_D_alpha_pop1, na.rm = T)[number_of_generations]
# final_res[i,5] <- rowMeans(mean_one_D_alpha_pop2 , na.rm = T)[number_of_generations]
# final_res[i,6] <- rowMeans(mean_two_D_alpha_pop1, na.rm = T)[number_of_generations]
# final_res[i,7] <- rowMeans(mean_two_D_alpha_pop2 , na.rm = T)[number_of_generations]
# final_res[i,8] <- rowMeans(mean_one_D_beta , na.rm = T)[number_of_generations]
# final_res[i,9] <- rowMeans(mean_two_D_beta , na.rm = T)[number_of_generations]
# final_res[i,10] <- rowMeans(mean_sha_pop1, na.rm = T)[number_of_generations]
# final_res[i,11] <- rowMeans(mean_sha_pop2 , na.rm = T)[number_of_generations]
# final_res[i,12] <- rowMeans(mean_pops , na.rm = T)[number_of_generations]
# final_res[i,13] <- rowMeans(mean_shua , na.rm = T)[number_of_generations]
# final_res[i,14] <- rowMeans(mean_expected_het_pop1, na.rm = T)[number_of_generations]
# final_res[i,15] <- rowMeans(mean_expected_het_pop2, na.rm = T)[number_of_generations]
# final_res[i,16] <- rowMeans(mean_het_pop1 , na.rm = T)[number_of_generations]
# final_res[i,17] <- rowMeans(mean_het_pop2, na.rm = T)[number_of_generations]
# final_res[i,18] <- rowMeans(mean_Ht , na.rm = T)[number_of_generations]
# final_res[i,19] <- rowMeans(mean_Fst , na.rm = T)[number_of_generations]
# final_res[i,20] <- rowMeans(mean_Fstp , na.rm = T)[number_of_generations]
# final_res[i,21] <- rowMeans(mean_Dest , na.rm = T)[number_of_generations]
# final_res[i,22] <- mean_deleterious_eliminated[number_of_generations]
# final_res[i,23] <- mean_mean_s_final[number_of_generations]
# final_res[i,24] <- mean_mean_h_final[number_of_generations]
# final_res[i,25] <- mean_g_load_a_final[number_of_generations]
# final_res[i,26] <- mean_g_load_m_final[number_of_generations]

final_res <- cbind(final_res,res_temp)
}
names(final_res_2) <- names(final_res)
final_res_2 <- rbind(final_res_2,final_res)
}

# colnames(final_res) <- c(
#   "vars",
#   "AFD",
#     "sha_diff",
#     "one_D_alpha_pop1",
#     "one_D_alpha_pop2",
#     "two_D_alpha_pop1",
#     "two_D_alpha_pop2",
#     "one_D_beta",
#     "two_D_beta",
#     "sha_pop1",
#     "sha_pop2",
#     "sha_pops",
#     "shua",
#     "expected_het_pop1",
#     "expected_het_pop2",
#     "het_pop1",
#     "het_pop2",
#     "Ht",
#     "Fst",
#     "Fstp",
#     "Dest",
#     "deleterious_eliminated",
#     "mean_s_final",
#     "mean_h_final",
#     "g_load_a_final", 
#     "g_load_m_final")
# 
 write.csv(final_res_2,file = "neutral_sim_fig_1_snps.csv")
