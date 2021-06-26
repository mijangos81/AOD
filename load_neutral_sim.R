# 
# folders <- dir(path=getwd(),pattern="^final_stats")
# path.folder_sim <- folders[1]    
# files_ave <- paste0(path.folder_sim,"/",dir(path.folder_sim,pattern = "final_stats_average"))

files_ave_exp <- paste0(getwd(),"/",dir(pattern = "neutral_expectation"))

number_loci <- length(neutral_loci_location)
number_of_stats_calculated <- length(by_gen_lists_a_b)
number_of_generations <- gen_number_dispersal
number_of_lines <- number_of_stats_calculated * number_of_generations
initial_breaks <- seq(1,number_of_lines,number_of_generations)
end_breaks <- c((initial_breaks-1)[-1],number_of_lines)

final_res_exp <- as.data.frame(matrix(nrow = length(files_ave_exp) , ncol = number_of_stats_calculated +1))

for(i in 1:length(files_ave_exp)){
  
  pop_history_exp <- read.table(files_ave_exp[i],header = F,sep=",",nrows = 36, colClasses = "character")
  vars_exp <- pop_history_exp[c(33),2]
  vars_2_exp <- paste0("offspring_variance_",vars_exp[1])
  
  # vars <- pop_history[c(3,19,31,24),2]
  # vars_2 <- paste0("adap_",vars[1],"_sel_",vars[2],"_exp_",vars[3],"_mut_",vars[4])
  
  mean_df_exp <- read.table(files_ave_exp[i],header=F,sep=",",skip= 36,fill=T)
  
  mean_AFD_exp <- mean_df_exp[initial_breaks[1]:end_breaks[1],1:number_loci]
  mean_sha_diff_exp <- mean_df_exp[initial_breaks[2]:end_breaks[2],1:number_loci]
  mean_one_D_alpha_pop1_exp  <-  mean_df_exp[initial_breaks[3]:end_breaks[3],1:number_loci]
  mean_one_D_alpha_pop2_exp <-  mean_df_exp[initial_breaks[4]:end_breaks[4],1:number_loci]
  mean_two_D_alpha_pop1_exp <-  mean_df_exp[initial_breaks[5]:end_breaks[5],1:number_loci]
  mean_two_D_alpha_pop2_exp <-  mean_df_exp[initial_breaks[6]:end_breaks[6],1:number_loci]
  mean_one_D_beta_exp <-  mean_df_exp[initial_breaks[7]:end_breaks[7],1:number_loci]
  mean_two_D_beta_exp <-  mean_df_exp[initial_breaks[8]:end_breaks[8],1:number_loci]
  mean_sha_pop1_exp <-  mean_df_exp[initial_breaks[9]:end_breaks[9],1:number_loci]
  mean_sha_pop2_exp <-  mean_df_exp[initial_breaks[10]:end_breaks[10],1:number_loci]
  mean_pops_exp <- mean_df_exp[initial_breaks[11]:end_breaks[11],1:number_loci]
  mean_shua_exp <-  mean_df_exp[initial_breaks[12]:end_breaks[12],1:number_loci]
  mean_expected_het_pop1_exp <-  mean_df_exp[initial_breaks[13]:end_breaks[13],1:number_loci]
  mean_expected_het_pop2_exp <-  mean_df_exp[initial_breaks[14]:end_breaks[14],1:number_loci]
  mean_het_pop1_exp <-  mean_df_exp[initial_breaks[15]:end_breaks[15],1:number_loci]
  mean_het_pop2_exp <-  mean_df_exp[initial_breaks[16]:end_breaks[16],1:number_loci]
  mean_Ht_exp <-  mean_df_exp[initial_breaks[17]:end_breaks[17],1:number_loci]
  mean_Fst_exp <-  mean_df_exp[initial_breaks[18]:end_breaks[18],1:number_loci]
  mean_Fstp_exp <-  mean_df_exp[initial_breaks[19]:end_breaks[19],1:number_loci]
  mean_Dest_exp <-  mean_df_exp[initial_breaks[20]:end_breaks[20],1:number_loci]
  mean_deleterious_eliminated_exp <-  mean_df_exp[initial_breaks[21]:end_breaks[21],1]
  mean_mean_s_final_exp <-  mean_df_exp[initial_breaks[22]:end_breaks[22],1]
  mean_mean_h_final_exp <-  mean_df_exp[initial_breaks[23]:end_breaks[23],1]
  mean_g_load_a_final_exp <-  mean_df_exp[initial_breaks[24]:end_breaks[24],1]
  mean_g_load_m_final_exp <-  mean_df_exp[initial_breaks[25]:end_breaks[25],1]
  
  final_res_exp[i,1] <- vars_2_exp
  final_res_exp[i,2] <- rowMeans(mean_AFD_exp, na.rm = T)[number_of_generations]
  final_res_exp[i,3] <- rowMeans(mean_sha_diff_exp, na.rm = T)[number_of_generations]
  final_res_exp[i,4] <- rowMeans(mean_one_D_alpha_pop1_exp, na.rm = T)[number_of_generations]
  final_res_exp[i,5] <- rowMeans(mean_one_D_alpha_pop2_exp , na.rm = T)[number_of_generations]
  final_res_exp[i,6] <- rowMeans(mean_two_D_alpha_pop1_exp, na.rm = T)[number_of_generations]
  final_res_exp[i,7] <- rowMeans(mean_two_D_alpha_pop2_exp , na.rm = T)[number_of_generations]
  final_res_exp[i,8] <- rowMeans(mean_one_D_beta_exp , na.rm = T)[number_of_generations]
  final_res_exp[i,9] <- rowMeans(mean_two_D_beta_exp , na.rm = T)[number_of_generations]
  final_res_exp[i,10] <- rowMeans(mean_sha_pop1_exp, na.rm = T)[number_of_generations]
  final_res_exp[i,11] <- rowMeans(mean_sha_pop2_exp , na.rm = T)[number_of_generations]
  final_res_exp[i,12] <- rowMeans(mean_pops_exp , na.rm = T)[number_of_generations]
  final_res_exp[i,13] <- rowMeans(mean_shua_exp , na.rm = T)[number_of_generations]
  final_res_exp[i,14] <- rowMeans(mean_expected_het_pop1_exp, na.rm = T)[number_of_generations]
  final_res_exp[i,15] <- rowMeans(mean_expected_het_pop2_exp, na.rm = T)[number_of_generations]
  final_res_exp[i,16] <- rowMeans(mean_het_pop1_exp , na.rm = T)[number_of_generations]
  final_res_exp[i,17] <- rowMeans(mean_het_pop2_exp, na.rm = T)[number_of_generations]
  final_res_exp[i,18] <- rowMeans(mean_Ht_exp , na.rm = T)[number_of_generations]
  final_res_exp[i,19] <- rowMeans(mean_Fst_exp , na.rm = T)[number_of_generations]
  final_res_exp[i,20] <- rowMeans(mean_Fstp_exp , na.rm = T)[number_of_generations]
  final_res_exp[i,21] <- rowMeans(mean_Dest_exp , na.rm = T)[number_of_generations]
  final_res_exp[i,22] <- mean_deleterious_eliminated_exp[number_of_generations]
  final_res_exp[i,23] <- mean_mean_s_final_exp[number_of_generations]
  final_res_exp[i,24] <- mean_mean_h_final_exp[number_of_generations]
  final_res_exp[i,25] <- mean_g_load_a_final_exp[number_of_generations]
  final_res_exp[i,26] <- mean_g_load_m_final_exp[number_of_generations]
}

colnames(final_res_exp) <- c(
  "vars_exp",
  "AFD_exp",
  "sha_diff_exp",
  "one_D_alpha_pop_exp1",
  "one_D_alpha_pop2_exp",
  "two_D_alpha_pop1_exp",
  "two_D_alpha_pop2_exp",
  "one_D_beta_exp",
  "two_D_beta_exp",
  "sha_pop1_exp",
  "sha_pop2_exp",
  "sha_pops_exp",
  "shua_exp",
  "expected_het_pop1_exp",
  "expected_het_pop2_exp",
  "het_pop1_exp",
  "het_pop2_exp",
  "Ht_exp",
  "Fst_exp",
  "Fstp_exp",
  "Dest_exp",
  "deleterious_eliminated_exp",
  "mean_s_final_exp",
  "mean_h_final_exp",
  "g_load_a_final_exp", 
  "g_load_m_final_exp")

#write.csv(final_res_exp,file = "final_res_exp_stats.csv")
