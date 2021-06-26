path.folder_sim <- getwd()
files_ave <- paste0(path.folder_sim,"/",dir(path.folder_sim,pattern = "^final_stats_sd"))
get_number_loci <-  read.table(files_ave[1],header=F,sep=",",skip= variables_number,fill=T)
number_loci <- ncol(get_number_loci)
variables_number <- 36
number_of_stats_calculated <- 25
get_number_generations <- read.table(files_ave[i],header = F,sep=",",nrows = variables_number, colClasses = "character")
number_of_generations <- as.numeric(get_number_generations[7,2])
number_of_lines <- number_of_stats_calculated * number_of_generations
initial_breaks <- seq(1,number_of_lines,number_of_generations)
end_breaks <- c((initial_breaks-1)[-1],number_of_lines)
final_res <- as.data.frame(matrix(nrow = length(files_ave) , ncol = number_of_stats_calculated +1))

for(i in 1:length(files_ave)){
pop_history <- read.table(files_ave[i],header = F,sep=",",nrows = variables_number, colClasses = "character")
vars <- pop_history[c(18,21,23,24),2]
vars_2 <- paste0("selection_",vars[1],"_logmean_",vars[2],"_dominance_",vars[3],"_mut_",vars[4])
sd_df <- read.table(files_ave[i],header=F,sep=",",skip= variables_number,fill=T)
sd_df[1,2]
sd_AFD <- sd_df[initial_breaks[1]:end_breaks[1],1:number_loci]
sd_sha_diff <- sd_df[initial_breaks[2]:end_breaks[2],1:number_loci]
sd_one_D_alpha_pop1  <-  sd_df[initial_breaks[3]:end_breaks[3],1:number_loci]
sd_one_D_alpha_pop2 <-  sd_df[initial_breaks[4]:end_breaks[4],1:number_loci]
sd_two_D_alpha_pop1 <-  sd_df[initial_breaks[5]:end_breaks[5],1:number_loci]
sd_two_D_alpha_pop2 <-  sd_df[initial_breaks[6]:end_breaks[6],1:number_loci]
sd_one_D_beta <-  sd_df[initial_breaks[7]:end_breaks[7],1:number_loci]
sd_two_D_beta <-  sd_df[initial_breaks[8]:end_breaks[8],1:number_loci]
sd_sha_pop1 <-  sd_df[initial_breaks[9]:end_breaks[9],1:number_loci]
sd_sha_pop2 <-  sd_df[initial_breaks[10]:end_breaks[10],1:number_loci]
sd_pops <- sd_df[initial_breaks[11]:end_breaks[11],1:number_loci]
sd_shua <-  sd_df[initial_breaks[12]:end_breaks[12],1:number_loci]
sd_expected_het_pop1 <-  sd_df[initial_breaks[13]:end_breaks[13],1:number_loci]
sd_expected_het_pop2 <-  sd_df[initial_breaks[14]:end_breaks[14],1:number_loci]
sd_het_pop1 <-  sd_df[initial_breaks[15]:end_breaks[15],1:number_loci]
sd_het_pop2 <-  sd_df[initial_breaks[16]:end_breaks[16],1:number_loci]
sd_Ht <-  sd_df[initial_breaks[17]:end_breaks[17],1:number_loci]
sd_Fst <-  sd_df[initial_breaks[18]:end_breaks[18],1:number_loci]
sd_Fstp <-  sd_df[initial_breaks[19]:end_breaks[19],1:number_loci]
sd_Dest <-  sd_df[initial_breaks[20]:end_breaks[20],1:number_loci]
sd_deleterious_eliminated <-  sd_df[initial_breaks[21]:end_breaks[21],1]
sd_sd_s_final <-  sd_df[initial_breaks[22]:end_breaks[22],1]
sd_sd_h_final <-  sd_df[initial_breaks[23]:end_breaks[23],1]
sd_g_load_a_final <-  sd_df[initial_breaks[24]:end_breaks[24],1]
sd_g_load_m_final <-  sd_df[initial_breaks[25]:end_breaks[25],1]

final_res[i,1] <- vars_2
final_res[i,2] <- rowMeans(sd_AFD, na.rm = T)[number_of_generations]
final_res[i,3] <- rowMeans(sd_sha_diff, na.rm = T)[number_of_generations]
final_res[i,4] <- rowMeans(sd_one_D_alpha_pop1, na.rm = T)[number_of_generations]
final_res[i,5] <- rowMeans(sd_one_D_alpha_pop2 , na.rm = T)[number_of_generations]
final_res[i,6] <- rowMeans(sd_two_D_alpha_pop1, na.rm = T)[number_of_generations]
final_res[i,7] <- rowMeans(sd_two_D_alpha_pop2 , na.rm = T)[number_of_generations]
final_res[i,8] <- rowMeans(sd_one_D_beta , na.rm = T)[number_of_generations]
final_res[i,9] <- rowMeans(sd_two_D_beta , na.rm = T)[number_of_generations]
final_res[i,10] <- rowMeans(sd_sha_pop1, na.rm = T)[number_of_generations]
final_res[i,11] <- rowMeans(sd_sha_pop2 , na.rm = T)[number_of_generations]
final_res[i,12] <- rowMeans(sd_pops , na.rm = T)[number_of_generations]
final_res[i,13] <- rowMeans(sd_shua , na.rm = T)[number_of_generations]
final_res[i,14] <- rowMeans(sd_expected_het_pop1, na.rm = T)[number_of_generations]
final_res[i,15] <- rowMeans(sd_expected_het_pop2, na.rm = T)[number_of_generations]
final_res[i,16] <- rowMeans(sd_het_pop1 , na.rm = T)[number_of_generations]
final_res[i,17] <- rowMeans(sd_het_pop2, na.rm = T)[number_of_generations]
final_res[i,18] <- rowMeans(sd_Ht , na.rm = T)[number_of_generations]
final_res[i,19] <- rowMeans(sd_Fst , na.rm = T)[number_of_generations]
final_res[i,20] <- rowMeans(sd_Fstp , na.rm = T)[number_of_generations]
final_res[i,21] <- rowMeans(sd_Dest , na.rm = T)[number_of_generations]
final_res[i,22] <- sd_deleterious_eliminated[number_of_generations]
final_res[i,23] <- sd_sd_s_final[number_of_generations]
final_res[i,24] <- sd_sd_h_final[number_of_generations]
final_res[i,25] <- sd_g_load_a_final[number_of_generations]
final_res[i,26] <- sd_g_load_m_final[number_of_generations]
}

colnames(final_res) <- c(
  "vars",
  "AFD",
    "sha_diff",
    "one_D_alpha_pop1",
    "one_D_alpha_pop2",
    "two_D_alpha_pop1",
    "two_D_alpha_pop2",
    "one_D_beta",
    "two_D_beta",
    "sha_pop1",
    "sha_pop2",
    "sha_pops",
    "shua",
    "expected_het_pop1",
    "expected_het_pop2",
    "het_pop1",
    "het_pop2",
    "Ht",
    "Fst",
    "Fstp",
    "Dest",
    "deleterious_eliminated",
    "mean_s_final",
    "mean_h_final",
    "g_load_a_final", 
    "g_load_m_final")

write.csv(final_res,file = "final_res_stats_sd.csv")
