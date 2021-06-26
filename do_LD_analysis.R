source('functions.R')
library(data.table)
library(fields)
library(rlist)
library(plyr)

# ld_max_pairwise <- 2000000  # maximun distance, in basepairs, at which pairwise LD should be calculated
# ld_resolution <- 200000 # resolution, in basepairs, at which LD should be measured
# region_size <- 200000 # the size of the window at which the LD statistics and number of NS are calculated
variables_number <- 36
# variables_number <- 37
folders <- dir(path=getwd(),pattern="^SIM")
# folders <- dir(path=getwd(),pattern="^neutral_chill")
# folders <- dir(path=getwd(),pattern="^AOD_fly")
# folders <- dir(path=getwd(),pattern="^first_sim")

# folders <- getwd()

for (i in 1:length(folders)){
path.folder_sim <- folders[i]    

# INPUT FILES
# LD analyses are done with offspring and not with parents and using pop1.   
# over 100 replications, LD patterns in pop1 and pop2 are the same.

#the following code lines are to bind together the LD matrices with their corresponding snps files because in some replicates either the matrix or the snps file is not created
files_LD_pop1 <- dir(path.folder_sim,pattern = "ld_pop1")
files_length_LD <- length(files_LD_pop1)

myvec_LD <- files_LD_pop1
myvec_LD_2 <- do.call(paste, c(read.table(text = myvec_LD, sep = "_")[8], sep = "_"))
myvec_LD_3 <- as.numeric(gsub(".csv","",myvec_LD_2))
myvec_LD_4 <- as.data.frame(cbind(1:files_length_LD,myvec_LD_3))
# myvec_LD_4 <- myvec_LD_4[order(myvec_LD_4$myvec_LD_3),]

files_snps_pop1 <- dir(path.folder_sim,pattern = "snps_pop1")
files_length_snps <- length(files_snps_pop1)

myvec_SNPS <- files_snps_pop1
myvec_SNPS_2 <- do.call(paste, c(read.table(text = myvec_SNPS, sep = "_")[8], sep = "_"))
myvec_SNPS_3 <- as.numeric(gsub(".csv","",myvec_SNPS_2))
myvec_SNPS_4 <- as.data.frame(cbind(1:files_length_snps,myvec_SNPS_3))
# myvec_SNPS_4 <- myvec_SNPS_4[order(myvec_SNPS_4$myvec_SNPS_3),]

if (files_length_LD<files_length_snps) {
  files_for_LD <- myvec_SNPS_4[myvec_SNPS_4$myvec_SNPS_3 %in% myvec_LD_4$myvec_LD_3,] 
  files_for_LD_2 <- cbind(files_for_LD,myvec_LD_4)
}
if (files_length_LD>files_length_snps) {
  files_for_LD <- myvec_LD_4[myvec_LD_4$myvec_LD_3 %in% myvec_SNPS_4$myvec_SNPS_3,] 
  files_for_LD_2 <- cbind(files_for_LD,myvec_SNPS_4)
  
}
if (files_length_LD==files_length_snps) {
  files_for_LD <- myvec_SNPS_4[myvec_SNPS_4$myvec_SNPS_3 %in% myvec_LD_4$myvec_LD_3,] 
  files_for_LD_2 <- cbind(files_for_LD,myvec_LD_4)
}

 # files_for_LD_2 <- files_for_LD_2[sample(1:nrow(files_for_LD_2),5),]

pop_history <- read.table(paste0(path.folder_sim,"/",files_LD_pop1[1]),header = F,sep=",",nrows = variables_number, colClasses = "character")
simulation_type <- pop_history[1,2]

# VARIABLES
# this is the width of the bins at which LD is measured for the pairwise LD plot 
ld_bins <- 1000
 minor <- 0.05
 # ld_max_pairwise <- 2000000  # maximun distance, in basepairs, at which pairwise LD should be calculated
 # ld_resolution <- 200000 # resolution, in basepairs, at which LD should be measured
 # region_size <- 200000 # the size of the window at which the LD statistics and number of NS are calculated
if(simulation_type == "general"){
chromosome_length <- as.numeric(pop_history[15,2]) * as.numeric(pop_history[34,2])# windows_gral * map_resolution
}
if(simulation_type == "fly"){
chromosome_length <- 23200000
  # chr_length
}

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
       "_B_", pop_history[21,2], # rate,
      "_C_", pop_history[23,2], # gamma_scale,
      "_D_", pop_history[24,2], # gamma_shape,
      "_F_", pop_history[18,2] #  mutation_rate,
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

# DATAFRAMES
link_iterations <- rep(list(as.data.frame(matrix(ncol=length(files_LD_pop1),nrow = chromosome_length/region_size))),ld_max_pairwise/ld_resolution)
number_nonsyn <- as.data.frame(matrix(ncol = length(files_LD_pop1), nrow = chromosome_length/region_size))
snpsum.col_mean <- rep(list(as.data.frame(matrix(ncol=20,nrow = chromosome_length/region_size))),length(files_snps_pop1))

window_size <- as.character(seq(ld_resolution,ld_max_pairwise,ld_resolution)/1000)
window_size <- paste0(window_size,"Kbp")
region_size_rownames <- as.character(seq(region_size,chromosome_length,region_size))

LD_final_pop1 <- do_LD_analysis(path_folder=path.folder_sim, list_ld=files_LD_pop1, list_snps= files_snps_pop1, list_files=files_for_LD_2)
link_iterations <- LD_final_pop1[[1]]
number_nonsyn <- LD_final_pop1[[2]]
snpsum.col_mean <- LD_final_pop1[[3]]

if (length(files_LD_pop1)>2) {
snpsum.col_mean_2 <- aaply(laply(snpsum.col_mean, as.matrix), c(2, 3), mean,na.rm=T)
snpsum.col_sd <- aaply(laply(snpsum.col_mean, as.matrix), c(2, 3), sd,na.rm=T)
rownames(snpsum.col_mean_2) <- region_size_rownames
rownames(snpsum.col_sd) <- region_size_rownames
}
if (length(files_LD_pop1)<2) {
snpsum.col_mean_2 <- as.data.frame(snpsum.col_mean)
snpsum.col_sd <- as.data.frame(snpsum.col_mean)
rownames(snpsum.col_mean_2) <- region_size_rownames
rownames(snpsum.col_sd) <- region_size_rownames
}

link_final <-  as.data.frame(sapply(link_iterations,rowMeans,na.rm=T))
colnames(link_final) <- window_size
rownames(link_final) <- region_size_rownames
link_col <- round(colMeans(link_final,na.rm = T),3)

number_nonsyn_b <- rowMeans(number_nonsyn, na.rm = T)
LD_final_df <- cbind(link_final,number_nonsyn_b) 

final_pairwise_ld <- as.data.frame(matrix(nrow = ld_max_pairwise/ld_bins,ncol = length(files_LD_pop1)))
break_bins <- c(seq(1,ld_max_pairwise,ld_bins),ld_max_pairwise)

for(i in 1:length(files_LD_pop1)){
ld_columns_b <- read.table(paste0(path.folder_sim,"/",files_LD_pop1[i]),sep = ",",skip = variables_number,row.names = 1)
# ld_columns_b <- read.table(files_LD_pop1[i],sep = ",",skip = variables_number,row.names = 1)
colnames(ld_columns_b) <- rownames(ld_columns_b)
ld_columns_b <- as.data.frame(as.table(as.matrix(ld_columns_b)))
ld_columns_b <- ld_columns_b[-ld_columns_b$Freq < 0,] #remove cases where LD was not calculated
ld_columns_b$Var1 <- as.numeric(as.character(ld_columns_b$Var1))
ld_columns_b$Var2 <- as.numeric(as.character(ld_columns_b$Var2))
#determine the distance at which LD was calculated
ld_columns_b$dis <- ld_columns_b$Var2 - ld_columns_b$Var1 
ld_columns_c <- ld_columns_b[order(ld_columns_b$dis),]
ld_columns_d <- ld_columns_c[which(ld_columns_c$dis<ld_max_pairwise),]
bins_ld <- stats.bin(ld_columns_d$dis,ld_columns_d$Freq,breaks = break_bins)
final_pairwise_ld[,i] <- unname(bins_ld$stats[2,])
}

# This is the final dataframe for LD per locus
LD_final_df <- cbind(rownames(link_final),link_final,number_nonsyn_b)
write.table(pop_history,file = paste0(getwd(),"/","final_LD_locus_",variables,".csv"),sep = ",",row.names = F,col.names = F)
suppressWarnings(write.table(LD_final_df,file = paste0(getwd(),"/","final_LD_locus_",variables,".csv"),sep = ",",row.names = F,append = T))

# This is the final dataframe for the SNP's stats average
write.table(pop_history,file = paste0(getwd(),"/","final_snps_stats_ave_",variables,".csv"),sep = ",",row.names = F,col.names = F)
suppressWarnings(write.table(snpsum.col_mean_2,file = paste0(getwd(),"/","final_snps_stats_ave_",variables,".csv"),sep = ",",row.names = F,append = T))

# This is the final dataframe for the SNP's stats standard deviation
write.table(pop_history,file = paste0(getwd(),"/","final_snps_stats_sd_",variables,".csv"),sep = ",",row.names = F,col.names = F)
suppressWarnings(write.table(snpsum.col_sd,file = paste0(getwd(),"/","final_snps_stats_sd_",variables,".csv"),sep = ",",row.names = F,append = T))

# This is the final dataframe for LD pairwise
final_pairwise_ld <- rowMeans(final_pairwise_ld,na.rm = T)
final_pairwise_ld <- as.data.frame(cbind(break_bins[-1],final_pairwise_ld))
colnames(final_pairwise_ld) <- c("distance","rsqr")
final_pairwise_ld <- final_pairwise_ld[complete.cases(final_pairwise_ld$rsqr),]
write.table(pop_history,file = paste0(getwd(),"/","final_LD_pairwise_",variables,".csv"),sep = ",",row.names = F,col.names = F)
suppressWarnings(write.table(final_pairwise_ld,file = paste0(getwd(),"/","final_LD_pairwise_",variables,".csv"),sep = ",",row.names = F,append = T))
}

