####### LOADING LIBRARIES #######
library(patchwork)
library(broom)
library(readxl)
library(reshape2)
library(readr)
library(bigsnpr)
library(fields)
library(snpStats)
library(ggthemes)
library(ggplot2)
library(plyr)
library(scales)
library(raster)
library(rasterVis)
library(viridis)
library(rgeos)
library(DescTools)
library(maptools)
library(zoo)
library(parallel)
library(ggnewscale)
source('gHap.R')
source('rotate_matrix_2.R')
source('functions.R')
#############################################################
################ VARIABLES ###############################
#############################################################
# maf_var <- 0.05 # Minimum Minor Allele Frequency (MAF) for a SNP to be kept
# geno_var <- 0.05 # Maximum proportion of missing values for a SNP to be kept
# mind_var <- 0.05 # Maximum proportion of missing values for a sample to be kept
# hwe_var <- 1e-20 # Filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold. 
chrom <-  "2L" # these are the chromosomes to analyse
map_resolution <- 100000 # this is the resolution of the recombination map
# these are the distance taht will be tested in the regression analyses
dist_cM <- seq(1,10,1) # distances in cM
dist_bp <- c(seq(200000,2000000,200000)) # distance in bp
# dist_cM <- seq(1,10,1) # distances in cM
# dist_bp <- c(seq(1000000,10000000,1000000)) # distance in bp
# PAIRWISE LD
ld_bins <- 1000 # this is the length , in base pairs, of the bins 
ld_max_pairwise <- 2000000  # maximun distance, in basepairs, at which pairwise LD should be calculated
ld_resolution <- 10000 # resolution, in basepairs, at which LD should be measured
resolution_regression_analyses <- 1000000
haploview_plot <- F
ld_pairwise_plot <- T
# recombination_plot
# non_synonymous
# r_squared_R1
# pairwise
#############################################################
################ LOADING DATA ###############################
#############################################################
# LODING DATA WITH PLINK
# snp_plinkQC(plink.path= "/Users/mijangos/Dropbox/PhD/R/parentage_platypus/plink", prefix.in= "/Users/mijangos/Dropbox/PhD/R/chillingham_data/Chillingham", file.type = "--file", prefix.out = "res_chill", maf = maf_var, geno = geno_var, mind = mind_var, hwe = hwe_var,extra.options = "--chr-set 29")
# chill_data_plink <- read.plink(bed= "res_chill.bed" , bim="res_chill.bim" , fam="res_chill.fam" )
# # the function read.plink does not read the map file properly, so it is attached separately
# map_plink <- read_table2("res_chill.bim", col_names = FALSE)
# colnames(map_plink) <- c("chromosome", "snp.name", "cM", "position", "allele.1", "allele.2")
# chill_data_plink$map <- map_plink
# genotype_plink <- chill_data_plink$genotypes # this is the input for pairwise LD analyses
# colnames(genotype_plink) <- chill_data_plink$map$position
# snpsum.col_plink<- col.summary(genotype_plink)
# chill_df_plink <- as.data.frame(cbind(snpsum.col_plink,chill_data_plink$map))

# LOADING DATA WITH read.pedfile
# chill_data_ped <- read.pedfile_b(paste0(getwd(),"/","Chillingham.ped"),sep = " ",show_warnings=T)
# map_ped <- read_table2("Chillingham.map", col_names = FALSE,cols(X1 = "c", X2 = "c", X3 = "d",X4="d"))
# colnames(map_ped) <- c("chromosome","snp.name","cM","position")
# chill_data_ped$map <- map_ped
# genotype_ped_filter <- chill_data_ped$genotypes
# # colnames(genotype_ped) <- chill_data_ped$map$position
# snpsum.col_ped_filter <- col.summary(genotype_ped_filter)
# filter_ped <- which(snpsum.col_ped_filter$MAF > 0.05 & snpsum.col_ped_filter$Call.rate > 0.9 & snpsum.col_ped_filter$Calls > 10)
# 
# genotype_chill_temp <- chill_data_ped$genotypes
# genotype_chill_temp <- genotype_chill_temp[,filter_ped]
# 
# map_chill_temp <- chill_data_ped$map
# map_chill_temp <- map_chill_temp[filter_ped,]
# 
# chill_data_ped_temp <- chill_data_ped
# chill_data_ped_temp$genotypes <- genotype_chill_temp
# chill_data_ped_temp$map <- map_chill_temp
#   
# genotype_ped <- chill_data_ped_temp$genotypes
# colnames(genotype_ped) <- chill_data_ped_temp$map$position
# snpsum.col_ped <- col.summary(genotype_ped)
# chill_df_ped <- as.data.frame(cbind(snpsum.col_ped,chill_data_ped_temp$map))

chill_df <- snp_final
# chill_df <- read.table("snps_pop1FLY_B_0.005_C_0.2_D_5e-05_F_TRUE_rep_1.csv",header=T,sep = ",",skip = variables_number)
chill_data <- genotype_pop1
genotype_chill <- genotype
  
# these are the recombination maps of all the chromosomes obtained from Ma, Li et al. (2016), Data from: Cattle sex-specific recombination and genetic control from a large pedigree analysis.
recom_chr_1_temp <- map

AFD_analyses <- data.frame(matrix(unlist(Fst), ncol=length(Fst), byrow=F))
heterozygosity_pop1 <- data.frame(matrix(unlist(het_pop1), ncol=length(het_pop1), byrow=F))
heterozygosity_pop2 <- data.frame(matrix(unlist(het_pop2), ncol=length(het_pop2), byrow=F))
# for (chrom in chromosome) {
  plot_het <- chill_df
  # [which(chill_df$chromosome==chrom),]

  if(ld_pairwise_plot==T){
#############################################################
################ PAIRWISE LD ANALYSES #######################
#############################################################
# use <-  chill_data$map$chromosome == chrom
genotype_temp <- genotype_chill
# [,use]
genotype <- genotype_temp[,order(as.numeric(colnames(genotype_temp)))]
snp_loc <- as.numeric(colnames(genotype))
#this is the mean distance between each snp which is used to determine the depth at which LD analyses are performed
mean_dis <- mean(diff(snp_loc))
ld_depth_b <- ceiling((ld_max_pairwise/mean_dis))
ld_snps <- ld(genotype,depth = ld_depth_b,stats = "R.squared") #function to calculate LD
ld_matrix <- as.matrix(ld_snps)
break_bins <- c(seq(1,ld_max_pairwise,ld_bins),ld_max_pairwise)
ld_columns_b <- ld_matrix
colnames(ld_columns_b) <- rownames(ld_columns_b)
ld_columns_b <- as.data.frame(as.table(as.matrix(ld_columns_b)))
ld_columns_b <- ld_columns_b[-ld_columns_b$Freq < 0,] #remove cases where LD was not calculated
ld_columns_b$Var1 <- as.numeric(as.character(ld_columns_b$Var1))
ld_columns_b$Var2 <- as.numeric(as.character(ld_columns_b$Var2))
ld_columns_b$dis <- ld_columns_b$Var2 - ld_columns_b$Var1 #determine the distance at which LD was calculated
ld_columns_c <- ld_columns_b[order(ld_columns_b$dis),]
ld_columns_d <- ld_columns_c[which(ld_columns_c$dis<ld_max_pairwise),]
bins_ld <- stats.bin(ld_columns_d$dis,ld_columns_d$Freq,breaks = break_bins)
final_pairwise_ld_temp <- unname(bins_ld$stats[2,])
# This is the final dataframe for LD pairwise
final_pairwise_ld <- as.data.frame(cbind(break_bins[-1],final_pairwise_ld_temp))
colnames(final_pairwise_ld) <- c("distance","rsqr")
final_pairwise_ld <- final_pairwise_ld[complete.cases(final_pairwise_ld$rsqr),]
# final_pairwise_ld <- final_pairwise_ld[sample(1:nrow(final_pairwise_ld),nrow(final_pairwise_ld)*0.1),]

ticks_lab_1 <- as.character(seq(0,ld_max_pairwise/1000000,0.5))
ticks_breaks_1 <- seq(0,ld_max_pairwise,500000)
print(

pairwise <- 
ggplot(final_pairwise_ld, aes(distance,rsqr)) +
  geom_line(size=0.5) +
  geom_smooth(method = "loess") +
  geom_hline(aes(yintercept=0.2,colour= "red"),size=1) +
  labs(x="Distance (Mbp)", y="r2 (LD)", title="Pairwise LD")+
  theme_bw(base_size = 12,base_family="Helvetica")+
  theme(legend.position = "none") +
  scale_x_continuous(breaks=ticks_breaks_1, labels=ticks_lab_1)

)
}
#############################################################
################ REGRESSION ANALYSES ##########################
############################################################
# ns_chr_2_chill_b <- read_csv(paste0("ns_chill_chr_ ", chrom, " .csv"))
ns_chr_1_temp <- chill_df
ns_chr_1_temp <- ns_chr_1_temp[which(!is.na(ns_chr_1_temp$z.HWE_pop1) & !is.na(ns_chr_1_temp$z.HWE_pop2)),]
ns_chr_1 <-    ns_chr_1_temp$loc_bp 
# ns_chr_1 <- reference$location
# ns_chr_1$midpoint <- (ns_chr_1$start + ns_chr_1$end)/2
# ns_chr_1 <- ns_chr_1$midpoint

  # chill_df[which(chill_df$MAF_pop2>0),"loc_bp"]
  # reference$location
  # chill_df$loc_bp
  # 
  # ns_chr_2_chill_b[,"chrom_start"]
chr_1 <- chill_df
# [which(chill_df$chromosome==chrom),]
chr_length <-  chromosome_length
  # unlist(unname(ns_chr_1[length(ns_chr_1)])) # this is the length of the chromosome based on the location of the last snp
recom_chr_1 <- recom_chr_1_temp
recom_maps <- list(recom_chr_1)
map_resolution <- 100000
i<-1
recom_maps_b <- NULL
# here the recombination rate is divided beween windows of 10,000 bp to obtain a higher resolution
for (map in recom_maps){
  temp_map <- NULL
  map <- as.data.frame(map)
  for (value in 1:nrow(map)) {
    temp <- as.data.frame(rep((as.numeric(map[value,1]) / 10), 10))
    distance <- map_resolution / nrow(temp)
    temp$loc_bp <- distance * (1:nrow(temp))
    if(value==1){
      temp$loc_bp <- temp$loc_bp
    }else{
      temp$loc_bp <- temp$loc_bp + ((value * map_resolution)-map_resolution)
    }
    temp_map <- rbind(temp_map, temp)
  }
  recom_maps_b[[i]] <- temp_map
  i<-i+1
}

r_map <-  recom_maps_b[[1]] 
colnames(r_map) <- c("cM","loc_bp")
r_map$cM <- r_map$cM * 100
  # recombination_map
# [which(recom_chr_1_temp$Chr==chrom),]
# here the resolution of the recombination map is calculated in windows of the size of the variable "map_resolution" (100 kbp) based on the recombination rate averaged over the two sexes
# map_cattle_binned <- stats.bin(recom_chr_1$locations_deleterious,recom_chr_1$c,breaks = seq(0,chr_length,map_resolution))
# map_cattle_binned_b <- unname(map_cattle_binned$stats[2,]* map_cattle_binned$stats[1,])
# map <- as.data.frame(map_cattle_binned_b)
# colnames(map) <- c("cM")
# map[is.na(map$cM),] <- 0
# recom_chr_1_b <- map$cM
# here the recombination rate is multiplied by 100 to convert it to centiMorgan. 1 centiMorgan = 0.01
# recom_chr_1_b <- recom_chr_1_b * 100
# recom_maps <- list(recom_chr_1_b)
# i<-1
# recom_maps_b <- NULL
# # here the recombination rate is divided beween windows of 10,000 bp to obtain a higher resolution
# for (map in recom_maps){
#   temp_map <- NULL
#   map <- as.data.frame(map)
#   for (value in 1:nrow(map)) {
#     temp <- as.data.frame(rep((as.numeric(map[value,1]) / 10), 10))
#     distance <- map_resolution / nrow(temp)
#     temp$loc_bp <- distance * (1:nrow(temp))
#     if(value==1){
#       temp$loc_bp <- temp$loc_bp
#     }else{
#       temp$loc_bp <- temp$loc_bp + ((value * map_resolution)-map_resolution)
#     }
#     temp_map <- rbind(temp_map, temp)
#   }
#   recom_maps_b[[i]] <- temp_map
#   i<-i+1
# }

# recom_maps_b[[1]]$cumsum <- cumsum(recom_maps_b[[1]][,1])
location_neutral_loci_analysis <- location_neutral_loci_analysis[order(location_neutral_loci_analysis)]

response_chr_1_temp <- as.data.frame(cbind(location_neutral_loci_analysis,unlist(AFD_analyses[number_of_generations,]),unlist(heterozygosity_pop1[number_of_generations,]),unlist(heterozygosity_pop2[number_of_generations,])))
colnames(response_chr_1_temp) <- c("loc_bp","P.AB_pop1","het_pop1","het_pop2")
response_chr_1_temp <- response_chr_1_temp[order(response_chr_1_temp$loc_bp),]
response_chr_1_temp <- response_chr_1_temp[which(response_chr_1_temp$het_pop1 > 0  & response_chr_1_temp$het_pop2 > 0), ]

brk_regression <- seq(0,chr_length,resolution_regression_analyses)
response_chr_1_temp_b <- as.numeric(stats.bin(response_chr_1_temp$loc_bp, response_chr_1_temp$P.AB_pop1, breaks = brk_regression)[[3]][2,])
# response_chr_1_temp_b <- as.numeric(stats.bin(response_chr_1_temp$loc_bp, response_chr_1_temp$het_pop2, breaks = brk_regression)[[3]][2,])


response_chr_1_temp_c <- as.data.frame(cbind(brk_regression[-1],response_chr_1_temp_b))
colnames(response_chr_1_temp_c) <- c("loc_bp","P.AB_pop1")
response_chr_1_temp_c <- response_chr_1_temp_c[!is.na(response_chr_1_temp_c$P.AB_pop1),]
# discarding those loci located outside of the range of the distances to be tested
keep_loci <- which(response_chr_1_temp_c$loc_bp > max(dist_bp) & response_chr_1_temp_c$loc_bp < (response_chr_1_temp_c[nrow(response_chr_1_temp_c),"loc_bp"] - max(dist_bp)))
response_chr_1_temp_c <- response_chr_1_temp_c[keep_loci,]
# brk_response <- seq(0,chr_length,1000000)
# response_chr_1_temp <- stats.bin(response_chr_1_temp$position, response_chr_1_temp$P.AB_pop1, breaks = brk_response)[[3]]
# # windows with low number of snps
# low_snps_response <-unique(which(response_chr_1_temp[1,]<=10))
#   response_chr_1_temp[2,low_snps_response] <- 0
#   response_chr_1_temp <- response_chr_1_temp[2,]
# # het_all_loci <- as.numeric(stats.bin(plot_het$position, plot_het$P.AB_pop1, breaks = brk)[[3]][2,])
# response_chr_1_temp[is.na(response_chr_1_temp)] <- 0
# response_chr_1_temp <- as.data.frame(cbind(unname(response_chr_1_temp),1:length(response_chr_1_temp)))
# colnames(response_chr_1_temp) <- c("P.AB_pop1","position")
# response_chr_1_temp$position <- response_chr_1_temp$position * 1000000
# response_chr_1_temp <- response_chr_1_temp[which(response_chr_1_temp$P.AB_pop1>0),]
response_chr_1 <- response_chr_1_temp_c
# response_chr_1 <- response_chr_1[which(response_chr_1$P.AB_pop1!="NaN"),]
# the locations of each locus in 10 kbp (the same as the resolution of the recombination map)
loc_chr_1 <- round(response_chr_1_temp_c$loc_bp/10000,0)
# ######  R1 #######
####LINKAGE
# link_final_regression <- link_final[keep_loci,]
# cm_R1 <- as.data.frame(matrix(nrow = length(dist_bp),ncol= length(loc_chr_1)))
# cm_R1 <- link_final_regression
# colnames(cm_R1) <- as.character(dist_bp/1000000)
# centiMorgans. determining how many cM there are in a given physical distance (bp)
cm_R1 <- as.data.frame(matrix(nrow = length(dist_bp),ncol= length(loc_chr_1)))
j <- 1
# here the distances to be tested (in Mbp) are converted to 10 Kbp (the same as the resolution of the recombination map)
dist_bp_b <- dist_bp / 10000
for (locus in loc_chr_1) {
  cm_temp_2 <-NULL
  for (d in dist_bp_b) {
    cm_left <- sum(r_map[locus:(locus - d),1])
    cm_right <- sum(r_map[locus:(locus + d),1]) - r_map[locus,1]
    cm_temp <- cm_left + cm_right
    cm_temp_2 <- rbind(cm_temp_2,cm_temp)
  }
  cm_R1[,j] <- cm_temp_2
  j <- j + 1
}
cm_R1 <- as.data.frame(t(cm_R1))
colnames(cm_R1) <- as.character(dist_bp/1000000)

#non-synonymous
number_ns_R1 <- as.data.frame(matrix(nrow = length(response_chr_1$loc_bp)))
number_ns_R1 <- number_ns_R1[,-1]
# here the distances to be tested (in Mbp) are converted to 10 Kbp (the same as the resolution of the recombination map)
dist_vector_c <- dist_bp 
for (dist_var in dist_vector_c) {
  distance_left_temp <- response_chr_1$loc_bp - dist_var
  distance_right_temp <- response_chr_1$loc_bp + dist_var
  distance_temp <- as.data.frame(cbind(distance_left_temp,distance_right_temp))
  ns_temp <- as.data.frame(matrix(nrow = length(response_chr_1$loc_bp)))
  for (r in 1:nrow(distance_temp)) {
    ns_temp[r,] <- length(which(ns_chr_1>distance_temp[r,1] & ns_chr_1<distance_temp[r,2]))
  }
  number_ns_R1 <- cbind(number_ns_R1,ns_temp)
}
colnames(number_ns_R1) <- as.character(dist_bp/1000000)


#### chr_1 R2 ####
# centiMorgans. determine for each locus the number of bp per cM
sum_tot_left_chr_1<- as.data.frame(matrix(nrow = length(loc_chr_1),ncol=2 ))
sum_tot_right_chr_1<- as.data.frame(matrix(nrow = length(loc_chr_1),ncol=2 ))
for (locus_2 in loc_chr_1) {
  sum_left_chr_1 <- NULL
  sum_right_chr_1 <- NULL
  for (dist_var_2 in dist_cM) {
    # accumulative addition of the number of cM from the location of each locus to the beginning of the chromosome (to the left)
    sum_left_chr_1 <- cumsum(r_map[locus_2:1,1])
    # getting all the rows that are farther than the genetic distance tested (dist_var_2) 
    n_rows_left_chr_1_temp <- which(sum_left_chr_1 > dist_var_2)
    # getting the first row at which the genetic distance tested is reached. every row represents 10 kbp
    n_rows_left_chr_1 <- min(n_rows_left_chr_1_temp)
    # placing the result in the correct location in the df
    sum_tot_left_chr_1[which(loc_chr_1==locus_2),which(dist_cM==dist_var_2)] <- n_rows_left_chr_1
    sum_right_chr_1 <- cumsum(r_map[locus_2:nrow(r_map),])
    n_rows_right_chr_1 <- min(which(sum_right_chr_1 > dist_var_2))
    sum_tot_right_chr_1[which(loc_chr_1==locus_2),which(dist_cM==dist_var_2)] <- n_rows_right_chr_1
  }
}
# multiplying by 10000 to obtain the real distance in bp 
sum_tot_left_chr_1 <- sum_tot_left_chr_1*10000
sum_tot_right_chr_1 <- sum_tot_right_chr_1*10000
colnames(sum_tot_left_chr_1) <- as.character(dist_cM)
colnames(sum_tot_right_chr_1) <- as.character(dist_cM)
# non_synonymous. obtain for each locus the number of ns occurring within the distances to be tested
loc_chr_1_ns <- response_chr_1$loc_bp
number_ns_R2_chr_1 <- as.data.frame(matrix(nrow = length(loc_chr_1_ns)))
number_ns_R2_chr_1 <- number_ns_R2_chr_1[,-1]
for (dist_var_3 in as.character(dist_cM)) {
location_left_temp_chr_1 <- loc_chr_1_ns - sum_tot_left_chr_1[dist_var_3]
location_right_temp_chr_1 <- loc_chr_1_ns + sum_tot_right_chr_1[dist_var_3]
location_temp_chr_1 <- as.data.frame(cbind(location_left_temp_chr_1,location_right_temp_chr_1))
number_temp_chr_1 <- as.data.frame(matrix(nrow = length(loc_chr_1_ns)))
for (r in 1:nrow(location_temp_chr_1)) {
  number_temp_chr_1[r,] <- length(which(ns_chr_1>location_temp_chr_1[r,1] & ns_chr_1<location_temp_chr_1[r,2]))
}
number_ns_R2_chr_1 <- cbind(number_ns_R2_chr_1,number_temp_chr_1)
}
colnames(number_ns_R2_chr_1) <- as.character(dist_cM)

ns_R1 <- number_ns_R1
ns_R2 <- number_ns_R2_chr_1

#### R1 ####

results_ns_R1 <- as.data.frame(matrix(nrow = 10))
# for (col_a in 1:ncol(response_var_low_snps)) {
  results_ns_R1[1,1] <- glance(lm(response_chr_1$P.AB_pop1~cm_R1[,1]+ns_R1[,1]))$r.squared
  results_ns_R1[2,1] <- glance(lm(response_chr_1$P.AB_pop1~cm_R1[,2]+ns_R1[,2]))$r.squared
  results_ns_R1[3,1] <- glance(lm(response_chr_1$P.AB_pop1~cm_R1[,3]+ns_R1[,3]))$r.squared
  results_ns_R1[4,1] <- glance(lm(response_chr_1$P.AB_pop1~cm_R1[,4]+ns_R1[,4]))$r.squared
  results_ns_R1[5,1] <- glance(lm(response_chr_1$P.AB_pop1~cm_R1[,5]+ns_R1[,5]))$r.squared
  results_ns_R1[6,1] <- glance(lm(response_chr_1$P.AB_pop1~cm_R1[,6]+ns_R1[,6]))$r.squared
  results_ns_R1[7,1] <- glance(lm(response_chr_1$P.AB_pop1~cm_R1[,7]+ns_R1[,7]))$r.squared
  results_ns_R1[8,1] <- glance(lm(response_chr_1$P.AB_pop1~cm_R1[,8]+ns_R1[,8]))$r.squared
  results_ns_R1[9,1] <- glance(lm(response_chr_1$P.AB_pop1~cm_R1[,9]+ns_R1[,9]))$r.squared
  results_ns_R1[10,1] <- glance(lm(response_chr_1$P.AB_pop1~cm_R1[,10]+ns_R1[,10]))$r.squared
# }

results_ns_R1$distance <- dist_bp/1000000
results_ns_R1$test <- "Non-synonymous"

 colnames(results_ns_R1) <-   c("het","distance","test")

#### R2####
results_ns_R2 <- as.data.frame(matrix(nrow = 10))

# for (col_a in 1:ncol(response_var_low_snps)) {
  results_ns_R2[1,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,1]))$p.value
  results_ns_R2[2,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,2]))$p.value
  results_ns_R2[3,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,3]))$p.value
  results_ns_R2[4,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,4]))$p.value
  results_ns_R2[5,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,5]))$p.value
  results_ns_R2[6,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,6]))$p.value
  results_ns_R2[7,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,7]))$p.value
  results_ns_R2[8,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,8]))$p.value
  results_ns_R2[9,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,9]))$p.value
  results_ns_R2[10,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,10]))$p.value
# }
  results_ns_R2$distance <- dist_cM
results_ns_R2$test <- "Non-synonymous"
colnames(results_ns_R2) <- c("het","distance","test")

results_ns_R2_rsqr <- as.data.frame(matrix(nrow = 10))
# for (col_a in 1:ncol(response_var_low_snps)) {
  results_ns_R2_rsqr[1,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,1]))$r.squared
  results_ns_R2_rsqr[2,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,2]))$r.squared
  results_ns_R2_rsqr[3,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,3]))$r.squared
  results_ns_R2_rsqr[4,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,4]))$r.squared
  results_ns_R2_rsqr[5,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,5]))$r.squared
  results_ns_R2_rsqr[6,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,6]))$r.squared
  results_ns_R2_rsqr[7,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,7]))$r.squared
  results_ns_R2_rsqr[8,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,8]))$r.squared
  results_ns_R2_rsqr[9,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,9]))$r.squared
  results_ns_R2_rsqr[10,1] <- glance(lm(response_chr_1$P.AB_pop1~ns_R2[,10]))$r.squared
# }
  results_ns_R2_rsqr$distance <- dist_cM
results_ns_R2_rsqr$test <- "Non-synonymous"
colnames(results_ns_R2_rsqr) <- c("het","distance","test")

plot_low_R2 <- melt(results_ns_R2_rsqr,id = c("distance","test"))
colnames(plot_low_R2) <- c("distance","test","Statistic","value")

plot_low_R1 <- melt(results_ns_R1,id = c("distance","test"))
colnames(plot_low_R1) <- c("distance","test","Statistic","value")


# inferred_distance <- (chr_1[nrow(chr_1),"loc_bp"] / (sum(recom_chr_1_b)/10) ) * 5
##########################################
 #######################################
# sig_ns_R2 <- which(results_ns_R2_rsqr$het==max(results_ns_R2_rsqr$het))
# # summary(lm(response_chr_1$P.AB_pop1~ns_R2[,sig_ns_R2]))
# plot_ns_R2 <- as.data.frame(cbind(response_chr_1$P.AB_pop1,ns_R2[,sig_ns_R2]))
# colnames(plot_ns_R2) <- c("Heterozygosity","Non_synonymous")
# print(
# ggplot(plot_ns_R2,aes(x=Heterozygosity,y=Non_synonymous))+
#   geom_point()+
#   geom_smooth(method='lm', formula= y~x)
# )
 # ggsave(paste("ns_R2_chr_",chrom,".pdf"),  width = 7, height =7, units = "in", dpi="retina", bg = "transparent" )
##########################################
 #######################################
 sig_R1 <- which(results_ns_R1$het==max(results_ns_R1$het))
sig_R1_distance <- results_ns_R1[sig_R1,"distance"]
  # summary(lm(response_chr_1$P.AB_pop1~cm_R1[,sig_R1]))
 plot_cm_R1 <- as.data.frame(cbind(response_chr_1$P.AB_pop1,cm_R1[,sig_R1]))
 colnames(plot_cm_R1) <- c("Fst","centiMorgans")

recombination_plot <- 
  
  ggplot(plot_cm_R1,aes(y=Fst,x=centiMorgans))+
  geom_point(color="deeppink",alpha=1,size=2)+
  geom_smooth(method='lm', formula= y~x)+
    labs(x = "centiMorgans",y="Fst")+
theme_bw(base_size = 12,base_family="Helvetica")+
        labs(title=paste("Distance",sig_R1_distance,"Mbp"))
  


  # ggsave(paste("cm_R1_chr_",chrom,".pdf"),  width = 7, height =7, units = "in", dpi="retina", bg = "transparent" )
##########################################
 #######################################
  # summary(lm(response_chr_1$P.AB_pop1~ns_R1[,sig_R1]))
  plot_ns_R1 <- as.data.frame(cbind(response_chr_1$P.AB_pop1,ns_R1[,sig_R1]))
colnames(plot_ns_R1) <- c("Fst","Non_synonymous")

non_synonymous <- 
ggplot(plot_ns_R1,aes(y=Fst,x=Non_synonymous))+
  geom_point(color="deeppink",alpha=1,size=2)+
  geom_smooth(method='lm', formula= y~x)+
theme_bw(base_size = 12,base_family="Helvetica")+
  labs(x = "Targets of selection",y="Fst")+
      labs(title=paste("Distance",sig_R1_distance,"Mbp"))

  # ggsave(paste("ns_R1_chr_",chrom,".pdf"),  width = 7, height =7, units = "in", dpi="retina", bg = "transparent" )

r_squared_R2 <-
  ggplot(plot_low_R2, aes(x=distance, y=value, colour=test)) +
    geom_line(size=1,color="deeppink") +
  geom_point(size=2,color="blue") +
    # facet_wrap(~Statistic) +
    # scale_y_continuous(name="Corrected p-value (holm)")+
    scale_y_continuous(name="R-squared")+
    # ,breaks=c(1,0.75,0.5,0.25,0.05)) +
    scale_x_continuous(name="Distance (cM)",breaks=results_ns_R2_rsqr$distance ) +
        # labs(subtitle=paste("SNP's in chromosome",chrom),title = "Regression using GENETIC distance")+ 
theme_bw(base_size = 12)+
    theme(legend.title=element_blank())
    # theme(legend.loc_bp="bottom")


sig_R2 <- which(results_ns_R2_rsqr$het ==max(results_ns_R2_rsqr$het))
sig_R2_distance <- results_ns_R2_rsqr[sig_R2,"distance"]
  # summary(lm(response_chr_1$P.AB_pop1~cm_R1[,sig_R1]))
 plot_ns_R2 <- as.data.frame(cbind(response_chr_1$P.AB_pop1,ns_R2[,sig_R1]))
 colnames(plot_ns_R2) <- c("Heterozygosity","Non_synonymous")

non_synonymous_R2 <- 
ggplot(plot_ns_R2,aes(y=Heterozygosity,x=Non_synonymous))+
  geom_point(color="deeppink",alpha=0.5)+
  geom_smooth(method='lm', formula= y~x)+
theme_bw(base_size = 12)+
      labs(title=paste("Distance",sig_R1_distance,"Mbp"))

r_squared_R1 <-
    ggplot(plot_low_R1, aes(x=distance, y=value, colour=test)) +
    geom_line(size=1,color="blue") +
  geom_point(size=2,color="deeppink") +
    # facet_wrap(~Statistic) +
    # scale_y_continuous(name="Corrected p-value (holm)")+
    labs(title="Sets of multiple regressions \nFst~cM+targets of selection")+
    scale_y_continuous(name="R-squared (regression)")+
    # ,breaks=c(1,0.75,0.5,0.25,0.05)) +
    scale_x_continuous(name="Distance (Mbp)",breaks=seq(0,2,0.4)) +
                         # results_ns_R2_rsqr$distance ) +
        # labs(subtitle=paste("SNP's in chromosome",chrom),title = "Regression using GENETIC distance")+ 
theme_bw(base_size = 12,base_family="Helvetica")+
    theme(legend.title=element_blank())

if(haploview_plot==T){

#############################################################
################ HAPLOVIEW ##########################
#############################################################
  snp_final <- chill_df
# [which(chill_df$chromosome==chrom),]
ld_matrix_2 <- ld_matrix
colnames(ld_matrix_2) <- rownames(ld_matrix_2)
rownames(ld_matrix_2) <- 1:nrow(ld_matrix_2)
colnames(ld_matrix_2) <- 1:ncol(ld_matrix_2)
ld_columns_2 <- as.data.frame(as.table(as.matrix(ld_matrix_2)))
ld_columns_2 <- ld_columns_2[-ld_columns_2$Freq < 0,] #remove cases where LD was not calculated
ld_columns_2$Var1 <- as.numeric(as.character(ld_columns_2$Var1))
ld_columns_2$Var2 <- as.numeric(as.character(ld_columns_2$Var2))
ld_columns_2 <- ld_columns_2[complete.cases(ld_columns_2),]
ld_columns_2$Freq <- round(ld_columns_2$Freq,1)
raster_haplo <- rasterFromXYZ(ld_columns_2)
polygon_haplo <- rasterToPolygons(raster_haplo, fun=NULL, n=4, na.rm=TRUE, digits=3, dissolve=T)
polygon_haplo <- elide(polygon_haplo,rotate=45)
polygon_haplo$id <- rownames(as.data.frame(polygon_haplo))
polygon_haplo.pts <- fortify(polygon_haplo, polygon_haplo="id") #this only has the coordinates
polygon_haplo.df <- merge(polygon_haplo.pts, polygon_haplo, by="id", type='left') # add the attributes back 
width_poly <-  round(max(polygon_haplo.df$long),0)
height_poly <-  max(polygon_haplo.df$lat)

reduce_factor <- nrow(snp_final)/width_poly
enlarge_factor <- width_poly/nrow(snp_final)
#this is to correct the cell size unit when the polygon figure is more or less than 1000 units 
correction_factor <- (width_poly/1000)

ld_matrix_3 <- ld_matrix
  # read.table(paste0(path.folder_sim,"/",files_LD_pop1[1]),sep = ",",skip = variables_number,row.names = 1)
ld_matrix_3 <- as.matrix(ld_matrix_3)
dimnames(ld_matrix_3)<-NULL 
matrix_rotate <- rotate.matrix_2(x=ld_matrix_3, angle=-45, method="simple")
matrix_rotate_2 <- matrix_rotate
matrix_rotate_2[matrix_rotate_2==0] <- NA
matrix_rotate_3 <- t(na.locf(t(matrix_rotate_2), fromLast = F, na.rm = T))

means_col <- apply(matrix_rotate_3,2,var,na.rm = T)
means_col[means_col==0] <- NA
means_col <- na.locf(means_col,fromLast =T)
means_col <- means_col * 2
df_col <- as.data.frame(matrix(ncol = 2 ,nrow = length(means_col))) 
df_col[,1] <- 1:nrow(df_col)
df_col[,2] <- means_col
df_col[,2] <-  rescale(df_col[,2],to = c(0,height_poly))
means_col_2 <- means_col

mean_column <- as.numeric(summary(means_col_2,na.rm=T)[1:6])
mean_column <- rescale(mean_column,to = c(0,height_poly))

# as the matrix was rotated the position of the snps is not correct anymore. the following code reassigns the snp position based on the second row from bottom to top of the rotated matrix. the first two and the last snps are removed to take in account the snps that are not present in the second row
 second_row_temp <- which(!is.na(matrix_rotate_3[,1])) -1
 second_row <-  second_row_temp[length(second_row_temp)]
 second_row_2 <- matrix_rotate_2[second_row,]
 second_row_3 <- second_row_2
 reassign_loc <- snp_final$loc_bp[3:nrow(snp_final)]
 reassign_loc <- reassign_loc[1:length(reassign_loc)-1]
 
 element_reassign_loc <- 1
 for(loc in 1:length(second_row_2)){
   if(is.na(second_row_2[loc])){
     next
   }else{
     second_row_3[loc] <- reassign_loc[element_reassign_loc]
      element_reassign_loc <-  element_reassign_loc + 1
   }
 }
 
 # putting back the first two snps and the last that were removed
 second_row_3[1:2] <- snp_final$loc_bp[1:2]
 second_row_3[length(second_row_3)] <- snp_final$loc_bp[length(snp_final$loc_bp)]
 #filling the NAs
 second_row_4 <- na.locf(second_row_3,fromLast=T, na.rm = F)
  second_row_4 <- na.locf(second_row_4,fromLast=F, na.rm = F)

 ld_threshold <- 0.1
  second_row_ver_2 <- na.locf(second_row_2,fromLast=T)
 first_row <- which(!is.na(matrix_rotate_3[,1])) 
  first_row_2 <- matrix_rotate_3[first_row,]
 haplo_loc_test <- first_row_2 >= ld_threshold
 haplo_loc_test <- c(haplo_loc_test,F)
start_haplo <- NULL
end_haplo <- NULL
 for(i in 1:(length(haplo_loc_test)-1)){
   if( haplo_loc_test[i] == haplo_loc_test[i+1]){
      next()
   }
  
   if( haplo_loc_test[i] == T ){
     end_haplo_temp <- i
     end_haplo <- c(end_haplo,end_haplo_temp )
   }
    if( haplo_loc_test[i] == F ){
     start_haplo_temp <- i
     start_haplo <- c(start_haplo,start_haplo_temp )
   }
 }
# #this is to not take in account those cases of a recent recombination rate, which causes to not deivide the complete haplotype under LD
# remove_recent_recombination <- which(diff(start_haplo)<15)
# remove_recent_recombination_2 <- c(remove_recent_recombination,remove_recent_recombination+1)
# start_haplo <- start_haplo[-remove_recent_recombination_2]
# end_haplo <- end_haplo[-remove_recent_recombination_2]
start_haplo <- c(1,start_haplo)

start_haplo_2 <-  unique(second_row_4[start_haplo])
end_haplo_2 <-  unique(second_row_4[end_haplo])
haplo_1_ver_2 <- as.data.frame(cbind(start_haplo_2,end_haplo_2))
haplo_1_ver_2$size <- (haplo_1_ver_2[,2]-haplo_1_ver_2[,1])
# haplo_1_ver_2 <- haplo_1_ver_2[which(haplo_1_ver_2$size>500000),]

n_snps <- as.matrix(haplo_1_ver_2[,1:2])
n_snps <- n_snps[which(!duplicated(n_snps[,1])),]
n_snps <- n_snps[which(!duplicated(n_snps[,2])),]
n_snps[,1] <- n_snps[,1]+1
n_snps <- n_snps[which(n_snps[,1]!=n_snps[,2]),]
# n_snps <- n_snps[1:(nrow(n_snps)-1),]
n_snps <- unique(n_snps[,])
  
df.4.cut <- as.data.frame(table(cut(snp_final$loc_bp, breaks=n_snps)),stringsAsFactors=F)
df.4.cut <- df.4.cut[which(df.4.cut$Freq>10),]
df.4.cut_3 <- gsub("[][()]", "", df.4.cut$Var1,",")
df.4.cut_3 <- strsplit(df.4.cut_3,",")
df.4.cut_4 <- lapply(df.4.cut_3, as.numeric)
df.4.cut_4 <- as.data.frame(laply(df.4.cut_4,rbind))
df.4.cut_4[,3] <- (df.4.cut_4[,2]-df.4.cut_4[,1])

# this is to calculate the real distance in bp of the polygon figure of LD

real_distance <- c(0,second_row_4)
real_distance_2 <- diff(real_distance)
real_distance_3 <- cumsum(real_distance_2)
real_distance_4 <- as.data.frame(cbind(1:length(real_distance_3),real_distance_3))

test_var <- unname(unlist(df.4.cut_4[,1:2]))

location_test <- lapply(test_var, findInterval,vec = as.numeric(paste(unlist(real_distance_4$real_distance_3))))
location_test_2 <-  unlist(location_test)
test_var_2 <- as.data.frame(cbind(1:length(location_test_2),location_test_2))

hap_blocks <- as.data.frame(cbind(paste0("CHR_",1:nrow(df.4.cut_4)),
                                  rep(1,times=nrow(df.4.cut_4)),
                                  df.4.cut_4[,1:3],
                                  df.4.cut$Freq))
colnames(hap_blocks) <- c("BLOCK", "CHR" ,  "BP1" ,  "BP2"   ,"SIZE",  "NSNP" )


locations_temp <- as.data.frame(cbind(hap_blocks$BP1,hap_blocks$BP2))
locations_temp_2 <- locations_temp
colnames(locations_temp_2) <- c("start","end")
locations_temp_2$start_ld_plot <- unlist(lapply(locations_temp_2$start, findInterval,vec = as.numeric(paste(unlist(real_distance_4$real_distance_3)))))
locations_temp_2$end_ld_plot<- unlist(lapply(locations_temp_2$end, findInterval,vec = as.numeric(paste(unlist(real_distance_4$real_distance_3)))))
locations_temp_2$midpoint <-  (locations_temp_2$start + locations_temp_2$end) /2
locations_temp_2$midpoint_ld_plot <-  (locations_temp_2$start_ld_plot + locations_temp_2$end_ld_plot) /2
locations_temp_2$labels <- paste0(
  as.character(round(locations_temp_2$start/1000000,0)),"-",as.character(round(locations_temp_2$end/1000000,0))  )

haplo_temp_a <- locations_temp_2[which(as.numeric(row.names(locations_temp_2))%%2==1),]
haplo_temp_b <- locations_temp_2[which(as.numeric(row.names(locations_temp_2))%%2==0),]

ticks_breaks <- c(locations_temp_2$start_ld_plot,locations_temp_2$end_ld_plot)
ticks_breaks <- ticks_breaks[order(ticks_breaks)]
ticks_lab <- c(locations_temp_2$start,locations_temp_2$end)
ticks_lab <- round(ticks_lab/1000000,0)
ticks_lab <- as.character(ticks_lab[order(ticks_lab)])

ticks_joint <- as.data.frame(cbind(ticks_breaks,ticks_lab))
ticks_joint <- ticks_joint[!duplicated(ticks_joint$ticks_lab),]
ticks_joint$ticks_lab <- as.character(ticks_joint$ticks_lab)
ticks_joint$ticks_breaks <- as.numeric(as.character(ticks_joint$ticks_breaks))
  
#this is SNP's heterozygosity to be calculated alone
snp_het_alone <- snp_final[,c("loc_bp","P.AB_pop1")]
snp_het_alone$loc_bp <- 1:nrow(snp_het_alone)
snp_het_alone[,1] <- rescale(snp_het_alone[,1],to = c(0,width_poly))
snp_het_alone[,2] <- snp_het_alone[,2] * height_poly

# snp_fst <- snp_final[,c("loc_bp","Fst")]
# snp_fst$loc_bp <- 1:nrow(snp_fst)
# snp_fst[,1] <- rescale(snp_fst[,1],to = c(0,width_poly))
# snp_fst[,2] <- snp_fst[,2] * height_poly
 
recom_haplo <- r_map
  # as.data.frame(cbind((1:nrow(r_map)*map_resolution),map$map))
# colnames(recom_haplo) <- c("V1","V2")
# location_recom <- unlist(lapply(recom_haplo$loc_bp, findInterval,vec = as.numeric(paste(unlist(real_distance_4$real_distance_3)))))
# recom_haplo$real_distance <- location_recom
recom_haplo$V2 <- rescale(recom_haplo$cM,to=c(0,height_poly))
recom_haplo <- unique(recom_haplo)

# rug_ns <-  ns_chr_1[sample(1:length(ns_chr_1),length(ns_chr_1)*0.01)]
# rug_ns <- as.data.frame(rescale(rug_ns,to=c(0,max(real_distance_4$V1))))
# colnames(rug_ns) <- "ns"
# colnames(plot_haplo) <- c("Var1", "Var2" ,"Variant")

 # cols <- c("SNP's"="lightgoldenrod","Heterozygosity"="deeppink","Msats expected"="cyan")

# p1 <- ggplot() +
#  geom_polygon(data=polygon_haplo.df, aes(long,lat,group=group, fill=Freq)) +
#     scale_fill_viridis(name="R square",option= "viridis") +
#   new_scale("fill") +
#   geom_vline(xintercept = test_var_2$location_test_2,color="red",size=1/2)+
#  # geom_rect(data = plot_haplo, aes(xmin=Var2, xmax=Var2+enlarge_factor, ymin=Var1, ymax=Var1+enlarge_factor,fill=Variant),color="black",size=1/10,inherit.aes = F) +
#   # scale_fill_manual(values = c("red", "blue","black"), labels = c("Deleterious", "Alternative"))+
# # geom_line(data=df_col,aes(x=df_col$V1,y=df_col$V2),color="orange",inherit.aes = F,size=1/2) +
# # geom_hline(yintercept = mean_column[2],color="red",size=1/2)+
# 
#  # geom_step(data=recom_haplo,aes(x=recom_haplo$V1,y=recom_haplo$V2),color="white",inherit.aes = F,size=1/2) +
# # labs(x=NULL, y="Haplotypes", title="Simulation with a chromosome of 500 cM",subtitle = paste(nrow(snp_final),"SNP's ","mean Fst = ",round(final_res[1,18],2)," mean He = ",round(final_res[1,15],2),"s = ", s_gral,"h = ", h_gral,"q = ",q_gral, " with ",loci_number, " initial loci" ))+
# # scale_x_continuous(breaks=ticks_breaks, labels=ticks_lab)+
# theme_tufte(base_family="Helvetica") +
#   theme(legend.loc_bp = "none",
#         axis.title.x=element_blank(),
#         axis.title.y=element_text(angle = 0,face = "bold",size = 12,hjust=0,vjust=0.5),
#         axis.ticks = element_blank(),
#         axis.ticks.y=element_blank() ,
#          axis.ticks.x=element_blank() ,
#          axis.text.y = element_blank(),
#         plot.margin=unit(c(0,1,-1,1), "in"),
#         # legend.justification =c(0,0),
#          # legend.margin=margin(0,0,0,0, unit = "in"),
#         # legend.box.margin=margin(-1,0,0,0, unit = "in"),
#          # legend.box.spacing= unit(0, "in"),
#         plot.title = element_text(face = "bold",size = 14,hjust=0,vjust=0.5),
#         axis.text.x = element_blank()) +
#   coord_fixed(ratio = 1/1)
colors_haploview <- c("Fst"="deeppink","Recombination"="blue","Deleterious_variants"="black","Haplotypes limits"="lightgoldenrod3")
labels_haplo <- as.character(1:nrow(locations_temp_2))


haploview <- ggplot() +
  geom_rect(aes(xmin=haplo_temp_a$start_ld_plot, xmax=haplo_temp_a$end_ld_plot, ymin=min(polygon_haplo.df$lat)-20, ymax=max(polygon_haplo.df$lat)+20),color="cornsilk3",fill="cornsilk3")+
    geom_rect(aes(xmin=haplo_temp_b$start_ld_plot, xmax=haplo_temp_b$end_ld_plot, ymin=min(polygon_haplo.df$lat)-20, ymax=max(polygon_haplo.df$lat)+20),color="cornsilk4",fill="cornsilk4")+
# geom_rect(aes(xmin=haplo_temp_b$V1, xmax=haplo_temp_b$V2, ymin=0, ymax=1),color="navajowhite",fill="navajowhite")+
  geom_polygon(data=polygon_haplo.df, aes(long,lat,group=group, fill=Freq)) +
  scale_fill_viridis(name="R squared (LD)",option= "viridis") +
  # geom_step(data=recom_haplo,aes(x=recom_haplo$V1,y=recom_haplo$V2,color="Recombination"),inherit.aes = F,size=0.5,alpha=0.8) +
     geom_vline(aes(xintercept =c(haplo_temp_b$start_ld_plot,haplo_temp_b$end_ld_plot),color="Haplotypes limits"),size=1/2)+

  geom_line(data=snp_het_alone,aes(x=loc_bp,y=P.AB_pop1,color= "Fst"),inherit.aes = F,size=0.5,alpha=0.8)+
  # geom_rug(data =rug_ns,aes(x=rug_ns$ns,y=2,color="Non_synonymous"),alpha = 1/4,sides="b",length = unit(0.05, "npc"))+
    annotate("text", x = locations_temp_2$midpoint_ld_plot, y = max(polygon_haplo.df$lat)+12, label = labels_haplo ,size=2) +
      annotate("text", x = locations_temp_2$midpoint_ld_plot, y = min(polygon_haplo.df$lat)-12, label = locations_temp_2$labels ,size=1) +
  labs( x= "Chromosome location (Mbp)", y="Het",title=paste("Chromosome",chrom,"-",nrow(plot_het),"SNP's"))+
  scale_x_continuous(breaks= ticks_joint$ticks_breaks, labels=ticks_joint$ticks_lab)+
  scale_colour_manual(name="",values=colors_haploview) +
  theme_tufte(base_family="Helvetica") +
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_text(hjust =0.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  coord_fixed(ratio = 1/1)
}
# print(haploview)
# p3 <- ggplot() +
#  geom_polygon(data=polygon_haplo.df, aes(long,lat,group=group, fill=Freq)) +
#   scale_fill_viridis(name="R square",option= "viridis") +
#   new_scale("fill") +
#  geom_rect(data = plot_haplo, aes(xmin=Var2, xmax=Var2+enlarge_factor, ymin=Var1, ymax=Var1+enlarge_factor,fill=Variant),color="black",size=1/10,inherit.aes = F) +
#    scale_fill_manual(values = c("red", "blue","black"), labels = c("Deleterious", "Alternative"))+
#   # new_scale("color") +
#     geom_col(data=fst_msats,aes(x=V1,y=V2),color="deeppink",inherit.aes = F,position = position_dodge2(10),size=0.75) +
#   geom_line(data=snp_fst, aes(x=position,y=Fst),color= "lightgoldenrodyellow",inherit.aes = F,size=1/2)+
#   geom_point(data=snp_fst,aes(x=position,y=Fst,color="SNP's"),inherit.aes = F,size =1/2,stroke=0,shape=16)+
#   geom_hline(yintercept = 0.6 * height_poly,color="cyan",size=1/2) +
#     scale_colour_manual(name="Markers",values=cols) +
#   labs(x= "Chromosome location (Mbp)", y="Fst", title="")+
#   scale_x_continuous(breaks=ticks_breaks, labels=ticks_lab)+
#   theme_tufte(base_family="Helvetica") +
#   theme(legend.position = "none",
#          # axis.title.x=element_blank(),
#          axis.title.y=element_text(angle = 0,face = "bold",size = 12,hjust=0.5,vjust=0.5),
#         axis.text.y=element_blank(),
#         axis.ticks = element_line(colour = "black", size = 0.5),
#         axis.ticks.y=element_blank() ,
#        plot.margin=unit(c(-1,1,1,1), "in"),
#       axis.title.x= element_text(hjust =0.5),
#         axis.text.x = element_text(face = "bold",size = 12)) +
#   coord_fixed(ratio = 1/1)

# p4 <- ggplot() +
# geom_polygon(data=polygon_haplo.df, aes(long,lat,group=group, fill=Freq)) +
# scale_fill_viridis(name="R square",option= "viridis") +
# new_scale("fill") +
# geom_rect(data = plot_haplo, aes(xmin=Var2, xmax=Var2+enlarge_factor, ymin=Var1, ymax=Var1+enlarge_factor,fill=Variant),color="black",size=1/30,inherit.aes = F) +
# scale_fill_manual(values = c("red", "blue"))+
# geom_line(data=HW_snp ,aes(x=HW_snp[,1],y=HW_snp[,2]),color="deeppink",inherit.aes = F,size=1/2) +
# geom_point(data=HW_snp ,aes(x=HW_snp[,1],y=HW_snp[,2]),color="deeppink4",inherit.aes = F,size = 1/2, stroke = 0, shape = 16) +
# labs(x="Chromosome location (Mbp)", y="h*s", title=NULL)+
# scale_x_continuous(breaks=ticks_breaks, labels=ticks_lab)+
# theme_tufte(base_family="Helvetica") +
# theme(legend.position = "none",
#          # axis.title.x=element_blank(),
#          axis.title.y=element_text(angle = 0,face = "bold",size = 12,hjust=0,vjust=0.5),
#         axis.text.y=element_blank(),
#         axis.ticks = element_line(colour = "black", size = 0.25),
#         axis.ticks.y=element_blank() ,
#        plot.margin=unit(c(-1,1,1,1), "in"),
#       axis.title.x= element_text(hjust =0.5),
#         axis.text.x = element_text(face = "bold",size = 12)) +
#   coord_fixed(ratio = 1/1)


# het_all_loci_temp <- stats.bin(plot_het$position, plot_het$P.AB_pop1, breaks = brk)[[3]]
# 
# # windows with low number of snps
# low_snps <-unique(which(het_all_loci_temp[1,]<=10))
#   het_all_loci_temp[2,low_snps] <- 0
#   het_all_loci <- het_all_loci_temp[2,]

het_all_loci <- as.data.frame(cbind(location_neutral_loci_analysis,unlist(AFD_analyses[number_of_generations,]),unlist(heterozygosity_pop1[number_of_generations,]),unlist(heterozygosity_pop2[number_of_generations,])))
colnames(het_all_loci) <- c("loc_bp","P.AB_pop1","het_pop1","het_pop2")
# het_all_loci <- het_all_loci[which(het_all_loci$het_pop1 > 0  & het_all_loci$het_pop2 > 0), ]
   # as.data.frame(cbind(location_neutral_loci_analysis,unlist(AFD_analyses[number_of_generations,])))
 # het_all_loci <- het_all_loci[which(het_all_loci$P.AB_pop1!="NA"),]
colnames(het_all_loci) <- c("loc_bp","P.AB_pop1")
# brk <- c(0,het_all_loci$loc_bp)
resolution_plot <- 200000
 brk <- seq(0,chr_length,resolution_plot)
het_all_loci <- as.numeric(stats.bin(het_all_loci$loc_bp, het_all_loci$P.AB_pop1, breaks = brk)[[3]][2,])
het_all_loci <- as.data.frame(cbind(brk[-1]/1000000,het_all_loci))
colnames(het_all_loci) <- c("loc_bp","P.AB_pop1")
# het_all_loci[is.na(het_all_loci)] <- 0
recom_all_loci <- stats.bin(r_map$loc_bp, r_map$cM, breaks = brk)
recom_all_loci_sum <- unlist(unname( ((recom_all_loci[[3]][2,] * recom_all_loci[[3]][1,]))  ))
# recom_all_loci_sum[is.na(recom_all_loci_sum)] <- 0
recom_all_loci_relative <- rescale(recom_all_loci_sum,to=c(0,max(het_all_loci$P.AB_pop1,na.rm = T)))

ns_all_loci <- as.data.frame(table(cut(ns_chr_1, breaks=brk)))
ns_all_loci  <- rescale(ns_all_loci$Freq,to=c(0,max(het_all_loci$P.AB_pop1,na.rm = T)))
 brk_2 <- brk[-1]

df_snps <- as.data.frame(cbind(brk_2/1000000,het_all_loci$P.AB_pop1,recom_all_loci_relative,ns_all_loci))
colnames(df_snps) <- c("location","het","recom","ns")

rescale_fun <- function(x){rescale(x,to=c(min(recom_all_loci_sum),max(recom_all_loci_sum)))}

het_plot <-
ggplot() +
# geom_rect(aes(xmin=haplo_temp_a$start/1000000, xmax=haplo_temp_a$end/1000000, ymin=0, ymax=max(het_all_loci$P.AB_pop1,na.rm = T)),color="cornsilk3",fill="cornsilk3")+
# geom_rect(aes(xmin=haplo_temp_b$start/1000000, xmax=haplo_temp_b$end/1000000, ymin=0, ymax=max(het_all_loci$P.AB_pop1,na.rm = T)),color="cornsilk4",fill="cornsilk4")+
# geom_area(aes(x=df_snps$location,y = df_snps$het,color="Heterozygosity"),fill="black",alpha=1/5)+
geom_col(aes(x=df_snps$location,y = df_snps$het,color="Fst"),fill="black",alpha=1/5)+
# geom_vline(xintercept = (locations_temp/1000000),color="red",size=1)+
geom_line(aes(x=df_snps$location,y = df_snps$recom,color="Recombination"),size=1) +
geom_line(aes(x=df_snps$location,y = df_snps$ns,color="Deleterious_variants"),size=1) +
  # annotate("text", x = locations_temp_2$midpoint/1000000, y =(max(het_all_loci$P.AB_pop1,na.rm = T)-0.05), label = labels_haplo ,size=3,color="black") +
scale_y_continuous(sec.axis = sec_axis(rescale_fun, name = "centiMorgans/Mbp")) +
scale_x_continuous(breaks = round(seq(0, max(df_snps$location), by = 1),1)) +
theme_bw(base_size = 12)+
# scale_colour_manual(name="",values=colors_haploview) +
labs(x="Chromosome location (Mbp)", y="Fst/Mpb", title="")+
  theme(legend.position = "bottom")


layout <- "
CCDD
EEFF
"
 #  print(
 #    (  haploview +het_plot+ pairwise + r_squared_R1 +recombination + targets) + 
 #  plot_layout(design = layout, heights= unit(c(1,1,1,1), c('in','null','null','null')))
 # )
   print(
    ( pairwise + r_squared_R1 +recombination_plot + non_synonymous) + 
  plot_layout(design = layout, heights= unit(c(1,1), c('null','null')))
 )
  
 # ggsave(paste0("chillingham_chr_",chrom,"_targets_c.pdf"),  width = 10, height =10, units = "in", dpi="retina", bg = "transparent" )
  ggsave(paste0("fly_sim_chr_",chrom,"_regression.pdf"),  width = 6, height =6, units = "in", dpi="retina", bg = "transparent" )





# layout <- "
# AAAAAA
# BBBBBB
# #CCDD#
# #EEFF#
# "
#   print(
#     (  haploview +het_plot+ pairwise + r_squared_R1 +recombination_plot + non_synonymous) + 
#   plot_layout(design = layout, heights= unit(c(1,1,1,1), c('in','null','null','null')))
#  )
#  
#  ggsave(paste0("chillingham_chr_",chrom,".pdf"),  width = 10, height =10, units = "in", dpi="retina", bg = "transparent" )

 # }


