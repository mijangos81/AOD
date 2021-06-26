####### LOADING LIBRARIES #######
library(patchwork)
library(broom)
library(reshape2)
library(readr)
library(fields)
library(ggplot2)
#############################################################
################ VARIABLES ###############################
############################################################
stat_to_analyse <- "He"# options "He" or "Fst"
variables_number <- 36
simulation_type_2 <- "fly" # simulation_type_2 options: "chillingham", "fly" or "general"
path.folder <- paste0(getwd(),"/AOD_fly_figure_SIM_FLY_B_0.005_C_0.2_D_5e-05_F_TRUE104704")
files_pop_history <- paste0(path.folder,"/",dir(path.folder,pattern = "^final_stats_average"))
pop_history <- read.table(files_pop_history[1],header = F,sep=",",nrows = variables_number, colClasses = "character")
chrom <-  pop_history[which(pop_history$V1=="chromosome_name"),2]
number_of_generations <- as.numeric(pop_history[which(pop_history$V1=="gen_number_dispersal"),2])
number_of_stats_calculated <- 25
map_resolution <- 100000 # this is the resolution of the recombination map
# these are the distance that will be tested in the regression analyses
dist_cM <- seq(1,10,1) # distances in cM
dist_bp <- c(seq(200000,2000000,200000)) # distance in bp
# PAIRWISE LD
ld_bins <- 1000 # this is the length , in base pairs, of the bins 
ld_max_pairwise <-  as.numeric(pop_history[which(pop_history$V1=="ld_max_pairwise"),2]) # maximun distance, in basepairs, at which pairwise LD should be calculated
ld_resolution <- 10000 # resolution, in basepairs, at which LD should be measured
resolution_regression_analyses <- 1000000
#############################################################
################ PAIRWISE LD ANALYSES #######################
#############################################################
files_LD <- paste0(path.folder,"/",dir(path.folder,pattern = "^final_LD_pairwise"))
final_pairwise_ld <-  read.table(files_LD[1],header=T,sep=",",skip= variables_number,fill=T)
final_pairwise_ld <- final_pairwise_ld[complete.cases(final_pairwise_ld$rsqr),]
if(nrow(final_pairwise_ld)>3000){
  final_pairwise_ld <- final_pairwise_ld[sample(1:nrow(final_pairwise_ld),2000),]
}
ticks_lab_1 <- as.character(seq(0,ld_max_pairwise/1000000,ld_max_pairwise/10000000))
ticks_breaks_1 <- seq(0,ld_max_pairwise,ld_max_pairwise/10)
print(
  pairwise <- 
ggplot(final_pairwise_ld, aes(distance,rsqr)) +
  geom_line(size=0.5) +
  geom_smooth(method = "loess",se=F) +
  geom_hline(aes(yintercept=0.2,colour= "red"),size=1) +
  labs(x="Distance (Mbp)", y="r2 (LD)", subtitle="Pairwise LD")+
  theme_bw(base_size = 10,base_family="Helvetica")+
  theme(legend.position = "none") +
  scale_x_continuous(breaks=ticks_breaks_1, labels=ticks_lab_1)
)
#############################################################
################ REGRESSION ANALYSES ##########################
############################################################
files_NS <- paste0(path.folder,"/",dir(path.folder,pattern = "^snps"))
ns_chr_1_temp <-  read.table(files_NS[1],header=T,sep=",",skip= variables_number,fill=T)
# ns_chr_1_temp <- ns_chr_1_temp[which(!is.na(ns_chr_1_temp$z.HWE_pop1) & !is.na(ns_chr_1_temp$z.HWE_pop2)),]
ns_chr_1_temp <- ns_chr_1_temp[which(!is.na(ns_chr_1_temp$z.HWE_pop1)),]
ns_chr_1 <-    ns_chr_1_temp$loc_bp 

if(simulation_type_2=="fly"){
map_temp <- read_csv("fly_recom_map.csv")
map_temp$Chr <- as.character(map_temp$Chr)
map_temp <- map_temp[which(map_temp$Chr==chrom),]
map_temp <- as.data.frame(map_temp$cM/1000)
colnames(map_temp) <- "cM"
map_temp[is.na(map_temp$cM),] <- 0
}
if(simulation_type_2=="chillingham"){
  RecRates_All_Chromosomes <- read_csv("chill_recom_map.csv") 
RecRates_All_Chromosomes <- as.data.frame(RecRates_All_Chromosomes)
RecRates_All_Chromosomes$Chr <- as.character(RecRates_All_Chromosomes$Chr)
RecRates_All_Chromosomes_chr <- RecRates_All_Chromosomes[which(RecRates_All_Chromosomes$Chr==chrom),]
chr_length <- RecRates_All_Chromosomes_chr[nrow(RecRates_All_Chromosomes_chr),"Location"]
map_cattle_binned <- stats.bin(RecRates_All_Chromosomes_chr$Location,RecRates_All_Chromosomes_chr$`mean sexes`,breaks = seq(0,chr_length,map_resolution))
map_cattle_binned_b <- unname(map_cattle_binned$stats[2,]* map_cattle_binned$stats[1,])
map_cattle_binned_b
map_temp <- as.data.frame(map_cattle_binned_b)
colnames(map_temp) <- "cM"
map_temp[is.na(map_temp$cM),] <- 0
}

chr_length <- (nrow(map_temp)+1) * map_resolution

recom_chr_1 <- map_temp
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

# location_msats_experiment <-  c(8325394,2373262,4781938,4960235,7040194,8348440,11015745,12507696, 13153885, 14615705, 14995570,20706003)
location_neutral_loci_analysis <- seq(map_resolution/2,(nrow(map_temp)*map_resolution),map_resolution)
 # location_neutral_loci_analysis <-  c(seq(map_resolution/2,(nrow(map_temp)*map_resolution),map_resolution),location_msats_experiment)
location_neutral_loci_analysis <- location_neutral_loci_analysis[order(location_neutral_loci_analysis)]

number_loci <- length(location_neutral_loci_analysis)
number_of_lines <- number_of_stats_calculated * number_of_generations
initial_breaks <- seq(1,number_of_lines,number_of_generations)
end_breaks <- c((initial_breaks-1)[-1],number_of_lines)

files_response <- paste0(path.folder,"/",dir(path.folder,pattern = "^final_stats_average"))
mean_df_response <- as.data.frame(read.table(files_response[1],header=F,sep=",",skip= variables_number,fill=T))
mean_AFD_response <- mean_df_response[initial_breaks[1]:end_breaks[1],1:number_loci]
mean_sha_diff_response <- mean_df_response[initial_breaks[2]:end_breaks[2],1:number_loci]
mean_het_pop1_response <-  mean_df_response[initial_breaks[15]:end_breaks[15],1:number_loci]
mean_het_pop2_response <-  mean_df_response[initial_breaks[16]:end_breaks[16],1:number_loci]
mean_Fst_response <-  mean_df_response[initial_breaks[18]:end_breaks[18],1:number_loci]
mean_Fstp_response <-  mean_df_response[initial_breaks[19]:end_breaks[19],1:number_loci]
mean_Dest_response <-  mean_df_response[initial_breaks[20]:end_breaks[20],1:number_loci]

response_chr_1_temp <- as.data.frame(cbind(location_neutral_loci_analysis,unlist(mean_Fstp_response[number_of_generations,]),unlist(mean_het_pop1_response[number_of_generations,]),unlist(mean_het_pop2_response[number_of_generations,])))
colnames(response_chr_1_temp) <- c("loc_bp","Fst","het_pop1","het_pop2")
response_chr_1_temp <- response_chr_1_temp[order(response_chr_1_temp$loc_bp),]
# response_chr_1_temp <- response_chr_1_temp[which(response_chr_1_temp$het_pop1 > 0  & response_chr_1_temp$het_pop2 > 0), ]
# response_chr_1_temp <- response_chr_1_temp[which(response_chr_1_temp$het_pop1 > 0), ]

brk_regression <- seq(0,chr_length,resolution_regression_analyses)

if(stat_to_analyse=="Fst"){
response_chr_1_temp_b <- as.numeric(stats.bin(response_chr_1_temp$loc_bp, response_chr_1_temp$Fst, breaks = brk_regression)[[3]][2,])
}
if(stat_to_analyse=="He"){
response_chr_1_temp_b <- as.numeric(stats.bin(response_chr_1_temp$loc_bp, response_chr_1_temp$het_pop1, breaks = brk_regression)[[3]][2,])
}
response_chr_1_temp_c <- as.data.frame(cbind(brk_regression[-1],response_chr_1_temp_b))
colnames(response_chr_1_temp_c) <- c("loc_bp","stat")
response_chr_1_temp_c <- response_chr_1_temp_c[!is.na(response_chr_1_temp_c$stat),]
# discarding those loci located outside of the range of the distances to be tested
keep_loci <- which(response_chr_1_temp_c$loc_bp > max(dist_bp) & response_chr_1_temp_c$loc_bp < (response_chr_1_temp_c[nrow(response_chr_1_temp_c),"loc_bp"] - max(dist_bp)))
response_chr_1_temp_c <- response_chr_1_temp_c[keep_loci,]
response_chr_1 <- response_chr_1_temp_c
loc_chr_1 <- round(response_chr_1_temp_c$loc_bp/10000,0)
# ######  R1 #######
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
  results_ns_R1[1,1] <- glance(lm(response_chr_1$stat~cm_R1[,1]+ns_R1[,1]))$r.squared
  results_ns_R1[2,1] <- glance(lm(response_chr_1$stat~cm_R1[,2]+ns_R1[,2]))$r.squared
  results_ns_R1[3,1] <- glance(lm(response_chr_1$stat~cm_R1[,3]+ns_R1[,3]))$r.squared
  results_ns_R1[4,1] <- glance(lm(response_chr_1$stat~cm_R1[,4]+ns_R1[,4]))$r.squared
  results_ns_R1[5,1] <- glance(lm(response_chr_1$stat~cm_R1[,5]+ns_R1[,5]))$r.squared
  results_ns_R1[6,1] <- glance(lm(response_chr_1$stat~cm_R1[,6]+ns_R1[,6]))$r.squared
  results_ns_R1[7,1] <- glance(lm(response_chr_1$stat~cm_R1[,7]+ns_R1[,7]))$r.squared
  results_ns_R1[8,1] <- glance(lm(response_chr_1$stat~cm_R1[,8]+ns_R1[,8]))$r.squared
  results_ns_R1[9,1] <- glance(lm(response_chr_1$stat~cm_R1[,9]+ns_R1[,9]))$r.squared
  results_ns_R1[10,1] <- glance(lm(response_chr_1$stat~cm_R1[,10]+ns_R1[,10]))$r.squared
results_ns_R1$distance <- dist_bp/1000000
results_ns_R1$test <- "Non-synonymous"
colnames(results_ns_R1) <-   c("het","distance","test")

#### R2####
results_ns_R2 <- as.data.frame(matrix(nrow = 10))
  results_ns_R2[1,1] <- glance(lm(response_chr_1$stat~ns_R2[,1]))$p.value
  results_ns_R2[2,1] <- glance(lm(response_chr_1$stat~ns_R2[,2]))$p.value
  results_ns_R2[3,1] <- glance(lm(response_chr_1$stat~ns_R2[,3]))$p.value
  results_ns_R2[4,1] <- glance(lm(response_chr_1$stat~ns_R2[,4]))$p.value
  results_ns_R2[5,1] <- glance(lm(response_chr_1$stat~ns_R2[,5]))$p.value
  results_ns_R2[6,1] <- glance(lm(response_chr_1$stat~ns_R2[,6]))$p.value
  results_ns_R2[7,1] <- glance(lm(response_chr_1$stat~ns_R2[,7]))$p.value
  results_ns_R2[8,1] <- glance(lm(response_chr_1$stat~ns_R2[,8]))$p.value
  results_ns_R2[9,1] <- glance(lm(response_chr_1$stat~ns_R2[,9]))$p.value
  results_ns_R2[10,1] <- glance(lm(response_chr_1$stat~ns_R2[,10]))$p.value
  results_ns_R2$distance <- dist_cM
results_ns_R2$test <- "Non-synonymous"
colnames(results_ns_R2) <- c("het","distance","test")

results_ns_R2_rsqr <- as.data.frame(matrix(nrow = 10))
  results_ns_R2_rsqr[1,1] <- glance(lm(response_chr_1$stat~ns_R2[,1]))$r.squared
  results_ns_R2_rsqr[2,1] <- glance(lm(response_chr_1$stat~ns_R2[,2]))$r.squared
  results_ns_R2_rsqr[3,1] <- glance(lm(response_chr_1$stat~ns_R2[,3]))$r.squared
  results_ns_R2_rsqr[4,1] <- glance(lm(response_chr_1$stat~ns_R2[,4]))$r.squared
  results_ns_R2_rsqr[5,1] <- glance(lm(response_chr_1$stat~ns_R2[,5]))$r.squared
  results_ns_R2_rsqr[6,1] <- glance(lm(response_chr_1$stat~ns_R2[,6]))$r.squared
  results_ns_R2_rsqr[7,1] <- glance(lm(response_chr_1$stat~ns_R2[,7]))$r.squared
  results_ns_R2_rsqr[8,1] <- glance(lm(response_chr_1$stat~ns_R2[,8]))$r.squared
  results_ns_R2_rsqr[9,1] <- glance(lm(response_chr_1$stat~ns_R2[,9]))$r.squared
  results_ns_R2_rsqr[10,1] <- glance(lm(response_chr_1$stat~ns_R2[,10]))$r.squared

  results_ns_R2_rsqr$distance <- dist_cM
results_ns_R2_rsqr$test <- "Non-synonymous"
colnames(results_ns_R2_rsqr) <- c("het","distance","test")

plot_low_R2 <- reshape2::melt(results_ns_R2_rsqr,id = c("distance","test"))
colnames(plot_low_R2) <- c("distance","test","Statistic","value")

plot_low_R1 <- reshape2::melt(results_ns_R1,id = c("distance","test"))
colnames(plot_low_R1) <- c("distance","test","Statistic","value")

sig_R1 <- which(results_ns_R1$het==max(results_ns_R1$het))
sig_R1_distance <- results_ns_R1[sig_R1,"distance"]
plot_cm_R1 <- as.data.frame(cbind(response_chr_1$stat,cm_R1[,sig_R1]))
colnames(plot_cm_R1) <- c("stat","centiMorgans")

recombination_plot <- 
ggplot(plot_cm_R1,aes(y=stat,x=centiMorgans))+
geom_smooth(method='lm', formula= y~x)+
geom_point(color="deeppink",alpha=1,size=2)+
labs(x = "centiMorgans",y=stat_to_analyse)+
theme_bw(base_size = 10,base_family="Helvetica")+
labs(subtitle=paste("Size genomic neighborhood =",sig_R1_distance,"Mbp"))

plot_ns_R1 <- as.data.frame(cbind(response_chr_1$stat,ns_R1[,sig_R1]))
colnames(plot_ns_R1) <- c("stat","Non_synonymous")

non_synonymous <- 
ggplot(plot_ns_R1,aes(y=stat,x=Non_synonymous))+
geom_point(color="deeppink",alpha=1,size=2)+
geom_smooth(method='lm', formula= y~x)+
theme_bw(base_size = 10,base_family="Helvetica")+
labs(x = "Targets of selection",y=stat_to_analyse)+
labs(subtitle=paste("Size genomic neighborhood =",sig_R1_distance,"Mbp"))

r_squared_R2 <-
ggplot(plot_low_R2, aes(x=distance, y=value, colour=test)) +
geom_line(size=1,color="deeppink") +
geom_point(size=2,color="blue") +
scale_y_continuous(name="R-squared")+
scale_x_continuous(name="Distance (cM)",breaks=results_ns_R2_rsqr$distance ) +
theme_bw(base_size = 10)+
theme(legend.title=element_blank())

sig_R2 <- which(results_ns_R2_rsqr$het ==max(results_ns_R2_rsqr$het))
sig_R2_distance <- results_ns_R2_rsqr[sig_R2,"distance"]
plot_ns_R2 <- as.data.frame(cbind(response_chr_1$stat,ns_R2[,sig_R1]))
colnames(plot_ns_R2) <- c("Heterozygosity","Non_synonymous")

r_squared_R1 <-
ggplot(plot_low_R1, aes(x=distance, y=value, colour=test)) +
geom_line(size=1,color="blue") +
geom_point(size=2,color="deeppink") +
scale_y_continuous(name="R-squared (regression)")+
scale_x_continuous(name="Size genomic neighbourhood (Mbp)",breaks=results_ns_R1$distance)  +
theme_bw(base_size = 10,base_family="Helvetica")+
labs(subtitle="Sets of multiple regressions\nHe~cM+targets of selection")+
theme(legend.title=element_blank())

layout <- "
CCDD
EEFF
"
   print(
    ( pairwise + r_squared_R1 +recombination_plot + non_synonymous) + 
  plot_layout(design = layout, heights= unit(c(1,1), c('null','null')))
 )

ggsave(paste0(simulation_type_2,"_chr_",chrom,"_regression.pdf"),  width = 6, height =6, units = "in", dpi="retina", bg = "transparent" )
