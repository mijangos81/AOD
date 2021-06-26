####### LOADING LIBRARIES #######
library(patchwork)
library(broom)
library(reshape2)
library(readr)
library(fields)
library(ggplot2)
#############################################################
################ VARIABLES ###############################
#############################################################
chromosome <- c("X")
  # c("2L","2R","3L","3R","X")
distance_reg <- 5
dispersal_regimes <- c("Low","Moderate","High")
response_var <- "Fst" #options : "Fst","Fstp","Dest", "het_mean"
stat_label <- "Fst"# options: "He lost from T0 to T2", "Fst"
marker <- "SNPs" #options: "SNPs", "microsatellites"
map_resolution <- 100000 # this is the resolution of the recombination map
# these are the distance that will be tested in the regression analyses
dist_cM <- seq(2.5,12.5,2.5) # distances in cM
dist_bp <- seq(1000000,distance_reg*1000000,1000000) # distance in bp
# these are the recombination maps of all the chromosomes obtained from Ma, Li et al. (2016), Data from: Cattle sex-specific recombination and genetic control from a large pedigree analysis.
recom_chr_1_temp <- read_csv("fly_recom_BDGP5.csv") 
targets_of_selection_temp <- read.csv("fly_targets_of_selection.csv")
targets_of_selection_temp$CUB_dmel <- (targets_of_selection_temp$ns - targets_of_selection_temp$s) 
targets_of_selection <- targets_of_selection_temp[,c("start","end", "chr_name","CUB_dmel")]
targets_of_selection$midpoint <- (targets_of_selection$start+targets_of_selection$end)/2
targets_of_selection$size <- targets_of_selection$end - targets_of_selection$start
targets_of_selection <- targets_of_selection[!duplicated(targets_of_selection$start),]
targets_of_selection <- targets_of_selection[!duplicated(targets_of_selection$end),]

stats_together <- read_csv("stats_all_together.csv")
stats_together$het_mean <-  (stats_together$het_pop1 + stats_together$het_pop2)/2
chr_lengths <- as.data.frame(rbind(c("2L",23513712),c("2R",25286936),c("3L",28110227),c("3R",32079331),c("X",	23542271)))
colnames(chr_lengths) <- c("chr","length")

for (dispersal in dispersal_regimes) {
  
final_response_R1 <- as.data.frame(matrix(ncol=3))
colnames(final_response_R1) <- c("position","loc_BDGP5","Fst")
final_targets_R1 <- as.data.frame(matrix(ncol=length(dist_bp)))
colnames(final_targets_R1)<- as.character(dist_bp/1000000)
final_recombination_R1 <- as.data.frame(matrix(ncol=length(dist_bp)))
colnames(final_recombination_R1)<- as.character(dist_bp/1000000)

final_response_R2 <- as.data.frame(matrix(ncol=3))
colnames(final_response_R2) <- c("position","loc_BDGP5","Fst")
final_targets_R2 <- as.data.frame(matrix(ncol=length(dist_cM)))
colnames(final_targets_R2)<- as.character(dist_cM)
final_recombination_R2 <- as.data.frame(matrix(ncol=length(dist_cM)))
colnames(final_recombination_R2)<- as.character(dist_cM)

    for(chrom in chromosome){
targets <- targets_of_selection[which(targets_of_selection$chr_name== chrom ),c("midpoint","CUB_dmel","chr_name")]
targets <-  targets[order(targets$midpoint),]
#############################################################
################ REGRESSION ANALYSES ##########################
#############################################################
chr_1 <-  stats_together[which(stats_together$chromosome==chrom & stats_together$dispersal==dispersal),]
# chr_1 <-  stats_together[which(stats_together$chromosome==chrom & stats_together$dispersal==dispersal & stats_together$marker==marker),]

colnames(chr_1)[4] <- "position"
map_temp <- recom_chr_1_temp[which(recom_chr_1_temp$Chr==chrom),]
map_temp <- as.data.frame(map_temp)
chr_length <- as.numeric(chr_lengths[which(chr_lengths$chr==chrom),"length"])
recom_maps <- list(map_temp)
map_resolution <- 100000
i<-1
recom_maps_b <- NULL
# here the recombination rate is divided between windows of 10,000 bp to obtain a higher resolution
for (map_res in recom_maps){
  temp_map <- NULL
  map_res <- as.data.frame(map_res)
  for (value in 1:nrow(map_res)) {
    temp <- as.data.frame(rep((as.numeric(map_res[value,"cM"]) / 10), 10))
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
r_map$cM <- r_map$cM /10
# discarding those loci located outside of the range of the distances to be tested
response_chr_1_temp <- chr_1[which(chr_1$position > max(dist_bp) & chr_1$position < (chr_length - max(dist_bp)) ), ]
response_chr_1_temp <- response_chr_1_temp[order(response_chr_1_temp$position),]
response_chr_1 <- response_chr_1_temp
response_chr_1_targets <- response_chr_1[,c("position","loc_BDGP5",response_var)]
colnames(response_chr_1_targets) <- c("position","loc_BDGP5","Fst")
loc_chr_1 <-  response_chr_1_targets$loc_BDGP5
# ######  R1 #######
# centiMorgans. determining how many cM there are in a given physical distance (bp)
cm_R1 <- as.data.frame(matrix(nrow = length(dist_bp),ncol= length(loc_chr_1)))
j <- 1
recom_maps_test <- as.data.frame(cbind(map_temp[,"cM"],map_temp[,"midpoint"]))
colnames(recom_maps_test) <- c("cM","position")
for (locus in loc_chr_1) {
  cm_temp_2 <-NULL
  for (d in dist_bp) {
    cm_temp <- sum(recom_maps_test[which(recom_maps_test$position>(locus-d) & recom_maps_test$position<(locus+d)),"cM"])
    cm_temp_2 <- rbind(cm_temp_2,cm_temp)
  }
  cm_R1[,j] <- cm_temp_2
  j <- j + 1
}
cm_R1 <- as.data.frame(t(cm_R1))
colnames(cm_R1) <- as.character(dist_bp/1000000)
##### targets R1######
number_targets_R1 <- as.data.frame(matrix(nrow = length(response_chr_1_targets$position)))
number_targets_R1 <- number_targets_R1[,-1]
dist_vector_c <- dist_bp 
for (dist_var in dist_vector_c) {
  distance_left_temp <- response_chr_1_targets$position - dist_var
  distance_right_temp <- response_chr_1_targets$position  + dist_var
  distance_temp <- as.data.frame(cbind(distance_left_temp,distance_right_temp))
  targets_temp <- as.data.frame(matrix(nrow = length(response_chr_1_targets$position)))
  for (r in 1:nrow(distance_temp)) {
    targets_temp_a <- targets[which(targets$midpoint>distance_temp[r,1] & targets$midpoint<distance_temp[r,2]),]
    targets_temp[r,] <- sum(targets_temp_a$CUB_dmel)
    }
  number_targets_R1 <- cbind(number_targets_R1,targets_temp)
}
colnames(number_targets_R1) <- as.character(dist_bp/1000000)
targets_R1 <- number_targets_R1
####### centiMorgans R2 ########
# centiMorgans. determine for each locus the number of bp per cM
sum_tot_left_chr_1<- as.data.frame(matrix(nrow = length(loc_chr_1),ncol=2 ))
sum_tot_right_chr_1<- as.data.frame(matrix(nrow = length(loc_chr_1),ncol=2 ))
loci_to_test <- round(loc_chr_1/10000,0)
for (locus_2 in loci_to_test) {
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
    sum_tot_left_chr_1[which(loci_to_test==locus_2),which(dist_cM==dist_var_2)] <- n_rows_left_chr_1
    sum_right_chr_1 <- cumsum(r_map[locus_2:nrow(r_map),])
    n_rows_right_chr_1 <- min(which(sum_right_chr_1 > dist_var_2))
    sum_tot_right_chr_1[which(loci_to_test==locus_2),which(dist_cM==dist_var_2)] <- n_rows_right_chr_1
  }
}
# multiplying by 10000 to obtain the real distance in bp
sum_tot_left_chr_1 <- sum_tot_left_chr_1*10000
sum_tot_right_chr_1 <- sum_tot_right_chr_1*10000
colnames(sum_tot_left_chr_1) <- as.character(dist_cM)
colnames(sum_tot_right_chr_1) <- as.character(dist_cM)
# non_synonymous. obtain for each locus the number of ns occurring within the distances to be tested
loc_chr_1_ns <- response_chr_1$position
number_ns_R2_chr_1 <- as.data.frame(matrix(nrow = length(loc_chr_1_ns)))
number_ns_R2_chr_1 <- number_ns_R2_chr_1[,-1]
for (dist_var_3 in as.character(dist_cM)) {
location_left_temp_chr_1 <- loc_chr_1_ns - sum_tot_left_chr_1[dist_var_3]
location_right_temp_chr_1 <- loc_chr_1_ns + sum_tot_right_chr_1[dist_var_3]
location_temp_chr_1 <- as.data.frame(cbind(location_left_temp_chr_1,location_right_temp_chr_1))
number_temp_chr_1 <- as.data.frame(matrix(nrow = length(loc_chr_1_ns)))
for (r in 1:nrow(location_temp_chr_1)) {
   number_temp_chr_1_temp <- targets[which(targets$midpoint>location_temp_chr_1[r,1] & targets$midpoint<location_temp_chr_1[r,2]),]
    number_temp_chr_1[r,] <- sum(number_temp_chr_1_temp$CUB_dmel)
}
number_ns_R2_chr_1 <- cbind(number_ns_R2_chr_1,number_temp_chr_1)
}
colnames(number_ns_R2_chr_1) <- as.character(dist_cM)
targets_R2 <- number_ns_R2_chr_1

  
final_response_R1 <- rbind(final_response_R1,response_chr_1_targets)
final_targets_R1 <- rbind(final_targets_R1, targets_R1)
final_recombination_R1 <- rbind(final_recombination_R1, cm_R1)

final_response_R2 <- rbind(final_response_R2,response_chr_1_targets)
final_targets_R2 <- rbind(final_targets_R2, targets_R2)
}
##########################
########## R1 ############
##########################
results_targets_R1 <-  as.data.frame(matrix(nrow = length(dist_bp)))
  results_targets_R1[1,1] <- glance(lm(final_response_R1$Fst~final_recombination_R1[,1]+final_targets_R1[,1]))$r.squared
  results_targets_R1[2,1] <- glance(lm(final_response_R1$Fst~final_recombination_R1[,2]+final_targets_R1[,2]))$r.squared
  results_targets_R1[3,1] <- glance(lm(final_response_R1$Fst~final_recombination_R1[,3]+final_targets_R1[,3]))$r.squared
  results_targets_R1[4,1] <- glance(lm(final_response_R1$Fst~final_recombination_R1[,4]+final_targets_R1[,4]))$r.squared
  results_targets_R1[5,1] <- glance(lm(final_response_R1$Fst~final_recombination_R1[,5]+final_targets_R1[,5]))$r.squared

results_targets_R1$distance <- dist_bp/1000000
results_targets_R1$test <- "Non-synonymous"
colnames(results_targets_R1) <-   c("het","distance","test")

results_targets_R1_pvalue <-  as.data.frame(matrix(nrow = length(dist_bp)))
  results_targets_R1_pvalue[1,1] <- glance(lm(final_response_R1$Fst~final_recombination_R1[,1]+final_targets_R1[,1]))$p.value
  results_targets_R1_pvalue[2,1] <- glance(lm(final_response_R1$Fst~final_recombination_R1[,2]+final_targets_R1[,2]))$p.value
  results_targets_R1_pvalue[3,1] <- glance(lm(final_response_R1$Fst~final_recombination_R1[,3]+final_targets_R1[,3]))$p.value
  results_targets_R1_pvalue[4,1] <- glance(lm(final_response_R1$Fst~final_recombination_R1[,4]+final_targets_R1[,4]))$p.value
  results_targets_R1_pvalue[5,1] <- glance(lm(final_response_R1$Fst~final_recombination_R1[,5]+final_targets_R1[,5]))$p.value

results_targets_R1$bonferroni <-  p.adjust(results_targets_R1_pvalue$V1,method="bonferroni") 
sig_R1 <- which(results_targets_R1$het==max(results_targets_R1$het))
sig_R1_distance <- results_targets_R1[sig_R1,"distance"]
regression_R1 <- as.data.frame(cbind(final_response_R1$Fst,final_targets_R1[,sig_R1],final_recombination_R1[,sig_R1]))
colnames(regression_R1) <- c("Fst","Targets","centiMorgans")

cM_plot_R1 <- 
ggplot(regression_R1)+
geom_smooth(aes(y=Fst,x=centiMorgans/10),color= "cyan2", fill="cyan2" ,method='lm', formula= y~x,size=2)+
geom_point(aes(y=Fst,x=centiMorgans/10),color="darkcyan",size=2)+
theme_bw(base_size = 10,base_family="Helvetica")+
labs(x = "centiMorgans (cM)",y= stat_label)+
labs(subtitle=paste("Size genomic neighborhood:\n",sig_R1_distance*2,"Mbp"))

targets_plot_R1 <- 
ggplot(regression_R1)+
geom_smooth(aes(y=Fst,x=Targets),color="chartreuse2",fill="chartreuse2",method='lm', formula= y~x,size=2)+
geom_point(aes(y=Fst,x=Targets),color="chartreuse4",size=2)+
theme_bw(base_size = 10,base_family="Helvetica")+
labs(x = "Targets of selection",y=stat_label)+
labs(subtitle=paste("Size genomic neighborhood:\n",sig_R1_distance*2,"Mbp"))

plot_R1 <- results_targets_R1[,c("distance","bonferroni","het")]
r_squared_R1 <-
ggplot(plot_R1) +
geom_line(aes(x=distance,y=het),size=1,color="blue") +
geom_point(aes(x=distance,y=het),size=2,color="gray30") +
geom_point(aes(x=distance[sig_R1], y=het[sig_R1]),size=3,color="deeppink") +
annotate("text",x=-Inf,y=Inf,hjust=-0.1,vjust=1,label=paste("adj.R^2 == ",round(plot_R1$het[sig_R1],3)),parse =TRUE)+
annotate("text",x=-Inf,y=Inf,hjust=-0.1,vjust=3,label=paste("P == ",round(plot_R1$bonferroni[sig_R1],3)),parse=TRUE)+
scale_y_continuous(expression(Adjusted~R^2))+
scale_x_continuous(name="Size genomic\nneighbourhood (Mbp)",breaks=seq(0,max(plot_R1$distance),1),labels = as.character(seq(0,max(plot_R1$distance),1)*2)) +
theme_bw(base_size = 10,base_family="Helvetica")+    
labs(subtitle = "Sets of multiple regression\nusing physical distance (bp)")+
theme(legend.title=element_blank())

print(r_squared_R1 + cM_plot_R1 + targets_plot_R1) 
ggsave(paste0("FLY_regression_",dispersal,"_",stat_label,"_",marker,"_R1.pdf"),  width =8, height =3, units = "in", dpi="retina", bg = "transparent" )
##########################
########## R2 ############
##########################
results_targets_R2 <-  as.data.frame(matrix(nrow = length(dist_cM)))
  results_targets_R2[1,1] <- glance(lm(final_response_R2$Fst~final_targets_R2[,1]))$r.squared
  results_targets_R2[2,1] <- glance(lm(final_response_R2$Fst~final_targets_R2[,2]))$r.squared
  results_targets_R2[3,1] <- glance(lm(final_response_R2$Fst~final_targets_R2[,3]))$r.squared
  results_targets_R2[4,1] <- glance(lm(final_response_R2$Fst~final_targets_R2[,4]))$r.squared
  results_targets_R2[5,1] <- glance(lm(final_response_R2$Fst~final_targets_R2[,5]))$r.squared

results_targets_R2$distance <- dist_cM
results_targets_R2$test <- "Non-synonymous"
colnames(results_targets_R2) <-   c("het","distance","test")

results_targets_R2_pvalue <-  as.data.frame(matrix(nrow = length(dist_cM)))
  results_targets_R2_pvalue[1,1] <- glance(lm(final_response_R2$Fst~final_targets_R2[,1]))$p.value
  results_targets_R2_pvalue[2,1] <- glance(lm(final_response_R2$Fst~final_targets_R2[,2]))$p.value
  results_targets_R2_pvalue[3,1] <- glance(lm(final_response_R2$Fst~final_targets_R2[,3]))$p.value
  results_targets_R2_pvalue[4,1] <- glance(lm(final_response_R2$Fst~final_targets_R2[,4]))$p.value
  results_targets_R2_pvalue[5,1] <- glance(lm(final_response_R2$Fst~final_targets_R2[,5]))$p.value
  
results_targets_R2$bonferroni <- p.adjust(results_targets_R2_pvalue$V1,method="bonferroni") 
sig_R2 <- which(results_targets_R2$het==max(results_targets_R2$het))
sig_R2_distance <- results_targets_R2[sig_R2,"distance"]
regression_R2 <- as.data.frame(cbind(final_response_R2$Fst,final_targets_R2[,sig_R2]))
colnames(regression_R2) <- c("Fst","Targets")

targets_plot_R2 <- 
ggplot(regression_R2)+
geom_smooth(aes(y=Fst,x=Targets),color="chartreuse2",fill="chartreuse2",method='lm', formula= y~x,size=2)+
geom_point(aes(y=Fst,x=Targets),color="chartreuse4",size=2)+
theme_bw(base_size = 10,base_family="Helvetica")+
labs(x = "Targets of selection",y=stat_label)+
labs(subtitle=paste("Size genomic neighborhood:\n",sig_R2_distance*2,"cM"))

plot_R2 <- results_targets_R2[,c("distance","bonferroni","het")]
r_squared_R2 <-
ggplot(plot_R2) +
geom_line(aes(x=distance, y=het),size=1,color="blue") +
geom_point(aes(x=distance, y=het),size=2,color="gray30") +
geom_point(aes(x=distance[sig_R2], y=het[sig_R2]),size=3,color="deeppink") +
annotate("text",x=-Inf,y=Inf,hjust=-0.1,vjust=1,label=paste("R^2 == ",round(plot_R2$het[sig_R2],3)),parse =TRUE)+
annotate("text",x=-Inf,y=Inf,hjust=-0.1,vjust=3,label=paste("P == ",round(plot_R2$bonferroni[sig_R2],3)),parse=TRUE)+
scale_y_continuous(expression(R^2))+
scale_x_continuous(name="Size genomic\nneighbourhood (cM)",breaks=seq(0,max(plot_R2$distance),2.5),labels = as.character(seq(0,max(plot_R2$distance),2.5)*2)) +
theme_bw(base_size = 10,base_family="Helvetica")+    
labs(subtitle = "Sets of simple regression\nusing genetic distance (cM)")+
theme(legend.title=element_blank())

print(r_squared_R2 + targets_plot_R2) 
ggsave(paste0("FLY_regression_",dispersal,"_",stat_label,"_",marker,"_R2.pdf"),  width = 5, height =3, units = "in", dpi="retina", bg = "transparent" )

}
