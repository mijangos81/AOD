library(fields)
library(scales)
library(ggplot2)
library(ggthemes)
library(readr)
chr_fly <- "2L"
map_resolution <- 100000
variables_number <- 36
map_fly <- read_csv("fly_recom_map.csv")
map_fly$Chr <- as.character(map_fly$Chr)
map_fly <- map_fly[which(map_fly$Chr==chr_fly),]
map_fly <- as.data.frame(map_fly$cM/1000)
colnames(map_fly) <- "cM"
map_fly[is.na(map_fly$cM),] <- 0
chr_length_fly <- (nrow(map_fly)+1) * map_resolution
number_loci <- 243
number_of_generations <- 34
number_of_stats_calculated <- 25
number_of_lines <- number_of_stats_calculated * number_of_generations
initial_breaks <- seq(1,number_of_lines,number_of_generations)
end_breaks <- c((initial_breaks-1)[-1],number_of_lines)

path.folder_sim_neutral <- paste0(getwd(),"/fig_sim_vs_real_data_neutral")
files_neutral <- paste0(path.folder_sim_neutral,"/",dir(path.folder_sim_neutral,pattern = "^final_stats_average"))
mean_df_neutral <- read.table(files_neutral[1],header=F,sep=",",skip= variables_number,fill=T)
mean_AFD_neutral <- mean_df_neutral[initial_breaks[1]:end_breaks[1],1:number_loci]
mean_sha_diff_neutral <- mean_df_neutral[initial_breaks[2]:end_breaks[2],1:number_loci]
mean_het_pop1_neutral <-  mean_df_neutral[initial_breaks[15]:end_breaks[15],1:number_loci]
mean_het_pop2_neutral <-  mean_df_neutral[initial_breaks[16]:end_breaks[16],1:number_loci]
mean_Fst_neutral <-  mean_df_neutral[initial_breaks[18]:end_breaks[18],1:number_loci]
mean_Fstp_neutral <-  mean_df_neutral[initial_breaks[19]:end_breaks[19],1:number_loci]
mean_Dest_neutral <-  mean_df_neutral[initial_breaks[20]:end_breaks[20],1:number_loci]

path.folder_sim_selection <- paste0(getwd(),"/fig_sim_vs_real_data_selection")
files_selection <- paste0(path.folder_sim_selection,"/",dir(path.folder_sim_selection,pattern = "^final_stats_average"))
mean_df_selection <- read.table(files_selection[11],header=F,sep=",",skip= variables_number,fill=T)
mean_AFD_selection <- mean_df_selection[initial_breaks[1]:end_breaks[1],1:number_loci]
mean_sha_diff_selection <- mean_df_selection[initial_breaks[2]:end_breaks[2],1:number_loci]
mean_het_pop1_selection <-  mean_df_selection[initial_breaks[15]:end_breaks[15],1:number_loci]
mean_het_pop2_selection <-  mean_df_selection[initial_breaks[16]:end_breaks[16],1:number_loci]
mean_Fst_selection <-  mean_df_selection[initial_breaks[18]:end_breaks[18],1:number_loci]
mean_Fstp_selection <-  mean_df_selection[initial_breaks[19]:end_breaks[19],1:number_loci]
mean_Dest_selection <-  mean_df_selection[initial_breaks[20]:end_breaks[20],1:number_loci]

resolution_df <- 100000
break_bins <- seq(0,chr_length_fly-resolution_df,resolution_df)
break_bins_2 <- break_bins[-1]

location_msats_experiment <-  c(8325394,2373262,4781938,4960235,7040194,8348440,11015745,12507696, 13153885, 14615705, 14995570,20706003)
location_neutral_loci_analysis <- c(seq(map_resolution/2,(nrow(map_fly)*map_resolution),map_resolution),location_msats_experiment)
loc_exp_loci <- location_neutral_loci_analysis
loc_exp_loci <- loc_exp_loci[order(loc_exp_loci)]
loc_exp_loci <- which(loc_exp_loci %% 10000 != 0)   

map_df_temp <- as.data.frame(cbind(seq(100000,chr_length_fly-100000,100000), map_fly))
colnames(map_df_temp) <- c("location","recombination")
map_df_temp_2 <- stats.bin(map_df_temp$location, map_df_temp$recombination, breaks = break_bins)
map_df <- unlist(unname(((map_df_temp_2[[3]][2,] * map_df_temp_2[[3]][1,]) * 100)))
map_df[is.na(map_df)] <- 0

targets <- read.csv("fly_targets_of_selection.csv") 
targets$chr_name <- as.character(targets$chr_name)
targets <- targets[which(targets$chr_name==chr_fly),]
targets <- targets[!duplicated(targets$start),]
targets <- targets[!duplicated(targets$end),]
targets_df <- targets
targets_df$midpoint <- (targets_df$start + targets_df$end)/2
targets_df$targets <- ceiling(targets_df$ns + (targets_df$s * targets_df$CUB_dmel))/10
targets_df_2 <- as.numeric(stats.bin(targets_df$midpoint, targets_df$targets, breaks = break_bins)[[3]][2,])
targets_df_2 <- stats.bin(targets_df$midpoint, targets_df$targets, breaks = break_bins)
targets_df_2 <- targets_df_2[[3]][2,] * targets_df_2[[3]][1,]
# targets_df_2[which(targets_df_2>7)] <- 7

# stat_df_neutral <- as.data.frame(cbind(location_neutral_loci_analysis,unlist(unname(mean_Fstp_neutral[34,]))))
stat_df_neutral_temp_pop1 <- unlist(unname(mean_het_pop1_neutral[34,]))
stat_df_neutral_temp_pop2 <- unlist(unname(mean_het_pop2_neutral[34,]))
stat_df_neutral_temp_mean <- as.data.frame(cbind(stat_df_neutral_temp_pop1,stat_df_neutral_temp_pop2))
stat_df_neutral_temp_mean <- rowMeans(stat_df_neutral_temp_mean)
stat_df_neutral <- as.data.frame(cbind(location_neutral_loci_analysis,stat_df_neutral_temp_mean))
colnames(stat_df_neutral) <- c("location","stat")

stat_df_neutral <- stat_df_neutral[-loc_exp_loci,]
stat_df_neutral <- as.numeric(stats.bin(stat_df_neutral$location,stat_df_neutral$stat, breaks = break_bins)[[3]][2,])

# stat_df_selection <- as.data.frame(cbind(location_neutral_loci_analysis,unlist(unname(mean_Fstp_selection[34,]))))
stat_df_selection_temp_pop1 <- unlist(unname(mean_het_pop1_selection[34,]))
stat_df_selection_temp_pop2 <- unlist(unname(mean_het_pop2_selection[34,]))
stat_df_selection_temp_mean <- as.data.frame(cbind(stat_df_selection_temp_pop1,stat_df_selection_temp_pop2))
stat_df_selection_temp_mean <- rowMeans(stat_df_selection_temp_mean)
stat_df_selection <- as.data.frame(cbind(location_neutral_loci_analysis,stat_df_selection_temp_mean))
colnames(stat_df_selection) <- c("location","stat")

stat_df_selection <- stat_df_selection[-loc_exp_loci,]
stat_df_selection <- as.numeric(stats.bin(stat_df_selection$location,stat_df_selection$stat, breaks = break_bins)[[3]][2,])

exp_loci <- read.csv("response_variables.csv")
exp_loci <- exp_loci[which(exp_loci$chromosome==chr_fly),]
# exp_loci <- exp_loci[,c("loc_BDGP6","Fstp")]
exp_loci <- exp_loci[,c("loc_BDGP6","het_mean")]

# map_df_2 <- rescale(map_df,to=c(0, max(stat_df_neutral)))
map_df_2 <- rescale(map_df,to=c(0, max(stat_df_selection,na.rm = T)))

test_df <- as.data.frame(cbind(break_bins_2,map_df_2,targets_df_2,stat_df_neutral,stat_df_selection))  
colnames(test_df) <- c("location","recombination","targets","stat_neutral","stat_selection")
# test_df$targets <- rescale(test_df$targets,to=c(0, max(stat_df_neutral)))
test_df$targets <- rescale(test_df$targets,to=c(0, max(stat_df_selection,na.rm = T)))

ticks_breaks <- seq(0,chr_length_fly,5000000)
ticks_lab <- ticks_breaks/1000000

rescale_fun <- function(x){rescale(x,to=c(0,max(map_df,na.rm = T)))}

print(
ggplot(test_df) +
geom_area(aes(x=location,y=stat_selection),fill="gray25",color="black",size=1,alpha=1) +
geom_area(aes(x=location,y=stat_neutral),fill="gray75",color="black",size=1,alpha=1) +

geom_col(data=exp_loci,aes(x=loc_BDGP6,y=het_mean),fill="deeppink",color="deeppink",size=2)+
# geom_col(data=exp_loci,aes(x=loc_BDGP6,y=Fstp),fill="deeppink",color="deeppink",size=2)+
geom_line(aes(x=location,y=targets),color="chartreuse3",size=2,alpha= 1)+
geom_line(aes(x=location,y=recombination),color="darkturquoise",size=2,alpha=1)+
theme_tufte(base_family="Helvetica",base_size=18) +
# scale_y_continuous(sec.axis = sec_axis(rescale_fun, name = "centiMorgans/Mbp",breaks=seq(0,max(map_df,na.rm = T),0.5),labels=as.character(seq(0,max(map_df,na.rm = T),0.5))),breaks=seq(0, max(stat_df_neutral),0.1),labels=as.character(seq(0, max(test_df$stat_neutral),0.1))) +
scale_y_continuous(sec.axis = sec_axis(rescale_fun, name = "centiMorgans/Mbp",breaks=seq(0,max(map_df,na.rm = T),0.5),labels=as.character(seq(0,max(map_df,na.rm = T),0.5))),breaks=seq(0, max(test_df$stat_selection,na.rm = T),0.1),labels=as.character(seq(0, max(test_df$stat_selection,na.rm = T),0.1))) +
xlab("Genome location (Mbp)") +
# ylab("Fst")+
ylab("He")+
scale_x_continuous(breaks=ticks_breaks, labels=ticks_lab) 
 )

# ggsave(paste("fst_real_vs_all_sim_fly_chr_",chr_fly,".pdf"),  width = 10, height =5, units = "in", dpi="retina", bg = "transparent" )

ggsave(paste("He_real_vs_all_sim_fly_chr_",chr_fly,".pdf"),  width = 10, height =5, units = "in", dpi="retina", bg = "transparent" )

#################################################################
########### SCHEMATIC FIGURE REGRESSION APPROACH ################
#################################################################

targets <- read.csv("fly_targets_of_selection.csv") 
targets$chr_name <- as.character(targets$chr_name)
targets <- targets[which(targets$chr_name==chr_fly),]
targets <- targets[!duplicated(targets$start),]
targets <- targets[!duplicated(targets$end),]
targets_df <- targets
targets_df$midpoint <- (targets_df$start + targets_df$end)/2
targets_df$targets <- ceiling(targets_df$ns + (targets_df$s * targets_df$CUB_dmel))
targets_df_2 <- stats.bin(targets_df$midpoint, targets_df$targets, breaks = break_bins)
targets_df_2 <- targets_df_2[[3]][2,] * targets_df_2[[3]][1,]
# targets_df_2[which(targets_df_2>7)] <- 7
ticks_breaks <- seq(0,chr_length_fly,1000000)
ticks_lab <- ticks_breaks/1000000
rescale_fun <- function(x){rescale(x,to=c(0,max(map_df,na.rm = T)))}
map_df_2 <- rescale(map_df,to=c(0, max(targets_df_2,na.rm = T)))
test_df <- as.data.frame(cbind(break_bins_2,map_df_2,targets_df_2,stat_df_neutral,stat_df_selection))  
colnames(test_df) <- c("location","recombination","targets","stat_neutral","stat_selection")

colors_plot <- c("Recombination rate"="darkturquoise","Targets of selection"="chartreuse3","Genotyped locus"="deeppink")

print(
ggplot(test_df[20:80,]) +
geom_line(aes(x=location,y=targets,color="Targets of selection"),size=2,alpha= 1)+
geom_line(aes(x=location,y=recombination,color="Recombination rate"),size=2,alpha=1)+
geom_vline(aes(xintercept=5000000,color="Genotyped locus"),size=3)+
theme_tufte(base_family="Helvetica",base_size=14) +
scale_y_continuous(sec.axis = sec_axis(rescale_fun, name = "Recombination rate (cM)",breaks=seq(0,max(map_df,na.rm = T),0.5),labels=as.character(seq(0,max(map_df,na.rm = T),0.5))),breaks=seq(0, max(test_df$targets,na.rm = T),500),labels=as.character(seq(0, max(test_df$targets,na.rm = T),500))) +
scale_colour_manual(name="",values=colors_plot) +
xlab("Genome location (Mbp)") +
ylab("Number of targets of selection")+
  theme(legend.position = "bottom")+
scale_x_continuous(breaks=ticks_breaks, labels=ticks_lab) 
 )

ggsave(paste("schema_bp_regression_approach_",chr_fly,".pdf"),  width = 9, height =5, units = "in", dpi="retina", bg = "transparent" )

map_df_temp <- as.data.frame(cbind(seq(100000,chr_length_fly-100000,100000), map_fly))
colnames(map_df_temp) <- c("location","recombination")
map_schema <- map_df_temp
map_schema$accum <- cumsum(map_schema[,"recombination"]) *100
test_df <- as.data.frame(cbind(map_schema$accum,targets_df_2))  
colnames(test_df) <- c("location","targets")


print(
ggplot(test_df[20:80,]) +
geom_line(aes(x=location,y=targets,color="Targets of selection"),size=2,alpha= 1)+
geom_vline(aes(xintercept=15.998,color="Genotyped locus"),size=3)+
theme_tufte(base_family="Helvetica",base_size=14) +
scale_y_continuous(breaks=seq(0, max(test_df$targets,na.rm = T),500),labels=as.character(seq(0, max(test_df$targets,na.rm = T),500))) +
theme(legend.position = "bottom")+
  scale_colour_manual(name="",values=colors_plot) +
xlab("Genome location (cM)") +
ylab("Number of targets of selection")
)

ggsave(paste("schema_cM_regression_approach_",chr_fly,".pdf"),  width = 9, height =5, units = "in", dpi="retina", bg = "transparent" )
