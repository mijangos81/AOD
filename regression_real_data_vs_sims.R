library(ggplot2)
library(ggthemes)
library(snpStats)
library(readr)
variables_number <- 36
number_of_stats_calculated <-25
chr_fly <- "2L"
number_loci <- 242
number_of_generations <- 34
number_of_lines <- number_of_stats_calculated * number_of_generations
initial_breaks <- seq(1,number_of_lines,number_of_generations)
end_breaks <- c((initial_breaks-1)[-1],number_of_lines)

path.folder_sim_neutral <- paste0(getwd(),"/fig_sim_vs_real_data_neutral")
files_neutral <- paste0(path.folder_sim_neutral,"/",dir(path.folder_sim_neutral,pattern = "^final_stats_average"))
mean_df_neutral <- read.table(files_neutral[1],header=F,sep=",",skip= variables_number,fill=T)
mean_AFD_neutral <- mean_df_neutral[initial_breaks[1]:end_breaks[1],1:number_loci]
mean_sha_diff_neutral <- mean_df_neutral[initial_breaks[2]:end_breaks[2],1:number_loci]
mean_expected_het_pop1_neutral <-  mean_df_neutral[initial_breaks[13]:end_breaks[13],1:number_loci]
mean_expected_het_pop2_neutral <-  mean_df_neutral[initial_breaks[14]:end_breaks[14],1:number_loci]
mean_het_pop1_neutral <-  mean_df_neutral[initial_breaks[15]:end_breaks[15],1:number_loci]
mean_het_pop2_neutral <-  mean_df_neutral[initial_breaks[16]:end_breaks[16],1:number_loci]
mean_Fst_neutral <-  mean_df_neutral[initial_breaks[18]:end_breaks[18],1:number_loci]
mean_Fstp_neutral <-  mean_df_neutral[initial_breaks[19]:end_breaks[19],1:number_loci]
mean_Dest_neutral <-  mean_df_neutral[initial_breaks[20]:end_breaks[20],1:number_loci]

path.folder_sim_selection <- paste0(getwd(),"/fig_sim_vs_real_data_selection")
files_selection <- paste0(path.folder_sim_selection,"/",dir(path.folder_sim_selection,pattern = "^final_stats_average"))
mean_df_selection <- read.table(files_selection[i],header=F,sep=",",skip= variables_number,fill=T)
mean_AFD_selection <- mean_df_selection[initial_breaks[1]:end_breaks[1],1:number_loci]
mean_sha_diff_selection <- mean_df_selection[initial_breaks[2]:end_breaks[2],1:number_loci]
mean_expected_het_pop1_selection <-  mean_df_selection[initial_breaks[13]:end_breaks[13],1:number_loci]
mean_expected_het_pop2_selection <-  mean_df_selection[initial_breaks[14]:end_breaks[14],1:number_loci]
mean_het_pop1_selection <-  mean_df_selection[initial_breaks[15]:end_breaks[15],1:number_loci]
mean_het_pop2_selection <-  mean_df_selection[initial_breaks[16]:end_breaks[16],1:number_loci]
mean_Fst_selection <-  mean_df_selection[initial_breaks[18]:end_breaks[18],1:number_loci]
mean_Fstp_selection <-  mean_df_selection[initial_breaks[19]:end_breaks[19],1:number_loci]
mean_Dest_selection <-  mean_df_selection[initial_breaks[20]:end_breaks[20],1:number_loci]

# loci_2L_experiment <- read.csv("stats_all_together.csv")
loci_2L_experiment <- read.csv("response_variables.csv")
# loci_2L_experiment <- loci_2L_experiment[which(loci_2L_experiment$marker=="msat"),]
# loci_2L_experiment <- loci_2L_experiment[which(loci_2L_experiment$marker=="snp"),]
# this snp by mistake is not included due to a mistake in setting the allele frequency in the final simulations
# loci_2L_experiment <- loci_2L_experiment[-which(loci_2L_experiment$loci=="Chr2-30"),]
# loci_2L_experiment <- loci_2L_experiment[complete.cases(loci_2L_experiment),]
loci_2L_experiment <- loci_2L_experiment[which(loci_2L_experiment$chromosome=="2L"),]
                                               # & loci_2L_experiment$marker=="snp" ),]
loci_2L_experiment <- loci_2L_experiment[order(loci_2L_experiment$loc_BDGP6),]
loci_2L_experiment$window_2 <-  loc_exp_loci - 2
  # loci_2L_experiment$window_1
loci_2L_experiment <- loci_2L_experiment[-2,]
  # loci_2L_experiment[-c(4,7,10,11),]

fst_loci_2L_sim_neutral <- unname(unlist(mean_Fst_neutral[34,loci_2L_experiment$window_2]))
fst_loci_2L_sim_selection <- unname(unlist(mean_Fst_selection[34,loci_2L_experiment$window_2]))

He_loci_2L_sim_neutral_pop1_T0 <- unname(unlist(mean_het_pop1_neutral[1,loci_2L_experiment$window_2]))
He_loci_2L_sim_neutral_pop2_T0 <- unname(unlist(mean_het_pop2_neutral[1,loci_2L_experiment$window_2]))

He_loci_2L_sim_neutral_pop1_T2 <- unname(unlist(mean_het_pop1_neutral[34,loci_2L_experiment$window_2]))
He_loci_2L_sim_neutral_pop2_T2 <- unname(unlist(mean_het_pop2_neutral[34,loci_2L_experiment$window_2]))

He_loci_2L_sim_neutral_pop1_expected <- unname(unlist(mean_expected_het_pop1_neutral[34,loci_2L_experiment$window_2]))
He_loci_2L_sim_neutral_pop2_expected <- unname(unlist(mean_expected_het_pop2_neutral[34,loci_2L_experiment$window_2]))

He_loci_2L_sim_neutral_pop1_bias <- (He_loci_2L_sim_neutral_pop1_T2 - He_loci_2L_sim_neutral_pop1_expected)
He_loci_2L_sim_neutral_pop2_bias <- (He_loci_2L_sim_neutral_pop2_T2 - He_loci_2L_sim_neutral_pop2_expected)

He_loci_2L_sim_neutral_pop1_T0_T2 <-  He_loci_2L_sim_neutral_pop1_T0 - He_loci_2L_sim_neutral_pop1_T2
He_loci_2L_sim_neutral_pop2_T0_T2 <-  He_loci_2L_sim_neutral_pop2_T0 - He_loci_2L_sim_neutral_pop2_T2

He_loci_2L_sim_neutral_pops_bias_mean <- rowMeans(as.data.frame(cbind(He_loci_2L_sim_neutral_pop1_bias,He_loci_2L_sim_neutral_pop2_bias)))
He_loci_2L_sim_neutral_pops_T0_T2_mean <- rowMeans(as.data.frame(cbind(He_loci_2L_sim_neutral_pop1_T0_T2,He_loci_2L_sim_neutral_pop2_T0_T2)))
He_loci_2L_sim_neutral_het_mean <- rowMeans(as.data.frame(cbind(He_loci_2L_sim_neutral_pop1_T2,He_loci_2L_sim_neutral_pop2_T2)))

He_loci_2L_sim_selection_pop1_T0 <- unname(unlist(mean_het_pop1_selection[1,loci_2L_experiment$window_2]))
He_loci_2L_sim_selection_pop2_T0 <- unname(unlist(mean_het_pop2_selection[1,loci_2L_experiment$window_2]))

He_loci_2L_sim_selection_pop1_T2 <- unname(unlist(mean_het_pop1_selection[34,loci_2L_experiment$window_2]))
He_loci_2L_sim_selection_pop2_T2 <- unname(unlist(mean_het_pop2_selection[34,loci_2L_experiment$window_2]))

He_loci_2L_sim_selection_pop1_expected <- unname(unlist(mean_expected_het_pop1_selection[34,loci_2L_experiment$window_2]))
He_loci_2L_sim_selection_pop2_expected <- unname(unlist(mean_expected_het_pop2_selection[34,loci_2L_experiment$window_2]))

He_loci_2L_sim_selection_pop1_bias <- (He_loci_2L_sim_selection_pop1_T2 - He_loci_2L_sim_selection_pop1_expected)
He_loci_2L_sim_selection_pop2_bias <- (He_loci_2L_sim_selection_pop2_T2 - He_loci_2L_sim_selection_pop2_expected)

He_loci_2L_sim_selection_pop1_T0_T2 <-  He_loci_2L_sim_selection_pop1_T0 - He_loci_2L_sim_selection_pop1_T2
He_loci_2L_sim_selection_pop2_T0_T2 <-  He_loci_2L_sim_selection_pop2_T0 - He_loci_2L_sim_selection_pop2_T2

He_loci_2L_sim_selection_pops_bias_mean <- rowMeans(as.data.frame(cbind(He_loci_2L_sim_selection_pop1_bias,He_loci_2L_sim_selection_pop2_bias)))
He_loci_2L_sim_selection_pops_T0_T2_mean <- rowMeans(as.data.frame(cbind(He_loci_2L_sim_selection_pop1_T0_T2,He_loci_2L_sim_selection_pop2_T0_T2)))
He_loci_2L_sim_selection_het_mean <- rowMeans(as.data.frame(cbind(He_loci_2L_sim_selection_pop1_T2,He_loci_2L_sim_selection_pop2_T2)))

plot_2L_loci_regression <- as.data.frame(cbind(loci_2L_experiment$Fst,fst_loci_2L_sim_neutral,fst_loci_2L_sim_selection))
# plot_2L_loci_regression <- as.data.frame(cbind(loci_2L_experiment$het_T0_T2,He_loci_2L_sim_neutral_pops_T0_T2_mean,He_loci_2L_sim_selection_pops_T0_T2_mean))
# # plot_2L_loci_regression <- as.data.frame(cbind(loci_2L_experiment$het_mean,He_loci_2L_sim_neutral_het_mean,He_loci_2L_sim_selection_het_mean))

colnames(plot_2L_loci_regression) <- c("experiment","neutral_sims","selection_sims")

print(
  ggplot(plot_2L_loci_regression)+
    geom_smooth(aes(y=experiment,x=selection_sims),fill="deeppink",color="gray20",method='lm', formula= y~x,size=2)+
  geom_point(aes(y=experiment,x=selection_sims),color="deeppink",alpha=1,size=2)+
  # geom_point(aes(y=experiment,x=neutral_sims),color="darkturquoise",alpha=1,size=2)+
  # geom_smooth(aes(y=experiment,x=neutral_sims),fill="darkturquoise",color="black",method='lm', formula= y~x)+
  labs(x = "",y="")+
theme_bw(base_size = 16,base_family="Helvetica")
)

 r_squared_AOD_fly <- summary(lm(plot_2L_loci_regression$experiment~plot_2L_loci_regression$selection_sims))
 print(r_squared_AOD_fly$r.squared  )

 ggsave(paste("Fst_real_vs_sim_AOD_fly_chr_",chr_fly,".pdf"),  width = 4, height =4, units = "in", dpi="retina", bg = "transparent" )
 
  # ggsave(paste("Fst_real_vs_sim_AOD_fly_chr_",chr_fly,".pdf"),  width = 4, height =4, units = "in", dpi="retina", bg = "transparent" )
 
print(
  ggplot(plot_2L_loci_regression)+
  geom_point(aes(y=experiment,x=neutral_sims),color="darkslategrey",alpha=1,size=2)+
  geom_smooth(aes(y=experiment,x=neutral_sims),fill="darkslategrey",color="gray20",method='lm', formula= y~x,size=2)+
  labs(x = "",y="")+
theme_bw(base_size = 16,base_family="Helvetica")
)
r_squared_neutral_fly <- summary(lm(plot_2L_loci_regression$experiment~plot_2L_loci_regression$neutral_sims))

 # ggsave(paste("Fst_real_vs_sim_neutral_fly_chr_",chr_fly,".pdf"),  width = 4, height =4, units = "in", dpi="retina", bg = "transparent" )
 
  ggsave(paste("He_real_vs_sim_neutral_fly_chr_",chr_fly,".pdf"),  width = 4, height =4, units = "in", dpi="retina", bg = "transparent" )
plot_sim <- cbind(1:34,mean_Fst_selection[,loci_2L_experiment$window_2])
colnames(plot_sim) <- as.character(1:ncol(plot_sim))
plot_sim_2 <- reshape2::melt(plot_sim,id="1")
plot_sim_2$variable <- as.numeric(as.character(plot_sim_2$variable))
plot_2L_loci_regression$selection_sims <- as.numeric(as.character(plot_2L_loci_regression$selection_sims))

  ggplot() +
          geom_line(data=plot_sim_2,aes(x=plot_sim_2$"1",y=value,group=variable,colour=rainbow(nrow(plot_sim_2)))) +
           geom_point(data=plot_2L_loci_regression,aes(x = 34,y = selection_sims),colour=rainbow(11))+
# scale_color_viridis(discrete = TRUE) +
    theme_bw(base_size = 14) +
    ylim(c(0.13,0.2))+
    xlim(c(25,34))+
         theme(legend.title=element_blank())+
         theme(legend.position="none") +
         xlab("GENERATIONS") +
         ylab("FSTp")
###################################
######## CHILLINGHAM #############
###################################
chromosome_name <- 18
chrom <- chromosome_name
map_resolution <- 100000
resolution_het_chill <- 1000000
brk_chill <- seq(0,chr_length_chill,resolution_het_chill)

RecRates_All_Chromosomes <- read_csv("chill_recom_map.csv") 
RecRates_All_Chromosomes <- as.data.frame(RecRates_All_Chromosomes)
RecRates_All_Chromosomes$Chr <- as.character(RecRates_All_Chromosomes$Chr)
RecRates_All_Chromosomes_chr <- RecRates_All_Chromosomes[which(RecRates_All_Chromosomes$Chr==chromosome_name),]
chr_length_chill <- RecRates_All_Chromosomes_chr[nrow(RecRates_All_Chromosomes_chr),"Location"]
map_cattle_binned <- stats.bin(RecRates_All_Chromosomes_chr$Location,RecRates_All_Chromosomes_chr$`mean sexes`,breaks = seq(0,chr_length_chill,map_resolution))
map_cattle_binned_b <- unname(map_cattle_binned$stats[2,]* map_cattle_binned$stats[1,])
map_chill <- as.data.frame(map_cattle_binned_b)
colnames(map_chill) <- "cM"
map_chill[is.na(map_chill$cM),] <- 0

path.folder_sim_neutral_chill <- paste0(getwd(),"/chillingham_chr_18_no_selection")
files_neutral_chill <- paste0(path.folder_sim_neutral_chill,"/",dir(path.folder_sim_neutral_chill,pattern = "^final_stats_average"))
mean_df_neutral_chill <- read.table(files_neutral_chill[1],header=F,sep=",",skip= variables_number,fill=T)

number_loci_chill <-ncol(mean_df_neutral_chill)
loc_snps <- seq(map_resolution/2,(nrow(map_chill)*map_resolution),map_resolution)
get_number_generations <- read.table(files_neutral_chill[1],header = F,sep=",",nrows = variables_number, colClasses = "character")
number_of_generations_chill <- as.numeric(get_number_generations[7,2])
number_of_lines <- number_of_stats_calculated * number_of_generations_chill
initial_breaks <- seq(1,number_of_lines,number_of_generations_chill)
end_breaks <- c((initial_breaks-1)[-1],number_of_lines)

mean_het_pop1_neutral_chill <-  mean_df_neutral_chill[initial_breaks[15]:end_breaks[15],1:number_loci_chill]
mean_het_pop2_neutral_chill <-  mean_df_neutral_chill[initial_breaks[16]:end_breaks[16],1:number_loci_chill]

loci_het_pop1_neutral_chill <- unname(unlist(mean_het_pop1_neutral_chill[number_of_generations_chill,]))
loci_het_pop2_neutral_chill <- unname(unlist(mean_het_pop2_neutral_chill[number_of_generations_chill,]))

het_pops_chill_neutral <- as.data.frame(cbind(loci_het_pop1_neutral_chill,loci_het_pop2_neutral_chill))
het_pops_chill_neutral$mean <- rowMeans(het_pops_chill_neutral)
het_pops_chill_neutral$loc <- seq(50000,chr_length_chill-map_resolution,map_resolution)
final_het_pops_chill_neutral <- as.numeric(stats.bin(het_pops_chill_neutral$loc, het_pops_chill_neutral$mean, breaks = brk_chill)[[3]][2,])

path.folder_sim_selection_chill <- paste0(getwd(),"/chillingham_chr_18_selection")
files_selection_chill <- paste0(path.folder_sim_selection_chill,"/",dir(path.folder_sim_selection_chill,pattern = "^final_stats_average"))
mean_df_selection_chill <- read.table(files_selection_chill[1],header=F,sep=",",skip= variables_number,fill=T)

mean_het_pop1_selection_chill <-  mean_df_selection_chill[initial_breaks[15]:end_breaks[15],1:number_loci_chill]
mean_het_pop2_selection_chill <-  mean_df_selection_chill[initial_breaks[16]:end_breaks[16],1:number_loci_chill]

loci_het_pop1_selection_chill <- unname(unlist(mean_het_pop2_selection_chill[number_of_generations_chill,]))
loci_het_pop2_selection_chill <- unname(unlist(mean_het_pop2_selection_chill[number_of_generations_chill,]))

het_pops_chill_selection <- as.data.frame(cbind(loci_het_pop1_selection_chill,loci_het_pop2_selection_chill))
het_pops_chill_selection$mean <- rowMeans(het_pops_chill_selection)
het_pops_chill_selection$loc <- seq(50000,chr_length_chill-map_resolution,map_resolution)
final_het_pops_chill_selection <- as.numeric(stats.bin(het_pops_chill_selection$loc, het_pops_chill_selection$mean, breaks = brk_chill)[[3]][2,])

chill_data_plink <- read.plink(bed= "res_chill.bed" , bim="res_chill.bim" , fam="res_chill.fam" )
# the function read.plink does not read the map file properly, so it is attached separately
map_plink <- read_table2("res_chill.bim", col_names = FALSE)
colnames(map_plink) <- c("chromosome", "snp.name", "cM", "position", "allele.1", "allele.2")
chill_data_plink$map <- map_plink
genotype_plink <- chill_data_plink$genotypes # this is the input for pairwise LD analyses
colnames(genotype_plink) <- chill_data_plink$map$position
snpsum.col_plink<- col.summary(genotype_plink)
chill_df_plink <- as.data.frame(cbind(snpsum.col_plink,chill_data_plink$map))
chill_df <- chill_df_plink
chr_1 <- chill_df[which(chill_df$chromosome==chrom),]
chr_1$exp_het <- chr_1$RAF * chr_1$MAF * 2

het_all_loci_chill <- as.numeric(stats.bin(chr_1$position, chr_1$exp_het, breaks = brk_chill)[[3]][2,])
het_all_loci_chill[is.na(het_all_loci_chill)] <- 0

plot_18_loci_regression <- as.data.frame(cbind(het_all_loci_chill,final_het_pops_chill_neutral,final_het_pops_chill_selection))
colnames(plot_18_loci_regression) <- c("real_data","neutral_chill_sims","selection_chill_sims")

print(
  ggplot(plot_18_loci_regression)+
  geom_point(aes(y=real_data,x=selection_chill_sims),color="deeppink",alpha=1/2,size=2)+
  geom_smooth(aes(y=real_data,x=selection_chill_sims),fill="deeppink",color="gray20",method='lm', formula= y~x,size=2)+
  labs(x = "",y="")+
theme_bw(base_size = 16,base_family="Helvetica")
)

r_squared_AOD_chill <- summary(lm(plot_18_loci_regression$real_data~plot_18_loci_regression$selection_chill_sims))

ggsave(paste("real_vs_sim_AOD_chill_chr_",chr_fly,".pdf"),  width = 4, height =4, units = "in", dpi="retina", bg = "transparent" )
 
print(
  ggplot(plot_18_loci_regression)+
   geom_point(aes(y=real_data,x=neutral_chill_sims),color="darkslategrey",alpha=1/2,size=2)+
  geom_smooth(aes(y=real_data,x=neutral_chill_sims),fill="darkslategrey",color="gray20",method='lm', formula= y~x,size=2)+
  labs(x = "",y="")+
theme_bw(base_size = 16,base_family="Helvetica")
)

r_squared_neutral_chill <- summary(lm(plot_18_loci_regression$real_data~plot_18_loci_regression$neutral_chill_sims))

ggsave(paste("real_vs_sim_neutral_chill_chr_",chr_fly,".pdf"),  width = 4, height =4, units = "in", dpi="retina", bg = "transparent" )
 

