library(fields)
library(snpStats)
library(ggthemes)
library(ggplot2)
library(scales)
library(readr)
library(stringi)
source("hierfstat.R")
chromosome <- 18 # these are the chromosomes to analyse
map_resolution <- 100000

RecRates_All_Chromosomes <- read_csv("chill_recom_map.csv") 
RecRates_All_Chromosomes <- as.data.frame(RecRates_All_Chromosomes)
RecRates_All_Chromosomes$Chr <- as.character(RecRates_All_Chromosomes$Chr)
RecRates_All_Chromosomes_chr <- RecRates_All_Chromosomes[which(RecRates_All_Chromosomes$Chr==chromosome),]
chr_length <- RecRates_All_Chromosomes_chr[nrow(RecRates_All_Chromosomes_chr),"Location"]

resolution_stat <- 1000000
break_bins <- seq(0,chr_length+resolution_stat,resolution_stat)
break_bins_2 <- break_bins[-length(break_bins)]

map_cattle_binned <- stats.bin(RecRates_All_Chromosomes_chr$Location,RecRates_All_Chromosomes_chr$`mean sexes`,breaks = seq(0,chr_length,map_resolution))
map_cattle_binned_b <- unname(map_cattle_binned$stats[2,]* map_cattle_binned$stats[1,])
map_cattle_binned_b
map <- as.data.frame(map_cattle_binned_b)
colnames(map) <- "cM"
map[is.na(map$cM),] <- 0
map_df_temp <- as.data.frame(cbind(seq(100000,chr_length,100000), map))
colnames(map_df_temp) <- c("location","recombination")
map_df_temp_2 <- stats.bin(map_df_temp$location, map_df_temp$recombination, breaks = break_bins)
map_df <- unlist(unname(((map_df_temp_2[[3]][2,] * map_df_temp_2[[3]][1,]) * 100)))
map_df[is.na(map_df)] <- 0
#############################################################
################ LOADING DATA ###############################
#############################################################
# LODING DATA WITH PLINK
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
plot_het <- chill_df[which(chill_df$chromosome==chromosome),]
  
use <-  chill_data_plink$map$chromosome == chromosome
genotype_temp <- genotype_plink[,use]
genotype <- genotype_temp[,order(as.numeric(colnames(genotype_temp)))]
df_he_temp <- genotype@.Data
df_he_temp <- as.data.frame(matrix(as.double(df_he_temp),nrow = nrow(df_he_temp),ncol = ncol(df_he_temp) ))
df_he_temp[df_he_temp==1] <- 11
df_he_temp[df_he_temp==2] <- 12
df_he_temp[df_he_temp==3] <- 22
df_he_temp <- as.data.frame(cbind(1,df_he_temp))
df_he_res_temp <- basic.stats(df_he_temp)
df_he <- unname(unlist(df_he_res_temp$Hs[,1]))
  
plot_het$exp_het <- df_he

folders <- paste0(getwd(),"/",dir(path=getwd(),pattern="^chillingham_chr_18_selection_exons"))

variables_number <- 36
number_of_stats_calculated <- 25
path.folder_stats_average <- folders[1]    
files_ave <-paste0(path.folder_stats_average,"/",dir(path.folder_stats_average,pattern = "generations"))
pops_df_pop1 <- as.data.frame(matrix(ncol=3))
colnames(pops_df_pop1) <- c("loc_msats","V2","V3")
pops_df_pop2 <- as.data.frame(matrix(ncol=3))
colnames(pops_df_pop2) <- c("loc_msats","V2","V3")

loc_msats <-  seq(map_resolution/2,(nrow(map)*map_resolution),map_resolution)
get_number_loci <-  read.table(files_ave[1],header=F,sep=",",skip= variables_number,fill=T)
number_loci <- ncol(get_number_loci)
get_number_generations <- read.table(files_ave[1],header = F,sep=",",nrows = variables_number, colClasses = "character")
number_of_generations <- as.numeric(get_number_generations[7,2])
number_of_lines <- number_of_stats_calculated * number_of_generations
initial_breaks <- seq(1,number_of_lines,number_of_generations)
end_breaks <- c((initial_breaks-1)[-1],number_of_lines)

for(pop in 1:length(files_ave) ){
      pops_df_temp <- read.table(files_ave[pop],header=F,sep=",",skip= variables_number,fill=T)
      pops_het_pop1 <- pops_df_temp[initial_breaks[15]:end_breaks[15],1:number_loci]
      pops_df_pop1 <- rbind(pops_df_pop1,as.data.frame(cbind(loc_msats,unname(unlist(pops_het_pop1[number_of_generations,])),pop+0)))
      pops_het_pop2 <-  pops_df_temp[initial_breaks[16]:end_breaks[16],1:number_loci]
      pops_df_pop2 <- rbind(pops_df_pop2,as.data.frame(cbind(loc_msats,unname(unlist(pops_het_pop2[number_of_generations,])),pop+length(files_ave))))
  }

pops_final <- rbind(pops_df_pop1,pops_df_pop2)
pops_final$V3 <- as.factor(pops_final$V3)

targets <- read_csv("chill_targets_of_selection_exons.csv") 
targets$chr_name <- as.character(targets$chr_name)
targets <- targets[which(targets$chr_name==chromosome),]
targets <- targets[!duplicated(targets$start),]
targets <- targets[!duplicated(targets$end),]
targets_df <- targets
targets_df$midpoint <- (targets_df$start + targets_df$end)/2
targets_df$targets <- targets_df$ns 
# / 50
  # ceiling(targets_df$ns + (targets_df$s * targets_df$CUB_dmel))/75
targets_df_2 <- as.numeric(stats.bin(targets_df$midpoint, targets_df$targets, breaks = break_bins)[[3]][2,])
targets_df_2[is.na(targets_df_2)] <- 0 
# targets_df_2[which(targets_df_2>4)] <- 4

targets_df_2 <- rescale(targets_df_2,to=c(0,max(pops_final$V2,na.rm = T)))
map_df_2 <- rescale(map_df,to=c(0,max(pops_final$V2,na.rm = T)))

test_df <- as.data.frame(cbind(break_bins_2,map_df_2,targets_df_2))  
colnames(test_df) <- c("location","recombination","targets")
ticks_breaks <- seq(0,chr_length,5000000)
ticks_lab <- ticks_breaks/1000000

rescale_fun <- function(x){rescale(x,to=c(0,max(map_df,na.rm = T)))}

#### plot simulations
print(
ggplot(data=pops_final,aes(x=loc_msats,y=V2,fill=V3)) +
geom_area(alpha=0.02,position ="identity" ) +
geom_line(data=test_df,aes(x=location,y=recombination),color="darkturquoise",size=1.5,inherit.aes=F)+
geom_line(data=test_df,aes(x=location,y=targets),color="chartreuse3",size=1.5,inherit.aes = F)+
geom_smooth(data=pops_final,aes(x=loc_msats,y=V2),color="deeppink",fill="deeppink",method="loess",span=1/400,size=2,se=F)+
theme_tufte(base_family="Helvetica",base_size=16) +
scale_fill_manual(values=rep("black",length(files_ave)*2))+
labs(x="Chromosome position (Mbp)", y="He", title=NULL)+ 
scale_x_continuous(breaks=ticks_breaks, labels=ticks_lab) +
scale_y_continuous(sec.axis = sec_axis(rescale_fun, name = "centiMorgans/Mbp",breaks=seq(0,max(map_df,na.rm = T),0.5),labels=as.character(seq(0,max(map_df,na.rm = T),0.5))),breaks=seq(0,max(pops_final$V2,na.rm = T),0.1),labels=as.character(seq(0,max(pops_final$V2,na.rm = T),0.1))) +
theme(legend.position =  "none")
)

 ggsave(paste0("chillingham_chr_",chromosome,"_selection.pdf"),  width = 10, height =3, units = "in", dpi="retina", bg = "transparent" )
# plot real data
resolution_het <- 50000
brk <- seq(0,chr_length,resolution_het)
het_all_loci <- as.numeric(stats.bin(plot_het$position, plot_het$exp_het, breaks = brk)[[3]][2,])
het_all_loci[is.na(het_all_loci)] <- 0
plot_chill <- as.data.frame(cbind(seq(0,chr_length,resolution_het)[-1],het_all_loci))

print(
ggplot(plot_chill,aes(x=V1,y = het_all_loci)) +
geom_area(fill ="black",position ="identity" ) +
geom_line(data=test_df,aes(x=location,y=recombination),color="darkturquoise",size=1.5,inherit.aes=F)+
geom_line(data=test_df,aes(x=location,y=targets),color="chartreuse3",size=1.5,inherit.aes = F)+
geom_smooth(data=plot_chill,aes(x=V1,y=het_all_loci),color="deeppink",fill="deeppink",method="loess",span=1/4,size=2,se=F)+
theme_tufte(base_family="Helvetica",base_size=16) +
labs(x="Chromosome position (Mbp)", y="He", title=NULL)+ 
scale_x_continuous(breaks=ticks_breaks, labels=ticks_lab) +
scale_y_continuous(sec.axis = sec_axis(rescale_fun, name = "centiMorgans/Mbp",breaks=seq(0,max(map_df,na.rm = T),0.5),labels=as.character(seq(0,max(map_df,na.rm = T),0.5))),breaks=seq(0,max(pops_final$V2,na.rm = T),0.1),labels=as.character(seq(0,max(pops_final$V2,na.rm = T),0.1))) +
theme(legend.position =  "none")
)

 ggsave(paste0("chillingham_chr_",chromosome,"_real_data.pdf"),  width = 10, height =3, units = "in", dpi="retina", bg = "transparent" )



