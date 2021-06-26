library(ggplot2)
library(ggthemes)
library(snpStats)
library(readr)
library(ggpointdensity)
library(viridis)
library(ggnewscale)
library(patchwork)


variables_number <- 36
number_of_stats_calculated <-25
number_loci <- 50
number_of_generations <- 62
number_of_lines <- number_of_stats_calculated * number_of_generations
initial_breaks <- seq(1,number_of_lines,number_of_generations)
end_breaks <- c((initial_breaks-1)[-1],number_of_lines)

path.folder_sim_neutral <- paste0(getwd(),"/plots_bill_test")
files_neutral <- paste0(path.folder_sim_neutral,"/",dir(path.folder_sim_neutral,pattern = "^generations"))
neutral_loci_fst <-  vector()
neutral_loci_he_pop1 <-   vector()
neutral_loci_he_pop2 <-  vector()
neutral_loci_shua <-  vector()
neutral_loci_sha_pop1 <-  vector()
neutral_loci_sha_pop2 <- vector()

for(df_gen in 1:length(files_neutral)){
mean_df_neutral <- read.table(files_neutral[df_gen],header=F,sep=",",skip= variables_number,fill=T)
mean_sha_pop1_neutral <-  mean_df_neutral[initial_breaks[9]:end_breaks[9],1:number_loci]
mean_sha_pop2_neutral <-  mean_df_neutral[initial_breaks[10]:end_breaks[10],1:number_loci]
mean_shua_neutral <-  mean_df_neutral[initial_breaks[12]:end_breaks[12],1:number_loci]
mean_het_pop1_neutral <-  mean_df_neutral[initial_breaks[15]:end_breaks[15],1:number_loci]
mean_het_pop2_neutral <-  mean_df_neutral[initial_breaks[16]:end_breaks[16],1:number_loci]
mean_Fst_neutral <-  mean_df_neutral[initial_breaks[18]:end_breaks[18],1:number_loci]
mean_Fstp_neutral <-  mean_df_neutral[initial_breaks[19]:end_breaks[19],1:number_loci]

neutral_loci_fst <-  c(neutral_loci_fst,unlist(unname(mean_Fstp_neutral[number_of_generations,])))
neutral_loci_he_pop1 <- c(neutral_loci_he_pop1, unlist(unname(mean_het_pop1_neutral[number_of_generations,])))
neutral_loci_he_pop2 <- c(neutral_loci_he_pop2, unlist(unname(mean_het_pop2_neutral[number_of_generations,])))

neutral_loci_shua <- c(neutral_loci_shua, unlist(unname(mean_shua_neutral[number_of_generations,])))
neutral_loci_sha_pop1 <- c(neutral_loci_sha_pop1,  unlist(unname(mean_sha_pop1_neutral[number_of_generations,])))
neutral_loci_sha_pop2 <-  c(neutral_loci_sha_pop2, unlist(unname(mean_sha_pop2_neutral[number_of_generations,])))

}
neutral_loci_fst_he_pop1 <- as.data.frame(cbind(neutral_loci_fst,neutral_loci_he_pop1))
colnames(neutral_loci_fst_he_pop1) <- c("fst","he")
neutral_loci_fst_he_pop2 <- as.data.frame(cbind(neutral_loci_fst,neutral_loci_he_pop2))
colnames(neutral_loci_fst_he_pop2) <- c("fst","he")
neutral_loci_fst_he <- rbind(neutral_loci_fst_he_pop1,neutral_loci_fst_he_pop2)
# neutral_loci_fst_he <- neutral_loci_fst_he[complete.cases(neutral_loci_fst_he),]
neutral_loci_shua_sha_pop1 <- as.data.frame(cbind(neutral_loci_shua,neutral_loci_sha_pop1))
colnames(neutral_loci_shua_sha_pop1) <- c("MI","shannon")
neutral_loci_shua_sha_pop2 <- as.data.frame(cbind(neutral_loci_shua,neutral_loci_sha_pop2))
colnames(neutral_loci_shua_sha_pop2) <- c("MI","shannon")
neutral_loci_shua_sha <- rbind(neutral_loci_shua_sha_pop1,neutral_loci_shua_sha_pop2)
# neutral_loci_shua_sha <- neutral_loci_shua_sha[complete.cases(neutral_loci_shua_sha),]

path.folder_sim_selection <- paste0(getwd(),"/plots_bill_test")
files_selection <- paste0(path.folder_sim_selection,"/",dir(path.folder_sim_selection,pattern = "^snps"))
selection_loci_fst <- vector()
selection_loci_exp_he_pop1 <- vector()
selection_loci_exp_he_pop2 <- vector()
selection_loci_shua <-  vector()
selection_loci_sha_pop1 <-  vector()
selection_loci_sha_pop2 <-  vector()

for(df_snps in 1:length(files_selection)){
mean_df_selection <- read.table(files_selection[df_snps],header=T,sep=",",skip= variables_number,fill=T)
mean_sha_pop1_selection <-  mean_df_selection$shannon_pop1
mean_sha_pop2_selection <-  mean_df_selection$shannon_pop2
mean_shua_selection <-  mean_df_selection$MI
mean_het_pop1_selection <-  mean_df_selection$P.AB_pop1
mean_het_pop2_selection <-  mean_df_selection$P.AB_pop2
mean_exp_het_pop1_selection <-  2*mean_df_selection$MAF_pop1*(1-mean_df_selection$MAF_pop1)
mean_exp_het_pop2_selection <-  2*mean_df_selection$MAF_pop2*(1-mean_df_selection$MAF_pop2)
mean_Fst_selection <-  mean_df_selection$Fst

selection_loci_fst <- c(selection_loci_fst, mean_Fst_selection)
selection_loci_exp_he_pop1 <- c(selection_loci_exp_he_pop1, mean_exp_het_pop1_selection)
selection_loci_exp_he_pop2 <-  c(selection_loci_exp_he_pop2, mean_exp_het_pop2_selection)

selection_loci_shua <-  c(selection_loci_shua,mean_shua_selection)
selection_loci_sha_pop1 <-  c(selection_loci_sha_pop1,mean_sha_pop1_selection)
selection_loci_sha_pop2 <-  c(selection_loci_sha_pop2,mean_sha_pop2_selection)

}

selection_loci_fst_he_pop1 <- as.data.frame(cbind(selection_loci_fst,selection_loci_exp_he_pop1))
colnames(selection_loci_fst_he_pop1) <- c("fst","he")
selection_loci_fst_he_pop2 <- as.data.frame(cbind(selection_loci_fst,selection_loci_exp_he_pop2))
colnames(selection_loci_fst_he_pop2) <- c("fst","he")
selection_loci_fst_he <- rbind(selection_loci_fst_he_pop1,selection_loci_fst_he_pop2)
# selection_loci_fst_he <- selection_loci_fst_he[complete.cases(selection_loci_fst_he),]
selection_loci_fst_he <- selection_loci_fst_he[sample(1:nrow(selection_loci_fst_he),nrow(neutral_loci_fst_he)),]

selection_loci_shua_sha_pop1 <- as.data.frame(cbind(selection_loci_shua,selection_loci_sha_pop1))
colnames(selection_loci_shua_sha_pop1) <- c("MI","shannon")
selection_loci_shua_sha_pop2 <- as.data.frame(cbind(selection_loci_shua,selection_loci_sha_pop2))
colnames(selection_loci_shua_sha_pop2) <- c("MI","shannon")
selection_loci_shua_sha <- rbind(selection_loci_shua_sha_pop1,selection_loci_shua_sha_pop2)
# selection_loci_shua_sha <- selection_loci_shua_sha[complete.cases(selection_loci_shua_sha),]
selection_loci_shua_sha <- selection_loci_shua_sha[sample(1:nrow(selection_loci_shua_sha),nrow(neutral_loci_shua_sha)),]

number_transfers <- 1
transfer_each_gen <- 10
population_size_dispersal <- 10
Ne_fst <- 10
Ne <- 15
dispersal_rate <- (number_transfers / transfer_each_gen) / (population_size_dispersal)
Fst_expected <- 1 / ((4 * Ne_fst * dispersal_rate) * ((2 / (2 - 1)) ^ 2) + 1)
shua_expected <- (0.22 / (sqrt(2 * Ne_fst * dispersal_rate))) - (0.69 / ((2*Ne_fst) * sqrt(dispersal_rate)))
rate_of_loss <- 1 - (1 / (2 * Ne))

cols_fst <- c("Expected Fst"="red","Trend line (loess)"="darkcyan")

neutral_fst <- 
print(
ggplot(data = neutral_loci_fst_he)+
geom_pointdensity(aes(he,fst),size=4)+
scale_color_viridis_c(option = "C",name="Loci\ndensity") +
new_scale_color() + 
scale_colour_manual(name="",values=cols_fst) +
geom_hline(aes(yintercept=Fst_expected,color="Expected Fst"),size=2)+
geom_smooth(aes(x=he,y=fst,color="Trend line (loess)"),fill="darkcyan",size=2,method = "loess")+
guides(color=guide_legend(override.aes=list(fill=NA))) +
labs(x = "Heterozygosity",y="Fst",title = "Neutral loci Fst/Heterozygosity")+
xlim(0,1)+
ylim(0,1)+
theme_bw(base_size = 16,base_family="Helvetica")
)

selection_fst <- 
print(
ggplot(data = selection_loci_fst_he)+
geom_pointdensity(aes(he,fst),size=4)+
scale_color_viridis_c(option = "C",name="Loci\ndensity") +
new_scale_color() + 
scale_colour_manual(name="",values=cols_fst) +
geom_hline(aes(yintercept=Fst_expected,color="Expected Fst"),size=2)+
geom_smooth(aes(x=he,y=fst,color="Trend line (loess)"),fill="darkcyan",size=2,method = "loess")+
guides(color=guide_legend(override.aes=list(fill=NA))) +
labs(x = "Heterozygosity",y="Fst",title="Loci under selection Fst/Heterozygosity")+
xlim(0,1)+
ylim(0,1)+
theme_bw(base_size = 16,base_family="Helvetica")
)

cols_shua <- c("Expected MI"="red","Trend line (loess)"="darkcyan")

neutral_shua <- 
print(
ggplot(data = neutral_loci_shua_sha)+
geom_pointdensity(aes(shannon,MI),size=4)+
scale_color_viridis_c(option = "C",name="Loci\ndensity") +
new_scale_color() + 
scale_colour_manual(name="",values=cols_shua) +
geom_hline(aes(yintercept=shua_expected,color="Expected MI"),size=2)+
geom_smooth(aes(x=shannon,y=MI,color="Trend line (loess)"),fill="darkcyan",size=2,method = "loess")+
guides(color=guide_legend(override.aes=list(fill=NA))) +
labs(x = "Shannon index",y="MI",title = "Neutral loci MI/Shannon")+
xlim(0,1)+
ylim(0,1)+
theme_bw(base_size = 16,base_family="Helvetica")
)

selection_shua <- 
print(
ggplot(data = selection_loci_shua_sha)+
geom_pointdensity(aes(shannon,MI),size=4)+
scale_color_viridis_c(option = "C",name="Loci\ndensity") +
new_scale_color() + 
scale_colour_manual(name="",values=cols_shua) +
geom_hline(aes(yintercept=shua_expected,color="Expected MI"),size=2)+
geom_smooth(aes(x=shannon,y=MI,color="Trend line (loess)"),fill="darkcyan",size=2,method = "loess")+
guides(color=guide_legend(override.aes=list(fill=NA))) +
labs(x = "Shannon index",y="MI",title = "Loci under selection MI/Shannon")+
xlim(0,1)+
ylim(0,1)+
theme_bw(base_size = 16,base_family="Helvetica")
)


 print(
     (neutral_fst + selection_fst)
     /
       (neutral_shua + selection_shua)  + plot_annotation(tag_levels = 'a',title = 'Simulations WITHOUT selection',theme = theme(plot.title = element_text(size = 30)))
 )

 # ggsave("beta_vs_alpha_neutral_2.pdf",  width = 12, height = 8, units = "in", dpi="retina" )


test_fst <-  neutral_loci_fst_he_pop1
test_fst$ratio <- (test_fst$fst*6)/mean(test_fst$fst,na.rm=T)
plot(density(test_fst$fst,na.rm = T))
plot(density(FST.empDF_res$beta,na.rm = T))

plot(test_fst$he, test_fst$fst)
plot(FST.empDF_res$alpha,FST.empDF_res$beta)

mean(test_fst$fst,na.rm=T)
mean(FST.empDF_res$beta,na.rm=T)

