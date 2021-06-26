library(adegenet)
library(ggpointdensity)
library(ggplot2)
library(ggforce)
library(ggalt)
library(concaveman)
library(ggthemes)
library(viridis)
library(ggnewscale)
library(data.table)
source("basic_stats_fix.R")

snps <- F
msats <-T

# which sampling time to analyse?
  T2 <-T
  T0 <- F
  
# which dispersal regime to analyse?
  dis <- "Low"   # options: "Hig","Med", "Low"
  
if(snps==T){# Import genepop file to genind object 
all_snps_auto<-read.genepop(file = "every_snp.gen")
}
if(msats==T){# Import genepop file to genind object 
all_snps_auto<-read.genepop(file = "all_msats_auto.gen")
}

# assign pop information to strata
strata(all_snps_auto) <- data.frame(pop = pop(all_snps_auto))
# split strata listing the name of each hierarchy level 
splitStrata(all_snps_auto)<-~time/l_pair/rep/dis/b_pair/pop1_2/p_pair
# Visualize the summary of all strata fields
# lapply(strata(all_snps_auto),summary)
if (T2==T) {
  # # separating T2
   test_to_keep<-all_snps_auto$strata$time==2
}
if (T0==T) {
# separating T0
  test_to_keep<-all_snps_auto$strata$time==0
}

t2_snps_auto<-all_snps_auto[test_to_keep]

test_to_keep_dispersal <- t2_snps_auto$strata$dis==dis
t2_snps_auto <- t2_snps_auto[test_to_keep_dispersal]

# separate populations based in the population pair
per_pop_snps<-seppop(t2_snps_auto,t2_snps_auto$strata$p_pair, drop=TRUE)

per_pop_snps_exp_temp <-  lapply(per_pop_snps,genind2hierfstat_fix)
per_pop_snps_exp_temp <- lapply(per_pop_snps_exp_temp,wc_2)
per_pop_snps_exp_res <- as.data.frame(rbindlist(per_pop_snps_exp_temp))

per_pop_msats_exp_temp <-  lapply(per_pop_snps,genind2hierfstat_fix)
per_pop_msats_exp_temp <- lapply(per_pop_msats_exp_temp,wc_2)
per_pop_msats_exp_res <- as.data.frame(rbindlist(per_pop_msats_exp_temp))

    neutral_fdist_ten_demes <- read.table("fdist_neutral_ten_demes.sims", header=TRUE)
    neutral_fdist_two_demes <- read.table("fdist_neutral_two_demes.sims", header=TRUE)

    path.folder_sim_US <- "fdist_neutral_msats"
    neutral_mine_msats_files_US <- dir(path.folder_sim_US,pattern = paste0("^outflank_snps"))
    neutral_mine_msats_temp_US <- lapply(paste0(path.folder_sim_US,"/",neutral_mine_msats_files_US),read.table,header=T,sep=",",skip=nrow(variables_file))
    neutral_mine_msats_res_US <- do.call(rbind, neutral_mine_msats_temp_US)
    neutral_mine_msats_res_US <- neutral_mine_msats_res_US[,c("He","FST")]
    
    path.folder_sim_US <- "fdist_selection_msats"
    selection_mine_msats_files_US <- dir(path.folder_sim_US,pattern = paste0("^outflank_snps"))
    selection_mine_msats_temp_US <- lapply(paste0(path.folder_sim_US,"/",selection_mine_msats_files_US),read.table,header=T,sep=",",skip=nrow(variables_file))
    selection_mine_msats_res_US <- do.call(rbind, selection_mine_msats_temp_US)
    selection_mine_msats_res_US <- selection_mine_msats_res_US[,c("He","FST")]
    
    path.folder_sim_US <- "fdist_neutral_snps"
    neutral_mine_snps_files_US <- dir(path.folder_sim_US,pattern = paste0("^outflank_snps"))
    neutral_mine_snps_temp_US <- lapply(paste0(path.folder_sim_US,"/",neutral_mine_snps_files_US),read.table,header=T,sep=",",skip=nrow(variables_file))
    neutral_mine_snps_res_US <- do.call(rbind, neutral_mine_snps_temp_US)
    neutral_mine_snps_res_US <- neutral_mine_snps_res_US[,c("He","FST")]
    
    path.folder_sim_US <- "fdist_selection_snps"
    selection_mine_snps_files_US <- dir(path.folder_sim_US,pattern = paste0("^outflank_snps"))
    selection_mine_snps_temp_US <- lapply(paste0(path.folder_sim_US,"/",selection_mine_snps_files_US),read.table,header=T,sep=",",skip=nrow(variables_file))
    selection_mine_snps_res_US <- do.call(rbind, selection_mine_snps_temp_US)
    selection_mine_snps_res_US <- selection_mine_snps_res_US[,c("He","FST")]
    
    path.folder_sim_NS <- "fdist_neutral_msats"
    neutral_mine_msats_files_NS <- dir(path.folder_sim_NS,pattern = paste0("^outflank_msats"))
    neutral_mine_msats_temp_NS <- lapply(paste0(path.folder_sim_NS,"/",neutral_mine_msats_files_NS),read.table,header=T,sep=",",skip=nrow(variables_file))
    neutral_mine_msats_res_NS <- do.call(rbind, neutral_mine_msats_temp_NS)
    neutral_mine_msats_res_NS <- neutral_mine_msats_res_NS[,c("He","FST")]
    
    path.folder_sim_NS <- "fdist_selection_msats"
    selection_mine_msats_files_NS <- dir(path.folder_sim_NS,pattern = paste0("^outflank_msats"))
    selection_mine_msats_temp_NS <- lapply(paste0(path.folder_sim_NS,"/",selection_mine_msats_files_NS),read.table,header=T,sep=",",skip=nrow(variables_file))
    selection_mine_msats_res_NS <- do.call(rbind, selection_mine_msats_temp_NS)
    selection_mine_msats_res_NS <- selection_mine_msats_res_NS[,c("He","FST")]
    
    path.folder_sim_NS <- "fdist_neutral_snps"
    neutral_mine_snps_files_NS <- dir(path.folder_sim_NS,pattern = paste0("^outflank_msats"))
    neutral_mine_snps_temp_NS <- lapply(paste0(path.folder_sim_NS,"/",neutral_mine_snps_files_NS),read.table,header=T,sep=",",skip=nrow(variables_file))
    neutral_mine_snps_res_NS <- do.call(rbind, neutral_mine_snps_temp_NS)
    neutral_mine_snps_res_NS <- neutral_mine_snps_res_NS[,c("He","FST")]
    
    path.folder_sim_NS <- "fdist_selection_snps"
    selection_mine_snps_files_NS <- dir(path.folder_sim_NS,pattern = paste0("^outflank_msats"))
    selection_mine_snps_temp_NS <- lapply(paste0(path.folder_sim_NS,"/",selection_mine_snps_files_NS),read.table,header=T,sep=",",skip=nrow(variables_file))
    selection_mine_snps_res_NS <- do.call(rbind, selection_mine_snps_temp_NS)
    selection_mine_snps_res_NS <- selection_mine_snps_res_NS[,c("He","FST")]
    
colors_plot <- c(
  "\nNeutral loci\n"="green",
  "\nMean FST\nneutral loci\n"="green3",
  "\nLoci under\nselection\n"="deepskyblue",
  "\nMean FST loci\nunder selection\n"="deepskyblue3",
  "\nExperiment\n"="deeppink",
  "\nMean FST\nexperiment\n"="deeppink3",
  "\nNeutral\nenvelope\n"="black")
################### Msats ####################    
betas.sim <- neutral_mine_msats_res_NS
betas.sim <- betas.sim[,c("He","FST")] 
colnames(betas.sim) <- c("alpha", "beta")
betas.sim <- betas.sim[complete.cases(betas.sim),]
maxalpha <- max(c(betas.sim$alpha))
minalpha <- min(c(betas.sim$alpha))
maxbeta <- max(c(betas.sim$beta))
minbeta <- min(c(betas.sim$beta))
alpha.cell <- 0.01
beta.cell <- 0.01
xseq<-seq(minalpha,maxalpha,alpha.cell)
yseq<-seq(minbeta,maxbeta,beta.cell )
envelope_temp <-  as.data.frame(get.CI.alpha(xseq=xseq, percentile=0.975, betas.sim=betas.sim))
envelope <- as.data.frame(rbind(cbind(xseq,envelope_temp$V1),cbind(xseq,envelope_temp$V2)))
envelope <- envelope[complete.cases(envelope$V2),]
envelope[which(envelope$V2<0),"V2"] <- 0
# envelope[which(envelope$xseq>0.47),"xseq"] <- 0.5
###################################################
##"Neutral simulations microsatellites (4 alleles) two demes"
###################################################
ggplot()+
  geom_point(data=neutral_mine_msats_res_NS,aes(x=He,y=FST,colour="\nNeutral loci\n"),size=4)+
  geom_point(data=neutral_mine_msats_res_US,aes(x=He,y=FST,colour="\nLoci under\nselection\n"),size=2)+
  geom_mark_hull(data=envelope,aes(x=xseq,y=V2,color="\nNeutral\nenvelope\n"),fill="black",size=1,alpha=1/7,concavity=2,expand=0,radius=0)+
  geom_point(data=per_pop_msats_exp_res,aes(x=He,y=FST,color="\nExperiment\n"),size=2)+
  geom_hline(data=neutral_mine_msats_res_NS,aes(yintercept = mean(neutral_mine_msats_res_NS$FST),color="\nMean FST loci\nunder selection\n"),size=1)+
  geom_hline(data=neutral_mine_msats_res_US,aes(yintercept = mean(neutral_mine_msats_res_US$FST),color="\nMean FST\nneutral loci\n"),size=1)+
  geom_hline(data=per_pop_msats_exp_res,aes(yintercept = mean(per_pop_msats_exp_res$FST),color="\nMean FST\nexperiment\n"),size=1)+
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_colour_manual(name="",values=colors_plot) +
  xlab("Total heterozygosity") + 
  ylab("FST")+
  ylim(c(0,1))+
  ggtitle("Neutral simulations microsatellites\n(4 alleles) in two populations")+
  theme_tufte(base_family="Helvetica")+
  theme(
     # legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA))

 ggsave("neutral_msats_two_demes.pdf",  width = 6, height = 6, units = "in", dpi="retina", bg = "transparent"  )
###################################################
##"AOD simulations microsatellites (4 alleles) two demes"
###################################################
# ggplot()+
#   geom_point(data=selection_mine_msats_res_NS,aes(x=He,y=FST,colour="\nNeutral loci\n"),size=4)+
#   geom_point(data=selection_mine_msats_res_US,aes(x=He,y=FST,colour="\nLoci under\nselection\n"),size=2)+
#   geom_mark_hull(data = envelope,concavity = 2,expand=0,radius=0,aes(x=xseq,y=V2,color="\nNeutral\nenvelope\n"),fill="black",size=1,alpha=1/7)+
#   geom_point(data=per_pop_msats_exp_res,aes(x=He,y=FST,color="\nExperiment\n"),size=2)+
#   geom_hline(data=selection_mine_msats_res_NS,aes(yintercept = mean(selection_mine_msats_res_NS$FST),color="\nMean FST loci\nunder selection\n"),size=1)+
#   geom_hline(data=selection_mine_msats_res_US,aes(yintercept = mean(selection_mine_msats_res_US$FST),color="\nMean FST\nneutral loci\n"),size=1)+
#   geom_hline(data=per_pop_msats_exp_res,aes(yintercept = mean(per_pop_msats_exp_res$FST),color="\nMean FST\nexperiment\n"),size=1)+
#   guides(color=guide_legend(override.aes=list(fill=NA))) +
#   scale_colour_manual(name="",values=colors_plot) +
#   xlab("Total heterozygosity") + 
#   ylab("FST")+
#   ylim(c(0,1))+
#   ggtitle("AOD simulations microsatellites\n(4 alleles) in two populations")+
#   theme_tufte(base_family="Helvetica")+
#   theme(
#      # legend.position="bottom",
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_rect(fill = "transparent"), # bg of the panel
#     plot.background = element_rect(fill = "transparent", color = NA))
# 
#  ggsave("AOD_msats_two_demes.pdf",  width = 6, height = 6, units = "in", dpi="retina", bg = "transparent"  )
 
 
 colors_plot_2 <- c(
  "\nNeutral loci\n"="green",
  "\nMean FST\nneutral loci\n"="green3",
  "\nLoci under\nselection\n"="deepskyblue",
  "\nExpected\nneutral FST\n"="black",
  "\nMean FST\n"="blue",
  "\nExperiment\n"="deeppink",
  "\nMean FST\nexperiment\n"="deeppink3",
  "\nNeutral\nenvelope\n"="green3")

ggplot()+
  # geom_pointdensity(data=selection_mine_snps_res_NS,aes(x=He,y=FST,colour="\nNeutral loci\n"),size=4)+
  #   geom_pointdensity(data=selection_mine_msats_res_NS,aes(x=He,y=FST),size=5)+
  #   scale_color_viridis(option ="B",name="Number of\nneighboring loci" )+
  # new_scale_color()+
  stat_density(data=selection_mine_msats_res_NS,aes(FST),adjust = 2,fill="blue",alpha=1/2)

  geom_point(data=selection_mine_msats_res_NS,aes(x=He,y=FST,color="\nSimulations\n"),size=5,alpha=1/10)+
  
  # geom_point(data=selection_mine_snps_res_US,aes(x=He,y=FST,colour="\nLoci under\nselection\n"),size=2)+
  # geom_point(data=per_pop_snps_exp_res,aes(x=He,y=FST,color="\nExperiment\n"),size=2)+
  geom_hline(aes(yintercept = mean(selection_mine_msats_res_NS$FST),color="\nMean FST\n"),size=2,key_glyph = "rect")+
  geom_hline(aes(yintercept = Fst_expected,color="\nExpected\nneutral FST\n"),size=2,linetype = "dashed")+
    # geom_hline(data=selection_mine_snps_res_US,aes(yintercept = mean(selection_mine_snps_res_US$FST),color="\nMean FST\nneutral loci\n"),size=1)+

  # geom_hline(data=per_pop_snps_exp_res,aes(yintercept = mean(per_pop_snps_exp_res$FST),color="\nMean FST\nexperiment\n"),size=1)+
    geom_mark_hull(data = envelope,concavity = 3,expand=0,radius=0,aes(x=xseq,y=V2,color="\nNeutral\nenvelope\n"),size=2)+
  # guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_colour_manual(name="",values=colors_plot_2) +
  # guides(color=guide_legend(ncol=2))+
  xlab("Total heterozygosity") + 
  ylab("FST")+
  ylim(c(0,1))+
  # ggtitle("AOD simulations SNPs\nin two populations")+
  theme_tufte(base_family="Helvetica")+
  theme(  
    # legend.position="bottom",
      text = element_text(size = 16)
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(),
        # panel.background = element_rect(fill = "transparent"), # bg of the panel
    # plot.background = element_rect(fill = "transparent", color = NA)
    )
  # theme(
  #    # legend.position="bottom",
  #       panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank(),
  #       panel.background = element_rect(fill = "transparent"), # bg of the panel
  #   plot.background = element_rect(fill = "transparent", color = NA))

 ggsave("AOD_msats_two_demes.pdf",  width = 6, height = 5, units = "in", dpi="retina", bg = "transparent"  )
 
 colors_plot_3 <- c(
  "\nNeutral loci\n"="green",
  "\nMean FST\nneutral loci\n"="green3",
  "\nLoci under\nselection\n"="deepskyblue",
  "\nExpected\nneutral FST\n"="black",
  "\nMean FST\n"="blue",
  "\nFly\nexperiment\n"="deeppink",
  "\nMean FST\nexperiment\n"="deeppink3",
  "\nNeutral\nenvelope\n"="green3")

ggplot()+
  # geom_pointdensity(data=selection_mine_snps_res_NS,aes(x=He,y=FST,colour="\nNeutral loci\n"),size=4)+
    # geom_pointdensity(data=selection_mine_snps_res_NS,aes(x=He,y=FST),size=5)+
    # scale_color_viridis(option ="B",name="Number of\nneighboring loci" )+
  # new_scale_color()+
  # geom_point(data=selection_mine_snps_res_US,aes(x=He,y=FST,colour="\nLoci under\nselection\n"),size=2)+
  geom_point(data=per_pop_msats_exp_res,aes(x=He,y=FST,color="\nFly\nexperiment\n"),size=5,alpha=1/2,key_glyph = draw_key_point)+
  geom_hline(aes(yintercept = mean(per_pop_msats_exp_res$FST),color="\nMean FST\n"),size=2,key_glyph = draw_key_smooth)+
  geom_hline(aes(yintercept = Fst_expected,color="\nExpected\nneutral FST\n"),size=2,linetype = "dashed")+
    # geom_hline(data=selection_mine_snps_res_US,aes(yintercept = mean(selection_mine_snps_res_US$FST),color="\nMean FST\nneutral loci\n"),size=1)+

  # geom_hline(data=per_pop_snps_exp_res,aes(yintercept = mean(per_pop_snps_exp_res$FST),color="\nMean FST\nexperiment\n"),size=1)+
    geom_mark_hull(data = envelope,concavity = 4,expand=0,radius=0,aes(x=xseq,y=V2,color="\nNeutral\nenvelope\n"),size=2,key_glyph = draw_key_smooth)+
  # guides(fill=guide_legend(override.aes=list(color=NA))) +
  scale_colour_manual(name="",values=colors_plot_3) +
  # guides(color=guide_legend(ncol=2))+
  xlab("Total heterozygosity") + 
  ylab("FST")+
  ylim(c(0,1))+
  # ggtitle("AOD simulations SNPs\nin two populations")+
  theme_tufte(base_family="Helvetica")+
  theme(  
    # legend.position="bottom",
      text = element_text(size = 16)
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(),
        # panel.background = element_rect(fill = "transparent"), # bg of the panel
    # plot.background = element_rect(fill = "transparent", color = NA)
    )
  # theme(
  #    # legend.position="bottom",
  #       panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank(),
  #       panel.background = element_rect(fill = "transparent"), # bg of the panel
  #   plot.background = element_rect(fill = "transparent", color = NA))

 ggsave("experiment_msats_two_demes.pdf",  width = 6, height = 5, units = "in", dpi="retina", bg = "transparent"  )
 
 
 
###################################################
###################################################
##################SNPS##########################
###################################################
betas.sim <- neutral_mine_snps_res_NS
betas.sim <- betas.sim[,c("He","FST")] 
colnames(betas.sim) <- c("alpha", "beta")
betas.sim <- betas.sim[complete.cases(betas.sim),]
maxalpha <- max(c(betas.sim$alpha))
minalpha <- min(c(betas.sim$alpha))
maxbeta <- max(c(betas.sim$beta))
minbeta <- min(c(betas.sim$beta))
alpha.cell <- 0.01
beta.cell <- 0.01
xseq<-seq(minalpha,maxalpha,alpha.cell)
yseq<-seq(minbeta,maxbeta,beta.cell )
envelope_temp <-  as.data.frame(get.CI.alpha(xseq=xseq, percentile=0.975, betas.sim=betas.sim))
envelope <- as.data.frame(rbind(cbind(xseq,envelope_temp$V1),cbind(xseq,envelope_temp$V2)))
envelope <- envelope[complete.cases(envelope$V2),]
envelope[which(envelope$V2<0),"V2"] <- 0
envelope[which(envelope$xseq>0.47),"xseq"] <- 0.5
###################################################
##"AOD simulations SNPs two demes"
###################################################
colors_plot_2 <- c(
  "\nNeutral loci\n"="green",
  "\nMean FST\nneutral loci\n"="green3",
  "\nLoci under\nselection\n"="deepskyblue",
  "\nExpected\nneutral FST\n"="black",
  "\nMean FST\n"="blue",
  "\nSimulations\n"="deeppink",
  "\nMean FST\nexperiment\n"="deeppink3",
  "\nNeutral\nenvelope\n"="deepskyblue")

ggplot()+
  # geom_pointdensity(data=selection_mine_snps_res_NS,aes(x=He,y=FST,colour="\nNeutral loci\n"),size=4)+
  #   geom_pointdensity(data=selection_mine_snps_res_NS,aes(x=He,y=FST),size=5)+
  #   scale_color_viridis(option ="B",name="Number of\nneighboring loci" )+
  # new_scale_color()+
  stat_density(data=selection_mine_snps_res_NS,aes(FST),adjust = 1/2)
  geom_point(data=selection_mine_snps_res_NS,aes(x=He,y=FST,color="\nSimulations\n"),size=5,alpha=1/10)+
  
  # geom_point(data=selection_mine_snps_res_US,aes(x=He,y=FST,colour="\nLoci under\nselection\n"),size=2)+
  # geom_point(data=per_pop_snps_exp_res,aes(x=He,y=FST,color="\nExperiment\n"),size=2)+
  geom_hline(aes(yintercept = mean(selection_mine_snps_res_NS$FST),color="\nMean FST\n"),size=2,key_glyph = "rect")+
  geom_hline(aes(yintercept = Fst_expected,color="\nExpected\nneutral FST\n"),size=2,linetype = "dashed")+
    # geom_hline(data=selection_mine_snps_res_US,aes(yintercept = mean(selection_mine_snps_res_US$FST),color="\nMean FST\nneutral loci\n"),size=1)+

  # geom_hline(data=per_pop_snps_exp_res,aes(yintercept = mean(per_pop_snps_exp_res$FST),color="\nMean FST\nexperiment\n"),size=1)+
    geom_mark_hull(data = envelope,concavity = 3,expand=0,radius=0,aes(x=xseq,y=V2,color="\nNeutral\nenvelope\n"),size=2)+
  # guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_colour_manual(name="",values=colors_plot_2) +
  # guides(color=guide_legend(ncol=2))+
  xlab("Total heterozygosity") + 
  ylab("FST")+
  ylim(c(0,1))+
  # ggtitle("AOD simulations SNPs\nin two populations")+
  theme_tufte(base_family="Helvetica")+
  theme(  
    # legend.position="bottom",
      text = element_text(size = 16)
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(),
        # panel.background = element_rect(fill = "transparent"), # bg of the panel
    # plot.background = element_rect(fill = "transparent", color = NA)
    )
  # theme(
  #    # legend.position="bottom",
  #       panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank(),
  #       panel.background = element_rect(fill = "transparent"), # bg of the panel
  #   plot.background = element_rect(fill = "transparent", color = NA))

 ggsave("AOD_snps_two_demes_2.pdf",  width = 6, height = 5, units = "in", dpi="retina", bg = "transparent"  )
 ###################################################
##"EXPERIMENT SNPs two demes"
###################################################
colors_plot_3 <- c(
  "\nNeutral loci\n"="green",
  "\nMean FST\nneutral loci\n"="green3",
  "\nLoci under\nselection\n"="deepskyblue",
  "\nExpected\nneutral FST\n"="black",
  "\nMean FST\n"="blue",
  "\nFly\nexperiment\n"="deeppink",
  "\nMean FST\nexperiment\n"="deeppink3",
  "\nNeutral\nenvelope\n"="deepskyblue")

ggplot()+
  # geom_pointdensity(data=selection_mine_snps_res_NS,aes(x=He,y=FST,colour="\nNeutral loci\n"),size=4)+
    # geom_pointdensity(data=selection_mine_snps_res_NS,aes(x=He,y=FST),size=5)+
    # scale_color_viridis(option ="B",name="Number of\nneighboring loci" )+
  # new_scale_color()+
  # geom_point(data=selection_mine_snps_res_US,aes(x=He,y=FST,colour="\nLoci under\nselection\n"),size=2)+
  geom_point(data=per_pop_snps_exp_res,aes(x=He,y=FST,color="\nFly\nexperiment\n"),size=5,alpha=1/3,key_glyph = draw_key_point)+
  geom_hline(data=selection_mine_snps_res_NS,aes(yintercept = mean(per_pop_snps_exp_res$FST),color="\nMean FST\n"),size=2,key_glyph = draw_key_smooth)+
  geom_hline(data=selection_mine_snps_res_US,aes(yintercept = Fst_expected,color="\nExpected\nneutral FST\n"),size=2,linetype = "dashed")+
    # geom_hline(data=selection_mine_snps_res_US,aes(yintercept = mean(selection_mine_snps_res_US$FST),color="\nMean FST\nneutral loci\n"),size=1)+

  # geom_hline(data=per_pop_snps_exp_res,aes(yintercept = mean(per_pop_snps_exp_res$FST),color="\nMean FST\nexperiment\n"),size=1)+
    geom_mark_hull(data = envelope,concavity = 3,expand=0,radius=0,aes(x=xseq,y=V2,color="\nNeutral\nenvelope\n"),size=2,key_glyph = draw_key_smooth)+
  # guides(fill=guide_legend(override.aes=list(color=NA))) +
  scale_colour_manual(name="",values=colors_plot_3) +
  # guides(color=guide_legend(ncol=2))+
  xlab("Total heterozygosity") + 
  ylab("FST")+
  ylim(c(0,1))+
  # ggtitle("AOD simulations SNPs\nin two populations")+
  theme_tufte(base_family="Helvetica")+
  theme(  
    # legend.position="bottom",
      text = element_text(size = 16)
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(),
        # panel.background = element_rect(fill = "transparent"), # bg of the panel
    # plot.background = element_rect(fill = "transparent", color = NA)
    )
  # theme(
  #    # legend.position="bottom",
  #       panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank(),
  #       panel.background = element_rect(fill = "transparent"), # bg of the panel
  #   plot.background = element_rect(fill = "transparent", color = NA))

 ggsave("experiment_two_demes.pdf",  width = 6, height = 5, units = "in", dpi="retina", bg = "transparent"  )
 ###################################################
##"Neutral simulations SNPs"
###################################################
ggplot()+
  geom_point(data=neutral_mine_snps_res_NS,aes(x=He,y=FST,colour="\nNeutral loci\n"),size=4)+
  geom_point(data=neutral_mine_snps_res_US,aes(x=He,y=FST,colour="\nLoci under\nselection\n"),size=2)+
  geom_mark_hull(data = envelope,concavity = 3,expand=0,radius=0,aes(x=xseq,y=V2,color="\nNeutral\nenvelope\n"),fill="black",size=1,alpha=1/7)+
  geom_point(data=per_pop_snps_exp_res,aes(x=He,y=FST,color="\nExperiment\n"),size=2)+
  geom_hline(data=neutral_mine_snps_res_NS,aes(yintercept = mean(neutral_mine_snps_res_NS$FST),color="\nMean FST loci\nunder selection\n"),size=1)+
  geom_hline(data=neutral_mine_snps_res_US,aes(yintercept = mean(neutral_mine_snps_res_US$FST),color="\nMean FST\nneutral loci\n"),size=1)+
  geom_hline(data=per_pop_snps_exp_res,aes(yintercept = mean(per_pop_snps_exp_res$FST),color="\nMean FST\nexperiment\n"),size=1)+
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_colour_manual(name="",values=colors_plot) +
  xlab("Total heterozygosity") + 
  ylab("FST")+
  ylim(c(0,1))+
  ggtitle("Neutral simulations SNPs\nin two populations")+
  theme_tufte(base_family="Helvetica")+
  theme(
     # legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA))

 ggsave("neutral_snps_two_demes.pdf",  width = 6, height = 6, units = "in", dpi="retina", bg = "transparent"  )
 
###################################################
##"Neutral simulations SNPs ten demes"
###################################################
betas.sim <- neutral_fdist_ten_demes
betas.sim <- betas.sim[,c("He","FST")] 
colnames(betas.sim) <- c("alpha", "beta")
betas.sim <- betas.sim[complete.cases(betas.sim),]
maxalpha <- max(c(betas.sim$alpha))
minalpha <- min(c(betas.sim$alpha))
maxbeta <- max(c(betas.sim$beta))
minbeta <- min(c(betas.sim$beta))
alpha.cell <- 0.01
beta.cell <- 0.01
xseq<-seq(minalpha,maxalpha,alpha.cell)
yseq<-seq(minbeta,maxbeta,beta.cell)
envelope_temp <-  as.data.frame(get.CI.alpha(xseq=xseq, percentile=0.975, betas.sim=betas.sim))
envelope <- as.data.frame(rbind(cbind(xseq,envelope_temp$V1),cbind(xseq,envelope_temp$V2)))
envelope <- envelope[complete.cases(envelope$V2),]
envelope[which(envelope$xseq>0.47),"xseq"] <- 0.5

ggplot()+
  geom_point(data=neutral_fdist_ten_demes,aes(x=He,y=FST,colour="\nNeutral loci\n"),size=4)+
  # geom_point(data=neutral_mine_snps_res_US,aes(x=He,y=FST,colour="\nLoci under\nselection\n"),size=1.5)+
  geom_mark_hull(data = envelope,concavity = 2,expand=0,radius=0,aes(x=xseq,y=V2,color="\nNeutral\nenvelope\n"),fill="black",size=1,alpha=1/7)+
  geom_point(data=per_pop_snps_exp_res,aes(x=He,y=FST,color="\nExperiment\n"),size=2)+
  # geom_hline(data=neutral_mine_msats_res_NS,aes(yintercept = mean(neutral_mine_msats_res_NS$FST),color="\nMean FST loci\nunder selection\n"),size=1)+
  geom_hline(data=neutral_fdist_ten_demes,aes(yintercept = mean(neutral_fdist_ten_demes$FST),color="\nMean FST\nneutral loci\n"),size=1)+
  geom_hline(data=per_pop_snps_exp_res,aes(yintercept = mean(per_pop_snps_exp_res$FST),color="\nMean FST\nexperiment\n"),size=1)+
  scale_colour_manual(name="",values=colors_plot) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  xlab("Total heterozygosity") + 
  ylab("FST")+
  ylim(c(0,1))+
  ggtitle("Neutral simulations SNPs\nin ten populations")+
  theme_tufte(base_family="Helvetica")+
  theme(
     # legend.position="bottom",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA))

 ggsave("neutral_snps_ten_demes.pdf",  width = 6, height = 6, units = "in", dpi="retina", bg = "transparent"  )
 
###################################################
###################################################
################## SNPS low dispersal ##########################
###################################################
number_transfers <- 1
transfer_each_gen <- 29
Ne_fst <- 14.3
population_size_dispersal <- 14
dispersal_rate <- (number_transfers / transfer_each_gen) / (population_size_dispersal)
Fst_expected <- 1 / ((4 * Ne_fst * dispersal_rate) * ((2 / (2 - 1)) ^ 2) + 1)
print(Fst_expected)

# high_dispersal_20_demes <- read.table("fdist_high_dispersal.sims", header=TRUE)
high_dispersal_20_demes <- read.table("fdist_high_dispersal_2.sims", header=TRUE)
betas.sim <- high_dispersal_20_demes
betas.sim <- betas.sim[,c("He","FST")] 
colnames(betas.sim) <- c("alpha", "beta")
betas.sim <- betas.sim[complete.cases(betas.sim),]
maxalpha <- max(c(betas.sim$alpha))
minalpha <- min(c(betas.sim$alpha))
maxbeta <- max(c(betas.sim$beta))
minbeta <- min(c(betas.sim$beta))
alpha.cell <- 0.01
beta.cell <- 0.01
xseq<-seq(minalpha,maxalpha,alpha.cell)
yseq<-seq(minbeta,maxbeta,beta.cell)
envelope_temp <-  as.data.frame(get.CI.alpha(xseq=xseq, percentile=0.975, betas.sim=betas.sim))
envelope <- as.data.frame(rbind(cbind(xseq,envelope_temp$V1),cbind(xseq,envelope_temp$V2)))
envelope_high_dispersal <- envelope[complete.cases(envelope$V2),]
envelope_high_dispersal[which(envelope_high_dispersal$xseq>0.48),"xseq"] <- 0.5
  
# low_dispersal_20_demes <- read.table("fdist_low_dispersal.sims", header=TRUE)
low_dispersal_20_demes <- read.table("fdist_low_dispersal_2.sims", header=TRUE)

betas.sim <- low_dispersal_20_demes
betas.sim <- betas.sim[,c("He","FST")] 
colnames(betas.sim) <- c("alpha", "beta")
betas.sim <- betas.sim[complete.cases(betas.sim),]
maxalpha <- max(c(betas.sim$alpha))
minalpha <- min(c(betas.sim$alpha))
maxbeta <- max(c(betas.sim$beta))
minbeta <- min(c(betas.sim$beta))
alpha.cell <- 0.01
beta.cell <- 0.01
xseq<-seq(minalpha,maxalpha,alpha.cell)
yseq<-seq(minbeta,maxbeta,beta.cell)
envelope_temp <-  as.data.frame(get.CI.alpha(xseq=xseq, percentile=0.975, betas.sim=betas.sim))
envelope <- as.data.frame(rbind(cbind(xseq,envelope_temp$V1),cbind(xseq,envelope_temp$V2)))
envelope_low_dispersal <- envelope[complete.cases(envelope$V2),]
envelope_low_dispersal[which(envelope_low_dispersal$xseq>0.46),"xseq"] <- 0.5

colors_fdist <- c(
  "\nFDist high\ndispersal\n"="deepskyblue",
  "\nFDist low\ndispersal\n"="green",
  "\nEnvelope based on\nexpected FST\n"="deepskyblue3",
  "\nEnvelope based on\nfly experiment FST\n"="green3",
  "\nExpected\nFST\n"="black",
  "\nSNPs\nfly experiment\n"="deeppink",
  "\nMean FST\nfly experiment\n"="deeppink4"
   )

 
ggplot()+
geom_point(data=low_dispersal_20_demes,aes(x=He,y=FST),colour="deepskyblue",size=2,alpha=1/6)+
geom_mark_hull(data=envelope_low_dispersal,concavity=2,expand=0,radius=0,aes(x=xseq,y=V2,color="\nEnvelope based on\nexpected FST\n"),fill="black",size=1,alpha=1/8)+
geom_point(data=high_dispersal_20_demes,aes(x=He,y=FST),colour="green",size=2,alpha=1/6,show.legend = F)+
geom_mark_hull(data=envelope_high_dispersal,concavity=2,expand=0,radius=0,aes(x=xseq,y=V2,color= "\nEnvelope based on\nfly experiment FST\n"),fill="black",size=1,alpha=1/8)+
geom_point(data=per_pop_snps_exp_res,aes(x=He,y=FST,color="\nSNPs\nfly experiment\n"),size=3,alpha=1/2)+
geom_hline(aes(yintercept=Fst_expected,color="\nExpected\nFST\n"),size=2,linetype = "dashed")+
geom_hline(data=per_pop_snps_exp_res,aes(yintercept = mean(per_pop_snps_exp_res$FST),color="\nMean FST\nfly experiment\n"),size=2)+
guides(color=guide_legend(override.aes=list(fill=NA))) +
scale_colour_manual(name="",values=colors_fdist) +
xlab("Total heterozygosity") + 
ylab("FST")+
# ylim(c(0,1))+
# ggtitle("FDist envelopes using two different\ndispersal rates in twenty populations")+
theme_tufte(base_family="Helvetica")+
theme( # legend.position="bottom",
   text = element_text(size = 16),
   panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA))

 ggsave("fdist_experiment_comparison.pdf",  width = 7, height = 6, units = "in", dpi="retina", bg = "transparent"  )
 
