# CIRCLE PLOT
library(scales)
source("circle_plot_script.R")

ld_max_pairwise <- 2000000  # maximun distance, in basepairs, at which pairwise LD should be calculated
ld_resolution <- 500000 # resolution, in basepairs, at which LD should be measured
region_size <- 500000 # the size of the window at which the LD statistics and number of NS are calculated
map_resolution <- 100000
gen_number_dispersal <- 34


map <- read_csv("fly_recom_map.csv")
map$Chr <- as.character(map$Chr)
map <- map[which(map$Chr==chromosome_name),]
map <- as.data.frame(map$cM/1000)
colnames(map) <- "cM"
map[is.na(map$cM),] <- 0
chromosome_length <- (nrow(map)+1) * map_resolution
chromosome_name <- "2L"

path.folder_sim <- paste0(getwd(),"/LD_analyses_fly_sims_neutral")
files_LD <- paste0(path.folder_sim,"/",dir(path.folder_sim,pattern = "^final_LD_locus"))
link_final <- read_csv(files_LD[1],skip = variables_number)
# link_final_plot <- as.data.frame(cbind(loc=as.numeric(row.names(link_final)), link_final))
link_final_plot <- as.data.frame(link_final[,1:(ncol(link_final)-1)])
colnames(link_final_plot) <- c("col",as.character(c(5:8)))

brk<- seq(0,chromosome_length,region_size)
sectors <- length(brk)-1
ld_window <- (ld_max_pairwise/ld_resolution)
ld_final <- NULL
for (ld in 2:(ld_window+1)) {
  ld_temp <- as.numeric(stats.bin(link_final_plot[,1], link_final_plot[,ld], breaks = brk)[[3]][2,])
  ld_final <- cbind(ld_final,ld_temp)
}
# ld_final <- link_final
# rsquare_ceiling<- 0.2
# ld_final[ld_final> rsquare_ceiling] <- rsquare_ceiling
# ld_final <- rescale(ld_final,to=c(0,2))
ld_final_b <- as.data.frame(ld_final)
ld_final_b <- ld_final_b[,c(ld_window:1)]

heatmap <- as.data.frame(matrix(nrow = sectors))
heatmap[,1] <- chromosome_name
heatmap_c <- cbind(heatmap,loc=as.numeric(1:sectors),ld_final_b)

segm <- as.data.frame(cbind(seg.name=chromosome_name,seg.Start=0:(sectors-1),seg.End=1:sectors,the.v="NA",NO="NA"))
# segm_b <- segAnglePo(segm,seg=chromosome_name,angle.start=0,angle.end=180)
segm_b <- segAnglePo(segm,seg=chromosome_name,angle.start=0,angle.end=360)

break_bins <- seq(0,chr_length,region_size)
posi<- seq(map_resolution/2,(nrow(map)*map_resolution),map_resolution)
  # seq(region_size/2,chromosome_length,region_size)/map_resolution
het_pop1 <- as.numeric(unlist(mean_het_pop1[gen_number_dispersal,]))
het_pop1 <- het_pop1[-loc_exp_loci]
NS_pop1_temp <- as.data.frame(cbind(posi,het_pop1))
NS_pop1_temp <- stats.bin(NS_pop1_temp$posi,NS_pop1_temp$het_pop1,breaks = break_bins)[[3]][2,]
NS_pop1 <- unname(unlist(NS_pop1_temp))
exp_het_pop1 <- as.numeric(unlist(mean_expected_het_pop1[gen_number_dispersal,]))
exp_het_pop1 <- exp_het_pop1[-loc_exp_loci]
NS_exp_pop1_temp <- as.data.frame(cbind(posi,exp_het_pop1))
NS_exp_pop1_temp <- stats.bin(NS_exp_pop1_temp$posi,NS_exp_pop1_temp$exp_het_pop1,breaks = break_bins)[[3]][2,]
NS_exp_pop1 <- unname(unlist(NS_exp_pop1_temp))
NS_pop1 <- (NS_pop1 - NS_exp_pop1) / NS_exp_pop1
# NS_pop1<- rescale(NS_pop1,to=c(0,2),)

het_pop2 <- as.numeric(unlist(mean_het_pop2[gen_number_dispersal,]))
het_pop2 <- het_pop2[-loc_exp_loci]
NS_pop2_temp <- as.data.frame(cbind(posi,het_pop2))
NS_pop2_temp <- stats.bin(NS_pop2_temp$posi,NS_pop2_temp$het_pop2,breaks = break_bins)[[3]][2,]
NS_pop2 <- unname(unlist(NS_pop2_temp))
exp_het_pop2 <- as.numeric(unlist(mean_expected_het_pop2[gen_number_dispersal,]))
exp_het_pop2 <- exp_het_pop2[-loc_exp_loci]
NS_exp_pop2_temp <- as.data.frame(cbind(posi,exp_het_pop2))
NS_exp_pop2_temp <- stats.bin(NS_exp_pop2_temp$posi,NS_exp_pop2_temp$exp_het_pop2,breaks = break_bins)[[3]][2,]
NS_exp_pop2 <- unname(unlist(NS_exp_pop2_temp))
NS_pop2 <- (NS_pop2 - NS_exp_pop2) / NS_exp_pop2
# NS_pop2<- rescale(NS_pop2,to=c(0,2),)

fst_ld <- as.numeric(unlist(mean_Fst[gen_number_dispersal,]))
fst_ld <- fst_ld[-loc_exp_loci]
NS_fst_temp <- as.data.frame(cbind(posi,fst_ld))
NS_fst_temp <- stats.bin(NS_fst_temp$posi,NS_fst_temp$fst_ld,breaks = break_bins)[[3]][2,]
NS_fst_temp <-unname(unlist(NS_fst_temp))
NS_fst_temp <- (NS_fst_temp - Fst_expected)/ Fst_expected
# NS_fst_temp <- NS_fst_temp * -1
# NS_fst_temp <- rescale(NS_fst_temp,to=c(0,2))

nonsyn <- link_final$number_nonsyn_b
nonsyn[nonsyn>100] <- 100
  # as.data.frame(cbind(loc=rep(1:sectors, each=region_size/map_resolution),number_nonsyn_b))
# nonsyn <- vaggregate(.value = nonsyn$number_nonsyn_b,.group=nonsyn$loc, .fun = sum)
# nonsyn <- rescale(nonsyn,to=c(0,2))

NS<- as.data.frame( matrix(nrow = sectors) )
NS[,1] <- chromosome_name
NS$position <-  1:sectors
NS$NS <- nonsyn
NS$fst <- NS_fst_temp
NS$het_pop2<- unlist(NS_pop2)
NS$het_pop1<- unlist(NS_pop1)
# NS$fst <- c(NS_fst_temp,NA)
# NS$het_pop2<- c(unlist(NS_pop2),NA)
# NS$het_pop1<- c(unlist(NS_pop1),NA)

reg_ld <- rowMeans(ld_final_b,na.rm = T)
regre <- cbind(NS,reg_ld)
regre_res <- summary(lm(regre$fst~regre$NS+regre$reg_ld))

grid <- as.data.frame( matrix(nrow = sectors) )
grid[,1] <- chromosome_name
grid$position <-  0.5:(sectors-0.5)
grid$grid <- 1
grid$arc <- 1

# par(mar=c(2,2,2,2))
# plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="",main="")
# circos(R=180,cir=segm_b,W=190,mapping=NS,col.v=3,type="ss",B=F,scale=F,cex=3)
# circos(R=350,cir=segm_b,W=60,mapping=heatmap_c,col.v=3,type="heatmap",B=F,col.bar=F,lwd = 32)
# circos(R=81,cir=segm_b,W=150,mapping=NS,col.v=4,type="heatmap",B=F,col.bar =T,col.bar.po = "topleft",lwd = 17)

heatmap_d <- link_final_plot
heatmap_d$col <- 1:nrow(link_final_plot)
heatmap_LD <- reshape2::melt(heatmap_d,id="col")
heatmap_LD <- heatmap_LD[order(heatmap_LD$variable),]

heatmap_stats_temp <- NS[,c(2,4:6)]
colnames(heatmap_stats_temp) <- c("position","1","2","3")
heatmap_stats <- reshape2::melt(heatmap_stats_temp,id="position")

heatmap_NS_temp <- NS[,2:3]
colnames(heatmap_NS_temp) <- c("position","4")
heatmap_NS <- reshape2::melt(heatmap_NS_temp,id="position")

# map_df_temp <- as.data.frame(cbind(seq(100000,chromosome_length-100000,100000), map))
# colnames(map_df_temp) <- c("location","recombination")
# map_df_temp_2 <- stats.bin(map_df_temp$location, map_df_temp$recombination, breaks = break_bins)
# map_df <- unlist(unname(((map_df_temp_2[[3]][2,] * map_df_temp_2[[3]][1,]) * 100)))
# map_df[is.na(map_df)] <- 0
# heatmap_map <- as.data.frame(cbind(1:nrow(link_final_plot),map_df))
# heatmap_map$map_df <- rescale()


library(viridis)
library(ggnewscale)
library(colorspace)
library(ggplot2)
library(ggthemes)
library(gplots)
heatmap_LD$variable <- as.numeric(as.character(heatmap_LD$variable))
heatmap_stats$variable <- as.numeric(as.character(heatmap_stats$variable))
heatmap_NS$variable <- as.numeric(as.character(heatmap_NS$variable))
# heatmap_NS$value <- heatmap_NS$value/10

ggplot()+
geom_tile(data=heatmap_stats,aes(x=position,y=variable,fill=value),color="black")+
scale_fill_continuous_diverging(palette = "Broc",name="Statistic \nbias",p1=2,p2=2)+
new_scale_fill() +
geom_raster(data=heatmap_LD,aes(x=col,y=variable,fill=value),interpolate = T)+
scale_fill_viridis(option = "D",name="r2 LD")+
 new_scale_fill() + 
# geom_point(data=heatmap_NS,aes(x=position,y=variable,size=value),inherit.aes = F,color="deeppink") +
 geom_raster(data=heatmap_NS,aes(x=position,y=variable,fill=value),interpolate = T) +
 scale_fill_continuous_sequential(palette = "Reds 3",name="Targets of \nselection") +
 # geom_line(data=heatmap_map,aes())+
  theme_void() +
 coord_fixed(ratio = 1/1)


  ggsave(paste0("fly_sim_LD.pdf"),  width = 12, height =6, units = "in", dpi="retina", bg = "transparent" )





