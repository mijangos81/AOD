# library(ggplot2)
 library(patchwork)
# library(viridis)
# library(ggalt)
library(parallel)
source('gHap.R')
# region_size <- 100000
# start_loc <- 1
# end_loc <- 10000000

offspring_ld_pop1 <- pop1
# offspring_ld_pop2 <- pop2

 offspring_ld_pop1$V1[offspring_ld_pop1$V1=="Male"]   <- 1
  offspring_ld_pop1$V1[offspring_ld_pop1$V1=="Female"] <- 2
  offspring_ld_pop1$V2 <- paste0("offspring_","pop1")
  offspring_ld_pop1$id <- paste0("pop1_",row.names(offspring_ld_pop1))
  
  # offspring_ld_pop2$V1[offspring_ld_pop2$V1=="Male"]   <- 1
  # offspring_ld_pop2$V1[offspring_ld_pop2$V1=="Female"] <- 2
  # offspring_ld_pop2$V2 <- paste0("offspring_","pop2")
  # offspring_ld_pop2$id <- paste0("pop2_",row.names(offspring_ld_pop2))
 # samples_haplotypes <-  rbind(offspring_ld_pop1,offspring_ld_pop2)
   samples_haplotypes <-  offspring_ld_pop1

    samples_haplo <- samples_haplotypes[,c("V2","id")]
    
  #  plink_temp_pop1 <- apply(offspring_ld_pop1,1,ped)
  # plink_temp_pop2 <- apply(offspring_ld_pop2,1,ped)
  # # 
  # plink_ped <- c(plink_temp_pop1,plink_temp_pop2)
  # 
  # # plink_ped <- apply(data,1,ped)
  # haploview <- gsub("a", "1", plink_ped) # converting allele names to numbers
  # haploview <- gsub("A", "2", haploview)
  # write.table(haploview,file = paste0(path.folder_sim,"/","haploview.ped"),quote = F,row.names = F,col.names = F)
  # snp_stats  <- read.pedfile_b(paste0(path.folder_sim,"/","haploview.ped"),sep = " ",snps = plink_map$V4,show_warnings=T)
  # genotype <- snp_stats$genotypes
  # snpsum.col <- col.summary(genotype)
  # ns<- snpsum.col[complete.cases(snpsum.col$z.HWE),]
  
  chromosomes_pop <- c(samples_haplotypes$V3,samples_haplotypes$V4)

chromosomes_pop <- gsub("[-^1-9]", "o", chromosomes_pop)
chromosomes_pop <- gsub("a", "0", chromosomes_pop) # converting allele names to numbers
chromosomes_pop <- gsub("A", "1", chromosomes_pop)
find_haplo <- strsplit(chromosomes_pop,split = "")
find_haplo_2 <-data.frame((sapply(find_haplo,c)), stringsAsFactors = F)
ref <- "A"
alt <- "a"
info_haplo <- as.data.frame(cbind(1,1:(nrow(recombination_map)-1),recombination_map$loc_bp[1:(nrow(recombination_map)-1)],ref,alt), stringsAsFactors = F)
colnames(info_haplo) <- c("chr","marker","pos","ref","alt")
find_haplo_3 <- cbind(info_haplo,find_haplo_2) 
snp_final_2 <- read.table(paste0(path.folder_sim,"/",files_snps_pop1[1]),header=T,sep = ",",skip = 36)

find_haplo_3$HW <- snp_final_2$z.HWE_pop1
find_haplo_3 <- find_haplo_3[complete.cases(find_haplo_3),]
find_haplo_3 <- find_haplo_3[,-ncol(find_haplo_3)]
markers_haplo <-  find_haplo_3[,1:5]
markers_haplo$pos <- as.numeric(as.character(markers_haplo$pos))

phase_haplo <- find_haplo_3[,6:ncol(find_haplo_3)]

phase_haplo_2 <- ghap.loadphase(
  samples.file=samples_haplo,
  markers.file=markers_haplo,
  phase.file=phase_haplo,
  verbose = TRUE
)

# blocks.mkr <- ghap.blockgen(phase_haplo_2, windowsize = 10, slide = 1, unit = "marker",nsnp = 4)

ghap.haplotyping(phase=phase_haplo_2, blocks=hap_blocks, outfile="res_haplo", freq = c(0,1), drop.minor = F, batchsize = 500, ncores = 1, verbose = TRUE)

haplo <- ghap.loadhaplo("res_haplo.hapsamples", "res_haplo.hapalleles", "res_haplo.hapgenotypes")
hapstats <- ghap.hapstats(haplo, ncores = 2)
blockstats <- ghap.blockstats(hapstats, ncores = 2)

location_hap <- lapply(hapstats$BP1, findInterval, 
                   vec = as.numeric(paste(unlist(snp_final$loc_bp))))
hapstats_2 <- hapstats[,c("BLOCK","BP1","ALLELE","FREQ")]
hapstats_2$BP1 <- c(location_hap)
hapstats_2$BLOCK <- as.factor(hapstats_2$BLOCK)

split_haplo <- split(hapstats_2,hapstats_2$BLOCK)
split_haplo <- lapply(split_haplo,function(x){x[order(x$FREQ,decreasing = T),]})
split_haplo2 <- lapply(split_haplo,"[",,c("BP1","ALLELE"))
position_hap <- lapply(split_haplo2,"[[",1,1)
position_hap <- floor(unlist(position_hap))

position_hap <-  floor(rescale(position_hap,to=c(1,floor(width_poly*0.964))))

split_haplo3 <- lapply(split_haplo2,function(x){strsplit(x$ALLELE,split = "") } )
split_haplo4 <- lapply(split_haplo3,sapply,rbind)
split_haplo5 <- lapply(split_haplo4,t)
for(i in 1:length(position_hap)){
  col_num <- ncol(split_haplo5[[i]]) -1
  row_num <- nrow(split_haplo5[[i]])
  colnames(split_haplo5[[i]]) <- position_hap[i] : (position_hap[i] + col_num)
  rownames(split_haplo5[[i]]) <-  1:row_num
}
split_haplo6 <- lapply(split_haplo5,function(x){as.data.frame(as.table(x))})

plot_haplo <- rbindlist(split_haplo6)
plot_haplo$Var1 <- as.numeric(as.character(plot_haplo$Var1))
plot_haplo$Var1 <- 1-(plot_haplo$Var1)
# plot_haplo$Var1 <- (plot_haplo$Var1)-3
plot_haplo$Var2 <- as.numeric(as.character(plot_haplo$Var2))
plot_haplo$Freq <- as.character(plot_haplo$Freq)
plot_haplo$Freq <- plot_haplo$Freq=="A"
plot_haplo$Freq <- as.numeric(plot_haplo$Freq)

# 
# 
# 
#   map2color <- function(x,pal,limits=NULL){
#      if(is.null(limits)) limits=range(x)
#      pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
#   }
# 
#   plot_haplo$Freq <- map2color(plot_haplo$Freq,pal= plasma(2))
#     
#   
#   plot_haplo_2 <- matrix(data = (x=as.vector(plot_haplo),pal ), ncol = ncol(plot_haplo), nrow = nrow(plot_haplo),byrow = F)
# 
# p2<-ggplot(data = plot_haplo, aes(x=Var2, y=Var1, fill=Freq,height=1,width=1)) + 
#   geom_tile() +
#    labs(x=NULL, y=NULL, title=NULL)+ 
#   #scale_x_continuous(breaks=ticks, labels=ticks_lab)+
#   theme_tufte(base_family="Helvetica") +
#   theme(legend.position = "none",
#         axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         #axis.ticks = element_line(colour = "black", size = 0.2),
#         axis.ticks.y=element_blank())+
#         #axis.text.x = element_text(face = "bold",size = 4)) +
#    coord_fixed(ratio = 1/1)
#   
# p1+p2 + plot_layout(ncol = 1, heights = c(1,1))
# 
# 
# plot.matrix
# as.data.frame(as.table(split_haplo5[[1]]))
# 
# 
# 
# 
# 
# 
# 
# blocks_haplo <- blockstats[,c(3,6)]
# 
# blocks_haplo_2 <- blocks_haplo[which(blocks_haplo$BP1 >start_loc & blocks_haplo$BP1 < end_loc),]
# blocks_haplo_2$N.ALLELES <- rescale(blocks_haplo_2$N.ALLELES,to=c(0,1))
# 
# test_haplo <- hapstats 
# 
# test_haplo$ALLELE <- gsub("A",1,test_haplo$ALLELE)
# test_haplo$ALLELE <- gsub("a",0,test_haplo$ALLELE)
# haplo_blocks <- split(test_haplo,test_haplo$BLOCK)
# 
# aod_results <-  as.data.frame(matrix(nrow = length(haplo_blocks) ,ncol = 2))
# for (i in 1:length(haplo_blocks)){
#       haplo_temp <- haplo_blocks[[i]]
#       haplo_temp_2 <- tail(haplo_temp,2)
#       haplo_temp_3 <- as.data.frame(strsplit(c(haplo_temp_2[1,"ALLELE"],haplo_temp_2[2,"ALLELE"]), split = ""),stringsAsFactors =F)
#       haplo_temp_3$test <- haplo_temp_3[,1]==haplo_temp_3[,2]
#       haplo_temp_3[,1] <- as.numeric(as.character(haplo_temp_3[,1]))
#       haplo_temp_3[,2] <- as.numeric(as.character(haplo_temp_3[,2]))
#       aod_test <- sum(apply(haplo_temp_3,1,prod)) 
#       # / nrow(haplo_temp_3[which(haplo_temp_3[,2]==1 & haplo_temp_3[,3]==T),])
#       aod_results[i,] <- cbind(haplo_temp$BP1[1],aod_test)
# }
# 
# 
# ld_res_pop2 <- LD_fun(data=offspring_pop2, name_pop="pop2", show_warnings=F)
#       ld_columns_pop2 <- ld_res_pop2[[1]]
#       snpsum.col_pop2 <- ld_res_pop2[[2]]
#       
# test_ld_pop1 <-  ld_columns_pop1[order(ld_columns_pop1$Var1),]
# test_ld_pop1 <- test_ld_pop1[which(test_ld_pop1$dis>10000 & test_ld_pop1$dis<20000),]
# test_ld_2_pop1 <- split(test_ld_pop1,test_ld_pop1$Var1)
# aod_ld_pop1  <-  as.data.frame(matrix(nrow = length(test_ld_2_pop1) ,ncol = 3))
# for (i in 1:length(test_ld_2_pop1)){
#       ld_temp_pop1  <- test_ld_2_pop1[[i]]
#       ld_temp_2_pop1 <- head(ld_temp_pop1,1)
#       aod_ld_pop1[i,] <- cbind(ld_temp_2_pop1$Var1,ld_temp_2_pop1$Freq,ld_temp_2_pop1$dis)
# }
# aod_ld_pop1$V2 <- 1-aod_ld_pop1$V2
# 
# test_ld_pop2 <-  ld_columns_pop2[order(ld_columns_pop2$Var1),]
# test_ld_pop2 <- test_ld_pop2[which(test_ld_pop2$dis>10000 & test_ld_pop2$dis<20000),]
# test_ld_2_pop2 <- split(test_ld_pop2,test_ld_pop2$Var1)
# aod_ld_pop2  <-  as.data.frame(matrix(nrow = length(test_ld_2_pop2) ,ncol = 3))
# for (i in 1:length(test_ld_2_pop2)){
#       ld_temp_pop2  <- test_ld_2_pop2[[i]]
#       ld_temp_2_pop2 <- head(ld_temp_pop2,1)
#       aod_ld_pop2[i,] <- cbind(ld_temp_2_pop2$Var1,ld_temp_2_pop2$Freq,ld_temp_2_pop2$dis)
# }
# aod_ld_pop2$V2 <- 1-aod_ld_pop2$V2
# 
# 
# 
# aod_res <- aod_results
# aod_res$V2 <- rescale(aod_res$V2,to= c(0,1))
# aod_res <- aod_res[order(aod_res$V1),]
# aod_res  <- aod_res[which(aod_res$V1 > start_loc & aod_res$V1 < end_loc),]