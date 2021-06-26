 library(readr)
 library(fields)
 NS_chr1_cattle <- read_csv("synonymous_cattle_chr_7.txt")
 NS_chr1_cattle <- as.data.frame(NS_chr1_cattle)
 colnames(NS_chr1_cattle) <- "NS"
 NS_chr1_cattle$NS <- NS_chr1_cattle[order(NS_chr1_cattle$NS),]
 NS_chr1_cattle$NS <- unique(NS_chr1_cattle$NS)
 chr_length <- NS_chr1_cattle[nrow(NS_chr1_cattle),1]
 resolution_bp <- 100000
# Calculating the number of non-synonymous mutations in windows of 100 Kbp
NS_chr1_cattle_binned <- as.data.frame(table(cut(NS_chr1_cattle$NS,breaks=seq(0,chr_length,resolution_bp))))[,2]

map_cattle <- read_csv("recom_map_chr_7_bos_taurus.csv")

map_cattle_binned <- stats.bin(map_cattle$Location,map_cattle$`mean sexes`,breaks = seq(0,chr_length,resolution_bp))
map_cattle_binned_b <- unname(map_cattle_binned$stats[2,]* map_cattle_binned$stats[1,])

map <- as.data.frame(cbind(map_cattle_binned_b,NS_chr1_cattle_binned))
colnames(map) <- c("cM","non_synonymous")
map[is.na(map$cM),] <- 0
write.csv(map,file = "map_chr_7.csv")                               
