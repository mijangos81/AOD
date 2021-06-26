library(readr)
chromosome <-c("2L","2R","3L","3R","X")
recom_chr_1_temp <- read_csv("fly_recom_map.csv") 
table_cM_bp <-  as.data.frame(matrix(ncol = 4))
colnames(table_cM_bp) <- c("Chromosome","Genome_size","centiMorgans","Mbp_cM")
for (chrom in chromosome) {
ns_chr_2_chill_b <- read_csv(paste0(getwd(),"/non_syn_&_syn_fly/","ns_", chrom, "_BDGP6.csv"))
ns_chr_1 <- ns_chr_2_chill_b[,"position start (bp)"]
recom_chr_1 <- recom_chr_1_temp[which(recom_chr_1_temp$Chr==chrom),]
table_cM_bp[chrom,"Chromosome"] <- chrom
table_cM_bp[chrom,"Genome_size"] <-  unlist(unname(ns_chr_1[nrow(ns_chr_1),]))
table_cM_bp[chrom,"centiMorgans"] <- sum(recom_chr_1$cM)/10
table_cM_bp[chrom,"Mbp_cM"] <- (sum(recom_chr_1$cM)/10) / (unlist(unname(ns_chr_1[nrow(ns_chr_1),]))/1000000)
}

