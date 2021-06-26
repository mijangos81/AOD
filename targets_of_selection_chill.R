###### R SETTINGS ######
library(coRdon)
library(refGenome)
library(readr)
library(stringr)
###### CODON USAGE BIAS DATAFRAME ######
dmel_CDS_flybase <- readSet(file="Bos_taurus.ARS-UCD1.2.cds.all.fa")# downloaded from ENSEMBL
# dmel_CDS_flybase <- readSet(file="dmel-all-CDS-r6.31.fasta")# downloaded from FLYBASE
codons_dmel <- codonTable(dmel_CDS_flybase)
# Calculation codon bias
CUB_dmel <- MCB(codons_dmel)
# extracting information from file 
chr_name <- str_match(codons_dmel@ID,"ARS-UCD1.2:(.*?):")[,2]
gene_id <- str_match(codons_dmel@ID,"gene:(.*?) gene_biotype")[,2]
transcript_id <- str_match(codons_dmel@ID,"(.*?).. cds")[,2]
CUB_df <-as.data.frame(cbind(CUB_dmel,chr_name,gene_id,transcript_id))
colnames(CUB_df) <- c("CUB_dmel","chr_name","gene_id","transcript_id")
CUB_df$CUB_dmel <- as.numeric(as.character(CUB_df$CUB_dmel))
CUB_df$chr_name <- as.character(CUB_df$chr_name)
CUB_df$gene_id <- as.character(CUB_df$gene_id)
CUB_df$transcript_id <- as.character(CUB_df$transcript_id)
###### GENE TRANSFER FORMAT (GTF) FILE (downloaded from ENSEMBL) #####
ens_df <- ensemblGenome()
read.gtf(ens_df,"Bos_taurus.ARS-UCD1.2.99.gtf")
GTF_df <- as.data.frame(ens_df@ev$gtf)
chr_name_list <- as.character(1:29)
final_df <- data.frame(Doubles=double(),Ints=integer(),Factors=factor(),Logicals=logical(),Characters=character(),stringsAsFactors=FALSE)

for(chr in 1:length(chr_name_list)){
chr_name <- chr_name_list[[chr]]
ns_BDGP6 <- read_csv(paste0("non_syn_&_syn_chill/ns_chill_chr_ ", chr_name, " .csv"))
ns_BDGP6 <- ns_BDGP6$chrom_start
s_BDGP6 <- read_csv(paste0("non_syn_&_syn_chill/s_chill_chr_ ", chr_name, " .csv"))
s_BDGP6 <- s_BDGP6$chrom_start
GTF <- GTF_df[which(GTF_df$seqid==chr_name),]
codon_bias <- CUB_df[which(CUB_df$chr_name==chr_name),]
GTF <- GTF[,c("feature","start","end","transcript_id","gene_id")]
GTF <- GTF[which(GTF$feature=="CDS"),]
GTF <- GTF[order(GTF$start),]
GTF$CDS_length <- GTF$end - GTF$start
GTF <- GTF[,c("start","end","CDS_length","transcript_id","gene_id")]

rownames(GTF) <- 1:nrow(GTF)
count_mutations <- as.data.frame(matrix(nrow =nrow(GTF),ncol = 2 ))
colnames(count_mutations) <- c("ns","s")
for (i in 1:nrow(GTF)){
  count_mutations[i,"ns"]  <- length(ns_BDGP6[ns_BDGP6 >= as.numeric(GTF[i,"start"]) & ns_BDGP6 <= as.numeric(GTF[i,"end"])])
  count_mutations[i,"s"]  <- length(s_BDGP6[s_BDGP6 >= as.numeric(GTF[i,"start"]) & s_BDGP6 <= as.numeric(GTF[i,"end"])])
}

merge_exons_temp <- as.data.frame(cbind(chr_name,GTF,count_mutations))
merge_exons_temp <- merge_exons_temp[order(merge_exons_temp$transcript_id),]
merge_exons_temp$chr_name <- as.character(merge_exons_temp$chr_name)
merge_exons_temp <- split(x=merge_exons_temp,f=merge_exons_temp$transcript_id)
merge_exons <- as.data.frame(matrix(nrow = length(merge_exons_temp),ncol = 10))
colnames(merge_exons) <- c("chr_name","start","end","transcript_size","CDS_length","number_exons","transcript_id","gene_id","ns","s")
 for(i in 1:length(merge_exons_temp)){
   merge_exons[i,"chr_name"] <- as.character(merge_exons_temp[[i]][1,"chr_name"])
   merge_exons[i,"start"] <- merge_exons_temp[[i]][1,"start"]
   merge_exons[i,"end"] <- merge_exons_temp[[i]][nrow(merge_exons_temp[[i]]),"end"]
   merge_exons[i,"transcript_size"] <-  merge_exons[i,"end"] - merge_exons[i,"start"]
   merge_exons[i,"CDS_length"] <- sum(merge_exons_temp[[i]][,"CDS_length"])
   merge_exons[i,"number_exons"] <- nrow(merge_exons_temp[[i]])
   merge_exons[i,"transcript_id"] <- merge_exons_temp[[i]][1,"transcript_id"]
   merge_exons[i,"gene_id"] <- merge_exons_temp[[i]][1,"gene_id"]
   merge_exons[i,"ns"] <- sum(merge_exons_temp[[i]][,"ns"])
   merge_exons[i,"s"] <- sum(merge_exons_temp[[i]][,"s"])
 }
merge_exons <- as.data.frame(merge_exons)
codon_bias <- codon_bias[order(codon_bias$transcript_id),]
codon_bias <- codon_bias[,c("CUB_dmel","transcript_id")]
final_df_temp <- merge(x=merge_exons, y=codon_bias, by="transcript_id", all=FALSE)
final_df <- rbind(final_df,final_df_temp)
}
write.csv(final_df,file = "chill_targets_of_selection.csv")

