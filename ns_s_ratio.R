library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(readxl)


RecRates_All_Chromosomes_100kb <- read_excel("~/OneDrive - UNSW/Flies DB/RecRates-All-Chromosomes-100kb.xlsx")

recom_X <- as.data.frame(na.omit(RecRates_All_Chromosomes_100kb[,3]))
recom_2L <- as.data.frame(na.omit(RecRates_All_Chromosomes_100kb[,9]))
recom_2R <- as.data.frame(na.omit(RecRates_All_Chromosomes_100kb[,15]))
recom_3L <- as.data.frame(na.omit(RecRates_All_Chromosomes_100kb[,21]))
recom_3R <- as.data.frame(na.omit(RecRates_All_Chromosomes_100kb[,27]))

ns_2L_BDGP6 <- read_csv("~/Dropbox/PhD/d_mel/ns_2L_BDGP6.csv")
ns_2L_BDGP6 <- ns_2L_BDGP6[order(ns_2L_BDGP6$`position start (bp)`),]
ns_2L_BDGP6 <- unname(unlist(ns_2L_BDGP6))

ns_2R_BDGP6 <- read_csv("~/Dropbox/PhD/d_mel/ns_2R_BDGP6.csv")

ns_3R_BDGP6 <- read_csv("~/Dropbox/PhD/d_mel/ns_3R_BDGP6.csv")
ns_3R_BDGP6 <- ns_3R_BDGP6[order(ns_3R_BDGP6$`position start (bp)`),]
ns_3R_BDGP6 <- unname(unlist(ns_3R_BDGP6))

ns_3L_BDGP6 <- read_csv("~/Dropbox/PhD/d_mel/ns_3L_BDGP6.csv")
ns_3L_BDGP6 <- ns_3L_BDGP6[order(ns_3L_BDGP6$`position start (bp)`),]
ns_3L_BDGP6 <- unname(unlist(ns_3L_BDGP6))

ns_X_BDGP6 <- read_csv("~/Dropbox/PhD/d_mel/ns_X_BDGP6.csv")
ns_X_BDGP6 <- ns_X_BDGP6[order(ns_X_BDGP6$`position start (bp)`),]
ns_X_BDGP6 <- unname(unlist(ns_X_BDGP6))

s_2L_BDGP6 <- read_csv("~/Dropbox/PhD/d_mel/s_2L_BDGP6.csv")
s_2L_BDGP6 <- s_2L_BDGP6[order(s_2L_BDGP6$`position start (bp)`),]
s_2L_BDGP6 <- unname(unlist(s_2L_BDGP6))

s_2R_BDGP6 <- read_csv("~/Dropbox/PhD/d_mel/s_2R_BDGP6.csv")

s_3R_BDGP6 <- read_csv("~/Dropbox/PhD/d_mel/s_3R_BDGP6.csv")
s_3R_BDGP6 <- s_3R_BDGP6[order(s_3R_BDGP6$`position start (bp)`),]
s_3R_BDGP6 <- unname(unlist(s_3R_BDGP6))

s_3L_BDGP6 <- read_csv("~/Dropbox/PhD/d_mel/s_3L_BDGP6.csv")
s_3L_BDGP6 <- s_3L_BDGP6[order(s_3L_BDGP6$`position start (bp)`),]
s_3L_BDGP6 <- unname(unlist(s_3L_BDGP6))

s_X_BDGP6 <- read_csv("~/Dropbox/PhD/d_mel/s_X_BDGP6.csv")
s_X_BDGP6 <- s_X_BDGP6[order(s_X_BDGP6$`position start (bp)`),]
s_X_BDGP6 <- unname(unlist(s_X_BDGP6))

genes_2L <- read_csv("~/Dropbox/PhD/d_mel/2L_coding_genes.csv")
genes_2L <- genes_2L[order(genes_2L$`Gene start (bp)`),]

genes_3L <- read_csv("~/Dropbox/PhD/d_mel/3L_coding_genes.csv")
genes_3L <- genes_3L[order(genes_3L$`Gene start (bp)`),]

genes_3R <- read_csv("~/Dropbox/PhD/d_mel/3R_coding_genes.csv")
genes_3R <- genes_3R[order(genes_3R$`Gene start (bp)`),]

genes_X <- read_csv("~/Dropbox/PhD/d_mel/X_coding_genes.csv")
genes_X <- genes_X[order(genes_X$`Gene start (bp)`),]

count_mutations_2L <- as.data.frame(matrix(nrow =nrow(genes_2L),ncol = 2 ))
for (i in 1:nrow(genes_2L)){
  count_mutations_2L[i,1]  <- length(ns_2L_BDGP6[ns_2L_BDGP6 >= as.numeric(genes_2L[i,1]) & ns_2L_BDGP6 <= as.numeric(genes_2L[i,2])])
  count_mutations_2L[i,2]  <- length(s_2L_BDGP6[between(s_2L_BDGP6, as.numeric(genes_2L[i,1]), as.numeric(genes_2L[i,2]))])
}
count_mutations_2L$ns_s <- count_mutations_2L$V1/count_mutations_2L$V2
count_mutations_2L$chr <- "2L"
count_mutations_2L <- count_mutations_2L[count_mutations_2L$V1>100,]
count_mutations_2L <- count_mutations_2L[count_mutations_2L$V2>100,]
count_mutations_2L <- count_mutations_2L[complete.cases(count_mutations_2L$ns_s),]

count_mutations_3R <- as.data.frame(matrix(nrow =nrow(genes_3R),ncol = 2 ))
for (i in 1:nrow(genes_3R)){
  count_mutations_3R[i,1]  <- length(ns_3R_BDGP6[ns_3R_BDGP6 >= as.numeric(genes_3R[i,1]) & ns_3R_BDGP6 <= as.numeric(genes_3R[i,2])])
  count_mutations_3R[i,2]  <- length(s_3R_BDGP6[between(s_3R_BDGP6, as.numeric(genes_3R[i,1]), as.numeric(genes_3R[i,2]))])
}
count_mutations_3R$ns_s <- count_mutations_3R$V1/count_mutations_3R$V2
count_mutations_3R$chr <- "3R"
count_mutations_3R <- count_mutations_3R[count_mutations_3R$V1>100,]
count_mutations_3R <- count_mutations_3R[count_mutations_3R$V2>100,]
count_mutations_3R <- count_mutations_3R[complete.cases(count_mutations_3R$ns_s),]


count_mutations_3L <- as.data.frame(matrix(nrow =nrow(genes_3L),ncol = 2 ))
for (i in 1:nrow(genes_3L)){
  count_mutations_3L[i,1]  <- length(ns_3L_BDGP6[ns_3L_BDGP6 >= as.numeric(genes_3L[i,1]) & ns_3L_BDGP6 <= as.numeric(genes_3L[i,2])])
  count_mutations_3L[i,2]  <- length(s_3L_BDGP6[between(s_3L_BDGP6, as.numeric(genes_3L[i,1]), as.numeric(genes_3L[i,2]))])
}
count_mutations_3L$ns_s <- count_mutations_3L$V1/count_mutations_3L$V2
count_mutations_3L$chr <- "3L"
count_mutations_3L <- count_mutations_3L[count_mutations_3L$V1>100,]
count_mutations_3L <- count_mutations_3L[count_mutations_3L$V2>100,]
count_mutations_3L <- count_mutations_3L[complete.cases(count_mutations_3L$ns_s),]

count_mutations_X <- as.data.frame(matrix(nrow =nrow(genes_X),ncol = 2 ))
for (i in 1:nrow(genes_X)){
  count_mutations_X[i,1]  <- length(ns_X_BDGP6[ns_X_BDGP6 >= as.numeric(genes_X[i,1]) & ns_X_BDGP6 <= as.numeric(genes_X[i,2])])
  count_mutations_X[i,2]  <- length(s_X_BDGP6[between(s_X_BDGP6, as.numeric(genes_X[i,1]), as.numeric(genes_X[i,2]))])
}
count_mutations_X$ns_s <- count_mutations_X$V1/count_mutations_X$V2
count_mutations_X$chr <- "X"
count_mutations_X <- count_mutations_X[count_mutations_X$V1>100,]
count_mutations_X <- count_mutations_X[count_mutations_X$V2>100,]
count_mutations_X <- count_mutations_X[complete.cases(count_mutations_X$ns_s),]

mutations_plot <- rbind(count_mutations_3L,count_mutations_2L,count_mutations_3R,count_mutations_X)

# colnames(mutations_plot) <- c("type","dispersal","T0","T1","T2")
mutations_plot_2 <- melt(mutations_plot,id="chr")
# colnames(mod_dis_fstp) <- c("Stock lines","dispersal","Sampling Time", "Fst")
# mod_dis_fstp$dispersal_f <- factor(mod_dis_fstp$dispersal, levels = c("High dispersal","Moderate dispersal","Low dispersal")) 
# mod_dis_fstp$`Sampling Time` <- factor(mod_dis_fstp$`Sampling Time`,levels = c("T0","T1","T2"))
# mod_dis_fstp$Fst <- as.numeric(mod_dis_fstp$Fst)

pal_2 <- c("deeppink","cyan","blue","green")


  ggplot(mutations_plot,aes(color=chr,fill=chr))  + 
  # facet_wrap(.~dispersal_f,scales = "free") + 
   # geom_point(aes(colour=chr),position=position_jitterdodge(dodge.width=0.8),show.legend = F) +
  geom_boxplot(aes(y=ns_s),alpha=0.25,position = position_dodge(width=0.8))+
    # geom_density(alpha=0.2)+
  scale_fill_manual(values = pal_2)+
  scale_colour_manual(values = pal_2)+
  theme_bw(base_family="Helvetica",base_size = 14) +
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank()) +
        labs( y="dN/dS ratio",title="")
  
library(fields)
X2L_linked_selection <- read_csv("~/Dropbox/PhD/d_mel/2L_linked_selection.csv")
X2L_linked_selection <- as.data.frame(X2L_linked_selection)
X2L_linked_selection_binned <- stats.bin(X2L_linked_selection$position,X2L_linked_selection$BS,breaks = seq(1,X2L_linked_selection[nrow(X2L_linked_selection),"position"],100000))
X2L_linked_selection_binned_2  <- unname(unlist(X2L_linked_selection_binned$stats[2,]))
X2L_linked_selection_plot <- as.data.frame(cbind("2L",X2L_linked_selection_binned_2))
colnames(X2L_linked_selection_plot) <- c("chr","BS")
X2L_linked_selection_plot$BS <- as.numeric(as.character(X2L_linked_selection_plot$BS))
X2L_linked_selection_plot$BS <- X2L_linked_selection_plot$BS/10

X3L_linked_selection <- read_csv("~/Dropbox/PhD/d_mel/3L_linked_selection.csv")
X3L_linked_selection <- as.data.frame(X3L_linked_selection)
X3L_linked_selection_binned <- stats.bin(X3L_linked_selection$position,X3L_linked_selection$BS,breaks = seq(1,X3L_linked_selection[nrow(X3L_linked_selection),"position"],100000))
X3L_linked_selection_binned_2  <- unname(unlist(X3L_linked_selection_binned$stats[2,]))
X3L_linked_selection_plot <- as.data.frame(cbind("3L",X3L_linked_selection_binned_2))
colnames(X3L_linked_selection_plot) <- c("chr","BS")
X3L_linked_selection_plot$BS <- as.numeric(as.character(X3L_linked_selection_plot$BS))
X3L_linked_selection_plot$BS <- X3L_linked_selection_plot$BS/10

X2R_linked_selection <- read_csv("~/Dropbox/PhD/d_mel/2R_linked_selection.csv")
X2R_linked_selection <- as.data.frame(X2R_linked_selection)
X2R_linked_selection_binned <- stats.bin(X2R_linked_selection$position,X2R_linked_selection$BS,breaks = seq(1,X2R_linked_selection[nrow(X2R_linked_selection),"position"],100000))
X2R_linked_selection_binned_2  <- unname(unlist(X2R_linked_selection_binned$stats[2,]))
X2R_linked_selection_plot <- as.data.frame(cbind("2R",X2R_linked_selection_binned_2))
colnames(X2R_linked_selection_plot) <- c("chr","BS")
X2R_linked_selection_plot$BS <- as.numeric(as.character(X2R_linked_selection_plot$BS))
X2R_linked_selection_plot$BS <- X2R_linked_selection_plot$BS/10

X3R_linked_selection <- read_csv("~/Dropbox/PhD/d_mel/3R_linked_selection.csv")
X3R_linked_selection <- as.data.frame(X3R_linked_selection)
X3R_linked_selection_binned <- stats.bin(X3R_linked_selection$position,X3R_linked_selection$BS,breaks = seq(1,X3R_linked_selection[nrow(X3R_linked_selection),"position"],100000))
X3R_linked_selection_binned_2  <- unname(unlist(X3R_linked_selection_binned$stats[2,]))
X3R_linked_selection_plot <- as.data.frame(cbind("3R",X3R_linked_selection_binned_2))
colnames(X3R_linked_selection_plot) <- c("chr","BS")
X3R_linked_selection_plot$BS <- as.numeric(as.character(X3R_linked_selection_plot$BS))
X3R_linked_selection_plot$BS <- X3R_linked_selection_plot$BS/10

linked_selection_plot <- rbind(X2L_linked_selection_plot,X2L_linked_selection_plot,X2R_linked_selection_plot,X3L_linked_selection_plot,X3R_linked_selection_plot)

  ggplot(linked_selection_plot,aes(color=chr,fill=chr))  + 
  # facet_wrap(.~dispersal_f,scales = "free") + 
   # geom_point(aes(colour=chr),position=position_jitterdodge(dodge.width=0.8),show.legend = F) +
  geom_boxplot(aes(y=BS),alpha=0.25,position = position_dodge(width=0.8))+
    # geom_density(alpha=0.2)+
  scale_fill_manual(values = pal_2)+
  scale_colour_manual(values = pal_2)+
  theme_bw(base_family="Helvetica",base_size = 14) +
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank()) +
        labs( y="Percentage",title="Percentage of loss of genetic variation\n due to background selection")
  
colnames(recom_2L) <- c("cM")
recom_2L$chr <- "2L"

colnames(recom_3L) <- c("cM")
recom_3L$chr <- "3L"

colnames(recom_3R) <- c("cM")
recom_3R$chr <- "3R"

colnames(recom_2R) <- c("cM")
recom_2R$chr <- "2R"

colnames(recom_X) <- c("cM")
recom_X$chr <- "X"

recombination_plot <- rbind(recom_2L,recom_2R,recom_3L,recom_3R,recom_X)

pal_3 <- c("deeppink","cyan","blue","green","brown")

  
ggplot(recombination_plot,aes(color=chr,fill=chr))  + 
   geom_boxplot(aes(y=cM),alpha=0.2)+
  scale_fill_manual(values = pal_3)+
  scale_colour_manual(values = pal_3)+
  theme_bw(base_family="Helvetica",base_size = 14) +
  theme(legend.position = "bottom",
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank()) +
        labs( y="dN/dS ratio",title="")



  
  
ns_2L_flybase_binned <- as.data.frame(table(cut(ns_2L_BDGP6,breaks=seq(0,ns_2L_BDGP6[length(ns_2L_BDGP6)],100000))))[,2]
ns_3L_flybase_binned <- as.data.frame(table(cut(ns_3L_BDGP6,breaks=seq(0,ns_3L_BDGP6[length(ns_3L_BDGP6)],100000))))[,2]
ns_3R_flybase_binned <- as.data.frame(table(cut(ns_3R_BDGP6,breaks=seq(0,ns_3R_BDGP6[length(ns_3R_BDGP6)],100000))))[,2]
ns_X_flybase_binned <- as.data.frame(table(cut(ns_X_BDGP6,breaks=seq(0,ns_X_BDGP6[length(ns_X_BDGP6)],100000))))[,2]

s_2L_flybase_binned <- as.data.frame(table(cut(s_2L_BDGP6,breaks=seq(0,s_2L_BDGP6[length(s_2L_BDGP6)],100000))))[,2]
s_3L_flybase_binned <- as.data.frame(table(cut(s_3L_BDGP6,breaks=seq(0,s_3L_BDGP6[length(s_3L_BDGP6)],100000))))[,2]
s_3R_flybase_binned <- as.data.frame(table(cut(s_3R_BDGP6,breaks=seq(0,s_3R_BDGP6[length(s_3R_BDGP6)],100000))))[,2]
s_X_flybase_binned <- as.data.frame(table(cut(s_X_BDGP6,breaks=seq(0,s_X_BDGP6[length(s_X_BDGP6)],100000))))[,2]


ratio_ns_cm_3R <- as.data.frame(cbind(recom_3R[,1],ns_3R_flybase_binned[41:length(ns_3R_flybase_binned)]))
ratio_ns_cm_3R$ns_cm <- ratio_ns_cm_3R$V2/ratio_ns_cm_3R$V1

