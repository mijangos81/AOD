library(ggplot2)
library(ggthemes)
library(viridis)
library(pracma)

# K_pops <- seq(2,10,1)
K_pops <- c(2,3,5,10,20)
M <- seq(0,1,0.001)
plot_fst <- as.data.frame(matrix(nrow = 1,ncol = 4))
colnames(plot_fst) <- c("M","Fst_max","HT","pops")
area_fst <- as.data.frame(matrix(nrow =length(K_pops),ncol = 2 ))
colnames(area_fst) <- c("pops","fst_area")
total_area_fst <- 1

for (K in K_pops) {
Fst_max <- (floor(K*M)+(K*M-floor(K*M))^2-K*M^2)/(K*M*(1-M))
HT <- (K*M-floor(K*M))^2
plot_fst_temp <- as.data.frame(cbind(M,Fst_max,HT))
plot_fst_temp$pops <- K
plot_fst_temp <- plot_fst_temp[complete.cases(plot_fst_temp$M),]
plot_fst_temp <- plot_fst_temp[complete.cases(plot_fst_temp$Fst_max),]
# area_fst_temp <- (trapz(plot_fst_temp$M, plot_fst_temp$Fst_max)*100)/total_area_fst
area_fst_temp <- (trapz(plot_fst_temp$HT, plot_fst_temp$Fst_max)*100)/total_area_fst
area_fst_temp_2 <- as.data.frame(cbind(K,area_fst_temp)) 
colnames(area_fst_temp_2) <- c("pops","fst_area")
area_fst <- as.data.frame(rbind(area_fst,area_fst_temp_2))
plot_fst <- as.data.frame(rbind(plot_fst,plot_fst_temp))
}
plot_fst$pops_order <- as.numeric(plot_fst$pops)
plot_fst$pops <- reorder(plot_fst$pops, plot_fst$pops_order, function(x) -max(x) )
# plot_fst$pops_order <- as.factor(plot_fst$pops)
plot_fst <- plot_fst[complete.cases(plot_fst$M),]

colors_fst_max <- c(
  "2"="deepskyblue",
  "3"="deepskyblue3",
  "5"="black",
  "10"="green",
  "100"="green3"
   )

# plot_fst <- plot_fst[which(plot_fst$pops=="2"),]
ggplot(plot_fst[plot_fst$pops==20,],aes(x= M, y= Fst_max,fill= pops,color= pops))+
  # ggplot(plot_fst,aes(x= HT, y= Fst_max,fill= pops,color= pops))+
 geom_area(position = "identity",alpha=0.25,size=2)+
  geom_vline(xintercept = 0.5,color="black",size=2)+
 # guides(color=guide_legend("fct_reorder(pops, M, .desc = T)"), color = FALSE)
  # guides(color=guide_legend(override.aes=list(fill=NA))) +
# scale_fill_manual(name="",values=colors_fst_max) +
scale_color_viridis(discrete = T) +
scale_fill_viridis(discrete = T,guide=F) +
xlab("Frequency of the most frequent allele") + 
ylab("FST")+
 # xlim(c(0,0.5))+
 # ggtitle("Ten populations")+
labs(color = "Number of\nsubpopulations")+
theme_tufte(base_family="Helvetica")+
theme(  legend.position="bottom",
      text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA))

 ggsave("fdist_max_20_multi.pdf",  width = 7, height = 5, units = "in", dpi="retina", bg = "transparent"  )
 
area_fst <- area_fst[complete.cases(area_fst$pops),]

ggplot(area_fst,aes(x= pops, y= fst_area))+
 geom_area(alpha=0.25,fill="black")+
  geom_line(alpha=0.75,size=2,color="black")+
   geom_point(color="deeppink",size=4)+
 # guides(color=guide_legend("fct_reorder(pops, M, .desc = T)"), color = FALSE)
  # guides(color=guide_legend(override.aes=list(fill=NA))) +
# scale_fill_manual(name="",values=colors_fst_max) +
# scale_color_viridis(discrete = T) +
# scale_fill_viridis(discrete = T,guide=F) +
    scale_x_continuous(breaks = seq(2, 10, 1)) +
  scale_y_continuous(breaks =seq(0, 100, 10))+
  # scale_x_discrete(limits = as.factor(area_fst$pops) )+
xlab("Number of subpopulations") + 
ylab("Percentage of possible FST values")+
# ylim(c(0,1))+
 # ggtitle("Ten populations")+
# labs(color = "Number of\nsubpopulations")+
  theme_gdocs(base_family="Helvetica")+
theme(  legend.position="bottom",
      text = element_text(size = 18),
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank(),
        # panel.background = element_rect(fill = "transparent"), # bg of the panel
    # plot.background = element_rect(fill = "transparent", color = NA)
    )

 ggsave("fst_area.pdf",  width = 6, height = 6, units = "in", dpi="retina", bg = "transparent"  )
 

