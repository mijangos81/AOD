library(ggplot2)
library(ggthemes)
library(viridis)
library(scales)
library(gtools)
library(readr)
library(png)
library(readxl)
library(wesanderson)
library(scico)
library(RColorBrewer)
library(stringr)
library(tidyr)
library(readr)

final_res_stats<- read_csv("final_res_stats_fly_fisrt_gen.csv")

final_res_stats <- final_res_stats[,c("x","y","load_difference")]
# final_res_stats <- final_res_stats[,c("x","y","Heterozygosity","het_mean_OBS")]
levels_plot <-12
final_res_stats$load_m <- round(final_res_stats$load_difference,2)
final_res_stats$groups <- quantcut(final_res_stats$load_m,q=levels_plot)
# final_res_stats$het_mean_OBS <- round(final_res_stats$het_mean_OBS,2)
# final_res_stats$Heterozygosity <- round(final_res_stats$Heterozygosity,2)
# final_res_stats$groups <- quantcut(final_res_stats$Heterozygosity,q=levels_plot)

# final_res_stats <- final_res_stats[,c("x","y","Fstp_BIAS","Fstp_OBS")]
# levels_plot <-8
# final_res_stats$Fstp_BIAS <- round(final_res_stats$Fstp_BIAS,2)
# final_res_stats$Fstp_OBS <- round(final_res_stats$Fstp_OBS,2)
# final_res_stats$groups <- quantcut(final_res_stats$Fstp_BIAS,q=levels_plot)

final_res_stats$groups <- gsub("[][()]", "", final_res_stats$groups ,",")
final_res_stats$groups <- gsub(",", " to ",final_res_stats$groups)

# final_res_stats[which(final_res_stats$`g_load_m_first`<0),4]<- paste0(as.character(min(final_res_stats$`g_load_m_first`))," to -0.01")
# final_res_stats[which(final_res_stats$`g_load_m_first`>0),4]<- paste0(" 0.01 to ",as.character(max(final_res_stats$`g_load_m_first`)))

final_res_stats$groups <- as.factor(final_res_stats$groups )
groups_n <- length(unique(final_res_stats$groups ))
levels_plot_2 <- levels(final_res_stats$groups)
levels_plot_2 <- as.character(levels_plot_2)
levels_plot_2 <- as.data.frame(levels_plot_2,stringsAsFactors=F)
levels_plot_2 <-separate(levels_plot_2,1,into=c("i","j","h"),sep = " ")
levels_plot_2 <-  levels_plot_2[,c("i","h")]
levels_plot_2 <- as.data.frame(apply(levels_plot_2,2,as.numeric))
levels_plot_2$sum <- levels_plot_2$i + levels_plot_2$i
levels_plot_2$levels <- levels(final_res_stats$groups)
levels_plot_2 <- levels_plot_2[order(levels_plot_2$sum),]
levels_plot_2 <- levels_plot_2$levels
final_res_stats$stat <- factor(final_res_stats$groups,levels = levels_plot_2)

# pal <- rev(wes_palette("Royal2", groups_n, type = "continuous"))
# pal <- cm.colors(groups_n)
pal <- rev(terrain.colors(groups_n))
# pal <- terrain.colors(groups_n)

 colors_3 <- pal

# colors_3 <- c("turquoise","slategrey",viridis::magma(groups_n-2,begin = 0.15,direction = 1))
 # colors_3 <- c("turquoise",viridis::magma(groups_n-1,begin = 0.15,direction = 1))
   # colors_3 <- c(viridis::plasma(groups_n,direction = 1))
  # colors_3 <- c(brewer.pal(groups_n, "Oranges"))

  # colors_3 <- c(viridis::magma(groups_n-2,begin = 0.15,direction = -1),"slategrey","turquoise")
# colors_3 <- c(viridis::magma(groups_n-1,begin = 0.15,direction = -1),"turquoise")

print(
ggplot(data = final_res_stats, aes(x=x,y=y,fill=stat))+
  geom_tile(colour="black" ,size=0.7) + 
    geom_text(aes(label = load_m),color="black", size=4,, fontface = "bold") +
  # geom_text(aes(label = het_mean_OBS),color="black", size=4,, fontface = "bold") +
    # geom_text(aes(label = Fstp_OBS),color="black", size=4,, fontface = "bold") +

     scale_fill_manual(values=setNames(colors_3,levels_plot_2),name="")+
   theme_tufte(base_family="Helvetica",base_size=20) +
   scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
  coord_equal() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.y=element_blank() ,
         axis.ticks.x=element_blank() ,
         axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position="bottom",
panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA))) # bg of the plot)+

   # ggsave("latin_fst.pdf",  width = 8, height =8, units = "in", dpi="retina", bg = "transparent" )
   ggsave("load_m_difference.pdf",  width = 8, height =8, units = "in", dpi="retina", bg = "transparent" )

   
