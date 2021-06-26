library(ggridges)
variables_number <- 37
number_of_stats_calculated <- 25
h <- "0.3"
path.folder_sim <- paste0(getwd(),"/","gral_sim_res_copy/")
files_ave <- paste0(path.folder_sim,"/",dir(path.folder_sim,pattern = "^final_stats_average"))
path.folder_sim_neutral <- paste0(getwd(),"/","neutral_gral/")
files_ave_neutral <- paste0(path.folder_sim_neutral,"/",dir(path.folder_sim,pattern = "^final_stats_average"))
get_number_loci <-  read.table(files_ave[1],header=F,sep=",",skip= variables_number,fill=T)
number_loci <- ncol(get_number_loci)

final_fst <- as.data.frame(matrix(ncol = 10))
colnames(final_fst) <- c("Fst","Fstp","Ne","c","h","s","mean_Fst","mean_Fstp","max_Fst","min_Fst")

false_positive <- as.data.frame(matrix(ncol = 5))
colnames(false_positive) <- c("fp","Ne","c","h","s")

for(i in 1:length(files_ave)){
  
get_number_generations <- read.table(files_ave[i],header = F,sep=",",nrows = variables_number, colClasses = "character")
number_of_generations <- as.numeric(get_number_generations[7,2])
number_of_lines <- number_of_stats_calculated * number_of_generations
initial_breaks <- seq(1,number_of_lines,number_of_generations)
end_breaks <- c((initial_breaks-1)[-1],number_of_lines)

pop_history <- read.table(files_ave[i],header = F,sep=",",nrows = variables_number, colClasses = "character")

vars <- pop_history[c(6,14,26,27),2]
vars[4] <- as.numeric(vars[4])
vars_2 <- paste0("final_stats_average_GRAL_","B",vars[1],"_C",vars[2],"_")

mean_df <- read.table(files_ave[i],header=F,sep=",",skip= variables_number,fill=T)

files_ave_neutral <- paste0(path.folder_sim_neutral,"/",dir(path.folder_sim_neutral,pattern = vars_2))
mean_df_neutral <- read.table(files_ave_neutral,header=F,sep=",",skip= variables_number,fill=T)

mean_Fst <-  mean_df[initial_breaks[18]:end_breaks[18],1:number_loci]
mean_Fstp <-  mean_df[initial_breaks[19]:end_breaks[19],1:number_loci]
temp_fst <- unname(unlist(mean_Fst[number_of_generations,]))
temp_fstp <- unname(unlist(mean_Fstp[number_of_generations,]))

mean_Fst_neutral <-  mean_df_neutral[initial_breaks[18]:end_breaks[18],1:number_loci]
mean_Fstp_neutral <-  mean_df_neutral[initial_breaks[19]:end_breaks[19],1:number_loci]
temp_fst_neutral <- unname(unlist(mean_Fst_neutral[number_of_generations,]))
temp_fstp_neutral <- unname(unlist(mean_Fstp_neutral[number_of_generations,]))


temp_fst_2_neutral <- as.data.frame(cbind(temp_fst_neutral,temp_fstp_neutral,vars[1],vars[2],h,"0",mean(temp_fst_neutral),mean(temp_fstp_neutral),max(temp_fst_neutral)*1.05,min(temp_fst_neutral)*0.95))
colnames(temp_fst_2_neutral) <- c("Fst","Fstp","Ne","c","h","s","mean_Fst","mean_Fstp","max_Fst","min_Fst")

temp_fst_2 <- as.data.frame(cbind(temp_fst,temp_fstp,vars[1],vars[2],vars[3],vars[4],mean(temp_fst),mean(temp_fstp),max(temp_fst)*1.05,min(temp_fst)*0.95))
colnames(temp_fst_2) <- c("Fst","Fstp","Ne","c","h","s","mean_Fst","mean_Fstp","max_Fst","min_Fst")

false_positive_temp <- density(temp_fst)
false_positive_y <- max(false_positive_temp$y)
false_positive_x <- mean(temp_fstp_neutral)
if(mean(temp_fst)<mean(temp_fst_neutral)){
seq_fp<- seq(from=mean(temp_fst) , to=mean(temp_fst_neutral),by=0.001)
}
if(mean(temp_fst)>mean(temp_fst_neutral)){
seq_fp<- mean(temp_fst_neutral)
}

false_positive_temp_2 <- as.data.frame(cbind(seq_fp ,vars[1],vars[2],vars[3],vars[4]))
colnames(false_positive_temp_2) <- c("fp","Ne","c","h","s")
# false_positive_temp_3 <- as.data.frame(cbind(mean(temp_fstp_neutral),vars[1],vars[2],vars[3],vars[4]))
# colnames(false_positive_temp_3) <- c("fp","Ne","c","h","s")
# false_positive_temp_4 <- rbind(false_positive_temp_2,false_positive_temp_3)

false_positive <- rbind(false_positive,false_positive_temp_2)

temp_fst_3 <- rbind(temp_fst_2,temp_fst_2_neutral)

final_fst <- rbind(final_fst,temp_fst_3)

}

final_fst$Fst <- as.numeric(final_fst$Fst)
final_fst$Fstp <- as.numeric(final_fst$Fstp)
final_fst$mean_Fst <- as.numeric(final_fst$mean_Fst)
final_fst$mean_Fstp <- as.numeric(final_fst$mean_Fstp)

plot_density_fst <- final_fst[which(final_fst$h==h),]
plot_density_fst <- plot_density_fst[-which(plot_density_fst$Ne==150),]
plot_density_fst$Ne <- as.numeric(plot_density_fst$Ne)
plot_density_fst[which(plot_density_fst$c=="0.25"),"c"] <- "12.5 cM"
plot_density_fst[which(plot_density_fst$c=="0.5"),"c"] <- "25 cM"
plot_density_fst[which(plot_density_fst$c=="1"),"c"] <- "50 cM"
plot_density_fst[which(plot_density_fst$c=="10"),"c"] <- "500 cM"

options(scipen = 999)
plot_density_fst$s <- as.numeric(plot_density_fst$s)
plot_density_fst$s <- as.factor(plot_density_fst$s)
plot_density_fst$max_Fst <- as.numeric(plot_density_fst$max_Fst)
plot_density_fst$min_Fst <- as.numeric(plot_density_fst$min_Fst)

true_positive<- plot_density_fst[which(plot_density_fst$s=="0"),]
true_positive<- true_positive[!duplicated(true_positive$mean_Fst),]
true_positive$max_Fst <- as.numeric(true_positive$max_Fst)

# false_positive <- plot_density_fst[which(plot_density_fst$s=="0.005"),]


true_positive[which(true_positive$Ne==10),"max_Fst"] <- as.numeric(max(true_positive[which(true_positive$Ne==10),"max_Fst"]))
 true_positive[which(true_positive$Ne==50),"max_Fst"] <- as.numeric(max(true_positive[which(true_positive$Ne==50),"max_Fst"]))
 true_positive[which(true_positive$Ne==100),"max_Fst"] <- as.numeric(max(true_positive[which(true_positive$Ne==100),"max_Fst"]))
 true_positive[which(true_positive$Ne==200),"max_Fst"] <- as.numeric(max(true_positive[which(true_positive$Ne==200),"max_Fst"]))
# 
# false_positive[which(false_positive$Ne==10),"min_Fst"] <- as.numeric(min(false_positive[which(false_positive$Ne==10),"min_Fst"]))
# false_positive[which(false_positive$Ne==50),"min_Fst"] <- as.numeric(min(false_positive[which(false_positive$Ne==50),"min_Fst"]))
#  false_positive[which(false_positive$Ne==100),"min_Fst"] <- as.numeric(min(false_positive[which(false_positive$Ne==100),"min_Fst"]))
#  false_positive[which(false_positive$Ne==200),"min_Fst"] <- as.numeric(min(false_positive[which(false_positive$Ne==200),"min_Fst"]))
# # 
# 
 false_positive$s <- as.factor( false_positive$s)
 false_positive$Ne <- as.numeric(false_positive$Ne)
 false_positive$fp <- as.numeric( false_positive$fp)
 
 false_positive <- false_positive[which(false_positive$h==h),]
false_positive <- false_positive[-which(false_positive$Ne==150),]
false_positive$Ne <- as.numeric(false_positive$Ne)
false_positive[which(false_positive$c=="0.25"),"c"] <- "12.5 cM"
false_positive[which(false_positive$c=="0.5"),"c"] <- "25 cM"
false_positive[which(false_positive$c=="1"),"c"] <- "50 cM"
false_positive[which(false_positive$c=="10"),"c"] <- "500 cM"

# false_positive$fp <- round(false_positive$fp,2)

ggplot() +
geom_rect(data=true_positive,aes(xmin=mean_Fst, xmax=max_Fst, ymin=1,ymax=Inf), alpha=0.3, fill="gray50",color=NA)+
 # geom_rect(data=false_positive,aes(xmin=mean_Fst, xmax=max_Fst, ymin=1,ymax=Inf), alpha=0.3, fill="red",color=NA)+
 # annotate("rect",xmin=false_positive$min_Fst, xmax=false_positive$max_Fst, ymin=10,ymax=1, alpha=0.3, fill="red",color=NA)+  
  # geom_density_ridges(data=false_positive,aes(x=fp,y=s),alpha=0.3,scale = 1,color="red",fill="red",stat="binline",size=7,rel_min_height = 0.0001)+
  geom_density_ridges(data=plot_density_fst,aes(x=Fst, color=s,fill=s,y=s),scale = 1)+
  geom_vline(data =true_positive,aes(xintercept =mean_Fst ),size=0.75, linetype="dashed" )+
  facet_grid(c~Ne, scales = "free") +
  # scale_fill_viridis_b()
   scale_fill_manual(values=c("black","darkgoldenrod2","forestgreen","blue","red4"))+
   scale_color_manual(values=c("black","darkgoldenrod2","forestgreen","blue","red4"))+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  # coord_cartesian(clip = "off") +
  theme_ridges(grid = F, center_axis_labels = TRUE,font_size = 18,font_family="Helvetica")+
   theme(legend.position = "none",strip.background = element_rect(colour="black",fill="white"))+
     labs(y = "Selection coefficient", x = "FST") 


 ggsave("fst_distributions.pdf",  width = 12, height = 8, units = "in", dpi="retina", bg = "transparent"  )

