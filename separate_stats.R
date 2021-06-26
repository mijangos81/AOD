variables_number <- 36
number_of_stats_calculated <- 25# 
 # folders <- dir(path=getwd(),pattern="^final_stats")
  # path.folder_sim <- paste0(getwd(),"/","results_neutral_simulations/")
    # path.folder_sim <- paste0(getwd(),"/","gral_sim_res_copy/")
     # path.folder_sim <- paste0(getwd(),"/","FLY/")
      # path.folder_sim <- paste0(getwd(),"/LD_analyses_fly_sims_neutral")
      path.folder_sim <- getwd()
       # path.folder_sim <- "/Users/mijangos/gral_sim_res"
files_ave <- paste0(path.folder_sim,"/",dir(path.folder_sim,pattern = "^final_stats_average"))
# files_ave <- paste0(path.folder_sim,"/",dir(path.folder_sim,pattern = "^generations"))

# files_ave <- paste0(getwd(),"/","results_neutral_simulations",dir(pattern = "final_stats_average"))
get_number_loci <-  read.table(files_ave[1],header=F,sep=",",skip= variables_number,fill=T)
number_loci <- ncol(get_number_loci)

final_res <- as.data.frame(matrix(nrow = length(files_ave) , ncol = number_of_stats_calculated +1))

for(i in 1:length(files_ave)){
  
get_number_generations <- read.table(files_ave[i],header = F,sep=",",nrows = variables_number, colClasses = "character")
number_of_generations <- as.numeric(get_number_generations[7,2])
number_of_lines <- number_of_stats_calculated * number_of_generations
initial_breaks <- seq(1,number_of_lines,number_of_generations)
end_breaks <- c((initial_breaks-1)[-1],number_of_lines)

pop_history <- read.table(files_ave[i],header = F,sep=",",nrows = variables_number, colClasses = "character")
# vars <- pop_history[c(6,14,26,27,18),2]
# vars_2 <- paste0("Ne_",vars[1],"_c_",vars[2],"_h_",vars[3],"_s_",vars[4],"_sel_",vars[5])

vars <- pop_history[c(17,14,26,27),2]
vars_2 <- paste0("Ne_",vars[1],"_c_",vars[2],"_h_",vars[3],"_s_",vars[4])

mean_df <- read.table(files_ave[i],header=F,sep=",",skip= variables_number,fill=T)

mean_AFD <- mean_df[initial_breaks[1]:end_breaks[1],1:number_loci]
mean_sha_diff <- mean_df[initial_breaks[2]:end_breaks[2],1:number_loci]
mean_one_D_alpha_pop1  <-  mean_df[initial_breaks[3]:end_breaks[3],1:number_loci]
mean_one_D_alpha_pop2 <-  mean_df[initial_breaks[4]:end_breaks[4],1:number_loci]
mean_two_D_alpha_pop1 <-  mean_df[initial_breaks[5]:end_breaks[5],1:number_loci]
mean_two_D_alpha_pop2 <-  mean_df[initial_breaks[6]:end_breaks[6],1:number_loci]
mean_one_D_beta <-  mean_df[initial_breaks[7]:end_breaks[7],1:number_loci]
mean_two_D_beta <-  mean_df[initial_breaks[8]:end_breaks[8],1:number_loci]
mean_sha_pop1 <-  mean_df[initial_breaks[9]:end_breaks[9],1:number_loci]
mean_sha_pop2 <-  mean_df[initial_breaks[10]:end_breaks[10],1:number_loci]
mean_pops <- mean_df[initial_breaks[11]:end_breaks[11],1:number_loci]
mean_shua <-  mean_df[initial_breaks[12]:end_breaks[12],1:number_loci]
mean_expected_het_pop1 <-  mean_df[initial_breaks[13]:end_breaks[13],1:number_loci]
mean_expected_het_pop2 <-  mean_df[initial_breaks[14]:end_breaks[14],1:number_loci]
mean_het_pop1 <-  mean_df[initial_breaks[15]:end_breaks[15],1:number_loci]
mean_het_pop2 <-  mean_df[initial_breaks[16]:end_breaks[16],1:number_loci]
mean_Ht <-  mean_df[initial_breaks[17]:end_breaks[17],1:number_loci]
mean_Fst <-  mean_df[initial_breaks[18]:end_breaks[18],1:number_loci]
mean_Fstp <-  mean_df[initial_breaks[19]:end_breaks[19],1:number_loci]
mean_Dest <-  mean_df[initial_breaks[20]:end_breaks[20],1:number_loci]
mean_deleterious_eliminated <-  mean_df[initial_breaks[21]:end_breaks[21],1]
mean_mean_s_final <-  mean_df[initial_breaks[22]:end_breaks[22],1]
mean_mean_h_final <-  mean_df[initial_breaks[23]:end_breaks[23],1]
mean_g_load_a_final <-  mean_df[initial_breaks[24]:end_breaks[24],1]
mean_g_load_m_final <-  mean_df[initial_breaks[25]:end_breaks[25],1]
generations <- as.data.frame(matrix(c(1:number_of_generations), ncol = 1,nrow = number_of_generations))
# generations <- as.data.frame(matrix(c(20:34), ncol = 1,nrow = 15))
# 
#  # print(
#       # ggplot(generations,aes(V1)) +
#       #    geom_line(aes(y = mean_shua[20:34,24], colour = "1"),size=1) +
#       #     # geom_line(aes(y = mean_het_pop2[,125], colour = "2"),size=1) +
#       #     # geom_line(aes(y = mean_shua[20:34,50], colour = "3"),size=1) +
#       #     geom_line(aes(y = mean_shua[20:34,70], colour = "4"),size=1) +
#       #     geom_line(aes(y = mean_shua[20:34,83], colour = "5"),size=1) +
#       #     # geom_line(aes(y = mean_shua[20:34,110], colour = "6"),size=1) +
#       #     geom_line(aes(y = mean_shua[20:34,132], colour = "7"),size=1) +
#       #    geom_line(aes(y = mean_shua[20:34,146], colour = "8"),size=1) +
#       #    # geom_line(aes(y = mean_expected_het_pop2[,84], colour = "2"),size=1) +
#       #    # geom_line(aes(y = mean_expected_het_pop2[,126], colour = "11"),size=1) +
#       #    # geom_line(aes(y = mean_expected_het_pop2[,151], colour = "9"),size=1) +
#       #    # geom_line(aes(y = mean_expected_het_pop2[,208], colour = "3"),size=1) +
#       #    # geom_point(aes(x = 1,y = 0.536,colour = "2"),size=3) +
#       #    # geom_point(aes(x = 1,y = 0.338,colour = "3"),size=3)+
#       #    # geom_point(aes(x = 1,y = 0.620,colour = "9"),size=3)+
#       #    # geom_point(aes(x = 1,y = 0.582,colour = "11"),size=3)+
#       #    geom_point(aes(x = 34,y = 0.053 ,colour = "1"),size=3) +
#       #    # geom_point(aes(x = 17,y = 0.338,colour = "2"),size=3)+
#       #    # geom_point(aes(x = 34,y = 0.05,colour = "3"),size=3)+
#       #    geom_point(aes(x = 34,y = 0.036,colour = "4"),size=3)+
#       #    geom_point(aes(x = 34,y = 0.063,colour = "5"),size=3) +
#       #     # geom_point(aes(x = 34,y = 0.07 ,colour = "6"),size=3)+
#       #    geom_point(aes(x = 34,y = 0.042,colour = "7"),size=3)+
#       #    geom_point(aes(x = 34,y = 0.046,colour = "8"),size=3)+
#       #    theme_bw(base_size = 14) +
#       #    theme(legend.title=element_blank())+
#       #    theme(legend.position="bottom")+
#       #    xlab("GENERATIONS") +
#       #    ylab("AFD")
#  # )
#  #
# 
#  # print(
#  #      ggplot(generations,aes(V1)) +
#  #         geom_line(aes(y = mean_deleterious_eliminated, colour = "1"),size=1) +
#  #         theme_bw(base_size = 14) +
#  #         theme(legend.title=element_blank())+
#  #         theme(legend.position="bottom")+
#  #         xlab("GENERATIONS") +
#  #         ylab("percent lost")
#  # )
#  #      
#  #              ggsave(paste("percent_lost",i,".pdf"),  width = 7, height =7, units = "in", dpi="retina", bg = "transparent" )
# 
# # print(
#        # ggplot(generations,aes(V1)) +
#        #   geom_line(aes(y = mean_het_pop2[,24], colour = "1"),size=1) +
#        #   # geom_line(aes(y = mean_het_pop2[,125], colour = "2"),size=1) +
#        #   # geom_line(aes(y = mean_AFD[20:34,51], colour = "3"),size=1) +
#        #   # geom_line(aes(y = mean_AFD[20:34,71], colour = "4"),size=1) +
#        #   # geom_line(aes(y = mean_AFD[20:34,81], colour = "5"),size=1) +
#        #    geom_line(aes(y = mean_het_pop2[,110], colour = "6"),size=1) +
#        #   # geom_line(aes(y = mean_AFD[20:34,133], colour = "7"),size=1) +
#        #   geom_line(aes(y = mean_het_pop2[,146], colour = "8"),size=1) +
#        #   # geom_line(aes(y = mean_expected_het_pop2[,84], colour = "2"),size=1) +
#        #   # geom_line(aes(y = mean_expected_het_pop2[,126], colour = "11"),size=1) +
#        #   # geom_line(aes(y = mean_expected_het_pop2[,151], colour = "9"),size=1) +
#        #   # geom_line(aes(y = mean_expected_het_pop2[,208], colour = "3"),size=1) +
#        #   # geom_point(aes(x = 1,y = 0.536,colour = "2"),size=3) +
#        #   # geom_point(aes(x = 1,y = 0.338,colour = "3"),size=3)+
#        #   # geom_point(aes(x = 1,y = 0.620,colour = "9"),size=3)+
#        #   # geom_point(aes(x = 1,y = 0.582,colour = "11"),size=3)+
#        #   geom_point(aes(x = 34,y = 0.295 ,colour = "1"),size=3) +
#        #   # geom_point(aes(x = 17,y = 0.338,colour = "2"),size=3)+
#        #   # geom_point(aes(x = 34,y = 0.155 + 0.15,colour = "3"),size=3)+
#        #   # geom_point(aes(x = 34,y = 0.191 + 0.15,colour = "4"),size=3)+
#        #   # geom_point(aes(x = 34,y = 0.247 + 0.15,colour = "5"),size=3) +
#        #    geom_point(aes(x = 34,y = 0.304 ,colour = "6"),size=3)+
#        #   # geom_point(aes(x = 34,y = 0.151 + 0.15,colour = "7"),size=3)+
#        #   geom_point(aes(x = 34,y = 0.405,colour = "8"),size=3)+
#        #   theme_bw(base_size = 14) +
#        #   theme(legend.title=element_blank())+
#        #   theme(legend.position="bottom")+
#        #   xlab("GENERATIONS") +
#        #   ylab("AFD")
# # )
# #        
# #         ggsave(paste("het_snps",i,".pdf"),  width = 7, height =7, units = "in", dpi="retina", bg = "transparent" )
# # 
# 
#       # 
#    
#   ggplot(generations, aes(V1)) +
#          geom_line(aes(y = mean_AFD[,83], colour = "2"),size=1) +
#          geom_line(aes(y = mean_AFD[,125], colour = "11"),size=1) +
#          geom_line(aes(y = mean_AFD[,150], colour = "9"),size=1) +
#          geom_line(aes(y = mean_AFD[,207], colour = "3"),size=1) +
#   #        geom_line(aes(y = mean_het_pop2[,83], colour = "2"),size=1) +
#   #        geom_line(aes(y = mean_het_pop2[,125], colour = "11"),size=1) +
#   #        geom_line(aes(y = mean_het_pop2[,150], colour = "9"),size=1) +
#   #        geom_line(aes(y = mean_het_pop2[,207], colour = "3"),size=1) +
#   #        geom_point(aes(x = 1,y = 0.536,colour = "2"),size=3) +
#   #        geom_point(aes(x = 1,y = 0.338,colour = "3"),size=3)+
#   #        geom_point(aes(x = 1,y = 0.620,colour = "9"),size=3)+
#   #        geom_point(aes(x = 1,y = 0.582,colour = "11"),size=3)+
#   #        geom_point(aes(x = 17,y = 0.504,colour = "2"),size=3) +
#   #        geom_point(aes(x = 17,y = 0.338,colour = "3"),size=3)+
#   #        geom_point(aes(x = 17,y = 0.521,colour = "9"),size=3)+
#   #        geom_point(aes(x = 17,y = 0.518,colour = "11"),size=3)+
#   #        geom_point(aes(x = 34,y = 0.423,colour = "2"),size=3) +
#   #        geom_point(aes(x = 34,y = 0.361,colour = "3"),size=3)+
#   #        geom_point(aes(x = 34,y = 0.414,colour = "9"),size=3)+
#   #        geom_point(aes(x = 34,y = 0.472,colour = "11"),size=3)+
#          geom_point(aes(x = 34,y = 0.321,colour = "2"),size=3) +
#          geom_point(aes(x = 34,y = 0.187,colour = "3"),size=3)+
#          geom_point(aes(x = 34,y = 0.290,colour = "9"),size=3)+
#          geom_point(aes(x = 34,y = 0.386,colour = "11"),size=3)+
#          theme_bw(base_size = 14) +
#          theme(legend.title=element_blank())+
#          theme(legend.position="bottom")+
#          xlab("GENERATIONS") +
#          ylab("HE")
# #   
#  # print(
#  #    ggplot(generations, aes(V1)) +
#  #         geom_line(aes(y = mean_shua[,83], colour = "2"),size=1) +
#  #         geom_line(aes(y = mean_shua[,125], colour = "11"),size=1) +
#  #         geom_line(aes(y = mean_shua[,150], colour = "9"),size=1) +
#  #         geom_line(aes(y = mean_shua[,207], colour = "3"),size=1) +
#  #         # geom_point(aes(x = 17,y = 0.1,colour = "2"),size=4) +
#  #         # geom_point(aes(x = 17,y = 0.102,colour = "3"),size=4)+
#  #         # geom_point(aes(x = 17,y = 0.130,colour = "9"),size=4)+
#  #         # geom_point(aes(x = 17,y = 0.170,colour = "11"),size=4)+
#  #         geom_point(aes(x = 34,y = 0.1,colour = "2"),size=4) +
#  #         geom_point(aes(x = 34,y = 0.06,colour = "3"),size=4)+
#  #         geom_point(aes(x = 34,y = 0.093,colour = "9"),size=4)+
#  #         geom_point(aes(x = 34,y = 0.172,colour = "11"),size=4)+
#  #         theme_bw(base_size = 14) +
#  #         theme(legend.title=element_blank())+
#  #         theme(legend.position="bottom")
#  # )
#     # 
#    
#    
#    
#     ggplot(generations, aes(V1)) +
#          geom_line(aes(y = mean_Fst[,83], colour = "2"),size=1) +
#          geom_line(aes(y = mean_Fst[,125], colour = "11"),size=1) +
#          geom_line(aes(y = mean_Fst[,150], colour = "9"),size=1) +
#          geom_line(aes(y = mean_Fst[,207], colour = "3"),size=1) +
#          # geom_point(aes(x = 17,y = 0.1,colour = "2"),size=4) +
#          # geom_point(aes(x = 17,y = 0.102,colour = "3"),size=4)+
#          # geom_point(aes(x = 17,y = 0.130,colour = "9"),size=4)+
#          # geom_point(aes(x = 17,y = 0.170,colour = "11"),size=4)+
#          geom_point(aes(x = 34,y = 0.180,colour = "2"),size=4) +
#          geom_point(aes(x = 34,y = 0.088,colour = "3"),size=4)+
#          geom_point(aes(x = 34,y = 0.158,colour = "9"),size=4)+
#          geom_point(aes(x = 34,y = 0.219,colour = "11"),size=4)+
#          theme_bw(base_size = 14) +
#          theme(legend.title=element_blank())+
#          theme(legend.position="bottom")
#     # 
#     # 
#     #   ggplot(generations, aes(V1)) +
#     #      geom_line(aes(y = mean_shua[,83], colour = "2"),size=1) +
#     #      geom_line(aes(y = mean_shua[,125], colour = "11"),size=1) +
#     #      geom_line(aes(y = mean_shua[,150], colour = "9"),size=1) +
#     #      geom_line(aes(y = mean_shua[,207], colour = "3"),size=1) +
#     #      # geom_point(aes(x = 17,y = 0.1,colour = "2"),size=4) +
#     #      # geom_point(aes(x = 17,y = 0.102,colour = "3"),size=4)+
#     #      # geom_point(aes(x = 17,y = 0.130,colour = "9"),size=4)+
#     #      # geom_point(aes(x = 17,y = 0.170,colour = "11"),size=4)+
#     #      geom_point(aes(x = 34,y = 0.1,colour = "2"),size=4) +
#     #      geom_point(aes(x = 34,y = 0.061,colour = "3"),size=4)+
#     #      geom_point(aes(x = 34,y = 0.093,colour = "9"),size=4)+
#     #      geom_point(aes(x = 34,y = 0.172,colour = "11"),size=4)+
#     #      theme_bw(base_size = 14) +
#     #      theme(legend.title=element_blank())+
#     #      theme(legend.position="bottom")
#     #   
#     #   
#     #    ggplot(generations, aes(V1)) +
#     #      geom_line(aes(y = mean_Dest[,83], colour = "2"),size=1) +
#     #      geom_line(aes(y = mean_Dest[,125], colour = "11"),size=1) +
#     #      geom_line(aes(y = mean_Dest[,150], colour = "9"),size=1) +
#     #      geom_line(aes(y = mean_Dest[,207], colour = "3"),size=1) +
#     #      # geom_point(aes(x = 17,y = 0.1,colour = "2"),size=4) +
#     #      # geom_point(aes(x = 17,y = 0.102,colour = "3"),size=4)+
#     #      # geom_point(aes(x = 17,y = 0.130,colour = "9"),size=4)+
#     #      # geom_point(aes(x = 17,y = 0.170,colour = "11"),size=4)+
#     #      geom_point(aes(x = 34,y = 0.161,colour = "2"),size=4) +
#     #      geom_point(aes(x = 34,y = 0.119,colour = "3"),size=4)+
#     #      geom_point(aes(x = 34,y = 0.156,colour = "9"),size=4)+
#     #      geom_point(aes(x = 34,y = 0.278,colour = "11"),size=4)+
#     #      theme_bw(base_size = 14) +
#     #      theme(legend.title=element_blank())+
#     #      theme(legend.position="bottom")
#     #    
#        # print(
#        #        ggplot(generations, aes(V1)) +
#        #   geom_line(aes(y = mean_AFD[,83], colour = "2"),size=1) +
#        #   geom_line(aes(y = mean_AFD[,125], colour = "11"),size=1) +
#        #   geom_line(aes(y = mean_AFD[,150], colour = "9"),size=1) +
#        #   geom_line(aes(y = mean_AFD[,207], colour = "3"),size=1) +
#        #   # geom_point(aes(x = 17,y = 0.1,colour = "2"),size=4) +
#        #   # geom_point(aes(x = 17,y = 0.102,colour = "3"),size=4)+
#        #   # geom_point(aes(x = 17,y = 0.130,colour = "9"),size=4)+
#        #   # geom_point(aes(x = 17,y = 0.170,colour = "11"),size=4)+
#        #   geom_point(aes(x = 34,y = 0.321,colour = "2"),size=4) +
#        #   geom_point(aes(x = 34,y = 0.187,colour = "3"),size=4)+
#        #   geom_point(aes(x = 34,y = 0.290,colour = "9"),size=4)+
#        #   geom_point(aes(x = 34,y = 0.386,colour = "11"),size=4)+
#        #   theme_bw(base_size = 14) +
#        #   theme(legend.title=element_blank())+
#        #   theme(legend.position="bottom")+
#        #   xlab("GENERATIONS") +
#        #   ylab("AFD")
#        # )
#        
#        #  print(
#        #        ggplot(generations, aes(V1)) +
#        #   geom_line(aes(y = mean_AFD[,24], colour = "1"),size=1) +
#        #   # geom_line(aes(y = mean_AFD[,82], colour = "2"),size=1) +
#        #   geom_line(aes(y = mean_AFD[,110], colour = "3"),size=1) +
#        #   geom_line(aes(y = mean_AFD[,132], colour = "4"),size=1) +
#        #   geom_line(aes(y = mean_AFD[,146], colour = "5"),size=1) +
#        #   # geom_point(aes(x = 17,y = 0.1,colour = "2"),size=4) +
#        #   # geom_point(aes(x = 17,y = 0.102,colour = "3"),size=4)+
#        #   # geom_point(aes(x = 17,y = 0.130,colour = "9"),size=4)+
#        #   # geom_point(aes(x = 17,y = 0.170,colour = "11"),size=4)+
#        #   geom_point(aes(x = 34,y = 0.194,colour = "1"),size=4) +
#        #   # geom_point(aes(x = 34,y = 0.247,colour = "2"),size=4)+
#        #   geom_point(aes(x = 34,y = 0.216,colour = "3"),size=4)+
#        #   geom_point(aes(x = 34,y = 0.151,colour = "4"),size=4)+
#        #   geom_point(aes(x = 34,y = 0.243,colour = "5"),size=4)+
#        #   theme_bw(base_size = 14) +
#        #   theme(legend.title=element_blank())+
#        #   theme(legend.position="bottom")+
#        #   xlab("GENERATIONS") +
#        #   ylab("AFD")
#        # )
# 
# 
# # print(
# #               ggplot(generations, aes(V1)) +
# #          geom_line(aes(y = mean_Dest[,24], colour = "1"),size=1) +
# #          # geom_line(aes(y = mean_AFD[,82], colour = "2"),size=1) +
# #          geom_line(aes(y = mean_Dest[,110], colour = "3"),size=1) +
# #          geom_line(aes(y = mean_Dest[,132], colour = "4"),size=1) +
# #          geom_line(aes(y = mean_Dest[,146], colour = "5"),size=1) +
# #          # geom_point(aes(x = 17,y = 0.1,colour = "2"),size=4) +
# #          # geom_point(aes(x = 17,y = 0.102,colour = "3"),size=4)+
# #          # geom_point(aes(x = 17,y = 0.130,colour = "9"),size=4)+
# #          # geom_point(aes(x = 17,y = 0.170,colour = "11"),size=4)+
# #          geom_point(aes(x = 34,y = 0.092,colour = "1"),size=4) +
# #          # geom_point(aes(x = 34,y = 0.247,colour = "2"),size=4)+
# #          geom_point(aes(x = 34,y = 0.119,colour = "3"),size=4)+
# #          geom_point(aes(x = 34,y = 0.059,colour = "4"),size=4)+
# #          geom_point(aes(x = 34,y = 0.117,colour = "5"),size=4)+
# #          theme_bw(base_size = 14) +
# #          theme(legend.title=element_blank())+
# #          theme(legend.position="bottom")+
# #          xlab("GENERATIONS") +
# #          ylab("Jost's D")
# #        )
# 
# 
#               
#               
#             ggplot(generations, aes(V1)) +
#   geom_line(aes(y = mean_sha_pop1[,83], colour = "2"),size=1) +
#   geom_line(aes(y = mean_sha_pop1[,125], colour = "11"),size=1) +
#   geom_line(aes(y = mean_sha_pop1[,150], colour = "9"),size=1) +
#   geom_line(aes(y = mean_sha_pop1[,207], colour = "3"),size=1) +
#   # geom_point(aes(x = 17,y = 0.1,colour = "2"),size=4) +
#   # geom_point(aes(x = 17,y = 0.102,colour = "3"),size=4)+
#   # geom_point(aes(x = 17,y = 0.130,colour = "9"),size=4)+
#   # geom_point(aes(x = 17,y = 0.170,colour = "11"),size=4)+
#   geom_point(aes(x = 34,y = 0.648,colour = "2"),size=4) +
#   geom_point(aes(x = 34,y = 0.568,colour = "3"),size=4)+
#   geom_point(aes(x = 34,y = 0.668,colour = "9"),size=4)+
#   geom_point(aes(x = 34,y = 0.788,colour = "11"),size=4)+
#   theme_bw(base_size = 14) +
#   theme(legend.title=element_blank())+
#   theme(legend.position="bottom")
#   +
#   xlab("GENERATIONS") +
#   ylab("FSTp")
#          # 
#   
#  ggplot(generations, aes(V1)) +
#          geom_line(aes(y = mean_Fst[,24], colour = "1"),size=1) +
#          geom_line(aes(y = mean_Fst[,50], colour = "3"),size=1) +
#          geom_line(aes(y = mean_Fst[,70], colour = "4"),size=1) +
#          geom_line(aes(y = mean_Fst[,82], colour = "5"),size=1) +
#          geom_line(aes(y = mean_Fst[,110], colour = "6"),size=1) +
#          geom_line(aes(y = mean_Fst[,132], colour = "7"),size=1) +
#          geom_line(aes(y = mean_Fst[,146], colour = "8"),size=1) +
#          geom_point(aes(x = 34,y = 0.160,colour = "1"),size=4) +
#          geom_point(aes(x = 34,y = 0.110,colour = "3"),size=4)+
#          geom_point(aes(x = 34,y = 0.100,colour = "4"),size=4)+
#          geom_point(aes(x = 34,y = 0.185,colour = "5"),size=4) +
#          geom_point(aes(x = 34,y = 0.168,colour = "6"),size=4)+
#          geom_point(aes(x = 34,y = 0.129,colour = "7"),size=4)+
#          geom_point(aes(x = 34,y = 0.133,colour = "8"),size=4)+
#          theme_bw(base_size = 14) +
#          theme(legend.title=element_blank())+
#          theme(legend.position="bottom") +
#          # xlab("GENERATIONS") +
#  ylab("FSTp")
#  )
# 
# ggplot(generations, aes(V1)) +
#         # geom_line(aes(y = mean_sha_diff[,24], colour = "1"),size=1) +
#         # geom_line(aes(y = mean_sha_diff[,50], colour = "3"),size=1) +
#          geom_line(aes(y = mean_sha_diff[,70], colour = "4"),size=1) +
#          geom_line(aes(y = mean_sha_diff[,83], colour = "5"),size=1) +
#         # geom_line(aes(y = mean_sha_diff[,110], colour = "6"),size=1) +
#          geom_line(aes(y = mean_sha_diff[,132], colour = "7"),size=1) +
#         geom_line(aes(y = mean_sha_diff[,146], colour = "8"),size=1) +
#         # geom_point(aes(x = 34,y = 0.077 ,colour = "1"),size=4) +
#         # geom_point(aes(x = 34,y = 0.072,colour = "3"),size=4)+
#          geom_point(aes(x = 34,y = 0.052,colour = "4"),size=4)+
#          geom_point(aes(x = 34,y = 0.091,colour = "5"),size=4) +
#         # geom_point(aes(x = 34,y = 0.101,colour = "6"),size=4)+
#         geom_point(aes(x = 34,y = 0.06,colour = "7"),size=4)+
#         geom_point(aes(x = 34,y = 0.066,colour = "8"),size=4)+
#         theme_bw(base_size = 14) +
#         theme(legend.title=element_blank())+
#         theme(legend.position="bottom") +
#         # xlab("GENERATIONS") +
# ylab("Shannon_diff")
# 

  # ggsave(paste("Dest_snps",i,".pdf"),  width = 7, height =7, units = "in", dpi="retina", bg = "transparent" )
# final_res[i,2] <- rowMeans(mean_AFD, na.rm = T)[number_of_generations]
first_gen <- 1
final_res[i,1] <- vars_2
final_res[i,2] <- rowMeans(mean_AFD, na.rm = T)[number_of_generations]
final_res[i,3] <- rowMeans(mean_sha_diff, na.rm = T)[number_of_generations]
final_res[i,4] <- rowMeans(mean_one_D_alpha_pop1, na.rm = T)[number_of_generations]
final_res[i,5] <- rowMeans(mean_one_D_alpha_pop2 , na.rm = T)[number_of_generations]
final_res[i,6] <- rowMeans(mean_two_D_alpha_pop1, na.rm = T)[number_of_generations]
final_res[i,7] <- rowMeans(mean_two_D_alpha_pop2 , na.rm = T)[number_of_generations]
final_res[i,8] <- rowMeans(mean_one_D_beta , na.rm = T)[number_of_generations]
final_res[i,9] <- rowMeans(mean_two_D_beta , na.rm = T)[number_of_generations]
final_res[i,10] <- rowMeans(mean_sha_pop1, na.rm = T)[number_of_generations]
final_res[i,11] <- rowMeans(mean_sha_pop2 , na.rm = T)[number_of_generations]
final_res[i,12] <- rowMeans(mean_pops , na.rm = T)[number_of_generations]
final_res[i,13] <- rowMeans(mean_shua , na.rm = T)[number_of_generations]
final_res[i,14] <- rowMeans(mean_expected_het_pop1, na.rm = T)[number_of_generations]
final_res[i,15] <- rowMeans(mean_expected_het_pop2, na.rm = T)[number_of_generations]
final_res[i,16] <- rowMeans(mean_het_pop1 , na.rm = T)[number_of_generations]
final_res[i,17] <- rowMeans(mean_het_pop2, na.rm = T)[number_of_generations]
final_res[i,18] <- rowMeans(mean_Ht , na.rm = T)[number_of_generations]
final_res[i,19] <- rowMeans(mean_Fst , na.rm = T)[number_of_generations]
final_res[i,20] <- rowMeans(mean_Fstp , na.rm = T)[number_of_generations]
final_res[i,21] <- rowMeans(mean_Dest , na.rm = T)[number_of_generations]
final_res[i,22] <- mean_deleterious_eliminated[number_of_generations]
final_res[i,23] <- mean_mean_s_final[number_of_generations]
final_res[i,24] <- mean_mean_h_final[number_of_generations]
final_res[i,25] <- mean_g_load_a_final[number_of_generations]
final_res[i,26] <- mean_g_load_m_final[number_of_generations]
}

colnames(final_res) <- c(
  "vars",
  "AFD",
    "sha_diff",
    "one_D_alpha_pop1",
    "one_D_alpha_pop2",
    "two_D_alpha_pop1",
    "two_D_alpha_pop2",
    "one_D_beta",
    "two_D_beta",
    "sha_pop1",
    "sha_pop2",
    "sha_pops",
    "shua",
    "expected_het_pop1",
    "expected_het_pop2",
    "het_pop1",
    "het_pop2",
    "Ht",
    "Fst",
    "Fstp",
    "Dest",
    "deleterious_eliminated",
    "mean_s_final",
    "mean_h_final",
    "g_load_a_final", 
    "g_load_m_final")

write.csv(final_res,file = "final_res_stats_fly_last_gen.csv")
