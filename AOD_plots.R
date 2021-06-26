##### LD PLOTS #####
# if (LD_analyses == TRUE) {
#     # PAIRWISE PLOT

final_pairwise_ld <-  read_csv("final_LD_pairwise_FLY_B_0.005_C_0.2_D_5e-05_F_TRUE.csv",skip = variables_number)
 p4 <- print(
   ggplot(final_pairwise_ld, aes(distance,rsqr)) +
    geom_line() +
    geom_smooth(method = "loess") +
    geom_hline(aes(yintercept=0.2,colour= "red"),size=1) +
  labs(x="Base pairs", y="R-squared", title="Pairwise LD")+
theme_tufte(base_family="Helvetica")+
theme(legend.position = "none"
      # , plot.margin=unit(c(-1,1,1,-1), "in")
      )
  )

# CIRCLE PLOT
# link_final_plot <- as.data.frame(cbind(loc=as.numeric(row.names(link_final)), link_final))
# brk<- seq(0,chromosome_length,region_size)
# sectors <- length(brk)-1
# ld_window <- (ld_max_pairwise/ld_resolution)
# ld_final <- NULL
# for (ld in 2:(ld_window+1)) {
#   ld_temp <- as.numeric(stats.bin(link_final_plot[,1], link_final_plot[,ld], breaks = brk)[[3]][2,])
#   ld_final <- cbind(ld_final,ld_temp)
# }
# rsquare_ceiling<- 0.2
# ld_final[ld_final> rsquare_ceiling] <- rsquare_ceiling
# ld_final <- rescale(ld_final,to=c(0,2))
# ld_final_b <- as.data.frame(ld_final)
# ld_final_b <- ld_final_b[,c(ld_window:1)]
# 
# heatmap <- as.data.frame(matrix(nrow = sectors))
# heatmap[,1] <- chromosome_name
# heatmap_c <- cbind(heatmap,loc=as.numeric(1:sectors),ld_final_b)
# 
# segm <- as.data.frame(cbind(seg.name=chromosome_name,seg.Start=0:(sectors-1),seg.End=1:sectors,the.v="NA",NO="NA"))
# segm_b <- segAnglePo(segm,seg=chromosome_name,angle.start=0,angle.end=180)
# 
# posi<- seq(region_size/2,chromosome_length,region_size)/map_resolution
# NS_pop1<- as.numeric(unlist(mean_het_pop1[gen_number_dispersal,]))
# NS_pop1<- NS_pop1[posi]
# NS_exp_pop1 <- as.numeric(unlist(mean_expected_het_pop1[gen_number_dispersal,]))
# NS_exp_pop1 <- NS_exp_pop1[posi]
# NS_pop1 <- (NS_pop1 - NS_exp_pop1) / NS_exp_pop1
# NS_pop1 <- rescale(NS_pop1,to=c(0,2))
# 
# NS_pop2 <- as.numeric(unlist(mean_het_pop2[gen_number_dispersal,]))
# NS_pop2 <- NS_pop2[posi]
# NS_exp_pop2 <- as.numeric(unlist(mean_expected_het_pop2[gen_number_dispersal,]))
# NS_exp_pop2 <- NS_exp_pop2[posi]
# NS_pop2 <- (NS_pop2 - NS_exp_pop2) / NS_exp_pop2
# NS_pop2<- rescale(NS_pop2,to=c(0,2),)
# 
# NS_fst_temp <- as.numeric(unlist(mean_Fst[gen_number_dispersal,]))
# NS_fst_temp <- NS_fst_temp[posi]
# NS_fst_temp <- (NS_fst_temp - Fst_expected)/ Fst_expected
# NS_fst_temp <- NS_fst_temp * -1
# NS_fst_temp <- rescale(NS_fst_temp,to=c(0,2))
# 
# nonsyn <- as.data.frame(cbind(loc=rep(1:sectors, each=region_size/map_resolution),number_nonsyn_b))
# nonsyn <- vaggregate(.value = nonsyn$number_nonsyn_b,.group=nonsyn$loc, .fun = sum)
# nonsyn <- rescale(nonsyn,to=c(0,2))
# 
# NS<- as.data.frame( matrix(nrow = sectors) )
# NS[,1] <- chromosome_name
# NS$position <-  1:sectors
# NS$NS <- nonsyn
# NS$fst <- NS_fst_temp
# NS$het_pop2<-unlist(NS_pop2)
# NS$het_pop1<-unlist(NS_pop1)
# 
# reg_ld <- rowMeans(ld_final_b,na.rm = T)
# regre <- cbind(NS,reg_ld)
# regre_res <- summary(lm(regre$fst~regre$NS+regre$reg_ld))
# 
# grid <- as.data.frame( matrix(nrow = sectors) )
# grid[,1] <- chromosome_name
# grid$position <-  0.5:(sectors-0.5)
# grid$grid <- 1
# grid$arc <- 1
# 
# par(mar=c(2,2,2,2))
# plot(c(1,800),c(1,800),type="n",axes=FALSE,xlab="",ylab="",main="")
# circos(R=180,cir=segm_b,W=190,mapping=NS,col.v=3,type="ss",B=F,scale=F,cex=3)
# circos(R=350,cir=segm_b,W=60,mapping=heatmap_c,col.v=3,type="heatmap",B=F,col.bar=F,lwd = 32)
# circos(R=81,cir=segm_b,W=150,mapping=NS,col.v=4,type="heatmap",B=F,col.bar =T,col.bar.po = "topleft",lwd = 17)
# }
##### GENERATIONS PLOTS ####
# generations <- as.data.frame(matrix(c(1:gen_number_dispersal), ncol = 1,nrow = gen_number_dispersal))
generations <- as.data.frame(matrix(c(1:number_generations), ncol = 1,nrow = number_generations))


# if(simulation_type=="fly"){
# print (
  ggplot(generations, aes(V1)) +
         geom_line(aes(y = mean_het_pop1[,83], colour = "2"),size=1) +
         geom_line(aes(y = mean_het_pop1[,125], colour = "11"),size=1) +
         geom_line(aes(y = mean_het_pop1[,150], colour = "9"),size=1) +
         geom_line(aes(y = mean_het_pop1[,207], colour = "3"),size=1) +
         # geom_line(aes(y = mean_expected_het_pop2[,84], colour = "2"),size=1) +
         # geom_line(aes(y = mean_expected_het_pop2[,126], colour = "11"),size=1) +
         # geom_line(aes(y = mean_expected_het_pop2[,151], colour = "9"),size=1) +
         # geom_line(aes(y = mean_expected_het_pop2[,208], colour = "3"),size=1) +
         geom_point(aes(x = 1,y = 0.536,colour = "2"),size=3) +
         geom_point(aes(x = 1,y = 0.338,colour = "3"),size=3)+
         geom_point(aes(x = 1,y = 0.620,colour = "9"),size=3)+
         geom_point(aes(x = 1,y = 0.582,colour = "11"),size=3)+
         geom_point(aes(x = 17,y = 0.504,colour = "2"),size=3) +
         geom_point(aes(x = 17,y = 0.338,colour = "3"),size=3)+
         geom_point(aes(x = 17,y = 0.521,colour = "9"),size=3)+
         geom_point(aes(x = 17,y = 0.518,colour = "11"),size=3)+
         geom_point(aes(x = 34,y = 0.423,colour = "2"),size=3) +
         geom_point(aes(x = 34,y = 0.361,colour = "3"),size=3)+
         geom_point(aes(x = 34,y = 0.414,colour = "9"),size=3)+
         geom_point(aes(x = 34,y = 0.472,colour = "11"),size=3)+
         theme_bw(base_size = 14) +
         theme(legend.title=element_blank())+
         theme(legend.position="bottom")+
         xlab("GENERATIONS") +
         ylab("HE")  
  
  
   ggplot(generations, aes(V1)) +
         geom_line(aes(y = mean_het_pop1[,24], colour = "1"),size=1) +
         # geom_line(aes(y = mean_het_pop2[,125], colour = "2"),size=1) +
         geom_line(aes(y = mean_het_pop1[,50], colour = "3"),size=1) +
         geom_line(aes(y = mean_het_pop1[,70], colour = "4"),size=1) +
         geom_line(aes(y = mean_het_pop1[,82], colour = "5"),size=1) +
         geom_line(aes(y = mean_het_pop1[,110], colour = "6"),size=1) +
         geom_line(aes(y = mean_het_pop1[,132], colour = "7"),size=1) +
         geom_line(aes(y = mean_het_pop1[,146], colour = "8"),size=1) +
         # geom_line(aes(y = mean_expected_het_pop2[,84], colour = "2"),size=1) +
         # geom_line(aes(y = mean_expected_het_pop2[,126], colour = "11"),size=1) +
         # geom_line(aes(y = mean_expected_het_pop2[,151], colour = "9"),size=1) +
         # geom_line(aes(y = mean_expected_het_pop2[,208], colour = "3"),size=1) +
         # geom_point(aes(x = 1,y = 0.536,colour = "2"),size=3) +
         # geom_point(aes(x = 1,y = 0.338,colour = "3"),size=3)+
         # geom_point(aes(x = 1,y = 0.620,colour = "9"),size=3)+
         # geom_point(aes(x = 1,y = 0.582,colour = "11"),size=3)+
         geom_point(aes(x = 34,y = 0.28,colour = "1"),size=3) +
         # geom_point(aes(x = 17,y = 0.338,colour = "2"),size=3)+
         geom_point(aes(x = 34,y = 0.24,colour = "3"),size=3)+
         geom_point(aes(x = 34,y = 0.44,colour = "4"),size=3)+
         geom_point(aes(x = 34,y = 0.38,colour = "5"),size=3) +
         geom_point(aes(x = 34,y = 0.34,colour = "6"),size=3)+
         geom_point(aes(x = 34,y = 0.21,colour = "7"),size=3)+
         geom_point(aes(x = 34,y = 0.40,colour = "8"),size=3)+
         theme_bw(base_size = 14) +
         theme(legend.title=element_blank())+
         theme(legend.position="bottom")+
         xlab("GENERATIONS") +
         ylab("HE")  
  
  
  
  
  
  
  
  
  
       # +
       #   ggtitle("ASSOCIATIVE OVERDOMINANCE",
       #           subtitle=paste(
       #             "G",gen_number_pre_adaptation,
       #             "P",population_size_pre_adaptation,
       #             "L",loci_number,
       #             "O",number_offspring,
       #             "M",msats)))

print (
  
  # location_msats_experiment <-  c(83, 125, 150, 207)
  # location_msats_experiment <-  c(88, 133, 161, 219)
  
msats_freq <- list(c(0.04,0.09,0.06,0.16,0.15,0.5) c(0.11,0.41,0.23,0.13,0.12), c(0.3,0.05,0.41,0.24), c(0.03,0.74,0.15,0.08))
  ggplot(generations, aes(V1)) +
         geom_line(aes(y = mean_Fst[,87], colour = "2"),size=1) +
         geom_line(aes(y = mean_Fst[,132], colour = "11"),size=1) +
         geom_line(aes(y = mean_Fst[,155], colour = "9"),size=1) +
         geom_line(aes(y = mean_Fst[,218], colour = "3"),size=1) +
         # geom_point(aes(x = 17,y = 0.1,colour = "2"),size=4) +
         # geom_point(aes(x = 17,y = 0.102,colour = "3"),size=4)+
         # geom_point(aes(x = 17,y = 0.130,colour = "9"),size=4)+
         # geom_point(aes(x = 17,y = 0.170,colour = "11"),size=4)+
         
    ### FSTp
    geom_point(aes(x = 34,y = 0.180,colour = "2"),size=4) +
         geom_point(aes(x = 34,y = 0.219,colour = "11"),size=4)+
         geom_point(aes(x = 34,y = 0.158,colour = "9"),size=4)+
         geom_point(aes(x = 34,y = 0.088,colour = "3"),size=4)+
    ## AFD
    # geom_point(aes(x = 34,y = 0.321,colour = "2"),size=4) +
    #    geom_point(aes(x = 34,y = 0.386,colour = "11"),size=4)+
    #    geom_point(aes(x = 34,y = 0.290,colour = "9"),size=4)+
    #    geom_point(aes(x = 34,y = 0.187,colour = "3"),size=4)+
    #### Ht
    # geom_point(aes(x = 34,y = 0.469,colour = "2"),size=4) +
    #      geom_point(aes(x = 34,y = 0.550,colour = "11"),size=4)+
    #      geom_point(aes(x = 34,y = 0.456,colour = "9"),size=4)+
    #      geom_point(aes(x = 34,y = 0.385,colour = "3"),size=4)+
        #### mean het
    # geom_point(aes(x = 34,y = 0.423,colour = "2"),size=4) +
    #      geom_point(aes(x = 34,y = 0.472,colour = "11"),size=4)+
    #      geom_point(aes(x = 34,y = 0.414,colour = "9"),size=4)+
    #      geom_point(aes(x = 34,y = 0.361,colour = "3"),size=4)+
         theme_bw(base_size = 14) +
         theme(legend.title=element_blank())+
         theme(legend.position="bottom") +
         xlab("GENERATIONS") +
         ylab("FSTp")
  
  ggplot(generations, aes(V1)) +
         geom_line(aes(y = mean_Fstp[,24], colour = "1"),size=1) +
         # geom_line(aes(y = mean_Fst[,48], colour = "2"),size=1) +
         geom_line(aes(y = mean_Fstp[,50], colour = "3"),size=1) +
         geom_line(aes(y = mean_Fstp[,70], colour = "4"),size=1) +
         geom_line(aes(y = mean_Fstp[,82], colour = "5"),size=1) +
         geom_line(aes(y = mean_Fstp[,110], colour = "6"),size=1) +
         geom_line(aes(y = mean_Fstp[,132], colour = "7"),size=1) +
         geom_line(aes(y = mean_Fstp[,146], colour = "8"),size=1) +
         # geom_point(aes(x = 17,y = 0.1,colour = "2"),size=4) +
         # geom_point(aes(x = 17,y = 0.102,colour = "3"),size=4)+
         # geom_point(aes(x = 17,y = 0.130,colour = "9"),size=4)+
         # geom_point(aes(x = 17,y = 0.170,colour = "11"),size=4)+
         # geom_point(aes(x = 34,y = 0.160,colour = "1"),size=4) +
         # # geom_point(aes(x = 34,y = 0.011,colour = "2"),size=4)+
         # geom_point(aes(x = 34,y = 0.110,colour = "3"),size=4)+
         # geom_point(aes(x = 34,y = 0.100,colour = "4"),size=4)+
         # geom_point(aes(x = 34,y = 0.185,colour = "5"),size=4) +
         # geom_point(aes(x = 34,y = 0.168,colour = "6"),size=4)+
         # geom_point(aes(x = 34,y = 0.129,colour = "7"),size=4)+
         # geom_point(aes(x = 34,y = 0.133,colour = "8"),size=4)+
    ##### this is fst
      # geom_point(aes(x = 34,y = 0.099,colour = "1"),size=4) +
      #    # geom_point(aes(x = 34,y = 0.011,colour = "2"),size=4)+
      #    geom_point(aes(x = 34,y = 0.079,colour = "3"),size=4)+
      #    geom_point(aes(x = 34,y = 0.058,colour = "4"),size=4)+
      #    geom_point(aes(x = 34,y = 0.111,colour = "5"),size=4) +
      #    geom_point(aes(x = 34,y = 0.115,colour = "6"),size=4)+
      #    geom_point(aes(x = 34,y = 0.074,colour = "7"),size=4)+
      #    geom_point(aes(x = 34,y = 0.078,colour = "8"),size=4)+
  ##### this is fstp
    geom_point(aes(x = 34,y = 0.160,colour = "1"),size=4) +
         # geom_point(aes(x = 34,y = 0.011,colour = "2"),size=4)+
         geom_point(aes(x = 34,y = 0.110,colour = "3"),size=4)+
         geom_point(aes(x = 34,y = 0.100,colour = "4"),size=4)+
         geom_point(aes(x = 34,y = 0.185,colour = "5"),size=4) +
         geom_point(aes(x = 34,y = 0.168,colour = "6"),size=4)+
         geom_point(aes(x = 34,y = 0.129,colour = "7"),size=4)+
         geom_point(aes(x = 34,y = 0.133,colour = "8"),size=4)+
         theme_bw(base_size = 14) +
         theme(legend.title=element_blank())+
         theme(legend.position="bottom") +
         xlab("GENERATIONS") +
         ylab("FSTp")
  
  
  
  
  
       # +
         # ggtitle("ASSOCIATIVE OVERDOMINANCE",
         #         subtitle=paste(
         #           # "G",gen_number_pre_adaptation,
         #           # "P",population_size_pre_adaptation,
         #           "L",loci_number,
         #           "O",number_offspring,
         #           "R",recombination)))
}
p6 <- print(ggplot(generations, aes(V1)) +
          geom_line(aes(y = rowMeans(mean_het_pop2), colour = "Simulations He"),size=1) +
            geom_line(aes(y = rowMeans(mean_expected_het_pop2), colour = "Expected He"),size=1) +
          # geom_line(aes(y = mean_df_2, colour = "Expected He"),size=1) +
        theme_bw(base_size = 18) +
         scale_fill_hc("darkunica")+
        labs(x="GENERATIONS", y="He", title=NULL)+ 
      theme(legend.title=element_blank())+
        theme(legend.position =  "bottom") +
        theme(legend.text=element_text(size=14)), plot.margin=unit(c(-1,1,1,-1), "in") )
  
p5 <-  print(ggplot(generations, aes(V1)) +
        geom_line(aes(y = rowMeans(mean_Fst,na.rm = T), colour = "Fst"),size=1) +
        # geom_line(aes(y = rowMeans(mean_Fstp,na.rm = T), colour = "Fstp"),size=1) +
        # geom_line(aes(y = Fst_expected, colour = "Fst expected"),size=1) +
        # geom_line(aes(y = 0.099, colour = "Fst expected"),size=1) +
        # geom_line(aes(y = 0.304, colour = "Fst expected"),size=1) +
        geom_line(aes(y = Fst_expected, colour = "Fst expected"),size=1) +
        theme_bw(base_size = 18) +
        scale_fill_hc("darkunica")+
        labs(x="GENERATIONS", y="Fst", title="Neutral simulations", subtitle ="Ne=50,dispersal rate=0.01" )+
        theme(legend.title=element_blank())+
        theme(legend.position = "bottom") +
        theme(legend.text=element_text(size=14)), plot.margin=unit(c(-1,-1,1,-1), "in") )
#   
#   print(ggplot(generations, aes(V1)) +
#         geom_line(aes(y = rowMeans(mean_shua), colour = "Shua"),size=1) +
#         geom_line(aes(y = shua_expected, colour = "Shua expected"),size=1) +
#         theme_bw(base_size = 18) +
#         scale_fill_hc("darkunica")+
#         xlab("GENERATIONS") +
#         ylab("Shua")+   
#         theme(legend.title=element_blank())+
#         theme(legend.position = c(0.7,0.25)) +
#         theme(legend.text=element_text(size=14)) )
# 
#  
# 
print(ggplot(generations, aes(V1)) +
        geom_line(aes(y = mean_g_load_a_final, colour = "Additive"),size=1) +
        geom_line(aes(y = mean_g_load_m_final, colour = "Multiplicative"),size=1) +
        theme_bw(base_size = 18) +
        scale_fill_hc("darkunica")+
        xlab("GENERATIONS") +
        ylab("Genetic load")+
        theme(legend.title=element_blank())+
        theme(legend.position = c(0.7,0.25)) +
        theme(legend.text=element_text(size=14)) )
# 
print(ggplot(generations, aes(V1)) +
        geom_line(aes(y = mean_deleterious_eliminated),size=1) +
        theme_bw(base_size = 18) +
        scale_fill_hc("darkunica")+
        xlab("GENERATIONS") +
        ylab("% SNP's lost")+
        theme(legend.title=element_blank())+
        theme(legend.position = c(0.7,0.25)) +
        theme(legend.text=element_text(size=14)) )
# 
print(ggplot(generations, aes(V1)) +
        geom_line(aes(y = mean_mean_s_final),size=1) +
        theme_bw(base_size = 18) +
        scale_fill_hc("darkunica")+
        xlab("GENERATIONS") +
        ylab("mean s")+
        theme(legend.title=element_blank())+
        theme(legend.position = c(0.7,0.25)) +
        theme(legend.text=element_text(size=14)) )

print(ggplot(generations, aes(V1)) +
        geom_line(aes(y = mean_mean_h_final),size=1) +
        theme_bw(base_size = 18) +
        scale_fill_hc("darkunica")+
        xlab("GENERATIONS") +
        ylab("mean h")+
        theme(legend.title=element_blank())+
        theme(legend.position = c(0.7,0.25)) +
        theme(legend.text=element_text(size=14)) )

##### SELECTION TESTS PLOTS ####
if(tests_selection==TRUE){
  if(B_B_test==T){print(plot_bayescan(B_B_res_1$bayescan))}
  if(B_B_test==T){print(plot_bayescan(B_B_res_2$bayescan))}

  if(W_L_test==T){print(OutFLANKResultsPlotter(OFoutput=W_L_res_a$results))
    print(OutFLANKResultsPlotter(OFoutput=W_L_res_b))
    }
  if(B_N_test==T){
    print(make.alpha.beta.plot(betas.sim=FST.simDF_res, All.Pvals=FST_fdist_res_b, truepos= which(FST_fdist_res_b$qval.tail<0.1), FODR.my.plot=F,title_plot= "FST vs He" ))
    print(make.alpha.beta.plot(betas.sim=Shua.simDF_res, All.Pvals=Shua_fdist_res_b, truepos= which(Shua_fdist_res_b$qval.tail<0.1), FODR.my.plot=F,title_plot= "Shua vs Shannon" ))
    print(make.alpha.beta.plot(betas.sim=Sha_diff.simDF_res, All.Pvals=Sha_diff_fdist_res_b, truepos= which(Sha_diff_fdist_res_b$qval.tail<0.1), FODR.my.plot=F,title_plot= "Sha_diff vs 1Dalpha" ))
    print(make.alpha.beta.plot(betas.sim=Sha_pops.simDF_res, All.Pvals=Sha_pops_fdist_res_b, truepos= which(Sha_pops_fdist_res_b$qval.tail<0.1), FODR.my.plot=F,title_plot= "Sha_pops vs Shannon" ))
  }
}
