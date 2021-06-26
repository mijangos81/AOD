library(ggplot2)
library(ggthemes)

hyp_sel_tests <- read.csv("selection_tests_res_chapter.csv")

directional <- hyp_sel_tests[,c("Fstp_BIAS","FDIST2_FST_directional","BAYESCAN_B_directional","OUTFLANK_A_directional")]
colnames(directional) <- c("Fst_bias", "FDist2","Bayescan","Outflank")

balancing <- hyp_sel_tests[,c("Fstp_BIAS","FDIST2_FST_balancing","BAYESCAN_B_balancing","OUTFLANK_A_balancing")]
colnames(balancing) <- c("Fst_bias", "FDist2","Bayescan","Outflank")


print(ggplot(directional) +
        geom_col(aes(y=FDist2 ,x = Fst_bias,fill = "FDist2"),width = 1/200) +
        geom_col(aes(y=Bayescan ,x = Fst_bias, fill = "Bayescan"),width = 1/200) +
        geom_col(aes(y=Outflank ,x = Fst_bias, fill = "Outflank"),width = 1/200) +
        theme_bw(base_size = 18) +
        scale_fill_hc("darkunica")+
        labs(y="Number of loci under directional selection", x="FST bias", title=NULL)+ 
        theme(legend.title=element_blank())+
        theme(legend.position =  "bottom") +
        theme(legend.text=element_text(size=14)))

print(ggplot(balancing) +
        geom_col(aes(y=FDist2 ,x = Fst_bias,fill = "FDist2"),width = 1/200) +
        geom_col(aes(y=Bayescan ,x = Fst_bias, fill = "Bayescan"),width = 1/200) +
        geom_col(aes(y=Outflank ,x = Fst_bias, fill = "Outflank"),width = 1/200) +
        theme_bw(base_size = 18) +
        scale_fill_hc("darkunica")+
        labs(y="Number of loci under balancing selection", x="FST bias", title=NULL)+ 
        theme(legend.title=element_blank())+
        theme(legend.position =  "bottom") +
        theme(legend.text=element_text(size=14)))
