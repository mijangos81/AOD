library(adegenet)
library(hierfstat)
library(dplyr)
library(rlist)
library(readxl)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(ggstance)
library(graph4lg)
library(data.table)
library(SuppDists)
library(qvalue)
library(ggpointdensity)
library(ggnewscale)
source('functions.R')
source('hierfstat.R')
source('Fdist_functions.R')
source('Bayescan_functions.R')
source('Outflank_functions.R')
source("basic_stats_fix.R")
source('snps_functions.R')
source("genind_to_genepop_2.R")

  B_N_test <- T # Beaumont and Nichols' test (FDIST2)
  B_B_test <- T # Balding and Beaumont's test (BAYESCAN)
  W_L_test <- T # Whitlock and Lotterhos' test (OUTFLANK)

 bayescan.path <- "/usr/local/bin/bayescan" # this is the location of the binary file of the program BAYESCAN
  number_cores <- 7 # this is the number of cores that the program BAYESCAN uses
  # Bayescan settings. these setting have been tested to ensure that convergence of the
  # Bayescan's RJ-MCMC algorithm has been reached. For this we used the function gelman.diag
  # from the R package coda.
  n_iter <- 3000 # number of outputted iterations
  thinning <- 10 # thinning interval size
  n_pilot <- 10 # number of pilot runs
  l_pilot <- 5000 # length of pilot runs
  burnin <- 5000 # Burn-in length
  # the width of the bins (number of observations) used to estimate the density of
  # the Johnson distribution.
  bins_Johnson_distribution <- 400

# loc_snps <- read_excel("loc_snps_past_version.xlsx")
snps <- F
msats <-T

# which sampling time to analyse?
  T2 <-T
  T0 <- F
  
# which dispersal regime to analyse?
  dis <- "Low"   # options: "Hig","Med", "Low"
  
if(snps==T){# Import genepop file to genind object 
all_snps_auto<-read.genepop(file = "every_snp.gen")
}
if(msats==T){# Import genepop file to genind object 
all_snps_auto<-read.genepop(file = "all_msats_auto.gen")
}

# assign pop information to strata
strata(all_snps_auto) <- data.frame(pop = pop(all_snps_auto))
# split strata listing the name of each hierarchy level 
splitStrata(all_snps_auto)<-~time/l_pair/rep/dis/b_pair/pop1_2/p_pair
# Visualize the summary of all strata fields
# lapply(strata(all_snps_auto),summary)
if (T2==T) {
  # # separating T2
   test_to_keep<-all_snps_auto$strata$time==2
}
if (T0==T) {
# separating T0
  test_to_keep<-all_snps_auto$strata$time==0
}

t2_snps_auto<-all_snps_auto[test_to_keep]

test_to_keep_dispersal <- t2_snps_auto$strata$dis==dis
t2_snps_auto <- t2_snps_auto[test_to_keep_dispersal]

# separate populations based in the population pair
per_pop_snps<-seppop(t2_snps_auto,t2_snps_auto$strata$p_pair, drop=TRUE)

  ##### Beaumont and Nichol's method (B_N_method or FDIST2) #####
  if(B_N_test==TRUE){
files_B_N_neutral <- dir(paste0(getwd(),"/","selection_sim_experiment"),pattern = "fdist")
B_N_files_neutral <- read.table(paste0(getwd(),"/","selection_sim_experiment","/",files_B_N_neutral),header=T,sep=",")
# removing the cases when one allele is fixed
B_N_files_neutral$het[B_N_files_neutral$het==1] <- NA
B_N_files_neutral <- B_N_files_neutral[complete.cases(B_N_files_neutral),]
B_N_files_neutral <- B_N_files_neutral[sample(x=1:nrow(B_N_files_neutral) ,size = 4500),]

FST.sim_high_dispersal_two_demes <- read.table("high_dispersal_two_demes.sims", header=TRUE)
FST.sim_low_dispersal_two_demes <- read.table("low_dispersal_two_demes.sims", header=TRUE)
FST.sim_high_dispersal_ten_demes<- read.table("high_dispersal_ten_demes.sims", header=TRUE)
FST.sim_low_dispersal_ten_demes <- read.table("low_dispersal_ten_demes.sims", header=TRUE)
# plot(FST.sim_high_dispersal$He,FST.sim_high_dispersal$FST,col="blue", cex=1)
# points(FST.sim_low_dispersal$He,FST.sim_low_dispersal$FST,col="deeppink", cex=0.25)
 FST.simDF_res <- FST.sim_high_dispersal_two_demes[,c("He","FST")]
# FST.simDF_res <- B_N_files_neutral[,c("het","FST")]
colnames(FST.simDF_res) <- c("alpha","beta")
# ranking by heterozygosity as in Beaumont & Nichols 1996
FST.simDF_res <- FST.simDF_res[order(FST.simDF_res$alpha),]
FST.simDF_res <- FST.simDF_res[complete.cases(FST.simDF_res),]
per_pop_snps_genepop <-  lapply(per_pop_snps,genind_to_genepop_2)
per_pop_snps_genepop <- lapply(per_pop_snps_genepop,function(x){
  x[x=="0000"]<- NA
  return(x)
  }) 
FST.experiment_res_temp <- lapply(per_pop_snps_genepop,CW1993.dataset, diploid=TRUE, ndig=2)

FST.experiment_res_temp_fdist <- mean(unlist(unname(lapply(FST.experiment_res_temp,function(x){x<- x[[1]]}))))

FST.experiment_res <- as.data.frame(rbindlist(lapply(FST.experiment_res_temp,function(x){x<- x[[2]][,]})))
# FST.experiment_res$selection <- T
# FST.experiment_res[as.numeric(neutral_loci_location),"selection"] <- F
# colnames(FST.experiment_res) <- c("alpha","beta","selection")
colnames(FST.experiment_res) <- c("alpha","beta","numer","locnames")

FST.experiment_res <- FST.experiment_res[complete.cases(FST.experiment_res),]

FST_fdist_res <- get.pval.dataset(beta.empDF=FST.experiment_res, beta.simDF=FST.simDF_res, beta_max=1.01)
# FST_fdist_res <- get.pval.dataset(FST.empDF=FST.experiment_res, FST.simDF=FST.simDF_res, He.bin=0.04,write.progress=TRUE)

FST_fdist_res_b <- correct.pval.dataframe(dataframe=FST_fdist_res, p.colName="p.val.cum")
FST_fdist_res_balancing <- FST_fdist_res_b[FST_fdist_res_b$tail=="L" & FST_fdist_res_b$FDR.05==TRUE, ]
# selection_tests[i,1] <- nrow(FST_fdist_res_balancing)
FST_fdist_res_directional <- FST_fdist_res_b[FST_fdist_res_b$tail=="R" & FST_fdist_res_b$FDR.05==TRUE, ]
# selection_tests[i,9] <- nrow(FST_fdist_res_directional)
print(make.alpha.beta.plot(betas.sim=FST.simDF_res, All.Pvals=FST_fdist_res_b, truepos= which(FST_fdist_res_b$Bonf==T), FODR.my.plot=T,title_plot= "FST vs He" ))
print(make.alpha.beta.plot(betas.sim=FST.simDF_res, All.Pvals=FST_fdist_res_b, truepos= which(FST_fdist_res_b$selection==T), FODR.my.plot=F,title_plot= "FST vs He" ))
}

  ##### Whitlock and Lotterhos' method (W_L_method or OUTFLANK)#####
  if (W_L_test==TRUE){
files_W_L_neutral <- dir(paste0(getwd(),"/","selection_sim_experiment"),pattern = "outflank_")

 outflank_neutral <- read.table(paste0(getwd(),"/","selection_sim_experiment","/",files_W_L_neutral),header=T,sep=",")
    Fstbar_neut <- outflank_neutral[1,1]
    df_infer <- outflank_neutral[2,1]
    outflank_msats <- 4
    if(outflank_msats>2){
      He_min <- 0.2
    }else{
      He_min <- 0.1
    }

per_pop_snps_outflank <-  lapply(per_pop_snps,genind2hierfstat_fix)
per_pop_snps_outflank_res_temp <- lapply(per_pop_snps_outflank,wc_2)
# per_pop_snps_outflank_res_temp <- lapply(per_pop_snps_outflank,wc)

per_pop_snps_outflank_res <- as.data.frame(rbindlist(per_pop_snps_outflank_res_temp))
outflank_input <- per_pop_snps_outflank_res
outflank_input$FSTNoCorr[outflank_input$FSTNoCorr==1] <- NA
outflank_input <- outflank_input[complete.cases(outflank_input$FSTNoCorr),]
outflank_input_a <- outflank_input[complete.cases(outflank_input$FST),]
outflank_input_b <- cbind(outflank_input_a,qvalues=NA, OutlierFlag=NA)
W_L_res_a <- OutFLANK(outflank_input_a,LeftTrimFraction = 0.05, RightTrimFraction = 0.05,Hmin = He_min, NumberOfSamples=2, qthreshold = 0.05)
W_L_res_b <- pOutlierFinderChiSqNoCorr(DataList=outflank_input_b, Fstbar=Fstbar_neut, dfInferred=df_infer, qthreshold = 0.05,Hmin = He_min)

OutFLANKResultsPlotter(W_L_res_a, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = He_min, binwidth = 0.05, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

OutFLANKResultsPlotter(W_L_res_b, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = He_min, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

}
  ##### Balding and Beaumont's method (B_B_method or BAYESCAN)#####
  if(B_B_test==TRUE){
per_pop_snps_bayescan <- lapply(per_pop_snps,genind2hierfstat_fix)
per_pop_snps_bayescan <- lapply(per_pop_snps_bayescan,function(x){x<-x[,-1]})
per_pop_snps_bayescan <- list.cbind(per_pop_snps_bayescan)
pops <- c(rep(1,nrow(per_pop_snps_bayescan)/2),rep(2,nrow(per_pop_snps_bayescan)/2))
per_pop_snps_bayescan <-  cbind(pops,per_pop_snps_bayescan)

files_bayescan_neutral <- dir(paste0(getwd(),"/","selection_sim_experiment"),pattern = "bayescan")
bayescan_neutral <- read.table(paste0(getwd(),"/","selection_sim_experiment","/",files_bayescan_neutral),header=T,sep=",")
per_pop_snps_bayescan_neutral <-  per_pop_snps_bayescan
  # cbind(per_pop_snps_bayescan,bayescan_neutral)
write.bayescan(per_pop_snps_bayescan_neutral,file_name = "bayescan_experiment_selection.txt")
# files_B_B <- dir(path.folder_sim,pattern = "^bayescan_experiment")
# Files should be in the working directory for bayescan to work.
working_dir <- getwd()
setwd(path.folder_sim)

B_B_res_1 <- run_bayescan(data = "bayescan_experiment_selection.txt",n=n_iter,thin=thinning,nbp=n_pilot,pilot=l_pilot,burn=burnin,pr_odds=10,parallel.core=number_cores,bayescan.path = bayescan.path,time_report = T)
# B_B_res_1_snps <- B_B_res_1
### BayeScan outputs a file with q-values, but not p-values.
# ### this function returns corrected p values for a Bayescan outfile
# # Function taken from Whitlock and Lotterhos 2014
B_B_res_1_b <- ReturnCorrectedPVal.BS(bayescan_ouput=B_B_res_1$bayescan)
B_B_res_1_balancing <- B_B_res_1_b[B_B_res_1_b$SELECTION == "balancing",]
# selection_tests[i,5] <- sum(B_B_res_1_balancing$FDR.05)
B_B_res_1_directional <- B_B_res_1_b[B_B_res_1_b$SELECTION == "diversifying",]
# selection_tests[i,13] <- sum(B_B_res_1_directional$FDR.05)
#DELETING BAYESCAN FOLDER
unlink("radiator_bayescan_*",recursive = T)
# After all the files have been analysed the working directory is set back
# to the original directory
setwd(working_dir)

print(plot_bayescan(B_B_res_1$bayescan))

print(plot_bayescan(B_B_res_1$bayescan))
}


FST.experiment_hierf <- lapply(per_pop_snps,genind2hierfstat_fix)
FST.experiment_alleles <- lapply(FST.experiment_hierf,allele.count)
FST.experiment_shannon_temp <- lapply(FST.experiment_alleles, function(x) {sha <- sapply(as.vector(x),shannon_pops)})

FST.experiment_shannon_temp_2 <- lapply(FST.experiment_shannon_temp,unlist)
FST.experiment_shannon_temp_2 <- unlist(FST.experiment_shannon_temp_2)
FST.experiment_shannon_temp_2 <- unname(FST.experiment_shannon_temp_2)

sha_experiment_pop1_temp <- lapply(FST.experiment_alleles, function(x) {sha <- unlist(sapply(x, shannon)[1, ])})
sha_experiment_pop1_temp_2 <- lapply(sha_experiment_pop1_temp,unlist)
sha_experiment_pop1_temp_2 <- unlist(sha_experiment_pop1_temp_2)
sha_experiment_pop1_temp_2 <- unname(sha_experiment_pop1_temp_2)

sha_experiment_pop2_temp <- lapply(FST.experiment_alleles, function(x) {sha <- unlist(sapply(x, shannon)[2, ])})
sha_experiment_pop2_temp_2 <- lapply(sha_experiment_pop2_temp,unlist)
sha_experiment_pop2_temp_2 <- unlist(sha_experiment_pop2_temp_2)
sha_experiment_pop2_temp_2 <- unname(sha_experiment_pop2_temp_2)


FST.experiment_shua_temp <- lapply(FST.experiment_alleles,function(x){shua <- sapply(as.vector(x),mutual_information)}) 

FST.experiment_shua_temp_2 <- lapply(FST.experiment_shua_temp,unlist)
FST.experiment_shua_temp_2 <- unlist(FST.experiment_shua_temp_2)
FST.experiment_shua_temp_2 <- unname(FST.experiment_shua_temp_2)

Shua.experiment_res <-  as.data.frame(cbind(FST.experiment_shannon_temp_2,FST.experiment_shua_temp_2))
colnames(Shua.experiment_res) <- c("alpha","beta")

fdist_selection <-  read.table(file = "fdist_selection.csv" ,header=T,sep=",",skip=variables_number)
fdist_neutral <- read.table(file = "fdist_neutral.csv" ,header=T,sep=",",skip=variables_number)

expected_fst <- mean(FST.experiment_res$beta/FST.experiment_res$alpha)

number_transfers <- 2
transfer_each_gen <- 3
dispersal_rate <- (number_transfers / transfer_each_gen) / (population_size_dispersal)
Fst_expected <- 1 / ((4 * Ne_fst * dispersal_rate) * ((2 / (2 - 1)) ^ 2) + 1)
print(Fst_expected)

FST.sim_high_dispersal_two_demes 
FST.sim_low_dispersal_two_demes 
FST.sim_high_dispersal_ten_demes
FST.sim_low_dispersal_ten_demes 

betas.sim <- FST.sim_low_dispersal_two_demes
colnames(betas.sim) <- c("numer","alpha", "beta")
All.Pvals <- FST.experiment_res

maxalpha <- max(c(betas.sim$alpha))
  minalpha <- min(c(betas.sim$alpha))
  maxbeta <- max(c(betas.sim$beta))
  minbeta <- min(c(betas.sim$beta))
  
  alpha.cell <- 0.02
  beta.cell <- 0.02
  
  xseq<-seq(minalpha,maxalpha,alpha.cell)
  yseq<-seq(minbeta,maxbeta,beta.cell )
envelope_temp <-  as.data.frame(get.CI.alpha(xseq=xseq, percentile=0.975, betas.sim=betas.sim))
envelope <- as.data.frame(rbind(cbind(xseq,envelope_temp$V1),cbind(xseq,envelope_temp$V2)))
# envelope[is.na(envelope$V2),"V2"] <- 1

fdist_sim <- betas.sim
colnames(fdist_sim) <- c("numer","He","FST")

ggplot(fdist_sim)+
  geom_pointdensity(aes(x=He,y=FST),size=3)+
  scale_color_viridis_c(option = "D",name="Loci\ndensity") +
  # geom_vline(xintercept = seq(0,1,0.05),color="black",alpha=0.75) +
  geom_path(data = envelope,aes(x=xseq,y=V2),color="red",size=1.5)+   
    geom_point(data = FST.experiment_res,aes(x=alpha,y=beta),color="deeppink",size=1.5,alpha=0.75)+
  geom_point(data = fdist_selection,aes(x=het,y=FST),color="deepskyblue",size=1.5,alpha=0.75)+
    geom_hline(yintercept = mean(fdist_sim$FST),color="cyan",size=1)+
ylim(c(0,1))+
 theme( panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA))

ggplot()+
  # geom_pointdensity(aes(x=He,y=FST),size=3)+
  # scale_color_viridis_c(option = "D",name="Loci\ndensity") +
  geom_point(data = Shua.experiment_res,aes(x=alpha,y=beta),color="deeppink",size=2,alpha=0.5)+
    geom_point(data = fdist_selection,aes(x=sha_pops,y=sha_pop2),color="deepskyblue",size=2,alpha=0.5)+
      geom_point(data = fdist_neutral,aes(x=sha_pops,y=shua),color="green",size=2,alpha=0.5)+

  # geom_vline(xintercept = seq(0,1,0.05),color="black",alpha=0.75) +
  # geom_path(data = envelope,aes(x=xseq,y=V2),color="red",size=1.5)+
    # geom_hline(yintercept = mean(fdist_sim$FST),color="cyan",size=1)+
ylim(c(0,1))+
 theme( panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA))
  

 ggsave("shannon_distribution.pdf",  width = 4, height = 4, units = "in", dpi="retina", bg = "transparent"  )


# new_scale_color() +
# # geom_line(aes(y=mean(s_per_bp),x=CUB_dmel,color="Mean"),size=1)+
# geom_smooth(aes(y=FST,x=He,color="Trend line\n(LOESS)"),fill="cadetblue1",method="loess",span=0.75,size=2)
path.folder_sim <- "SIM_GRAL_B_50_C_10_D_0.2_F_0.005180543"
     files_B_N_neutral <- dir(path.folder_sim,pattern = paste0("^fdist_"))
    B_N_files_neutral <- lapply(paste0(path.folder_sim,"/",files_B_N_neutral),read.table,header=T,sep=",",skip=nrow(variables_file))
    B_N_files_neutral_b <- do.call(rbind, B_N_files_neutral)
    
    path.folder_sim <- "SIM_GRAL_B_14_C_1_D_0.2_F_0.005152105"
     files_outflank_neutral_temp <- dir(path.folder_sim,pattern = paste0("^outflank"))
    files_outflank_neutral <- lapply(paste0(path.folder_sim,"/",files_outflank_neutral_temp),read.table,header=T,sep=",",skip=nrow(variables_file))
    files_outflank_neutral_res <- do.call(rbind, files_outflank_neutral)
    
      
    neutral_fdist <- read.table("fdist_neutral.sims", header=TRUE)
    
    fst_mine <- B_N_files_neutral_b[complete.cases(B_N_files_neutral_b$FST),]
    fst_mine <- fst_mine[,c("het", "FST")] 
    fst_mine$sim <- "mine"
    fst_fdist <- neutral_fdist[sample(1:2500,nrow(fst_mine)),]
    fst_fdist <- fst_fdist[,c("He","FST")]
    fst_fdist$sim <- "fdist"
    colnames(fst_fdist) <- colnames(fst_mine)
    
    plot_density_fdist <- as.data.frame(rbind(fst_mine,fst_fdist))
    
  ggplot(plot_density_fdist, aes(x=FST, color=sim)) +
  geom_density()
  
  fst_sim <- unlist(unname(lapply(Fstp,function(x){x[[100]]})))
  he_sim <- unlist(unname(lapply(Ht,function(x){x[[100]]})))
  
  plot(he_sim,fst_sim)
  
  
