  list_a <- list(
    AFD,
    sha_diff,
    one_D_alpha_pop1,
    one_D_alpha_pop2,
    two_D_alpha_pop1,
    two_D_alpha_pop2,
    one_D_beta,
    two_D_beta,
    sha_pop1,
    sha_pop2,
    sha_pops,
    shua,
    expected_het_pop1,
    expected_het_pop2,
    het_pop1,
    het_pop2,
    Ht,
    Fst,
    Fstp,
    Dest
     )
     
    list_b <- list(
    deleterious_eliminated,
    mean_s_final,
    mean_h_final,
    g_load_a_final,
    g_load_m_final
    )
    
    names_list <- c(
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
    "g_load_m_final"
    )

##### FDIST AND OUTFLANK NEUTRAL SIMULATIONS ######
# neutral simulations for the selection tests are produced for every combination 
# of variables that change population history
if (neutral_simulations==T){
  if(simulation_type=="fly"){
    variables_to_write <- variables_file_b[c(1:13,31)]
    # variables_to_write_b is used to save the file that will be read by the selection tests
    variables_to_write_b <- paste(variables_to_write,collapse = "_")
    # variables_to_write_c is used to save the files that will be read by the neutral_simulations.R script
    variables_to_write_c <- paste0(variables_to_write_b, "_rep_", iteration,".csv")
  }
  if(simulation_type=="general"){
    variables_to_write <- variables_file_b[c(1,3:15)]
    # variables_to_write_b is used to save the file that will be read by the selection tests
    variables_to_write_b <- paste(variables_to_write,collapse = "_")
    # variables_to_write_c is used to save the files that will be read by the neutral_simulations.R script
    variables_to_write_c <- paste0(variables_to_write_b, "_rep_", iteration,".csv")
  }
  # if(testing_sel_tests==F){
    if(B_N_test==T){
    for (value in 1:length(loci_location)) {
        alleles_B_N[,value] <- c(
        paste0("0",alleles_res[[value]][, 1],"0",alleles_res[[value]][, 2]),
        paste0("0",alleles_res[[value]][, 3],"0",alleles_res[[value]][, 4])
        )
    }
    pops_B_N <- cbind(hierf[[1]][,1],alleles_B_N)
    FST.empDF_res <- CW1993.dataset(data1=pops_B_N, diploid=TRUE, ndig=2)[[2]][,1:2]
    input_fdist[,1] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(sha_pop1[[x]])[gen_dispersal]))
    input_fdist[,2] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(sha_pop2[[x]])[gen_dispersal]))
    input_fdist[,3] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(one_D_alpha_pop1[[x]])[gen_dispersal]))
    input_fdist[,4] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(one_D_alpha_pop2[[x]])[gen_dispersal]))
    input_fdist[,5] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(shua[[x]])[gen_dispersal]))
    input_fdist[,6] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(sha_diff[[x]])[gen_dispersal]))
    input_fdist[,7] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(sha_pops[[x]])[gen_dispersal]))
    input_fdist[,8:9] <- FST.empDF_res
    colnames(input_fdist) <- c("sha_pop1","sha_pop2","one_D_alpha_pop1","one_D_alpha_pop2","shua","sha_diff","sha_pops","het","FST")
    write.table(variables_file,file = paste0(path.folder_sim,"/","fdist_neutral_",variables_to_write_c),sep = ",",row.names = F,col.names = F)
    suppressWarnings(write.table(input_fdist,file = paste0(path.folder_sim,"/","fdist_neutral_",variables_to_write_c),sep = ",",row.names = F,append = T))
    }
    if(W_L_test==T){
        for (value in 1:length(loci_location)) {
            hierf_W_L[,value+1] <- as.numeric(c(
                paste0(alleles_res[[value]][, 1], alleles_res[[value]][, 2]),
                paste0(alleles_res[[value]][, 3], alleles_res[[value]][, 4])
            ))
            }
            input_outflank <- wc_2(hierf_W_L)
            write.table(variables_file,file = paste0(path.folder_sim,"/","outflank_msats_",variables_to_write_c),sep = ",",row.names = F,col.names = F)
            suppressWarnings(write.table(input_outflank,file = paste0(path.folder_sim,"/","outflank_msats_",variables_to_write_c),sep = ",",row.names = F,append = T))
     }
  # }
    
 # if(testing_sel_tests==T){      
ind_pop1 <- nrow(pop1)
ind_pop2 <- nrow(pop2)
df_sel_tests <- rbind(pop1,pop2)

  df_sel_tests$V1[df_sel_tests$V1=="Male"]   <- 1
  df_sel_tests$V1[df_sel_tests$V1=="Female"] <- 2
  df_sel_tests[1:ind_pop1,2] <- "pop1"
  df_sel_tests[(ind_pop1+1):nrow(df_sel_tests),2] <- "pop2"
  df_sel_tests$id <- paste0(df_sel_tests[,2],"_",1:(ind_pop1+ind_pop2))
  plink_ped <- apply(df_sel_tests,1,ped)
  haploview <- gsub("a", "1", plink_ped) # converting allele names to numbers
  haploview <- gsub("A", "2", haploview)
  write.table(haploview,file = paste0(path.folder_sim,"/","haploview.ped"),quote = F,row.names = F,col.names = F)
  snp_stats  <- read.pedfile_b(paste0(path.folder_sim,"/","haploview.ped"),sep = " ",snps = plink_map$V4)
  genotype_pops <- snp_stats$genotypes
  plink_map_2 <- plink_map
  colnames(plink_map_2) <- c("chr_name","row_name","loc_cM","loc_bp")
  #This is mutual information foe each SNP
  MI_snps <- genotype_pops@.Data
  MI_snps <- as.data.frame(matrix(as.double(MI_snps),nrow = nrow(MI_snps),ncol = ncol(MI_snps) ))
  MI_snps[MI_snps==1] <- 11
  MI_snps[MI_snps==2] <- 12
  MI_snps[MI_snps==3] <- 22
  pop_names <- as.data.frame(c(rep(1, population_size_dispersal),rep(2, population_size_dispersal)))
  MI_snps_2 <- as.data.frame(cbind(pop_names,MI_snps))
  
# input_fdist <- CW1993.dataset(data1=MI_snps_2, diploid=TRUE, ndig=1)[[2]][,1:2]
# colnames(input_fdist) <- c("het","FST")
# write.table(variables_file,file = paste0(path.folder_sim,"/","fdist_neutral_",variables_to_write_c),sep = ",",row.names = F,col.names = F)
# suppressWarnings(write.table(input_fdist,file = paste0(path.folder_sim,"/","fdist_neutral_",variables_to_write_c),sep = ",",row.names = F,append = T))

input_outflank <- wc_2(MI_snps_2)
write.table(variables_file,file = paste0(path.folder_sim,"/","outflank_snps_",variables_to_write_c),sep = ",",row.names = F,col.names = F)
suppressWarnings(write.table(input_outflank,file = paste0(path.folder_sim,"/","outflank_snps_",variables_to_write_c),sep = ",",row.names = F,append = T))
    # }
    
}
    
    
    # else{ 
    ##### NORMAL SIMULATIONS #####
      # list_a <- lapply(list_a, function(x)x[-1] )
    by_gen_list_a <- lapply(list_a, function(x)suppressMessages(bind_cols(x)))
    by_gen_lists_a_b <- c(by_gen_list_a,list_b)
    by_gen_lists_a_b <- lapply(by_gen_lists_a_b,function(x) data.frame(matrix(unlist(x), ncol=length(x), byrow=F)))
    suppressWarnings(write.table(variables_file,file = paste0(path.folder_sim,"/","generations",variables,"_rep_", iteration,".csv"),sep = ",",row.names = F,col.names = F))
    lapply(by_gen_lists_a_b,write.table,file=paste0(path.folder_sim,"/","generations",variables,"_rep_", iteration,".csv"),append=T, quote =F, sep= ",",row.names=F,col.names=F)
    ##### SELECTION TESTS #####
    if (neutral_simulations==F){
    if (tests_selection == T){
        ##### FDIST INPUT FILE  #####
        if (B_N_test==T){
            for (value in 1:length(loci_location)) {
                alleles_B_N[,value] <- c(
                paste0("0",alleles_res[[value]][, 1],"0",alleles_res[[value]][, 2]),
                paste0("0",alleles_res[[value]][, 3],"0",alleles_res[[value]][, 4])
                )
            }
            pops_B_N <- cbind(hierf[[1]][,1],alleles_B_N)
            FST.empDF_res <- CW1993.dataset(data1=pops_B_N, diploid=TRUE, ndig=2)[[2]][,1:2]
            input_fdist[,1] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(sha_pop1[[x]])[gen_dispersal]))
            input_fdist[,2] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(sha_pop2[[x]])[gen_dispersal]))
            input_fdist[,3] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(one_D_alpha_pop1[[x]])[gen_dispersal]))
            input_fdist[,4] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(one_D_alpha_pop2[[x]])[gen_dispersal]))
            input_fdist[,5] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(shua[[x]])[gen_dispersal]))
            input_fdist[,6] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(sha_diff[[x]])[gen_dispersal]))
            input_fdist[,7] <- unlist(lapply(1:length(neutral_loci_location),function(x) unlist(sha_pops[[x]])[gen_dispersal]))
            input_fdist[,8:9] <- FST.empDF_res
            colnames(input_fdist) <- c("sha_pop1","sha_pop2","one_D_alpha_pop1","one_D_alpha_pop2","shua","sha_diff","sha_pops","het","FST")
            write.table(variables_file,file = paste0(path.folder_sim,"/","fdist_selection",variables,"_rep_", iteration,".csv"),sep = ",",row.names = F,col.names = F)
            suppressWarnings(write.table(input_fdist,file = paste0(path.folder_sim,"/","fdist_selection",variables,"_rep_", iteration,".csv"),sep = ",",row.names = F,append = T))
        }
        ##### BAYESCAN INPUT FILE  #####
        if (B_B_test==T){
            for (value in 1:length(loci_location)){
                alleles_B_B[,value] <- as.numeric(c(
                paste0(alleles_res[[value]][, 1], alleles_res[[value]][, 2]),
                paste0(alleles_res[[value]][, 3], alleles_res[[value]][, 4])
                ))
            } 
            input_bayescan <- alleles_B_B
            # input_bayescan <- cbind(hierf[[1]][,1],alleles_B_B)
            write.csv(input_bayescan,file = paste0(path.folder_sim,"/","bayescan",variables,"_rep_", iteration,".csv"))
            # write.bayescan(input_bayescan,file_name = paste0(path.folder_sim,"/","bayescan",variables,"_rep_", iteration,".txt"))
        }
        ##### OUTFLANK INPUT FILE  #####
        if(W_L_test==T){
      for (value in 1:length(loci_location)) {
            hierf_W_L[,value+1] <- as.numeric(c(
                paste0(alleles_res[[value]][, 1], alleles_res[[value]][, 2]),
                paste0(alleles_res[[value]][, 3], alleles_res[[value]][, 4])
            ))
            }
            input_outflank <- wc_2(hierf_W_L)
            write.table(variables_file,file = paste0(path.folder_sim,"/","outflank",variables,"_rep_", iteration,".csv"),sep = ",",row.names = F,col.names = F)
            suppressWarnings(write.table(input_outflank,file = paste0(path.folder_sim,"/","outflank",variables,"_rep_", iteration,".csv"),sep = ",",row.names = F,append = T))
        }
    }
    }
    ##### LD ANALYSES ######
    if(LD_analyses==T){
    # LD analyses are done with offspring and not with parents and using pop1. 
    # over 100 replications, LD patterns in pop1 and pop2 are the same.
        write.table(variables_file,file = paste0(path.folder_sim,"/","snps_pop1",variables,"_rep_", iteration,".csv"),sep = ",",row.names = F,col.names = F)
        suppressWarnings(write.table(snpsum.col_pop1,file = paste0(path.folder_sim,"/","snps_pop1",variables,"_rep_", iteration,".csv"),sep = ",",row.names = F,append = T))
        write.table(variables_file,file = paste0(path.folder_sim,"/","ld_pop1",variables,"_rep_", iteration,".csv"),sep = ",",row.names = F,col.names = F)
        suppressWarnings(write.table(ld_columns_pop1,file = paste0(path.folder_sim,"/","ld_pop1",variables,"_rep_", iteration,".csv"),sep = ",",row.names = T,col.names = F,append = T))
        
        #  write.table(variables_file,file = paste0(path.folder_sim,"/","snps_pop2",variables,"_rep_", iteration,".csv"),sep = ",",row.names = F,col.names = F)
        # suppressWarnings(write.table(snpsum.col_pop2,file = paste0(path.folder_sim,"/","snps_pop2",variables,"_rep_", iteration,".csv"),sep = ",",row.names = F,append = T))
        # write.table(variables_file,file = paste0(path.folder_sim,"/","ld_pop2",variables,"_rep_", iteration,".csv"),sep = ",",row.names = F,col.names = F)
        # suppressWarnings(write.table(ld_columns_pop2,file = paste0(path.folder_sim,"/","ld_pop2",variables,"_rep_", iteration,".csv"),sep = ",",row.names = F,append = T))
    }
# } 
