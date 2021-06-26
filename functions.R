ped <- function(df_ped){
  if(msats==2){ #if the neutral locus is biallelic, it is not deleted from the ped file 
    chromosome1 <- df_ped[3]
    chromosome2 <- df_ped[4]
  }else{ #if the neutral locus has more than 2 alleles it is deleted from the ped file 
    chromosome1 <- gsub("[-^1-9]", "0", df_ped[3])
    chromosome2 <- gsub("[-^1-9]", "0", df_ped[4])
  }
  split_seqs <- strsplit(c(chromosome1, chromosome2), split = "")
  genotypes <- as.data.frame(matrix(nrow =loci_number,ncol =2 ))
  genotypes$V1 <-split_seqs[[1]]
  genotypes$V2 <-split_seqs[[2]]
  genotypes_final <- paste0(paste(genotypes[,1],genotypes[,2]),collapse = " ")
  return(paste(df_ped[2],paste0("O",df_ped["id"],collapse = ""),df_ped[5],df_ped[6],df_ped[1],1,genotypes_final))
}

delta <- function(a, b, c) {
  # # Constructing delta
  b ^ 2 - 4 * a * c
}

q_equilibrium <- function(a, b, c) {
  # Constructing Quadratic Formula
  x_1 = (-b + sqrt(delta(a, b, c))) / (2 * a)
  return(x_1)
}

migration <- function(population1, population2, generation) {
  if (generation != 1 & generation %% transfer_each_gen == 0) {
    if(number_transfers == 1){
         if (maletran) {
      malepoptran <- sample(c(1:(population_size / 2)), size = 1)
      temppop1 <- population1[malepoptran,]
      temppop2 <- population2[malepoptran,]
      population2[malepoptran,] <- temppop1
      population1[malepoptran,] <- temppop2
    }
    if (femaletran) {
      fempoptran <-
        sample(c(((
          population_size / 2
        ) + 1):population_size), size = 1)
      temppop1 <- population1[fempoptran,]
      temppop2 <- population2[fempoptran,]
      population2[fempoptran,] <- temppop1
      population1[fempoptran,] <- temppop2
    }
    }
    if (number_transfers >= 2) {
         size_malepoptran <- ceiling(number_transfers/2)
    size_femalepoptran <- floor(number_transfers/2)
    if (maletran) {
      malepoptran <- sample(c(1:(population_size / 2)), size = size_malepoptran)
      temppop1 <- population1[malepoptran,]
      temppop2 <- population2[malepoptran,]
      population2[malepoptran,] <- temppop1
      population1[malepoptran,] <- temppop2
    }
    if (femaletran) {
      fempoptran <- sample(c(((population_size / 2) + 1):population_size), size = size_femalepoptran)
      temppop1 <- population1[fempoptran,]
      temppop2 <- population2[fempoptran,]
      population2[fempoptran,] <- temppop1
      population1[fempoptran,] <- temppop2
    }
    }
    # if necessary flip transfer
    if (number_transfers == 1) {
      maletran <- !maletran
      femaletran <- !femaletran
    }
  }
  return (list(population1, population2, maletran, femaletran))
}

selection_fun <- function(offspring,reference_pop) {
  offspring$fitness <- apply(offspring, 1, fitness, ref = reference_pop)
  if (natural_selection_model=="absolute"){
  offspring$random_deviate <- runif(nrow(offspring),min = 0, max=genetic_load)
  offspring$alive <- offspring$fitness > offspring$random_deviate 
  offspring <- offspring[which(offspring$alive == TRUE), ]
   }
   if (natural_selection_model=="relative"){
      fitnes_proportion <- sum(offspring$fitness)
      offspring$relative_fitness <- offspring$fitness / fitnes_proportion
      offspring$relative_fitness[offspring$relative_fitness<0] <- 0 
     }
  return(offspring)
}

# this is the function to calculate fitness
fitness <- function(df_fitness, ref) {
  #remove neutral loci from chromosomes
  chromosome1 <- gsub("[-^1-9]", "0", df_fitness[3])
  chromosome2 <- gsub("[-^1-9]", "0", df_fitness[4])
  split_seqs <- strsplit(c(chromosome1, chromosome2), split = "")
  ref$hom <- (split_seqs[[1]] == split_seqs[[2]] & split_seqs[[1]] == "a")
  ref$het_chr1 <- (split_seqs[[1]] == "a" & split_seqs[[2]] == "A")
  ref$het_chr2 <- (split_seqs[[1]] == "A" & split_seqs[[2]] == "a")
  ref$hom <- as.numeric(ref$hom)
  sum_hom <- sum(as.numeric(ref$hom))
  ref$het_chr1 <- as.numeric(ref$het_chr1)
  ref$het_chr2 <- as.numeric(ref$het_chr2)
  ref$het_sel_chr1 <- ref$h * ref$s * ref$het_chr1
  ref$het_sel_chr2 <- ref$h * ref$s * ref$het_chr2
  ref$hom_sel <- ref$s * ref$hom 
    ref$tot_sel <- ref$het_sel_chr1 + ref$het_sel_chr2 + ref$hom_sel
   ref$fitness <- 1 - (ref$tot_sel) 
  # keeping only loci with NS
  deleterious <- ref[ref$fitness<1,]
  if (fitness_for_repulsion == TRUE){
    # resetting row.names
    row.names(deleterious) <- 1:nrow(deleterious)
    # splitting the chromososme in genes of size gene_size, and getting the indices of the NS within genes
    genes <- split(rownames(deleterious),cut(deleterious$location,breaks = seq(1,tail(ref$location,1),gene_size)))
    # getting the genes with only one NS 
    one_NS <- which(lapply(genes,length)==1)
    # getting the genes with more than one NS
    more_than_one_NS <- which(lapply(genes,length)>1)
    # getting the indices of the NS in genes with just one NS
    alone_NS <- as.numeric(unlist(genes[one_NS]))
    # fitnesses of alone_NS do not change if they are alone
    fitness_alone_NS <- deleterious[alone_NS,"fitness"]
    if(length(more_than_one_NS)>1){
    # genes with more than one NS
    accompanied_NS <- genes[more_than_one_NS]
    accompanied_NS <- lapply(accompanied_NS,as.numeric)
    # getting the NS of the genes with more than one NS
    accompanied_NS_df <- lapply(accompanied_NS, function(x){deleterious[x,]})
    df_hom <- accompanied_NS_df[sapply(accompanied_NS_df, 
                                       function(x) any(sum(x$hom) > 0))]  
    df_repulsion <- accompanied_NS_df[sapply(accompanied_NS_df, 
                                             function(x) any(sum(x$het_chr1) > 0 
                                                             & sum(x$het_chr2) > 0 
                                                             & sum(x$hom) == 0))] 
    df_no_repulsion_chr1 <- accompanied_NS_df[sapply(accompanied_NS_df, 
                                                     function(x) any(sum(x$het_chr1) > 0
                                                                    & sum(x$het_chr2) == 0 
                                                                    & sum(x$hom) == 0))] 
    df_no_repulsion_chr2 <- accompanied_NS_df[sapply(accompanied_NS_df, 
                                                     function(x) any(sum(x$het_chr2) > 0
                                                                    & sum(x$het_chr1) == 0 
                                                                    & sum(x$hom) == 0))] 
    hom_fitness <- 1-unlist(lapply(df_hom,"[",,"s"))
    no_repulsion_chr1_fitness <- 1-unlist(lapply(df_no_repulsion_chr1,"[",,"tot_sel"))
    no_repulsion_chr2_fitness <- 1-unlist(lapply(df_no_repulsion_chr2,"[",,"tot_sel"))
    repulsion_fitness <- 1-unlist(lapply(df_repulsion,"[",,"s"))
    net_fitness <- prod(c(fitness_alone_NS,hom_fitness,no_repulsion_chr1_fitness,no_repulsion_chr2_fitness,repulsion_fitness))
    }else{
      net_fitness <- prod(deleterious$fitness)
    }
    }else{
      net_fitness <- prod(deleterious$fitness)
  }
  return(net_fitness)
}

shannon <- function(allele_matrix) {
  column_a <- allele_matrix[, 1]
  column_a <- column_a[column_a > 0]
  U_a <- sum(column_a)
  pi_a <- column_a / U_a
  est_a <- -sum(pi_a * log(pi_a))
  column_b <- allele_matrix[, 2]
  column_b <- column_b[column_b > 0]
  U_b <- sum(column_b)
  pi_b <- column_b / U_b
  est_b <- -sum(pi_b * log(pi_b))
  return(list(est_a, est_b))
}

shannon_pops <- function(alleles) {
  vector_alleles <- alleles
  vector_alleles <- vector_alleles[vector_alleles > 0]
  U <- sum(vector_alleles)
  pi <- vector_alleles / U
  est <- -sum(pi * log(pi))
  return(est)
}

mutual_information <- function(mat) {
  EstMLEFun <- function(mat) {
    # MLE
    entropyFun <- function(p) {
      p <- p[p > 0]
      out <- -sum(p * log(p))
      return(out)
    }
    n <- sum(mat)
    prob.hat <- mat / n
    px.hat <- apply(prob.hat, 1, sum)
    py.hat <- apply(prob.hat, 2, sum)
    I.hat <- entropyFun(px.hat) + entropyFun(py.hat) - entropyFun(prob.hat)
    # MLE of Mutual Information!
    return(I.hat)
  }
  mydata <- as.matrix(mat)
  est <- EstMLEFun(mydata)
  return(est)
}

# this is to calculate the frequency of deleterious alleles
split_seq <- function(df_split_seq) {
  if (msats == 2) {
    #if the neutral locus is biallelic, it is not deleted from the ped file
    chromosome1 <- df_split_seq[3]
    chromosome2 <- df_split_seq[4]
  }else{
    #if the neutral locus has more than 2 alleles it is set as missing in the ped file
    chromosome1 <- gsub("[-^1-9]", "0", df_split_seq[3])
    chromosome2 <- gsub("[-^1-9]", "0", df_split_seq[4])
  }
  split_seqs <- strsplit(c(chromosome1, chromosome2), split = "")
  genotypes <- as.data.frame(matrix(nrow = loci_number, ncol = 2))
  genotypes$V1 <- split_seqs[[1]]
  genotypes$V2 <- split_seqs[[2]]
  return(genotypes)
}

# this code returns the same LD output as the software haploview
# this is the input file for the LD analyses
LD_fun <- function(df_LD_fun ,show_warnings=T){
  df_LD_fun$V1[df_LD_fun$V1=="Male"]   <- 1
  df_LD_fun$V1[df_LD_fun$V1=="Female"] <- 2
  df_LD_fun[1:ind_pop1,2] <- "pop1"
  df_LD_fun[(ind_pop1+1):nrow(df_LD_fun),2] <- "pop2"
  df_LD_fun$id <- paste0(df_LD_fun[,2],"_",1:(ind_pop1+ind_pop2))
  plink_ped <- apply(df_LD_fun,1,ped)
  haploview <- gsub("a", "1", plink_ped) # converting allele names to numbers
  haploview <- gsub("A", "2", haploview)
  write.table(haploview,file = paste0(path.folder_sim,"/","haploview.ped"),quote = F,row.names = F,col.names = F)
  snp_stats  <- read.pedfile_b(paste0(path.folder_sim,"/","haploview.ped"),sep = " ",snps = plink_map$V4)
  genotype_pops <- snp_stats$genotypes
  genotype_pop1 <- genotype_pops[1:ind_pop1,]
  genotype_pop2 <- genotype_pops[(ind_pop1+1):nrow(genotype_pops),]
  snpsum.col_pop1 <- col.summary(genotype_pop1)
  colnames(snpsum.col_pop1) <-  paste0(colnames(snpsum.col_pop1),"_pop1")
  snpsum.col_pop2 <- col.summary(genotype_pop2)
  colnames(snpsum.col_pop2) <-  paste0(colnames(snpsum.col_pop2),"_pop2")
  snpFst <- Fst(genotype_pops,group=substr(rownames(genotype_pops),2,5))[1]
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
  allele_matrix_2 <- allele.count(MI_snps_2)
  MI <- as.numeric(paste(lapply(allele_matrix_2, mutual_information))) 
  AFD_snps <- as.numeric(paste(lapply(allele_matrix_2, AFD_fun))) 
 #This is Shannon index for each SNP
  sha_snps_pops <- lapply(allele_matrix_2,shannon)
  shannon_pop1 <- unname(unlist(lapply(sha_snps_pops,"[",1)))
  shannon_pop2 <- unname(unlist(lapply(sha_snps_pops,"[",2)))
  snp_final <- cbind(plink_map_2,snpsum.col_pop1[,4:9],snpsum.col_pop2[,4:9],snpFst,reference[,c("q","h","s")],MI,shannon_pop1,shannon_pop2,AFD_snps)
  # minor allele frequency (MAF) threshold
   minor <- 0.05 
   # Filter on MAF
  use <- with(snp_final, (!is.na(z.HWE_pop1) & MAF_pop1 >= minor)) 
  if(sum(use)<=9){return(list(NA,NA))}
  if(sum(use)>=10){
  # Remove NA's
  use[is.na(use)] <- FALSE 
  snp_final <- snp_final[use,]
  # Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
  genotype <- genotype_pop1[,use]
  # genotype_sample <- sample(1: ncol(genotype),size = ncol(genotype)/2)
  # genotype <- genotype[,genotype_sample]
  snp_loc  <- plink_map[use,] # this is the location of each snp in cM and in bp
  # snp_loc <- snp_loc[genotype_sample,]
  snp_loc$V2 <- as.numeric(snp_loc$V2)
  #this is the mean distance between each snp which is used to determine the depth at which 
  # LD analyses are performed 
  mean_dis <- mean(diff(snp_loc$V4))
  # ld_depth_b <- ceiling((ld_max_pairwise/mean_dis)*2)
  ld_depth_b <- ceiling((ld_max_pairwise/mean_dis))
  ld_snps <- ld(genotype,depth = ld_depth_b,stats = "R.squared") #function to calculate LD
  ld_matrix <- as.matrix(ld_snps)
  return(list(ld_matrix,snp_final))
  }
}

do_LD_analysis <- function(path_folder,list_ld,list_snps,list_files){
for( i in 1:nrow(list_files)){
    ld_columns <- read.table(paste0(path_folder,"/",list_ld[list_files[i,3]]),sep = ",",skip = variables_number,row.names = 1)
    colnames(ld_columns) <- rownames(ld_columns)
    snpsum.col <- read.table(paste0(path_folder,"/",list_snps[list_files[i,1]]),header=T,sep = ",",skip = variables_number)
    # if(is.na(snpsum.col)){next()}
    ld_columns <- as.data.frame(as.table(as.matrix(ld_columns)))
    ld_columns <- ld_columns[-ld_columns$Freq < 0,] #remove cases where LD was not calculated
    ld_columns$Var1 <- as.numeric(as.character(ld_columns$Var1))
    ld_columns$Var2 <- as.numeric(as.character(ld_columns$Var2))
    #determine the distance at which LD was calculated
    ld_columns$dis <- ld_columns$Var2 - ld_columns$Var1 
    #remove pairwise LD results that were calculated at larger distances than the required in the 
    # settings and then filtering and rearraging dataframes to match each other and then merge them
  test_linkageb <-ld_columns[which(ld_columns$dis<=ld_max_pairwise),] 
  ldtb <- data.table( test_linkageb , key = "Var1" )
  ldtc <- data.table( test_linkageb , key = "Var2" )
  use <- with(snpsum.col, (!is.na(MAF_pop1)))
  # Remove NA's
  use[is.na(use)] <- FALSE 
  # Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
  snp_loc  <- snpsum.col[use,c("chr_name","row_name","loc_cM","loc_bp")] # this is the location of each snp in cM and in bp
  dtb <- data.table(snp_loc, key = "loc_bp" )
  t_locationb <- ldtb[dtb, "row_name", nomatch=0]
  t_locationc <- ldtc[dtb, "row_name", nomatch=0]
  test_linkageb <- test_linkageb[order(test_linkageb$Var1),]
  test_linkageb$locus_a <- t_locationb
  test_linkageb <- test_linkageb[order(test_linkageb$Var2),]
  test_linkageb$locus_b <- t_locationc
  test_linkageb <- test_linkageb[order(test_linkageb$Var1),]
  test_linkage <- cbind(test_linkageb[,5],test_linkageb[,3], test_linkageb[,4], test_linkageb[,1])
  colnames(test_linkage) <- c("L1","r^2","Dist","location")
  test_split <- split(test_linkage,test_linkage$L1) #separating snps into individual dataframes
  test_split <- lapply(test_split,"[",,-1)
  # function to calculate the mean r^2 by locus in windows of the length specified in the settings
  func_rsqr <- function(dfx){stats.bin(dfx$Dist, dfx$`r^2`,breaks = seq(0,ld_max_pairwise,ld_resolution))}
  # this function is just to obtain the same output as the above funtion so the two outputs match each 
  # other and therefeore possible to merge them
  func_loc  <- function(dfx){stats.bin(dfx$Dist, dfx$location,breaks = seq(0,ld_max_pairwise,ld_resolution))}
  #applying functions and merging dataframes
  t_rsqr <- lapply(test_split,func_rsqr)
  t_rsqr_b <-lapply(t_rsqr,"[[",3)
  t_rsqr_c <- lapply(t_rsqr_b,"[",2,)
  t_rsqr_d <- lapply(t_rsqr_c, as.data.frame)
  t_loc <- lapply(test_split,func_loc)
  t_loc_b <- lapply(t_loc,"[[",3)
  t_loc_c <- lapply(t_loc_b,"[",2,)
  t_loc_d <- lapply(t_loc_c, as.data.frame)
  linkage <- Map(cbind,t_loc_d, t_rsqr_d)
  linkage_b <- rbindlist(linkage)
  linkage_b <- as.data.frame(linkage_b)
  # LD has been now calculated at the same maximum distance and at the same resolution for each snp. 
  # A vector (based on the maximum distance of LD and its resolution) is used to separate LD 
  # resolutions into bins so LD can be calculated in each genomic region
  linkage_b$factor <- 1:(ld_max_pairwise/ld_resolution)
  linkage_b <- linkage_b[complete.cases(linkage_b),]
  colnames(linkage_b) <- c("loc","rsqr","factor")
  linkage_c<-linkage_b[order(linkage_b$factor),]
  linkage_c$factor<- factor(linkage_c$factor)
  linkage_c<- as.data.frame(linkage_c)
  # separating the groups of LD resolutions into individual dataframes
  linkage_split<-split(linkage_c,linkage_c$factor) 
  linkage_split<-lapply(linkage_split,as.data.frame)
  # function to calculate the mean r^2 of each genomic region in windows of length of the region_size
  fun_rsqr<- function(df){stats.bin(df$loc, df$rsqr,  breaks = seq(0,chromosome_length,region_size))}
  fun_loc<- function(df){stats.bin(df$loc, df$factor,  breaks = seq(0,chromosome_length,region_size))}
  link_rsqr<-lapply(linkage_split,fun_rsqr)
  link_rsqr_b<-lapply(link_rsqr,"[[",3)
  link_rsqr_c<- lapply(link_rsqr_b,"[",2,)
  link_rsqr_d<-lapply(link_rsqr_c, as.data.frame)
  link <- list.cbind(link_rsqr_d)
  # this is the number of non-synonymous mutations at the end of the simulations
  use_b <- with(snpsum.col,(RAF_pop1>0)) 
  # Remove NA's
  use_b[is.na(use_b)] <- FALSE 
  # Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
  NS_sim <- snpsum.col[use_b,"loc_bp"] 
  number_nonsyn[,i] <- as.data.frame(table(cut(NS_sim, breaks=seq(0,chromosome_length,region_size))))[,2]
  # this is to  calculate the means of the SNP's statistics in bins of size of the region_size 
  snpsum.col_2 <- snpsum.col[use_b,] 
   snpsum.col_2 <- split(snpsum.col_2,cut(snpsum.col_2$loc_bp,breaks=seq(0,chromosome_length,region_size)))
    snpsum.col_2 <- lapply(snpsum.col_2,"[",,5:24)
  snpsum.col_2 <- lapply(snpsum.col_2,colMeans,na.rm=T)
  snpsum.col_2 <- list.rbind(snpsum.col_2)
  snpsum.col_mean[[i]] <- snpsum.col_2
  for (ld_distance in 1: ncol(link) ) {
    link_iterations[[ld_distance]][,i] <- link[,ld_distance]
  }
}
   return(list(link_iterations,number_nonsyn,snpsum.col_mean))
}

read.pedfile_b <- function (file, n, snps, which, split = "\t| +", sep = ".", na.strings = "0", lex.order = FALSE,show_warnings=T){
    r0 <- as.raw(0)
    r1 <- as.raw(1)
    r2 <- as.raw(2)
    r3 <- as.raw(3)
    con <- gzfile(file)
    open(con)
    if (missing(n)) {
        n <- 0
        repeat {
            line <- readLines(con, n = 1)
            if (length(line) == 0) 
                break
            n <- n + 1
        }
        if (n == 0) 
            stop("Nothing read")
        seek(con, 0)
    }
    gen <- missing(snps)
    map <- NULL
    if (!gen) {
        m <- length(snps)
        if (m == 1) {
            map <- read.table(snps, comment.char = "")
            m <- nrow(map)
            if (missing(which)) {
                which <- 1
                repeat {
                  snps <- map[, which]
                  if (!any(duplicated(snps))) 
                    break
                  if (which == ncol(map)) 
                    stop("No unambiguous snp names found on file")
                  which <- which + 1
                }
            }
            else {
                snps <- map[, which]
            }
        }
    }
    else {
        line <- readLines(con, n = 1)
        fields <- strsplit(line, split)[[1]]
        nf <- length(fields)
        if (nf%%2 != 0) 
            stop("Odd number of fields")
        m <- (nf - 6)/2
        seek(con, 0)
    }
    nf <- 6 + 2 * m
    result <- matrix(raw(n * m), nrow = n)
    ped <- character(n)
    mem <- character(n)
    pa <- character(n)
    ma <- character(n)
    sex <- numeric(n)
    aff <- numeric(n)
    rownms <- character(n)
    a1 <- a2 <- rep(NA, m)
    a1m <- a2m <- rep(TRUE, m)
    mallelic <- rep(FALSE, m)
    for (i in 1:n) {
        line <- readLines(con, n = 1)
        fields <- strsplit(line, "\t| +")[[1]]
        to.na <- fields %in% na.strings
        fields[to.na] <- NA
        ped[i] <- fields[1]
        mem[i] <- fields[2]
        pa[i] <- fields[3]
        ma[i] <- fields[4]
        sex[i] <- as.numeric(fields[5])
        aff[i] <- as.numeric(fields[6])
        alleles <- matrix(fields[7:nf], byrow = TRUE, ncol = 2)
        one <- two <- rep(FALSE, m)
        for (k in 1:2) {
            ak <- alleles[, k]
            akm <- is.na(ak)
            br1 <- !akm & a1m
            a1[br1] <- ak[br1]
            a1m[br1] <- FALSE
            br2 <- !akm & (a1 == ak)
            one[br2] <- TRUE
            br3 <- !akm & !a1m & (a1 != ak)
            br4 <- br3 & a2m
            a2[br4] <- ak[br4]
            a2m[br4] <- FALSE
            br5 <- br3 & (a2 == ak)
            two[br5] <- TRUE
            mallelic <- mallelic | !(akm | one | two)
        }
        gt <- rep(r0, m)
        gt[one & !two] <- r1
        gt[one & two] <- r2
        gt[two & !one] <- r3
        result[i, ] <- gt
    }
    close(con)
    if (any(a1m & show_warnings==T)){ 
        warning("no data for ", sum(a1m), " loci")
    }
    mono <- (a2m & !a1m)
    if (any(mono & show_warnings==T)){ 
         warning(sum(mono), " loci were monomorphic")
    }
    if (any(mallelic & show_warnings==T)) {
        result[, mallelic] <- r0
        warning(sum(mallelic), " loci were multi-allelic --- set to NA")
    }
    if (gen) 
        snps <- paste("locus", 1:m, sep = sep)
    if (any(duplicated(ped))) {
        if (any(duplicated(mem))) {
            rnames <- paste(ped, mem, sep = sep)
            if (any(duplicated(rnames))) 
                stop("could not create unique subject identifiers")
        }
        else rnames <- mem
    }
    else rnames <- ped
    dimnames(result) <- list(rnames, snps)
    result <- new("SnpMatrix", result)
    if (lex.order) {
        swa <- (!(is.na(a1) | is.na(a2)) & (a1 > a2))
        switch.alleles(result, swa)
        a1n <- a1
        a1n[swa] <- a2[swa]
        a2[swa] <- a1[swa]
        a1 <- a1n
    }
    fam <- data.frame(row.names = rnames, pedigree = ped, member = mem, 
        father = pa, mother = ma, sex = sex, affected = aff)
    if (is.null(map)) 
        map <- data.frame(row.names = snps, snp.name = snps, 
            allele.1 = a1, allele.2 = a2)
    else {
        map$allele.1 <- a1
        map$allele.2 <- a2
        names(map)[which] <- "snp.names"
    }
    list(genotypes = result, fam = fam, map = map)
}

reproduction <- function(pop, pop_number) {
  parents_matrix <- as.data.frame(matrix(nrow = population_size / 2, ncol = 2))
  # select males and females that are mating with each other - males column 1 females column 2
  # 20% of males do not mate, 30% mate once, 30% mate twice and 20% mate three times according to
  # Markow & Sawka 1992
  male_pool <- sample(rownames(pop[1:(population_size / 2), ]), size = population_size / 2)
  no_mating <- sample(male_pool, size = (population_size / 2) * 0.2)
  one_mating <- sample(male_pool[!male_pool %in% no_mating], size = (population_size / 2) * 0.3)
  two_matings <- sample(male_pool[!male_pool %in% c(no_mating, one_mating)], size = (population_size /2) * 0.3)
  three_matings <- sample(male_pool[!male_pool %in% c(no_mating, one_mating, two_matings)], size = (population_size / 2) * 0.2)
  mating_males <- c(one_mating, rep(two_matings, times = 2),rep(three_matings, times = 3))
  parents_matrix[, 1] <- sample(mating_males, size = population_size / 2)
  parents_matrix[, 2] <- sample(rownames(pop[((population_size / 2) + 1):population_size, ]), size =population_size / 2)
  offspring <- NULL
  for (parent in 1:dim(parents_matrix)[1]) {
    pairing_offspring <- rnbinom(1, size = variance_offspring, mu = number_offspring)
    offspring_temp <- as.data.frame(matrix(nrow = pairing_offspring, ncol = 6))
    if (pairing_offspring < 1) {
      next
    }
    offspring_temp[, 1] <- sample(c("Male", "Female"), size = pairing_offspring, replace = TRUE) # sex
    offspring_temp[, 2] <- pop_number # source population
    male_chromosomes <- list(pop[parents_matrix[parent, 1], 3], pop[parents_matrix[parent, 1], 4])
    female_chromosomes <- list(pop[parents_matrix[parent, 2], 3], pop[parents_matrix[parent, 2], 4])
    for (offs in 1:pairing_offspring) {
        males_recom_events <- rpois(1,recom_event)
        females_recom_events <- rpois(1,recom_event)
        #recombination in males
        if (recombination == TRUE & recombination_males == TRUE & males_recom_events > 1){
          for (event in males_recom_events){male_chromosomes<-recomb(male_chromosomes[[1]],male_chromosomes[[2]])}
          offspring_temp[offs, 3] <- male_chromosomes[[sample(c(1, 2), 1)]]
        }else{offspring_temp[offs, 3] <- male_chromosomes[[sample(c(1, 2), 1)]]}
      #recombination in females
        if (recombination == TRUE & females_recom_events > 1){
          for (event in females_recom_events){female_chromosomes<-recomb(female_chromosomes[[1]],female_chromosomes[[2]])}
      offspring_temp[offs, 4] <- female_chromosomes[[sample(c(1, 2), 1)]]
        }else{offspring_temp[offs, 4] <- female_chromosomes[[sample(c(1, 2), 1)]]}
        }
    offspring_temp[, 5] <- parents_matrix[parent, 1] #id father
    offspring_temp[, 6] <- parents_matrix[parent, 2] #id mother
    offspring <- rbind(offspring, offspring_temp)
    }
  return(offspring)
}

# this is the function for recombination
recomb <- function(r_chromosome1,r_chromosome2,r_map = recombination_map,loci = loci_number){
    chiasma <- as.numeric(sample(row.names(r_map),size = 1,prob = r_map[, "c"]))
    if (chiasma < (loci + 1)) {
      split_seqs <- strsplit(c(r_chromosome1, r_chromosome2), split = "")
      r_chr1 <- paste0(c(split_seqs[[1]][1:chiasma], split_seqs[[2]][(chiasma + 1):loci]), collapse = "")
      r_chr2 <- paste0(c(split_seqs[[2]][1:chiasma], split_seqs[[1]][(chiasma + 1):loci]), collapse = "")
      return(list(r_chr1, r_chr2))
    } else{
      return(list(r_chromosome1, r_chromosome2))
    }
}

initialise <- function(pop_number) {
  pop <- as.data.frame(matrix(ncol = 4, nrow = population_size))
  pop[, 1] <- rep(c("Male", "Female"), each = population_size / 2)
  pop[, 2] <- pop_number # second column stores population
  for (individual_pop in 1:population_size) {
    chromosome1 <- paste0(mapply(sample_alleles, q = reference$q,  USE.NAMES = F),collapse = "")
    chromosome2 <- paste0(mapply(sample_alleles, q = reference$q,  USE.NAMES = F),collapse = "")
    for (element in neutral_loci_location) {
      substr(chromosome1, start = element, stop = element) <-
        sample(as.character(c(1:msats)),size = 1,prob = c(rep(1 / msats, msats)))
      substr(chromosome2, start = element, stop = element) <-
        sample(as.character(c(1:msats)), size = 1, prob = c(rep(1 / msats, msats)))
    }
    if(experiment_freq == TRUE & simulation_type=="fly"){
      counter_msats <- 1
    for (element in 1:length(loc_exp_loci_2)) {
      substr(chromosome1, start = as.numeric(loc_exp_loci_2[element]), stop = as.numeric(loc_exp_loci_2[element])) <- sample(as.character(c(1:length(msats_freq[[counter_msats]]))),size = 1,prob = msats_freq[[counter_msats]])
      substr(chromosome2, start = as.numeric(loc_exp_loci_2[element]), stop = as.numeric(loc_exp_loci_2[element])) <- sample(as.character(c(1:length(msats_freq[[counter_msats]]))),size = 1,prob = msats_freq[[counter_msats]])
      counter_msats <- counter_msats + 1
      }
    }
    pop[individual_pop, 3] <- chromosome1
    pop[individual_pop, 4] <- chromosome2
  }
  return(pop)
}


sample_alleles <- function(alleles = c("a", "A"), q) {
  sample(alleles, size = 1, prob = c(q, 1 - q))
}

# allele frequency difference (AFD) from Berner 2019
AFD_fun <- function(allele_matrix) {
    total_count <-  colSums(allele_matrix,na.rm = T)
    allele_freq <- allele_matrix/total_count
    AFD <- sum(abs(allele_freq[,1] - allele_freq[,2]))/2
  return(AFD)
}

reproduction_2 <- function(pop, pop_number) {
  parents_matrix <- as.data.frame(matrix(nrow = population_size / 2, ncol = 2))
  # select males and females that are mating with each other - males column 1 females column 2
  # 20% of males do not mate, 30% mate once, 30% mate twice and 20% mate three times according to
  # Markow & Sawka 1992
  male_pool <- sample(rownames(pop[1:(population_size / 2), ]), size = population_size / 2)
  # no_mating <- sample(male_pool, size = (population_size / 2) * 0.2)
  # one_mating <- sample(male_pool[!male_pool %in% no_mating], size = (population_size / 2) * 0.3)
  # two_matings <- sample(male_pool[!male_pool %in% c(no_mating, one_mating)], size = (population_size /2) * 0.3)
  # three_matings <- sample(male_pool[!male_pool %in% c(no_mating, one_mating, two_matings)], size = (population_size / 2) * 0.2)
  # mating_males <- c(one_mating, rep(two_matings, times = 2),rep(three_matings, times = 3))
  parents_matrix[, 1] <- sample(male_pool, size = population_size / 2)
  parents_matrix[, 2] <- sample(rownames(pop[((population_size / 2) + 1):population_size, ]), size =population_size / 2)
  offspring <- NULL
  for (parent in 1:dim(parents_matrix)[1]) {
    pairing_offspring <-  rnbinom(1, size = variance_offspring, mu = number_offspring)
    offspring_temp <- as.data.frame(matrix(nrow = pairing_offspring, ncol = 6))
    if (pairing_offspring < 1) {
      next
    }
    offspring_temp[, 1] <- sample(c("Male", "Female"), size = pairing_offspring, replace = TRUE) # sex
    offspring_temp[, 2] <- pop_number # source population
    male_chromosomes <- list(pop[parents_matrix[parent, 1], 3], pop[parents_matrix[parent, 1], 4])
    female_chromosomes <- list(pop[parents_matrix[parent, 2], 3], pop[parents_matrix[parent, 2], 4])
    for (offs in 1:pairing_offspring) {
        males_recom_events <- rpois(1,recom_event)
        females_recom_events <- rpois(1,recom_event)
        #recombination in males
        if (recombination == TRUE & recombination_males == TRUE & males_recom_events > 1){
          for (event in males_recom_events){male_chromosomes<-recomb(male_chromosomes[[1]],male_chromosomes[[2]])}
          offspring_temp[offs, 3] <- male_chromosomes[[sample(c(1, 2), 1)]]
        }else{offspring_temp[offs, 3] <- male_chromosomes[[sample(c(1, 2), 1)]]}
      #recombination in females
        if (recombination == TRUE & females_recom_events > 1){
          for (event in females_recom_events){female_chromosomes<-recomb(female_chromosomes[[1]],female_chromosomes[[2]])}
      offspring_temp[offs, 4] <- female_chromosomes[[sample(c(1, 2), 1)]]
        }else{offspring_temp[offs, 4] <- female_chromosomes[[sample(c(1, 2), 1)]]}
        }
    offspring_temp[, 5] <- parents_matrix[parent, 1] #id father
    offspring_temp[, 6] <- parents_matrix[parent, 2] #id mother
    offspring <- rbind(offspring, offspring_temp)
    }
  return(offspring)
}


