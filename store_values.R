# Statistics are only calculated for the dispersal phase for computing efficiency. To analyze the 
# pre-adaptation phase, the simulations should be run using the dispersal phase using the parameters 
# of interest of the pre-adaptation phase. 
if (store_values == TRUE) {
  alleles_res <- list()
  for (locus in loci_location) {
    alleles <- sapply(c(pop1[, 3:4], pop2[, 3:4]), substr, locus, locus)
    alleles_res <- append(alleles_res, list(alleles))
  }
   # for (value in 2:length(loci_location)) {
       for (value in 1:length(loci_location)) {
    hierf[[value]][, gen_dispersal + 1] <- as.numeric(c(
        paste0(alleles_res[[value]][, 1], alleles_res[[value]][, 2]),
        paste0(alleles_res[[value]][, 3], alleles_res[[value]][, 4])
      ))
  }
  
  # this is to calculate the number of deleterious loci segregating in the population every generation
  freq <- apply(pop1, 1, split_seq)
  freq <- do.call(cbind, freq)
  lost_deleterious  <- freq[-as.numeric(neutral_loci_location), ]
  lost_deleterious[lost_deleterious == "a"] <- 1
  lost_deleterious[lost_deleterious == "A"] <- 0
  lost_deleterious[] <- lapply(lost_deleterious, as.numeric)
  lost_deleterious$fixation <- rowSums(lost_deleterious)
  q_deleterious <- lost_deleterious$fixation/(population_size*2)
  p_deleterious <- 1-q_deleterious
  lost_deleterious$he <- 2*p_deleterious*q_deleterious
  deleterious_eliminated[gen_dispersal, ] <- (nrow(lost_deleterious[lost_deleterious$fixation == 0, ]) / loci_number) * 100
  t1 <- reference[-as.numeric(neutral_loci_location), ]
  t2 <- cbind(t1[, 2:4], lost_deleterious$fixation)
  t3 <- t2[t2$`lost_deleterious$fixation` > 0, ]
  mean_h_gen <- mean(t3$h)
  mean_h_final[gen_dispersal, ] <- mean_h_gen
  mean_s_gen <- mean(t3$s)
  mean_s_final[gen_dispersal, ] <- mean_s_gen
  # this is to calculate the genetic load based on the equation:  2pqhs + q^s
  q_gen <- t3[, 4] / (population_size * 2)
  s_gen <- t3$s
  h_gen <- t3$h
  # this is genetic load based on additive selection
  g_load_a_gen <- sum((2 * (1 - q_gen) * q_gen * h_gen * s_gen) + ((q_gen ^ 2) * s_gen))
  g_load_a_final[gen_dispersal, ] <- g_load_a_gen
  # this is genetic load based on multiplicative selection
  g_load_m_gen <- 1 - prod(1 - ((2 * (1 - q_gen) * q_gen * h_gen * s_gen) + ((q_gen ^ 2) * s_gen)))
  g_load_m_final[gen_dispersal, ] <- g_load_m_gen
  
  if (generation == number_generations) {
    allele_matrix <- lapply(hierf, allele.count)
    results_generation <- lapply(hierf, basic.stats)
    ##### SUMMARY STATISTICS ####
    source('summary_stats.R')
    ##### LD ANALYSES #####
    # LD analyses are done with offspring and not with parents and using pop1. 
    # over 100 replications, LD patterns in pop1 and pop2 are the same.
    # LD_analyses <- sample(c(T,F,F,F),1)
    if(LD_analyses==T){
      ind_pop1<- nrow(pop1)
      ind_pop2<- nrow(pop2)
      ld_res_pop1 <- LD_fun(df_LD_fun=rbind(pop1,pop2), show_warnings=F)
      ld_columns_pop1 <- ld_res_pop1[[1]]
      snpsum.col_pop1 <- ld_res_pop1[[2]]
    }
  }
}