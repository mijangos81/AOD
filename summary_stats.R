# for (loc in 2:length(location_neutral_loci_analysis)) {
   for (loc in 1:length(location_neutral_loci_analysis)) {
  AFD[[loc]] <- as.numeric(paste(lapply(allele_matrix[[loc]], AFD_fun)))
  sha_pop1[loc] <- c(lapply(allele_matrix[loc], function(x) {sha <- unlist(sapply(x, shannon)[1, ])}))
  sha_pop2[loc] <- c(lapply(allele_matrix[loc], function(x) {sha <- unlist(sapply(x, shannon)[2, ])}))
  sha_pops[loc] <- c(lapply(allele_matrix[loc], function(x) {sha <- sapply(as.vector(x),shannon_pops)}))
  # Mutual information
  shua[[loc]] <- as.numeric(paste(lapply(allele_matrix[[loc]], mutual_information)))
  # get the first value of Shannon index to calculate the rate of loss of shannon index
  first_sha_pop1[loc] <- sha_pop1[loc][[1]][1]
  first_sha_pop2[loc] <- sha_pop1[loc][[1]][1]
  expected_sha_pop1[[loc]] <- sapply(first_sha_pop1[[loc]], "*", (rate_of_loss ^ c(1:gen_number_dispersal - 1)))
  expected_sha_pop2[[loc]] <- sapply(first_sha_pop2[loc], "*", (rate_of_loss ^ c(1:gen_number_dispersal - 1)))
  # get the first value of He to calculate the rate of loss of heterozygosity
  first_het_pop1[loc] <- lapply(results_generation[loc], "[[", "Hs")[[1]][1, 1]
  first_het_pop2[loc] <- lapply(results_generation[loc], "[[", "Hs")[[1]][1, 2]
  expected_het_pop1[[loc]] <- sapply(first_het_pop1[[loc]], "*", (rate_of_loss ^ c(1:gen_number_dispersal - 1)))
  expected_het_pop2[[loc]] <- sapply(first_het_pop2[[loc]], "*", (rate_of_loss ^ c(1:gen_number_dispersal - 1)))
  het_pop1[[loc]] <- lapply(results_generation[loc], "[[", "Hs")[[1]][, 1]
  het_pop2[[loc]] <- lapply(results_generation[loc], "[[", "Hs")[[1]][, 2]
  # Heterozygosity averaged over the two populations
  het_ave[[loc]]  <- lapply(results_generation[loc], "[[", "perloc")[[1]][, 2]
  # The overall gene diversity
  Ht[[loc]] <- lapply(results_generation[loc], "[[", "perloc")[[1]][, 3]
  Fst[[loc]]  <- lapply(results_generation[loc], "[[", "perloc")[[1]][, 7]
  Fstp[[loc]] <- lapply(results_generation[loc], "[[", "perloc")[[1]][, 8]
  Dest[[loc]] <- lapply(results_generation[loc], "[[", "perloc")[[1]][, 10]
  # Shannon differentiation Sherwin et al. 2017
  sha_diff[[loc]] <- shua[[loc]]/log(2)
  # q-profile of ‘effective-number’ diversity metrics "D"
  one_D_beta[[loc]] <- exp(shua[[loc]])
  one_D_alpha_pop1[[loc]]<- exp(sha_pop1[[loc]])
  one_D_alpha_pop2[[loc]] <- exp(sha_pop2[[loc]])
  two_D_alpha_pop1[[loc]] <- 1/(1-het_pop1[[loc]])
  two_D_alpha_pop2[[loc]]<- 1/(1-het_pop2[[loc]])
  two_D_gamma[[loc]] <- 1/(1-Ht[[loc]])
  two_D_alpha_average[[loc]] <- 1/(1-het_ave[[loc]])
  two_D_beta[[loc]] <- two_D_gamma[[loc]] / two_D_alpha_average[[loc]]
}




