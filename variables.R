number_iterations <- 10
# if(gen_number_dispersal==34){
# number_transfers <- 1
# transfer_each_gen <- 8
# }
# if(gen_number_dispersal==26){
# number_transfers <- 1
# transfer_each_gen <- 2
# }
# if(gen_number_dispersal==12){
# number_transfers <- 2
# transfer_each_gen <- 1
# }
number_transfers <- 1
transfer_each_gen <- 29
gen_number_dispersal <- 100
same_line <- T # populations should be sampled form the same founding population?
######## VARIABLES FINAL FLY SIMULATIONS #######
# this is used just to determine q
mutation_rate <- 5*10^-5
  # 5*10^-5 # this is for fly 
  # 1*10^-3 # this is for chillingham
log_mean <-  0.002
  # 0.02 # this is for chillingham
  # 0.002 # this is for fly 
dominance_mean <- 0.2
selection <- T
#############
log_sd <- 4
targets_factor <- 1
  # 0.04 # this is for chillingham
  # 1 # this is for fly 
# OVERALL VARIABLES
simulation_type <- "fly" # simulation_type options: "general" or "fly"
simulation_type_2 <- "fly" # simulation_type_2 options: "chillingham", "fly" or "general"
msats <- 4 # number of alleles for each neutral locus
chromosome_name <- "2L" # options: 1 to 29 for chillingham simulations and 2L, 2R, 3L, 3R and X for fly simulations
# POP HISTORY VARIABLES
dispersal_pre_adaptation <- F # dispersal occurs during pre_adaptation period?
dispersal_dispersal <- T # dispersal occurs during dispersal phase period?
recombination <- T
recombination_males <- F
pre_adaptation <- F # pre_adaptation period?
number_offspring <- 10
population_size_pre_adaptation <- 1000 # population size must be even
population_size_dispersal <- 50 # population size must be even
Ne <- 16
Ne_fst <- 14
gen_number_pre_adaptation <- 50
# FLY SIMULATIONS VARIABLES
# this option has to be updated in the initialize function before using it 
experiment_freq <- T
experiment_loci <- T
# these are the frequencies of all the markers of the experiment in chromosome 2L

location_msats_experiment <- c(2373262,4960235,7040194,8325394,8348440,11015745,12507696,13153885,14615705,14995570,20706003)
msats_freq <- list(c(0.43,0.57),c(0.38,0.62),c(0.4,0.6),c(0.04,0.09,0.06,0.16,0.15,0.5),c(0.21,0.79),c(0.43,0.57),c(0.11,0.41,0.23,0.13,0.12),c(0.66,0.34),c(0.53,0.47),c(0.3,0.05,0.41,0.24),c(0.03,0.74,0.15,0.08))
# GENERAL SIMULATIONS VARIABLES
testing_sel_tests <- F
loci_number_to_simulate <- 2000 #total number of loci to simulate in the general simulations
h_gral <- 0.2
s_gral <- 0.005
c_gral <- 0.05
q_gral <- 0.1465
windows_gral <- 250 # number of windows of size map_resolution. windows_gral * map_resolution = chromosome length
# SELECTION VARIABLES
natural_selection_model <- "relative" # natural_selection_model options: "absolute" or "relative"
genetic_load <- 0.8
# distribution of selection coefficients
gamma_scale <- 0.03 # scale of the gamma distribution to sample sel. coeff.
gamma_shape <- 0.35 # shape of the gamma distribution to sample sel. coeff.
intercept <- 0.5 # intercept for the dominance equation
rate <- 500 # this the rate from the dominance equation
fitness_for_repulsion <- F
gene_size <- 10000
# LD ANALYSES VARIABLES
LD_analyses <- F
ld_max_pairwise <- 1000000  # maximun distance, in basepairs, at which pairwise LD should be calculated
ld_resolution <- 500000 # resolution, in basepairs, at which LD should be measured
region_size <- 1000000 # the size of the window at which the LD statistics and number of NS are calculated
# SELECTION TESTS VARIABLES
neutral_simulations <- F
tests_selection <- F
B_N_test <- F # Beaumont and Nichols' test (FDIST2)
B_B_test <- F # Balding and Beaumont's test (BAYESCAN)
W_L_test <- F # Whitlock and Lotterhos' test (OUTFLANK)
# CONSTANT VARIABLES
map_resolution <- 100000
variables_number <- 36 # these are the number of variables of the model
if(simulation_type_2 == "chillingham"  & simulation_type=="fly"){
variance_offspring <-  1000000 # variance in family size
###### input recombination map 
RecRates_All_Chromosomes <- read_csv("chill_recom_map.csv") 
RecRates_All_Chromosomes <- as.data.frame(RecRates_All_Chromosomes)
RecRates_All_Chromosomes$Chr <- as.character(RecRates_All_Chromosomes$Chr)
RecRates_All_Chromosomes_chr <- RecRates_All_Chromosomes[which(RecRates_All_Chromosomes$Chr==chromosome_name),]
chr_length <- RecRates_All_Chromosomes_chr[nrow(RecRates_All_Chromosomes_chr),"Location"]
map_cattle_binned <- stats.bin(RecRates_All_Chromosomes_chr$Location,RecRates_All_Chromosomes_chr$`mean sexes`,breaks = seq(0,chr_length,map_resolution))
map_cattle_binned_b <- unname(map_cattle_binned$stats[2,]* map_cattle_binned$stats[1,])
map_cattle_binned_b
map <- as.data.frame(map_cattle_binned_b)
colnames(map) <- "cM"
map[is.na(map$cM),] <- 0
targets <-  read_csv("chill_targets_of_selection.csv") 
  # read_csv("chill_targets_of_selection_exons.csv") 
  # read_csv("chill_targets_of_selection.csv") 
targets$chr_name <- as.character(targets$chr_name)
targets <- targets[which(targets$chr_name==chromosome_name),]
targets <- targets[!duplicated(targets$start),]
targets <- targets[!duplicated(targets$end),]
}
if(simulation_type_2 == "fly" & simulation_type=="fly"){
variance_offspring <-  0.4 # variance in family size
###### input recombination map 
map <- read_csv("fly_recom_map.csv")
map$Chr <- as.character(map$Chr)
map <- map[which(map$Chr==chromosome_name),]
map <- as.data.frame(map$cM/1000)
colnames(map) <- "cM"
map[is.na(map$cM),] <- 0
chr_length <- (nrow(map)+1) * map_resolution
targets <- read.csv("fly_targets_of_selection_exons.csv") 
targets$chr_name <- as.character(targets$chr_name)
targets <- targets[which(targets$chr_name==chromosome_name),]
targets <- targets[!duplicated(targets$start),]
targets <- targets[!duplicated(targets$end),]
}
if(simulation_type_2=="general"){
  variance_offspring <-  1000000 # variance in family size
}

# targets_tests <- targets
# targets_tests <- targets_tests[!duplicated(targets_tests$start),]
# targets_tests <- targets_tests[!duplicated(targets_tests$end),]
# 
# # 
#  targets_tests$midpoint <- (targets_tests$start + targets_tests$end)/2
#  # targets$targets <- (targets$s * targets$CUB_dmel) + targets$ns
#   targets_tests$targets <- targets_tests$ns
# 
#  break_bins <-  seq(0,chr_length,1000000)
# 
# test_bins <- stats.bin(x= targets_tests$midpoint,y= targets_tests$targets,breaks=break_bins)
# test_bins_2 <- test_bins$stats[2,] * test_bins$stats[1,]
# # 
# # ns_test <- read_csv("ns_2L_BDGP6.csv")
# str(ns_test)
# ns_test_2<-as.data.frame(cbind(ns_test,1))
# colnames(ns_test_2) <- c("ns","test")
# test_bins_ns <- stats.bin(x= ns_test_2$ns,y= ns_test_2$test,breaks=break_bins)
# test_bins_ns_2 <- test_bins_ns$stats[1,]
# 
