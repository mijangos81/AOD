
if(pre_adaptation == FALSE){gen_number_pre_adaptation <- 0}
if(neutral_simulations==T){selection <- F}
# This is the total number of generations
number_generations <- gen_number_pre_adaptation + gen_number_dispersal
location_neutral_loci_simulations <-  seq(map_resolution/2,(nrow(map_input)*map_resolution),map_resolution)
chromosome_length <- 23200000
recom_event <- ceiling(sum(recombination_map[,"c"]))
# loci_number_to_simulate <- 46340
loc_exp_loci <- c( 25,  50  ,53,  74  ,88,  89, 117, 133, 141, 156, 161, 219)

neutral_loci_location <-  as.character(location)
# neutral_loci_location <- neutral_loci_location[-1]

loci_location <- location

location_neutral_loci_analysis <- reference[unlist(loci_location),"location"]

freq_deleterious <- reference[-as.numeric(neutral_loci_location), ]
freq_deleterious_b <- mean(2 * (freq_deleterious$q) * (1 - freq_deleterious$q))
density_mutations_per_cm <- (freq_deleterious_b * nrow(freq_deleterious))
# / (reference[loci_number, "accum"] * 100number_of_generations)
loci_number <- nrow(recombination_map)-1

if(LD_analyses==T){
  # one is subtracted from the recombination map to account for the last row that was added in the 
  # recombination map to avoid that the recombination function crashes
  loci_number <- nrow(recombination_map)-1
  plink_map <- as.data.frame(matrix(nrow = loci_number,ncol = 4 ))
  plink_map[,1] <- chromosome_name
  plink_map[,2] <- rownames(recombination_map[-(loci_number+1),])
  plink_map[,3] <- recombination_map[-(loci_number+1),"accum"]
  plink_map[,4] <- recombination_map[-(loci_number+1),"locations_deleterious"]
  
}
 # MIGRATION VARIABLES
  # pick which sex is going to be transferred first
  if (number_transfers >= 2) {
    maletran <- TRUE
    femaletran <- TRUE
  } else if (number_transfers == 1) {
    maletran <- TRUE
    femaletran <- FALSE
  }

dispersal_rate <- (number_transfers / transfer_each_gen) / (population_size_dispersal)
Fst_expected <- 1 / ((4 * Ne_fst * dispersal_rate) * ((2 / (2 - 1)) ^ 2) + 1)
shua_expected <- (0.22 / (sqrt(2 * Ne_fst * dispersal_rate))) - (0.69 / ((2*Ne_fst) * sqrt(dispersal_rate)))
rate_of_loss <- 1 - (1 / (2 * Ne))

##### VARIABLES TO WRITE TO EACH FILE #####
variables_file_a <- c(
  "simulation_type",
"chromosome_name",
"pre_adaptation",
"population_size_pre_adaptation",
"gen_number_pre_adaptation",
"population_size_dispersal",
"gen_number_dispersal",
"number_transfers",
"transfer_each_gen",
"same_line",
"msats",
"recombination_occurring",
"recombination_males",
"c_gral",
"windows_gral",
"dispersal_dispersal",
"Ne",
"selection",
"natural_selection_model",
"intercept",
"rate",
"gamma_scale",
"gamma_shape",
"mutation_rate",
"loci_number",
"h_gral",
"s_gral",
"q_gral",
"fitness_for_repulsion",
# "gene_size",
"experiment_freq",
"number_offspring",
"variance_offspring",
"map_resolution",
"ld_max_pairwise",
"density_mutations_per_cm",  
"number_iterations"
)
variables_file_b <- c(
simulation_type,
chromosome_name,
pre_adaptation,
population_size_pre_adaptation,
gen_number_pre_adaptation,
population_size_dispersal,
gen_number_dispersal,
number_transfers,
transfer_each_gen,
same_line,
msats,
recombination_occurring,
recombination_males,
c_gral,
windows_gral,
dispersal_dispersal,
Ne,
selection,
natural_selection_model,
intercept,
rate,
gamma_scale,
gamma_shape,
mutation_rate,
loci_number,
h_gral,
s_gral,
q_gral, 
fitness_for_repulsion,
# gene_size,
experiment_freq,
number_offspring,
variance_offspring,
map_resolution,
ld_max_pairwise,
density_mutations_per_cm,
 number_iterations
)
variables_file <- cbind(variables_file_a,variables_file_b)

 if(simulation_type=="fly"){
    variables <-  paste0(
      "FLY",
       "_B", rate,
      "_C",   gamma_scale,
      "_D",  gamma_shape,
      "_F",  mutation_rate
      # "_G", fdist_pop_history[11,2], #  msats,
      # "_H", fdist_pop_history[3,2], #  pre_adaptation,
      # "_I", fdist_pop_history[22,2], #  gamma_scale,
      # "_J", fdist_pop_history[23,2], #  gamma_shape,
      # "_K", fdist_pop_history[21,2], #  rate,
      # # "_L", fdist_pop_history[20,2], #  intercept,
      # "_M", fdist_pop_history[24,2], #  mutation_rate,
      # # "_N", fdist_pop_history[27,2], # s_gral
      # "_O", fdist_pop_history[26,2], # h_gral,
      # "_P", fdist_pop_history[28,2], #q_gral,
      # "_Q", fdist_pop_history[14,2] #c_gral
      # "_R",fdist_pop_history[15,2], # windows_gral,
      # "_S",fdist_pop_history[31,2] # experiment_freq
      # "_T",fdist_pop_history[29,2] # fitness_for_repulsion,
      # "_rep_", iteration,
      # ".csv"
    )
  }

    # CREATING THE FOLDER TO SAVE RESULTS EACH TIME THE PROGRAM IS RUN
if(neutral_simulations==T){
  file.date <- format(Sys.time(), "%H%M%S")
path.folder_sim <- stri_join(getwd(), "/","SIM_Neutral_", variables,file.date, sep = "")
dir.create(file.path(path.folder_sim))
}
if(neutral_simulations==F){
file.date <- format(Sys.time(), "%H%M%S")
path.folder_sim <- stri_join(getwd(), "/","SIM_", variables,file.date, sep = "")
dir.create(file.path(path.folder_sim))
}

