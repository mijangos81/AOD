##### R SETTINGS#####
rm(list=ls())
starting_pops <- F
 # library(tictoc)
# seed <- as.numeric(format(Sys.time(), "%H%M%S"))
# set.seed(42)
# libraries for simulations
library(dplyr)
library(fields)
library(stringi)
library(readr)
library(plyr)
# libraries for selection tests
library(scales)
library(qvalue)
library(SuppDists)
library(stringi)
library(tibble)
library(rlang)
library(purrr)
# libraries for LD analysis
library(data.table)
library(rlist)
library(snpStats)
# # libraries for plots
library(ggplot2)
library(ggthemes)
##### FUNCTIONS#####
source('functions.R')
source('hierfstat.R')
source('Fdist_functions.R')
source('Bayescan_functions.R')
source('Outflank_functions.R')

if(starting_pops==T){
##### LOADING STARTING POPULATIONS, REFERENCE AND RECOMBINATION MAP DATA SETS #######
load("start_pop1.Rdata")
load("start_pop2.Rdata")
load("reference.Rdata")
load("recombination_map.Rdata")
load("location.Rdata")
}

 # for(gen_number_dispersal in c(12,26,34)){
  # for(selection in c(T,F)){
# tic()
 # for(mutation_rate in c(5*10^-4,5*10^-5,5*10^-6)){
  # for(same_line in c(FALSE,TRUE)){
     # message("dispersal_",c(gen_number_dispersal,"_",as.character(same_line)))
##### SIMULATIONS VARIABLES######
source('variables.R')
##### REFERENCE VALUES#####
source('reference.R')
##### DATAFRAMES TO STORE VALUES#######
source('dataframes.R')
##### START ITERATION LOOP #####
for (iteration in 1:number_iterations) {
    if (iteration %% 1 == 0) {message("iteration = ", iteration)}
  
  if(starting_pops==F){
    ##### VARIABLES PRE_ADAPTATION PHASE #######
    if (pre_adaptation == TRUE) {
        population_size <- population_size_pre_adaptation
        dispersal <- dispersal_pre_adaptation
        store_values <- FALSE
    }else{
        population_size <- population_size_dispersal
        dispersal <- dispersal_dispersal
        store_values <- TRUE
    }
  }
  
    if(starting_pops==T){
        population_size <- population_size_dispersal
        dispersal <- dispersal_dispersal
        store_values <- TRUE
    }
  
      if(starting_pops==F){
    ##### INITIALISE POPS #####
    pop1 <- initialise(pop_number = 1)
    pop2 <- initialise(pop_number = 2)
      }
  
  if(starting_pops==T){
  ##### INITIALISE POPS #####
      source("sample_starting_population.R")
    }
  
    ##### START GENERATION LOOP ######
    for (generation in 1:number_generations) {
            if(generation %% 10 == 0){message("generation = ", generation)}
        ##### VARIABLES DISPERSAL PHASE ######
        if (dispersal == TRUE) {
            res <- migration(population1=pop1,population2=pop2,generation=generation)
            pop1 <- res[[1]]
            pop1$V2 <- 1
            pop2 <- res[[2]]
            pop2$V2 <- 2
            maletran <- res[[3]]
            femaletran <- res[[4]]
        }
        if (generation == (gen_number_pre_adaptation + 1)) {
            population_size <- population_size_dispersal
            dispersal <- dispersal_dispersal
            store_values <- TRUE
            gen_dispersal <- 0 # counter to store values every generation
            if (pre_adaptation == TRUE) {
                if (same_line == TRUE) {
                    # pop1_temp is used because pop1 is used to sample pop2
                    pop1_temp <- rbind(pop1[sample(which(pop1$V1 == "Male"), size = population_size / 2), ],
                    pop1[sample(which(pop1$V1 == "Female"), size = population_size / 2), ])
                    pop1$V2 <- 1
                    pop2 <- rbind(pop1[sample(which(pop1$V1 == "Male"), size = population_size / 2), ],
                    pop1[sample(which(pop1$V1 == "Female"), size = population_size / 2), ])
                    pop2$V2 <- 2
                    pop1 <- pop1_temp
                } 
              if (same_line == FALSE) {
                    pop1 <- rbind(pop1[sample(which(pop1$V1 == "Male"), size = population_size / 2), ],
                    pop1[sample(which(pop1$V1 == "Female"), size = population_size /2), ])
                    pop1$V2 <- 1
                    pop2 <- rbind(pop2[sample(which(pop2$V1 == "Male"), size = population_size / 2), ],
                    pop2[sample(which(pop2$V1 == "Female"), size = population_size /  2), ])
                    pop2$V2 <- 2
                }
            }
        }
        # generation counter
        if (store_values == TRUE) {gen_dispersal <- gen_dispersal + 1}
        ##### REPRODUCTION#########
      if(simulation_type_2=="fly"){
          offspring_pop1 <- reproduction(pop = pop1, pop_number = 1)
          offspring_pop2 <- reproduction(pop = pop2, pop_number = 2)
      }
      if(simulation_type_2=="chillingham"){
        offspring_pop1 <- reproduction_2(pop = pop1, pop_number = 1)
        offspring_pop2 <- reproduction_2(pop = pop2, pop_number = 2)
      }
      if(simulation_type_2=="general"){
        offspring_pop1 <- reproduction_2(pop = pop1, pop_number = 1)
        offspring_pop2 <- reproduction_2(pop = pop2, pop_number = 2)
        }
        ##### SELECTION####
        if (selection == TRUE) {
            offspring_pop1 <- selection_fun(offspring=offspring_pop1,reference_pop = reference)
            offspring_pop2 <- selection_fun(offspring=offspring_pop2,reference_pop = reference)
        }
        ##### SAMPLING NEXT GENERATION #########
        source('sample_next_gen.R')
        ##### STORE VALUES########
        source('store_values.R')
    }
    ##### SAVING RESULTS ########
    # message("SAVING RESULTS")
    source('save_script.R')
}
    
 # }
 # }
  # toc()
            # }}}}
##### PRODUCING NEUTRAL SIMULATIONS FOR FDIST AND OUTFLANK #####
# if(neutral_simulations==T){
#   message("PRODUCING NEUTRAL SIMULATIONS FOR FDIST AND OUTFLANK")
 # source('neutral_simulations.R')
#   }
# #### LD ANALYSES #####
# if(LD_analyses==T){
   # message("LD ANALYSES")
  # source('do_LD_analysis.R')
#    }
# #### SELECTION TESTS ######
# if(tests_selection==T){
#   message("SELECTION TESTS")
    # source('selection_tests.R')
#   }
# ### STATS AVERAGES ACROSS REPLICATES #####
# message("STATS AVERAGES ACROSS REPLICATES")
    # source('stats_average.R')
 # #### STATS AVERAGES ACROSS REPLICATES #####
#  message("SEPARATE STATS")
      # source('separate_stats.R') 
# source("separate_stats_sd.R")
# ##### PLOTS #####
# if(plots == TRUE){
#   message("PLOTS")
   # source('AOD_plots.R')
# }
# #### PLOT HAPLOTYPES #####
# # message("PLOTTING HAPLOTYPES")
 # source('regions_AOD.R')
# }

# }


