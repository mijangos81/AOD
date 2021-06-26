if (length(which(offspring_pop1$V1 == "Male")) < population_size / 2 |
    length(which(offspring_pop1$V1 == "Female")) < population_size / 2) {
  print("EXTINCT")
  break()
}
if (length(which(offspring_pop2$V1 == "Male")) < population_size / 2 |
    length(which(offspring_pop2$V1 == "Female")) < population_size /2) {
  print("EXTINCT")
  break()
}
if (selection==FALSE){
  pop1 <- 
    rbind(offspring_pop1[sample(which(offspring_pop1$V1 == "Male"), size =  population_size / 2), ],
          offspring_pop1[sample(which(offspring_pop1$V1 == "Female"), size = population_size / 2), ])
  pop2 <- 
    rbind(offspring_pop2[sample(which(offspring_pop2$V1 == "Male"), size = population_size / 2), ],
          offspring_pop2[sample(which(offspring_pop2$V1 == "Female"), size = population_size / 2), ])
}
if (selection==TRUE & natural_selection_model=="absolute"){
  pop1 <-
    rbind(offspring_pop1[sample(which(offspring_pop1$V1 == "Male"), size = population_size / 2), ],
          offspring_pop1[sample(which(offspring_pop1$V1 == "Female"), size = population_size / 2), ])
  pop2 <-
    rbind(offspring_pop2[sample(which(offspring_pop2$V1 == "Male"), size = population_size / 2), ],
          offspring_pop2[sample(which(offspring_pop2$V1 == "Female"), size = population_size / 2), ])
}
if (selection== TRUE & natural_selection_model=="relative"){
  # We modeled selection as Lesecque et al. 2012: offspring are randomly selected to become parents of
  # the next generation in proportion to their relative fitnesses, for example, if we had four individuals
  # with fitnesses (W) of 0.1, 0.2, 0.3, and 0.2 the first individual would be selected on average
  # 0.1/(0.1+0.2+0.3+0.2)=0.125 of the time to become parent of the next generation.
  # The vector of probabilities used in sample is multiplied by two because in the selection function
  # (selection_fun), the proportional relative fitness was calculated for all offspring together, and
  # below the males and females are separeted in groups, with the objective that exactly the parents of the
  # next generation are half males and half females
  males_pop1 <- offspring_pop1[which(offspring_pop1$V1 == "Male"), ]
  females_pop1 <- offspring_pop1[which(offspring_pop1$V1 == "Female"), ]
  pop1 <- rbind(
    males_pop1[sample(row.names(males_pop1), size = (population_size / 2), prob= (males_pop1$relative_fitness * 2)),],
    females_pop1[sample(row.names(females_pop1), size = (population_size / 2), prob= (females_pop1$relative_fitness * 2)),]
  )
  males_pop2 <- offspring_pop2[which(offspring_pop2$V1 == "Male"), ]
  females_pop2 <- offspring_pop2[which(offspring_pop2$V1 == "Female"), ]
  pop2 <- rbind(
    males_pop2[sample(row.names(males_pop2), size = (population_size / 2), prob= (males_pop2$relative_fitness * 2)),],
    females_pop2[sample(row.names(females_pop2), size = (population_size / 2), prob= (females_pop2$relative_fitness * 2)),]
  )
}