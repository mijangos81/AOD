
  first_het_pop1 <- mean_het_pop1[1,]
  first_het_pop2 <- mean_het_pop2[1,]
  
  Ne_100 <- 100
  rate_of_loss_100 <- 1 - (1 / (2 * Ne_100))
  expected_het_pop1_100 <- sapply(first_het_pop1, "*", (rate_of_loss_100 ^ c(1:gen_number_dispersal - 1)))
  expected_het_pop2_100 <- sapply(first_het_pop2, "*", (rate_of_loss_100 ^ c(1:gen_number_dispersal - 1)))
  
   Ne_125 <- 125
  rate_of_loss_125 <- 1 - (1 / (2 * Ne_125))
  expected_het_pop1_125 <- sapply(first_het_pop1, "*", (rate_of_loss_125 ^ c(1:gen_number_dispersal - 1)))
  expected_het_pop2_125 <- sapply(first_het_pop2, "*", (rate_of_loss_125 ^ c(1:gen_number_dispersal - 1)))
   
  Ne_150 <- 150
  rate_of_loss_150 <- 1 - (1 / (2 * Ne_150))
  expected_het_pop1_150 <- sapply(first_het_pop1, "*", (rate_of_loss_150 ^ c(1:gen_number_dispersal - 1)))
  expected_het_pop2_150 <- sapply(first_het_pop2, "*", (rate_of_loss_150 ^ c(1:gen_number_dispersal - 1)))
  
   Ne_175 <- 175
  rate_of_loss_175 <- 1 - (1 / (2 * Ne_175))
  expected_het_pop1_175 <- sapply(first_het_pop1, "*", (rate_of_loss_175 ^ c(1:gen_number_dispersal - 1)))
  expected_het_pop2_175 <- sapply(first_het_pop2, "*", (rate_of_loss_175 ^ c(1:gen_number_dispersal - 1)))
  

p6 <- print(ggplot(generations, aes(V1)) +
          geom_line(aes(y = rowMeans(mean_het_pop2), colour = "Simulations He"),size=2) +
            # geom_line(aes(y = rowMeans(mean_expected_het_pop2), colour = "Ne exp"),size=1/2) +
            geom_line(aes(y = rowMeans(expected_het_pop2_100), colour = "Ne 100"),size=1/2) +
            geom_line(aes(y = rowMeans(expected_het_pop2_125), colour = "Ne 125"),size=1/2) +
            geom_line(aes(y = rowMeans(expected_het_pop2_150), colour = "Ne 150"),size=1/2) +
            geom_line(aes(y = rowMeans(expected_het_pop2_175), colour = "Ne 175"),size=1/2) +
        theme_bw(base_size = 18) +
         scale_fill_hc("darkunica")+
        labs(x="GENERATIONS", y="He", title=NULL)+ 
      theme(legend.title=element_blank())+
        theme(legend.position =  "bottom") +
        theme(legend.text=element_text(size=14)), plot.margin=unit(c(-1,1,1,-1), "in") )
