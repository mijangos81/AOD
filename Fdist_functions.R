
getAlleles <- function(column, diploid, ndig){
### This function gets alleles from a vector of diploid genotypes in FSTAT format
### diploid is a logical TRUE or FALSE
### ndig is the number of digits in the coding for alleles
### returns a vector of alleles  
  
  if (diploid==TRUE){
			if (ndig==1){
				b <- cbind(substr(as.character(column), start=1,stop=1), substr(as.character(column), start=2,stop=2))
				}
			if (ndig==2){
				b <- cbind(substr(as.character(column), 1,2), substr(as.character(column), 3,4))
				}	
			 if (ndig==3){
				b <- cbind(substr(as.character(column), 1,3), substr(as.character(column), 4,6))
				}
		}	
		if (diploid==FALSE){	
				b <- column												 
		}	
	return(b)
}

getAllCounts <- function(pops, column, diploid, ndig){
### This function returns allele counts in each population
### pops is a vector assigning each individual to each population
### column is a vector of diploid genotypes
### diploid is a logical
### ndig is the number of digits used to code the genotypes
	a<- as.vector(getAlleles(column, diploid, ndig))
	if (diploid==TRUE){b<-tapply(a, list( c(pops,pops), a), length)}else{
	b<-tapply(a, list(pops,a), length)}
	b[which(is.na(b))] <- 0
	return(b)
	}

###############################################################################
CW1993.Beta.Eq <- function(AllCounts){
##########################################
### CW93 FST for infinite sample of allele freqs.
### Ref:Cockerham, C. C., and B. S. Weir. 1993. 
###### Estimation of gene flow from F-statistics. Evolution 855-863.
### This function takes a matrix of allele counts, 
### where rows are populations and columns are alleles.
### This function was developed from equations in Cockerham & Weir 1993
### Note that it assumes equal sample sizes
###########################################    
	# AllCounts is a matrix of allele counts at a locus in each population
  # Populations in rows and alleles are in columns
			pops <- rownames(AllCounts)
			numpops <- dim(AllCounts)[1]
			numalleles <- dim(AllCounts)[2]	
			sample.sizes <- rowSums(AllCounts)
      
      if(length(unique(sample.sizes))!=1){
        print("Error in CW1993.Beta.Eq(): sample sizes are not equal")
      }
			
		###### Uncorrected by sample size #####
		#######################################
			a2 <- AllCounts
			p <- AllCounts/sample.sizes
			X <- sum(p^2)
			Y <- sum(colSums(p)^2)
			M <- sample.sizes[1]	
			F0 <- (M*X-numpops)/((M-1)*numpops) 
			F1 <- (Y-X)/(numpops*(numpops-1))
			He <- 1-F1
			FST <- (F0-F1)/(1-F1)
			#FST <- 1 - (1-F0)/(1-F1); an alternative way of writing the previous line
			return(list(numer=F0-F1, He=He, FST=FST))
        # to calculate FST overall loci, sum(numer)/sum(He)
        # because the mean of a ratio is equal to the ratio of the sum(numerator)/sum(denominator)
			}

###############################################################################
CW1993.Beta.Beaumont <- function(AllCounts){
##########################################
### CW93 FST as implemented with sample size correction in FDIST2 
### Ref:Cockerham, C. C., and B. S. Weir. 1993. 
###### Estimation of gene flow from F-statistics. Evolution 855-863.
### Ref: Beaumont, M. A., and R. A. Nichols. 1996. 
###### Evaluating loci for use in the genetic analysis of 
###### population structure. P Roy Soc Lond B Bio 263:1619-1626.
### This function takes a matrix of allele counts, 
### where rows are populations and columns are alleles.
### This function was developed from directly from FDIST2 code
###########################################   
  #AllCounts is a matrix of allele counts at a locus in each population
  #Populations in rows and alleles are in columns

		pops <- rownames(AllCounts)
		numpops <- dim(AllCounts)[1]
		numalleles <- dim(AllCounts)[2]
		sample.sizes <- rowSums(AllCounts)
			
		a2 <- AllCounts
		p <- AllCounts/sample.sizes
		x0 <- sum((rowSums(AllCounts^2)-sample.sizes)/(sample.sizes*(sample.sizes-1)))
				# This is the line that differs from the Beta published in CW1993
				# In CW1993, equal sample sizes are assumed
				# Here, Beaumont applies a sample size correction to each population 
        # sample and then sums over all pops
			
		# The code below is included because it was copied from my.thetacal, 
    # but gives the same x0 as the line above
				#x0<-0
				#for (j in 1:numpops){
				#	x2 <- 0;
				#	for (i in 1:numalleles){
				#		x2<- x2 + AllCounts[j,i]*AllCounts[j,i]
						#print(x2)
				#	}
				#	x0 <- as.numeric(x0 + (x2-sample.sizes[j])/(sample.sizes[j]*(sample.sizes[j]-1)))
				#}
			
			# This code was copied from my.thetacal, 
      # but gives het1 that is the same as (1-F1) in WC1993
			  yy=0
			  for (j in 1:(numpops-1)){
				  for (k in (j+1):numpops){
					  y1=0;
					  for (i in 1:numalleles){
						  y1 <- y1 + AllCounts[j,i]*AllCounts[k,i]
						    #print(y1)
					   }
					  yy <- yy + y1/(sample.sizes[j]*sample.sizes[k])
				  }
			  }
			
			q2<- x0/numpops
			q3 <- 2*yy/(numpops*(numpops-1))
			
			het0 <- 1-q2
			het1 <- 1-q3 #this is what is output as heterozygosity by fdist2
			FST <- (q2-q3)/(1-q3)
  		#FST <- 1-het0/het1 #alternative way of writing    
			return(list(numer=(q2-q3), He=het1, FST=FST))
}

###############################################################################
WCtheta.FST.diploids.book <- function(AllCounts){
##########################################
### WC-Theta FST as implemented (but not used) in FDIST2 
### This function takes a matrix of allele counts, 
### where rows are populations and columns are alleles.
### This function was developed from directly from FDIST2 code
###### Matches "thetacal" function which is included, 
###### but not implemented, with fdist2 code
### nc is the term for the sample size correction
### Estimator taken from page 150 of Weir's book:
### B. S. Weir (1990).  Methods for discrete population genetic data.  
###### Sunderland, MA. Sinauer Publishers.
###########################################   
			pops <- rownames(AllCounts)
			numpops <- dim(AllCounts)[1]
			numalleles <- dim(AllCounts)[2]
			
			sample.sizes <- rowSums(AllCounts)
			a2 <- AllCounts
			p <- AllCounts/sample.sizes
			
		#######################################
			XX <- sum(a2*a2/sample.sizes)
			nbar <- mean(sample.sizes)
			
			q2 <- (XX-numpops)/((nbar-1)*numpops) 
			nc = 1.0/(numpops - 1.0)*(sum(sample.sizes)- sum(sample.sizes^2)/sum(sample.sizes));
			yy <- sum(colSums(a2)^2)
			q3 <- 1/(numpops*(numpops-1)*nbar*nc)*(yy - nbar*(nc-1.0)/(nbar-1.0)*XX) + 
        (nbar-nc)/(nc*(nbar-1.0))*(1.0-1.0/(numpops - 1.0)*XX);
			a0 <- 1-q2
			a1 <- 1-q3
			He <- a1 
			num <- a1-a0
			#FST <- 1 - a0/a1 alternative way of same equation
			FST <- num/He2
			
			return(list(numer=num, He=He, FST=FST))
}

###############################################################################
WCtheta.FST.Haploids.2Alleles<-function(AllCounts){
##########################################
### WC-Theta FST for finite sample of haploid allele freqs
### From page 145-148 of Weir's book
### B. S. Weir (1990).  Methods for discrete population genetic data.  
###### Sunderland, MA. Sinauer Publishers.
### This function is used to calculate FST from the haploid IM simulations
###########################################  
	# Input a matrix of the counts of each allele (columns) in each population (rows)
	if(ncol(AllCounts)!=2){print("Error in WC.FST.Haploids.2Alleles: 
                              only 2 allelic states allowed" )}	
  
	n.pops<-nrow(AllCounts)
	counts1 <- AllCounts[,1]
	sample.sizes <- rowSums(AllCounts)
	n.ave <- mean(as.numeric(sample.sizes))
	n.c = (n.pops*n.ave - sum(sample.sizes^2)/(n.pops*n.ave))/(n.pops-1)	
	p.freqs = counts1/sample.sizes
	p.ave = sum(sample.sizes*p.freqs)/(n.ave*n.pops)
		#note: this differs slightly from mean(p.freqs) in R
	He2 <- 2*p.ave*(1-p.ave)
			
	s2 = sum(sample.sizes*(p.freqs - p.ave)^2)/((n.pops-1)*n.ave)
		#note: this differs slightly from var(p.freqs) in R
  
	T1 <- s2 - 1/(n.ave-1)*(p.ave*(1-p.ave) -(s2*(  n.pops-1)/  n.pops))
	T2 <- (n.c - 1)*(p.ave*(1-p.ave))/(n.ave-1) + 
            (1 + (  n.pops-1)*(n.ave-n.c)/(n.ave-1))*s2/  n.pops

  FST <- T1/T2 
	return(list(numer=T1, denom=T2, He=He2, FST=FST))
}

###############################################################################
WCtheta.FST.dataset <- function(data1, diploid, ndig){
##########################################
### WC-theta FST for finite sample of diploid allele freqs
### return FST and He estimate for each locus in dataset
### Based on Weir and Cockerham's 1984 paper
### Using wc() function from package hierfstat
###########################################

	FST.list <- wc(data1, diploid=diploid)
	FST <- as.numeric(FST.list$per.loc$FST)
  writeLines("Calculating He and FST (Weir & Cockerham 1984) for each locus...")

  numer <- FST.list$sigma.loc$lsiga
  denom <- FST.list$sigma.loc$lsiga + FST.list$sigma.loc$lsigb +
            FST.list$sigma.loc$lsigw
  
	FST.overall <- as.numeric(FST.list$FST)
  FIS.overall <- as.numeric(FST.list$FIS)
  He <- as.numeric(getHe.dataset(data1, diploid, ndig))
	locnames <- colnames(data1)[-1]

	return(list(FST.overall=FST.overall, FIS.overall=FIS.overall, 
              FST.emp =data.frame(locnames =locnames, He=He, FST=FST, 
              numer=numer, denom=denom)))
}

###############################################################################
mean.of.ratios <- function(numer, denom){
  #numer and denom are each vectors
  sum(numer)/sum(denom)
}

###############################################################################
CW1993.dataset <- function(data1, diploid, ndig){
########################   
### getCW1993 from a dataframe in FSTAT format
### df is the dataframe
### diploid is TRUE for diploids and FALSE for haploids
### ndig is the number of digits used to code genotypes
########################
  pops <- data1[,1]
  nloc <- ncol(data1)-1
	locnames <- colnames(data1)[-1]

	F0F1 <- matrix( unlist(apply(data1[,-1], MARGIN=2,get.CW93.locus, diploid,
                               ndig, pops)), 
                  ncol= 3, byrow=TRUE)
	colnames(F0F1) <- c("Numer", "He", "Fst")
  F0F1[which(is.na(F0F1))]<-NA
  numer<-F0F1 [,1] 
  He <- F0F1[,2]
  FST <- F0F1[,3]
  FST.overall <- mean.of.ratios(numer,He)
  return(list(FST.overall=FST.overall, FST.emp =data.frame( He=He, FST=FST, 
              numer=numer, locnames =locnames)))
  # in the below original version of the code the denom variable to return is not 
  # declared within the function, so it was deleted
  # return(list(FST.overall=FST.overall, FST.emp =data.frame( He=He, FST=FST, 
  #             numer=numer, denom=denom, locnames =locnames)))
}

get.CW93.locus <- function(vector, diploid, ndig, pops){
### Internal function to CW1993.dataset
    	countMat <- getAllCounts(pops, vector, diploid, ndig)
      return(CW1993.Beta.Beaumont(countMat))
}      

get.pval <- function(beta.obs, alpha.obs, beta.sim, alpha.sim,beta_max){
########################   
##### Get p value for a single observed FST, based on simulated values
########################
	# FST.obs is a single value for the observed locus
	# He.obs is a single value for the observed locus
	# FST.sim is the vector of simulated FSTs
	# He.sim is the vector of simulated He values
	# He.bin is the bin around the observed He for which the p-value will be estimated
    # default is 0.04 following Excoffier et al (2009)
### Get the list of FST values that are within the Heterozygosity-bin of the observed He
  # location of the observed He within the vector He.sim
 	  loc_alpha.obs <- findInterval(x=alpha.obs,vec = alpha.sim)
 	  length_alpha.sim <- length(alpha.sim)
  if (loc_alpha.obs < bins_Johnson_distribution){
    #if the observed alpha is close to the lower bound, then select the values of the first bin
     alpha.ind <- 1:(bins_Johnson_distribution*2)
  }else if (loc_alpha.obs > (length_alpha.sim-bins_Johnson_distribution)){
     #if the observed alpha is close to the upper bound, then select the values of the last bin
 	     alpha.ind <- (length_alpha.sim-(bins_Johnson_distribution*2)):length_alpha.sim
 	  }else{
 	    alpha.ind <- (loc_alpha.obs-bins_Johnson_distribution):(loc_alpha.obs+bins_Johnson_distribution)
 	  }
 	  # The below code is the used in Lotterhos & Whitlock. The upper version of the code was use 
 	  # instead to specify the number of points that each bin has, as described in Beaumont & Nichols 1996
#   MinHe <- min(He.sim)
# 	MaxHe <- max(He.sim)
# 	if ((He.obs-MinHe)<He.bin/2){ 
#     #if the observed He is near the lower bound, then arrange the bin
# 			He.ind <- which(He.ind==He.obs which(He.sim<(MinHe+He.bin)  & He.sim>MinHe)
# 
# 	  			He.ind <- which(He.sim<(MinHe+He.bin)  & He.sim>MinHe)
# 	}else if((MaxHe-He.obs)<He.bin/2){
#     #if the observed He is near the upper bound, then arrange the bin
# 			He.ind <- which(He.sim<(MaxHe)  & He.sim>(MaxHe-He.bin))
# 	}else{
# 			He.ind <- which(He.sim<(He.obs+He.bin/2) & He.sim>(He.obs-He.bin/2))
# 	}	
	beta.alpha <- beta.sim[alpha.ind]	
### Fit the simulated betas to a Johnson distribution and get the density
	minbeta <- min(c(min(beta.sim), beta.obs))
	# beta_max take in account the difference in the maximumn value of Shua and FST
	beta <- seq(minbeta,beta_max,by=0.0001)
		# just to make sure "1" is included in sequence
	
	# beta.dens <- dJohnsonSU(x=beta, params=eJohnsonSU(beta.alpha)[1:4])
		beta.dens <- dJohnson(x=beta, parms=JohnsonFit(beta.alpha))
	beta.dens[which(beta.dens=="NaN")]=0
		# 0 probabilities are undefined in the function
### Normalize the density and get the p-value for the observed beta	
	beta.dens2 <- beta.dens/sum(beta.dens)
	beta.ind <- max(which(beta<=beta.obs))
	p.val<- sum(beta.dens2[1:beta.ind])
	return(p.val)
}
 
get.pval.dataset <- function(beta.empDF, beta.simDF,beta_max){
##########################################
### Applies the get.pval function to all loci in the dataset
### beta.empDF is a dataframe with a column for the observed He and beta for each locus
    # it is OK for this dataframe to have other columns, P-values will just be appended
### beta.simDF is a dataframe with a column for simulated He and beta for several thousand loci
###### note that these two dataframes should not have the same number of rows
### write.progress will write a line every 1000 loci
###########################################  
 
	beta.empVect <- beta.empDF$beta
	alpha.empVect <- beta.empDF$alpha
	
	beta.simVect <- beta.simDF$beta
	alpha.simVect <- beta.simDF$alpha
	
	ntimes <- length(beta.empVect)
	p.val <- rep(NA, ntimes)
	for (i in 1:ntimes){
			p.val[i] <- get.pval(beta.obs=beta.empVect[i], alpha.obs=alpha.empVect[i],  
                           beta.sim=beta.simVect, alpha.sim=alpha.simVect,beta_max = beta_max)
	}
	out1 <- as.data.frame(cbind(beta.empDF, p.val.cum=p.val))
	return(out1)
}

correct.pval.dataframe <- function(dataframe, p.colName="p.val.cum", p.colNum){
##########################################
### Correct a list of cumulative p-values to indicate tails,
###### and indicate significance at the Bonferroni, FDR=0.05, and FDR=0.01 levels
###### Depends on Storey's qvalue.R function
### infilepath is path to the file, if it is not yet loaded as a dataframe
### assumes the column in the dataframe is named p.val.cum
    ### if the column has another name, please specify name and number of the column 
    ### the p-value column should be based on cumulative probabilities 
    ###	(i.e. starting at 0 in the left tail and ending at 1 in the right tail)
### In the output, "L" indicates left-tail, and "R" indicates right-tail	
### if write.outfile==TRUE
    ### will write to a file in outfilepath, but with the extension ".Cpval"
    ### if outfilepath is NOT specified, will write to infilepath with the extension ".Cpval"
    ### if neither outfile path or infilepath is specified, will give an error

  if (p.colName!="p.val.cum"){
    p.val <- assign(p.colName, dataframe[,p.colNum])
  }else{
    p.val <- dataframe$p.val.cum
  }
		
	L.p <- p.val
	R.p <- 1-L.p
	num.obs <- length(L.p)
	Bonf.p <- 0.05/num.obs
	
  ### Get qvalues for left hand side of distribution
 # If you do not have enough p-values such that you do not have coverage in the interval (0,1) then you should set qvalue(p, pi0=1) since there will not be enough information for a reliable pi0 estimate (to get a reliable pi0 estimate you usually need a few hundred p-values). Setting pi0=1 is equivalent to the Benjimini-Hochberg procedure.

	   q2 <- qvalue(L.p, pi0=1)
  	L.q <- q2$qvalues
  		
  ### Get qvalues for right hand side of distribution	
  	q3 <- qvalue(R.p, pi0=1)
  	R.q <- q3$qvalues
  			
  ### Convert p-values and q-values to the right and left sides of the distribution
  	p.val.tail <- L.p
  	p.val.tail[L.p>0.5] <- R.p[L.p>0.5] 
  	
  	qval <- L.q
  	maxq <- min(q2$pi0, q3$pi0) 
      #if pi0 is different for the left and right sides, take the minimum
  	qval[L.q<maxq] <- L.q[L.q<maxq]
  	qval[R.q<maxq] <- R.q[R.q<maxq]
  		
  	Tail <- L.p
  	Tail[L.p>0.5] <- "R"
  	Tail[L.p<=0.5] <- "L"
  	tail <- as.factor(Tail)
  		
  	Bonf <- rep(FALSE, num.obs)
  	Bonf[p.val.tail<Bonf.p] <- TRUE
  
  	FDR.01 <- rep(FALSE, num.obs)
  	FDR.01[qval<0.01] <-  TRUE
  	
  	FDR.05 <- rep(FALSE, num.obs) 
  	FDR.05[qval<0.05] <-  TRUE
  	
  	out <- data.frame(dataframe, tail=tail, p.val.tail=p.val.tail, qval.tail=qval, 
                      Bonf=Bonf, FDR.01=FDR.01, FDR.05=FDR.05)

	return(out)
}

getCI <- function(vector, percentile){
  a<-sort(vector)
  index2=round(length(a)*percentile,0)
  index1=round(length(a)*(1-percentile),0)
  return(c(a[index1], a[index2]))
}

get.CI.alpha <- function(xseq, percentile, betas.sim) {
  CI.alpha<-  matrix(NA, length(xseq),2)
  for (i in 1:(length(xseq)-1)){
     indexes <- which(xseq[i]<betas.sim$alpha & xseq[i+1]>=betas.sim$alpha)
     if(length(indexes>0)){
     CI.alpha[i,]<-getCI(betas.sim$beta[indexes], percentile)
     }
  }
  return(CI.alpha)
}

makeCIplot <- function(percentile,c, betas.sim,xseq=xseq){
  CI.alpha<-get.CI.alpha(xseq, percentile, betas.sim)
  points(xseq, CI.alpha[,1], type="l", col=c, lwd=1)
  points(xseq, CI.alpha[,2], type="l", col=c, lwd=1)
}

make.alpha.beta.plot <- function(betas.sim, All.Pvals, truepos, FODR.my.plot,title_plot){
  maxalpha <- max(c(betas.sim$alpha, All.Pvals$alpha))
  minalpha <- min(c(betas.sim$alpha, All.Pvals$alpha))
  maxbeta <- max(c(betas.sim$beta, All.Pvals$beta))
  minbeta <- min(c(betas.sim$beta, All.Pvals$beta))
  
  alpha.cell <- 0.02
  beta.cell <- 0.02
  
  xseq<-seq(minalpha,maxalpha,alpha.cell)
  yseq<-seq(minbeta,maxbeta,beta.cell )
  
  plot(NULL, NULL, xlim=c(minalpha,maxalpha), ylim=c(min(yseq), max(yseq)), 
       bty="l", las=1, 
       ylab="BETA", xlab="ALPHA", main= title_plot)
  
  points(All.Pvals$alpha, All.Pvals$beta, col=alpha("darkolivegreen",0.5), pch=19, cex=1)

  if(FODR.my.plot==TRUE){
     # makeCIplot(0.75,"red", betas.sim,xseq=xseq)
    # makeCIplot(0.8,"orange", betas.sim,xseq=xseq)
    # makeCIplot(0.85,"yellow", betas.sim,xseq=xseq)
     makeCIplot(0.975,"green", betas.sim,xseq=xseq)
    # makeCIplot(0.95,"blue", betas.sim,xseq=xseq)
    # makeCIplot(0.99,"purple", betas.sim,xseq=xseq)
     points(All.Pvals$alpha[truepos], All.Pvals$beta[truepos], col="deeppink", pch=19, cex=0.5)
  }else{
    makeCIplot(0.95, "firebrick1", betas.sim,xseq=xseq)
    # makeCIplot(0.999,"blue", betas.sim,xseq=xseq)
     points(All.Pvals$alpha[truepos], All.Pvals$beta[truepos], col="deeppink", pch=19, cex=0.5)
    
  }
}
