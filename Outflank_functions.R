OutFLANK <- function(FstDataFrame, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples, qthreshold=0.05){
  #Setting up necessary columns in dataframe
  Fstdata= outputDFStarterNoCorr(FstDataFrame,Hmin)
  # making working dataframe with real Fst (no NAs), storing NAs to add back later
  # This also removes loci with He values lower than Hmin from the working data frame
  nonkeepers = which((is.na(Fstdata$FSTNoCorr))|(Fstdata$He<Hmin))
  if(length(nonkeepers)>0){ 
      workingDataFrame = Fstdata[-nonkeepers,]
  }else{
    workingDataFrame = Fstdata
  }
  storedDataFrameNA = Fstdata[nonkeepers,]
  #Finding upper and lower bounds for trimming (eliminating NAs, but not negative FSTs)
  sortedDataFrame=workingDataFrame[order(workingDataFrame$FSTNoCorr),]
  NLociTotal=length(sortedDataFrame$FSTNoCorr)
  SmallestKeeper=ceiling(NLociTotal*LeftTrimFraction)
  LargestKeeper=floor(NLociTotal*(1-RightTrimFraction))
  LowTrimPoint=sortedDataFrame$FSTNoCorr[[SmallestKeeper]]
  HighTrimPoint=sortedDataFrame$FSTNoCorr[[LargestKeeper]]
  if(LowTrimPoint<0) {writeLines("ERROR: The smallest FST in the trimmed set must be > 0. Please use a larger LeftTrimFraction."); return()}
  if(HighTrimPoint>=1) {writeLines("ERROR: The largest FST in the trimmed set must be < 1. Please use a larger RightTrimFraction."); return()}
  #finding dfInferred and Fstbar iteratively  
  putativeNeutralListTemp=ifelse(workingDataFrame$FSTNoCorr>0,TRUE,FALSE)
  oldOutlierFlag=rep(FALSE,NLociTotal)
  #Note: All negative FST loci are marked as putative outliers, which will need
  #to be tested with the coalescent model later. In the meantime, they are
  #removed so as to not confuse the likelihood function.
  keepGoing=TRUE
  count = 0
  #writeLines(paste(mean(workingDataFrame$FSTNoCorr[putativeNeutralListTemp])))
  while(keepGoing){
    count=count+1
    if(count>19) {
      keepGoing=FALSE  
      writeLines("Exceeded iteration maximum.") ###Try with increased maximum value for count two lines above.
    }
    FstbarNoCorrTemp=fstBarCalculatorNoCorr(workingDataFrame[putativeNeutralListTemp,])  
    dfInferredTemp=EffectiveNumberSamplesMLE(FstVect=workingDataFrame$FSTNoCorr[putativeNeutralListTemp],Fstbar=FstbarNoCorrTemp,NumberOfSamples=NumberOfSamples,SmallestFstInTrimmedList=LowTrimPoint,LargestFstInTrimmedList=HighTrimPoint)
    workingDataFrame = pOutlierFinderChiSqNoCorr(workingDataFrame,FstbarNoCorrTemp,dfInferredTemp,qthreshold, Hmin)
    #### mark all negative FSTs as outliers if lowest nonneg FST is outlier
    #### (because negative Fst estimates can't be evaluated through the
    #### chi-square approach on their own)
    if(any(workingDataFrame$OutlierFlag[workingDataFrame$FSTNoCorr<LowTrimPoint])) workingDataFrame$OutlierFlag[workingDataFrame$FSTNoCorr<0]=TRUE
    ####Any loci previously marked as $OutlierFlag=TRUE remain so, even if the new iteration doesn't flag them as outliers
    #     workingDataFrame$OutlierFlag=!as.logical((!workingDataFrame$OutlierFlag)*(!oldOutlierFlag))
    #Resetting neutral list, and checking whether the outlier list has stabilized
    putativeNeutralListTemp=ifelse((!workingDataFrame$OutlierFlag),TRUE,FALSE)
    if(sum(putativeNeutralListTemp)==0) {writeLines("No loci in neutral list..."); return("FAIL")}
    if(identical(oldOutlierFlag,workingDataFrame$OutlierFlag)) keepGoing=FALSE
    ######if all in trimmed get IDed as outlier - return to user with warning
    if(all(workingDataFrame$OutlierFlag[workingDataFrame$FSTNoCorr<LowTrimPoint])){
      writeLines("All loci with Fst below the lower (lefthand) trim point were marked as outliers. Re-run with larger LeftTrimFraction or smaller qthreshold.")
        #MODIFIED
      # return(0)
      return(NA)  
    }
    if(all(workingDataFrame$OutlierFlag[workingDataFrame$FSTNoCorr>HighTrimPoint])){
      writeLines("All loci with Fst above the upper (righthand) trim point were marked as outliers. Re-run with smaller RightTrimFraction or smaller qthreshold.")
      #MODIFIED
      # return(0)
      return(NA)

    }
    oldOutlierFlag=workingDataFrame$OutlierFlag
    #writeLines(paste(as.character(count),"   ",as.character(sum(putativeNeutralListTemp))))
  }
  if(count>19) writeLines("Loop iteration limit exceeded.")
  numberLowFstOutliers=sum(workingDataFrame$OutlierFlag[(workingDataFrame$FSTNoCorr<LowTrimPoint)])
  numberHighFstOutliers=sum(workingDataFrame$OutlierFlag[(workingDataFrame$FSTNoCorr>HighTrimPoint)])
  FSTbar=fstBarCalculator(workingDataFrame[putativeNeutralListTemp,])  
  #merge NA list back to working list, and sort by original order	
  resultsDataFrame=rbind(workingDataFrame,storedDataFrameNA)
  resultsDataFrame=resultsDataFrame[order(resultsDataFrame$indexOrder),]
  #return new dataframe
  list(FSTbar=FSTbar,FSTNoCorrbar=FstbarNoCorrTemp,dfInferred=dfInferredTemp,numberLowFstOutliers=numberLowFstOutliers,numberHighFstOutliers=numberHighFstOutliers,results=resultsDataFrame)
}

outputDFStarterNoCorr <- function(FstDataFrame,Hmin=0.1) {
  #This will take a given dataframe with $LocusName, $FST,$He, $T1,  $T2, etc. and 
  #    initialize $indexOrder,$GoodH,$OutlierFlag (to 0), and $q (to 1).
  len=length(FstDataFrame$FSTNoCorr)
  indexOrder=seq(1,len)
  GoodH=ifelse(FstDataFrame$He<Hmin,"lowH","goodH")
  OutlierFlag=ifelse(is.na(FstDataFrame$FSTNoCorr),NA,FALSE)
  # MODIFIED
  OutlierFlag_left=ifelse(is.na(FstDataFrame$FSTNoCorr),NA,FALSE)
  qvalues=rep(NA,len)
  # MODIFIED
  qvalues_left=rep(NA,len)
  pvalues=rep(NA,len)
  pvaluesRightTail=rep(NA,len)
  # MODIFiED
   pvaluesLeftTail=rep(NA,len)
  # MODIFiED
  cbind(FstDataFrame, indexOrder, GoodH, qvalues,qvalues_left,pvalues,pvaluesRightTail,pvaluesLeftTail,OutlierFlag,OutlierFlag_left)
}

#' Calculates q-values for test of neutrality for a list of loci, using input of an inferred degrees of freedom for the chi-square and mean Neutral FST
pOutlierFinderChiSqNoCorr <- function(DataList, Fstbar, dfInferred, qthreshold=0.05, Hmin=0.1){
  #Finds outliers based on chi-squared distribution
  #Takes given values of dfInferred and Fstbar, and returns a list of p-values and q-values for all loci based on chi-square.
  #Assumes that the DataList input has a column called $FSTNoCorr and that empty columns exist for $qvalues and $OutlierFlag 
  #Divide DataList into 3 lists:  DataListGood has $FST>0; DataListNeg has cases where $FST <=0; and
  #   DataListNA has cases where $FST is NA.
  #DataListNeg is necessary to keep separate here because these cases do not have meaningful results with the chi-square approach;
  #   however, they do carry information.
  keepers = which((DataList$FSTNoCorr > 0) & (DataList$He >= Hmin))
  DataListGood = DataList[keepers,]
  DataListOthers = DataList[-keepers,]
  numOthers = length(DataListOthers[,1])
  #Putting NAs in the results columns for all loci that don't meet Hmin or positive Fst criteria
  DataListOthers$pvalues = rep(NA,numOthers)
  DataListOthers$pvaluesRightTail = rep(NA,numOthers)
  # MODIFIED
  DataListOthers$pvaluesLeftTail = rep(NA,numOthers)
  DataListOthers$qvalues = rep(NA,numOthers)
  # MODIFIED
  DataListOthers$qvalues_left = rep(NA,numOthers)
  DataListOthers$OutlierFlag = rep(NA,numOthers)
  # MODIFIED
  DataListOthers$OutlierFlag_left = rep(NA,numOthers)
  #Calculating p values and q-values for loci with high enough He and positive Fst
  pList = pTwoSidedFromChiSq(DataListGood$FSTNoCorr*(dfInferred)/Fstbar,dfInferred)
  pListRightTail = 1-pchisq(DataListGood$FSTNoCorr*(dfInferred)/Fstbar,dfInferred)
  # MODIFIED
  pListLeftTail = pchisq(DataListGood$FSTNoCorr*(dfInferred)/Fstbar,dfInferred)
  # MODIFIED, the function crashes if the pListLeftTail does not have a p-value > 0.95
  qtemp <- qvalue(c(pListRightTail,0.95),fdr.level=qthreshold,pi0.method="bootstrap")
  qtemp$qvalues <-  qtemp$qvalues[1:length(pListRightTail)]
  qtemp$significant <- qtemp$significant[1:length(pListRightTail)]
  # MODIFIED, the function crashes if the pListLeftTail does not have a p-value > 0.95
  qtemp_left <- qvalue(c(pListLeftTail,0.95),fdr.level=qthreshold,pi0.method="bootstrap")
  qtemp_left$qvalues <-  qtemp_left$qvalues[1:length(pListLeftTail)]
  qtemp_left$significant <- qtemp_left$significant[1:length(pListLeftTail)]
  #Note:  Using the bootstrap method here seems OK, but if this causes problems remove the pi0.method="bootstrap" in the previous line to revert to the default.
  DataListGood$pvalues = pList
  DataListGood$pvaluesRightTail = pListRightTail
  # MODIFIED
  DataListGood$pvaluesLeftTail = pListLeftTail
  DataListGood$qvalues = qtemp$qvalues
  # MODIFIED
  DataListGood$qvalues_left = qtemp_left$qvalues
  DataListGood$OutlierFlag = qtemp$significant
  # MODIFIED
  DataListGood$OutlierFlag_left = qtemp_left$significant
  #Combining the good and bad loci back and sorting
  resultsDataFrame = rbind(DataListGood,DataListOthers) 
  #resultsDataFrame=resultsDataFrame[order(resultsDataFrame$indexOrder),]
}

#' Calculates q-values for test of neutrality for a list of loci, using input of an inferred degrees of freedom for the chi-square and mean Neutral FST, and returns the results in the same row order as the input
pOutlierFinderInOrder <- function(DataList, Fstbar, dfInferred, qthreshold=0.05, Hmin=0.1){
  #Assign a temporary index to each row
  len = length(DataList$FSTNoCorr)
  indexOrderTEMP = seq(1,len)
  DataListTEMP = cbind(DataList, indexOrderTEMP)
  #Calculate p and q values using pOutlierFinderChiSqNoCorr
  resultsDataFrame = pOutlierFinderChiSqNoCorr(DataListTEMP, Fstbar, dfInferred, qthreshold, Hmin)
  #Sort to index and delete temporary index
  resultsDataFrame=resultsDataFrame[order(resultsDataFrame$indexOrderTEMP),]
  within(resultsDataFrame, rm(indexOrderTEMP))
}

pTwoSidedFromChiSq <- function(x,df){
  #Takes a value x, finds the two-sided p-value for comparison to a chi-square distribution with df degrees of freedom.
  pOneSided=pchisq(x,df)
  ifelse(pOneSided>.5,(1-pOneSided)*2,pOneSided*2)
}

####### Likelihood functions needed by OutFLANK
EffectiveNumberSamplesMLE=function(FstVect, Fstbar, NumberOfSamples, SmallestFstInTrimmedList, LargestFstInTrimmedList){
  #This function should find the maximum likelihood value 
  #of the effective number of samples, for a given list of
  #Fst values.
  #The FstVect should already have been purged of NaN values and of loci with 
  #too low heterozygosity or MAF. 
  sortedFst=FstVect[order(FstVect)]  
  #The Minimum Fst considered in the trimmed data is the larger of the amount
  #specified by the user or the mean FSt over 100. This is to prevent extremely
  #small Fsts from causing estimation errors (Especially when R interprets a
  #small Fst as FSt=0.)
  LowTrimPoint=max(Fstbar/100,SmallestFstInTrimmedList)
  trimmedFstVect =FstVect[which((FstVect>=LowTrimPoint)&(FstVect<=LargestFstInTrimmedList))]
  trimmedFstArray=as.array(trimmedFstVect)
  localNLLAllData=function(dfInferred){
    localNLLOneLocus=function(Fst){
      negLLdfFstTrim(Fst,dfInferred,Fstbar,LowTrimPoint,LargestFstInTrimmedList)
    }
    sum(localNLLOneLocus(trimmedFstVect))
  }
  optim(NumberOfSamples, localNLLAllData, lower=2, method="L-BFGS-B")$par
}

IncompleteGammaFunction=function(a, z) {
  #equivalence to Mathematica Gamma[a,z] according to 
  #   http://r.789695.n4.nabble.com/Incomplete-Gamma-function-td833545.html
  pgamma(z,a,lower=FALSE)*gamma(a)
}

negLLdfFstTrim=function(Fst, dfInferred, Fstbar, LowTrimPoint, HighTrimPoint){
  #Fst is the Fst from a locus, and dfInferred is the candidate value for the
  #degrees of freedom for the chi-squared distribution of neutral Fst, and
  #Fstbar is the mean Fst of all neutral loci (sequentially inferred from
  #non-outlier loci) LowTrimPoint and HighTrimPoint are the values of the lowest
  #and highest Fst values allowed to be included in the Fst list.
  #Finds contribution to the negative log likelihood of a given locus' Fst for a
  #given dfInferred #CHECKED AGAINST MATHEMATICA DERIVATION##
  df=dfInferred
  1/(2*Fstbar)*(df * Fst +df * Fstbar * log(2) - df * Fstbar *log(df)-(df-2)*Fstbar * log(Fst)+df * Fstbar * log(Fstbar) + 2*Fstbar * log(-IncompleteGammaFunction(df/2,df*HighTrimPoint/(2*Fstbar))+IncompleteGammaFunction(df/2,df*LowTrimPoint/(2*Fstbar))))
}

#' Calculates FST both with and without correction for local sample sizes, for diploid biallelic data. Based on Weir and Cockerham (1984)
WC_FST_Diploids_2Alleles<-function(Sample_Mat){
  ##Calculate both Fst and Fst NoCorr at the same time, from WC84
  #Sample Mat has three columns (homo_p,m heterozygotes, and homo_q) and a row for each population
  sample_sizes = rowSums(Sample_Mat)
  n_ave = mean(sample_sizes)
  n_pops = nrow(Sample_Mat) #r
  r = n_pops
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)
  p_freqs = (Sample_Mat[,1] + Sample_Mat[,2]/2) /sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
  if(s2==0){return(0); break}  
  h_freqs = Sample_Mat[,2]/sample_sizes
  h_ave = sum(sample_sizes*h_freqs)/(n_ave*n_pops)
  a <- n_ave/n_c*(s2 - 1/(n_ave-1)*(p_ave*(1-p_ave)-((r-1)/r)*s2-(1/4)*h_ave))
  b <- n_ave/(n_ave-1)*(p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave - 1)/(4*n_ave)*h_ave)
  c <- 1/2*h_ave
  aNoCorr <- n_ave/n_c*(s2)
  bNoCorr <- (p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave)/(4*n_ave)*h_ave)
  cNoCorr <- 1/2*h_ave
  He <- 1-sum(p_ave^2, (1-p_ave)^2)
  FST <- a/(a+b+c) 
  FSTNoCorr = aNoCorr/(aNoCorr+bNoCorr+cNoCorr)
  return(list(He=He,FST=FST, T1=a, T2=(a+b+c),FSTNoCorr=FSTNoCorr, T1NoCorr=aNoCorr, T2NoCorr=(aNoCorr+bNoCorr+cNoCorr),meanAlleleFreq = p_ave))
}

####### Create input file for OutFLANK
MakeDiploidFSTMat = function(SNPmat,locusNames,popNames){
  # SNPmat is a matrix with individuals in rows and snps in columns
  # 0, 1, or 2 represent the number of copies of the focal allele, and 9 is for missing data
  # locusNames is a character vector of names of each SNP
  # popNames is a character vector with the population identifier for each individual 
  locusname <- unlist(locusNames)
  popname <- unlist(popNames)
  ### Check that SNPmat has appropriate values (0, 1, 2, or 9, only)
  snplevs <- levels(as.factor(unlist(SNPmat)))
  ls <- paste(snplevs, collapse="")
  if(ls!="012" & ls!="0129"){print("Error: Your snp matrix does not have 0,1, and 2"); break}
  ### Checking that locusNames and popNames have the same lengths as the columns and rows of SNPmat
  if(dim(SNPmat)[1]!=length(popname) ){
    print("Error: your population names do not match your SNP matrix")
    break}
  if(dim(SNPmat)[2]!=length(locusname)){
    print("Error:  your locus names do not match your SNP matrix")
    break}
  # writeLines("Calculating FSTs, may take a few minutes...")
  nloci <- length(locusname)
  FSTmat <- matrix(NA, nrow=nloci, ncol=8)
  for (i in 1:nloci){
    FSTmat[i,]=unlist(getFSTs_diploids(popname,SNPmat[,i]))
    if (i%%10000==0){print(paste(i, "done of", nloci))}
  }
  outTemp=as.data.frame(FSTmat)
  outTemp = cbind(locusname,outTemp)
  colnames(outTemp)= c("LocusName","He", "FST", "T1", "T2", "FSTNoCorr", "T1NoCorr", "T2NoCorr", "meanAlleleFreq")
  return (outTemp)
}

### Calculates FST etc. from a single locus from a column of individual data
getFSTs_diploids <- function(popNameList, SNPDataColumn){  
  #eliminating the missing data for this locus
  popnames=unlist(as.character(popNameList))
  popNameTemp=popnames[which(SNPDataColumn!=9)]
  snpDataTemp=SNPDataColumn[SNPDataColumn!=9]
  
  HetCounts <- tapply(snpDataTemp, list(popNameTemp,snpDataTemp), length)
  HetCounts[is.na(HetCounts)] = 0
  
  #Case: all individuals are genetically identical at this locus
  if(dim(HetCounts)[2]==1){
    return (list(He=NA,FST=NA, T1=NA, T2=NA,FSTNoCorr=NA, T1NoCorr=NA, T2NoCorr=NA,meanAlleleFreq = NA))
  }
  
  if(dim(HetCounts)[2]==2){
    if(paste(colnames(HetCounts),collapse="")=="01"){HetCounts=cbind(HetCounts,"2"=0)}
    if(paste(colnames(HetCounts),collapse="")=="12"){HetCounts=cbind("0"=0,HetCounts)} 
    if(paste(colnames(HetCounts),collapse="")=="02"){HetCounts=cbind(HetCounts[,1],"1"=0, HetCounts[,2])}
  }
  
  out = WC_FST_Diploids_2Alleles(HetCounts)	
  return(out)
}

fstBarCalculatorNoCorr <- function(DataList){
  #Calculates mean FstNoCorr from the dataframe, using sum(T1NoCorr) / sum(T2NoCorr) as the estimate of mean Fst.
  #Uses only data for which qvalues > qthreshold (i.e. $OutlierFlag==FALSE)
  #Does not internally screen for low MAF or low He values (but that can be added by only sending the
  #  high MAF rows to this function)
  sum(DataList$T1NoCorr[which(!DataList$OutlierFlag)])/sum(DataList$T2NoCorr[which(!DataList$OutlierFlag)])
}

fstBarCalculator <- function(DataList){
  #Calculates mean Fst from the dataframe, using sum(T1) / sum(T2) as the estimate of mean Fst.
  #Uses only data for which qvalues > qthreshold (i.e. $OutlierFlag==FALSE)
  #Does not internally screen for low MAF or low He values (but that can be added by only sending the
  #  high MAF rows to this function)
  sum(DataList$T1[which(!DataList$OutlierFlag)])/sum(DataList$T2[which(!DataList$OutlierFlag)])
}

#' Calculates FST without correction for local sample sizes, for diploid biallelic data. This is necessary for using OutFLANK, which depends on these uncorrected values for reliable function. (Otherwise, sampling corrections can sometimes cause negative estimates of FST.)
WC_FST_FiniteSample_Diploids_2Alleles_NoCorr <- function(Sample_Mat){
  #Sample Mat has three columns (homo_p,m heterozygotes, and homo_q) and a row for each population
  sample_sizes = rowSums(Sample_Mat)
  n_ave = mean(sample_sizes)
  n_pops = nrow(Sample_Mat) #r
  r = n_pops
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)
  p_freqs = (Sample_Mat[,1] + Sample_Mat[,2]/2) /sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
  if(s2==0){return(1); break}  
  h_freqs = Sample_Mat[,2]/sample_sizes
  h_ave = sum(sample_sizes*h_freqs)/(n_ave*n_pops)
  a <- n_ave/n_c*(s2)
  b <- (p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave)/(4*n_ave)*h_ave)
  c <- 1/2*h_ave
  He <- 1-sum(p_ave^2, (1-p_ave)^2)
  FST <- a/(a+b+c) 
  return(list(He=He,FSTNoCorr=FST, T1NoCorr=a, T2NoCorr=(a+b+c)))
}

#' Calculates FST with correction for local sample sizes, for diploid biallelic data. 
WC_FST_FiniteSample_Diploids_2Alleles<-function(Sample_Mat){
  #Sample Mat has three columns (homo_p,m heterozygotes, and homo_q) and a row for each population
  sample_sizes = rowSums(Sample_Mat)
  n_ave = mean(sample_sizes)
  n_pops = nrow(Sample_Mat) #r
  r = n_pops
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)
  p_freqs = (Sample_Mat[,1] + Sample_Mat[,2]/2) /sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
  if(s2==0){return(1); break}	
  h_freqs = Sample_Mat[,2]/sample_sizes
  h_ave = sum(sample_sizes*h_freqs)/(n_ave*n_pops)
  a <- n_ave/n_c*(s2 - 1/(n_ave-1)*(p_ave*(1-p_ave)-((r-1)/r)*s2-(1/4)*h_ave))
  b <- n_ave/(n_ave-1)*(p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave - 1)/(4*n_ave)*h_ave)
  c <- 1/2*h_ave
  He <- 1-sum(p_ave^2, (1-p_ave)^2)
  FST <- a/(a+b+c) 
  return(list(He=He,FST=FST, T1=a, T2=(a+b+c)))
}

# this function is modified from the function wc from the package hierfstat
# this function calculates the same statistics as the original function of
# outflank, but this function calculates the statistics also for multiallelic 
# data. the variance components a, b and c in outflank are calculated for one 
# allele. in a biallelic locus, if a, b and c are calculated for each allele 
# and then summed (as recommended by Weir and Cockerham 1984 page 1363), the FST
# calculated is the same in either case (using one allele or both of them).
wc_2 <- function (ndat,diploid=TRUE) {
    pop <- ndat[, 1]
    ni <- length(pop)
    dat <- ndat
    loc.names <- names(dat)[-1]
    n <- t(ind.count(dat))
    nt <- apply(n, 1, sum, na.rm = TRUE)
    # if all the alleles are monomorphic return NA
    dum_2 <- allele.count(df_allele.count=cbind(rep(1, ni), dat[, -1]))
    if(length(dum_2)==1){return(NA)}
    if(length(dum_2)>1){
    alploc <- nb.alleles(df_nb.alleles=cbind(rep(1, ni), dat[, -1]))
    np <- dim(n)[2]
    npl <- apply(n, 1, tempfun <- function(x) sum(!is.na(x)))
    nl <- dim(n)[1]
    p <- pop.freq(dat, diploid)
    pb <- pop.freq(cbind(rep(1, length(pop)), dat[, -1]), diploid)
    n <- matrix(unlist(n), ncol = np)
    nal <- n[rep(1:nl, alploc), ]
    nc <- (nt - apply(n^2, 1, sum, na.rm = TRUE)/nt)/(npl - 1)
    ntal <- rep(nt, alploc)
    ncal <- rep(nc, alploc)
    p <- matrix(unlist(lapply(p, t)), ncol = np, byrow = TRUE)
    pb <- matrix(unlist(pb), ncol = 1)
    dum <- getal.b(dat[, -1])
    all.loc <- apply(dum, 2, function(y) as.numeric(dimnames(table(y))[[1]]))
    hetpl <- apply(dum, 2, function(z) {
            lapply(as.numeric(dimnames(table(z))[[1]]), function(y) apply(z == 
                y, 1, function(x) xor(x[1], x[2])))
        })
    mho <- lapply(hetpl, function(x) matrix(unlist(lapply(x, 
            function(y) tapply(y, pop, sum, na.rm = TRUE))), 
            ncol = np))
    mho <- matrix(unlist(mho), ncol = np, byrow = TRUE)
    mhom <- (2 * nal * p - mho)/2
    SSG <- apply(nal * p - mhom, 1, sum, na.rm = TRUE)
    dum <- nal * (p - 2 * p^2) + mhom
    SSi <- rowSums(dum, na.rm = TRUE)
    dum1 <- nal * (sweep(p, 1, pb))^2
    SSP <- 2 * rowSums(dum1, na.rm = TRUE)
    ntalb <- rep(npl, alploc)
    MSG <- SSG/ntal
    MSP <- SSP/(ntalb  - 1)
    Fxy <- function(x) x[1]/sum(x, na.rm = TRUE)
    loc <- rep(1:nl, alploc)

    MSI <- SSi /(ntal - ntalb)
    sigw <- MSG
    sigb <- 0.5 * (MSI - MSG)
    siga <- 1/2/ncal * (MSP - MSI)
    FST.pal <- apply(cbind(siga, sigb, sigw), 1, Fxy)
    lsiga <- tapply(siga, loc, sum, na.rm = TRUE)
    lsigb <- tapply(sigb, loc, sum, na.rm = TRUE)
    lsigw <- tapply(sigw, loc, sum, na.rm = TRUE)
    lFST <- apply(cbind(lsiga, lsigb, lsigw), 1, Fxy)
    
    MSINoCorr <- SSi /(ntal )
    sigwNoCorr <- MSG
    sigbNoCorr <- 0.5 * (MSINoCorr - MSG)
    sigaNoCorr <- 1/2/ncal  * MSP  
    FST.palNoCorr <- apply(cbind(sigaNoCorr, sigbNoCorr, sigwNoCorr), 1, Fxy)
    lsigaNoCorr <- tapply(sigaNoCorr, loc, sum, na.rm = TRUE)
    lsigbNoCorr <- tapply(sigbNoCorr, loc, sum, na.rm = TRUE)
    lsigwNoCorr <- tapply(sigwNoCorr, loc, sum, na.rm = TRUE)
    lFSTNoCorr <- apply(cbind(lsigaNoCorr, lsigbNoCorr, lsigwNoCorr), 1, Fxy)
    
    # this calculates the pooled expected heterozygosity
    lhet <- tapply(pb, loc, function(x){1-(sum(x^2, na.rm = TRUE))})
    # this is the mean allele frequency of the first allele, for both biallelic 
    # and multiallelic loci. However, the allele frequency is not used by outflank
    lmean_freq <- tapply(pb, loc,"[",1)
    
    res <- data.frame(LocusName=loc.names,
                 He=lhet,
                 FST= lFST,
                 T1= lsiga,
                 T2= lsiga+lsigb+lsigw ,
                 FSTNoCorr = lFSTNoCorr,
                 T1NoCorr = lsigaNoCorr,
                 T2NoCorr = lsigaNoCorr+lsigbNoCorr+lsigwNoCorr,
                 meanAlleleFreq=lmean_freq)
    res <- res[complete.cases(res$FSTNoCorr),]
    res[res$FSTNoCorr==0,2:9]<- 0
      return(res)
}
}

####  Plotting functions fr Fst distributions after OutFLANK
##### OutFLANKResultsPlotter #### 


#'This function takes the output of OutFLANK as
#'input with the OFoutput parameter.  It plots a histogram of the FST (by
#'default, the uncorrected FSTs used by OutFLANK) of loci and overlays the
#'inferred null histogram.
#'
#'#'@title Plot the distributions of Fst from OutFLANK output
#'
#'
#'@param OFoutput The output of the function OutFLANK() 
#
#' @param withOutliers Determines whether the loci marked as outliers (with $OutlierFlag) are included in the histogram.
#' 
#' @param NoCorr Plots the distribution of FSTNoCorr when TRUE. Recommended, because this is the data used by OutFLANK to infer the distribution.
#' 
#' @param Hmin The minimum heterozygosity required before including a locus in the plot.
#' 
#' @param binwidth The width of bins in the histogram.
#' 
#' @param Zoom If Zoom is set to TRUE, then the graph will zoom in on the right tail of the distribution (based on argument RightZoomFraction)
#' 
#' @param RightZoomFraction Used when Zoom = TRUE. Defines the proportion of the distribution to plot.
#' 
#' @param titletext Allows a test string to be printed as a title on the graph
#' 
#' @return
#' 
#' The function returns a plot, containing a histogram of Fst with the inferred neutral distribution superimposed.
#' 
#' See the read me file at github.
#'@export

OutFLANKResultsPlotter = function(OFoutput,withOutliers = TRUE, NoCorr= TRUE, Hmin=0.1, binwidth=0.005, Zoom = FALSE,RightZoomFraction = 0.05,titletext=NULL){
  # MODIFIED
  data=OFoutput[which(OFoutput$He>Hmin),]
  if(NoCorr) {
    flist=data$FSTNoCorr
    fbar=sum(data$T1NoCorr)/sum(data$T2NoCorr)
    titletext= paste(c(titletext,"Fst without sample size correction"))
  }
  
  if(!NoCorr) {
    flist=data$FST
    fbar=OFoutput$FSTbar
    
    titletext= paste(c(titletext,"Fst with sample size correction"))
  }
  
  flist = flist[which(!is.na(flist))]
  keeperlist=which(!data$OutlierFlag)
  
  
  if(!withOutliers) flist = flist[keeperlist]
  
  if(Zoom) {FstDistPlotterZoom(df = OFoutput$dfInferred, FSTlist  = flist,  FSTbar = fbar, binwidth,titletext, RightZoomFraction)} else {
    # MODIFIED
    FstDistPlotter(df = 2, FSTlist = flist,  FSTbar = fbar, binwidth, titletext = titletext)}
  
}


################################

FstDistPlotter = function(df, FSTlist, FSTbar, binwidth=0.005,titletext=NULL){
  xPlotUpperBound=ceiling(max(FSTlist)*100)/100
  breakslist=seq(0,xPlotUpperBound+binwidth,by=binwidth)
  breaks = length(breakslist)
  
  x = breakslist
  y=rep(0,length(x))
  for(i in 1:breaks) y[i] = pchisq(((i-.5)*binwidth)/FSTbar*df , df=df) - pchisq((((i-1.5)*binwidth))/FSTbar*df , df=df)
  y=length(FSTlist)*y
  
  hist(FSTlist,col="darkgoldenrod1", breaks=breakslist, prob=F, xlab="Fst",  main=titletext)

  lines(x,y,col="darkblue", lwd=3)
}

###  FstDistPlotterZoom  #####
#This is a function that plots the right tail of the distribution.  

FstDistPlotterZoom = function(df, FSTlist,  FSTbar, binwidth = 0.005,titletext = NULL, RightZoomFraction = 0.1){
  
  FSTlistNoNA=FSTlist[which(!is.na(FSTlist))]
  
  xPlotUpperBound=ceiling(max(FSTlistNoNA)*100)/100
  xPlotLowerBound=floor(as.numeric(quantile(FSTlistNoNA, prob = 1 - RightZoomFraction, na.rm=TRUE)) * 100) / 100
  flist=FSTlistNoNA[which(FSTlistNoNA>xPlotLowerBound)]
  
  
  breakslist=seq(xPlotLowerBound,xPlotUpperBound,by=binwidth)
  breaks = length(breakslist)
  
  x = breakslist
  
  y=rep(0,length(x))
  for(i in 1:breaks) y[i] = pchisq((xPlotLowerBound + (i-.5)*binwidth)/FSTbar*df , df=df) - pchisq((xPlotLowerBound+ (i-1.5)*binwidth)/FSTbar*df , df=df)
  
  y=length(FSTlistNoNA)*y
  
  hist(flist,col="darkgoldenrod1", breaks=breakslist, prob=F, xlab="Fst",  main=titletext)
 
  lines(x,y,col="darkblue", lwd=3)
 
}

FstDistPlotterAddBadCurve = function(df, FSTlist,  FSTbar, binwidth = 0.005, RightZoomFraction = 0.99){
  
  FSTlistNoNA=FSTlist[which(!is.na(FSTlist))]
  
  xPlotUpperBound=ceiling(max(FSTlistNoNA)*100)/100
  xPlotLowerBound=floor(as.numeric(quantile(FSTlistNoNA, prob = 1 - RightZoomFraction, na.rm=TRUE)) * 100) / 100
  
  
  breakslist=seq(xPlotLowerBound,xPlotUpperBound,by=binwidth)
  breaks = length(breakslist)
  
  x = breakslist
  
  y=rep(0,length(x))
  for(i in 1:breaks) y[i] = pchisq((xPlotLowerBound + (i-.5)*binwidth)/FSTbar*df , df=df) - pchisq((xPlotLowerBound+ (i-1.5)*binwidth)/FSTbar*df , df=df)
  
  y=length(FSTlistNoNA)*y
  
  
  lines(x,y,col="red", lwd=3)
  
}


#  OutFLANKBadCurvePlotter draws a curve based on the same Fstbar but with some differnet degrees of freedom
OutFLANKBadCurvePlotter = function(badDF,OFoutput,withOutliers = TRUE, NoCorr= TRUE, Hmin=0.1, binwidth=0.005, Zoom = FALSE,RightZoomFraction = 0.99,titletext=NULL){
  data=OFoutput$results[which(OFoutput$results$He>Hmin),]
  if(NoCorr) {
    flist=data$FSTNoCorr
    fbar=sum(data$T1NoCorr)/sum(data$T2NoCorr)
  }
  
  if(!NoCorr) {
    flist=data$FST
    fbar=OFoutput$FSTbar
    
  }
  
  flist = flist[which(!is.na(flist))]
  keeperlist=which(!data$OutlierFlag)
  
  
  if(!withOutliers) flist = flist[keeperlist]
  
  FstDistPlotterAddBadCurve(badDF, FSTlist  = flist,  FSTbar = fbar, binwidth,RightZoomFraction)
  
}
