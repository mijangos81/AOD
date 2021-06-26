library(raster)
library(rasterVis)
library(viridis)
library(rgeos)
library(DescTools)
library(maptools)
library(zoo)
library(patchwork)
library(parallel)
library(ggnewscale)
source('gHap.R')
source('rotate_matrix_2.R')

chromosome_length <- 23100000

# folders_hap <- dir(path=getwd(),pattern="^SIM")
# path.folder_sim <- folders_hap[1] 
files_LD_pop1 <- dir(path.folder_sim,pattern = "ld_pop1")
files_snps_pop1 <- dir(path.folder_sim,pattern = "snps_pop1")

snp_final <- read.table(paste0(path.folder_sim,"/",files_snps_pop1[1]),header=T,sep = ",",skip = variables_number)
#this is the dataframe with the SNPs values 
#minor <- 0.05 
# Filter on MAF
use <- with(snp_final, (!is.na(z.HWE_pop1) 
                        # & MAF_pop1 >= minor
                        )) 
# Remove NA's
use[is.na(use)] <- FALSE 
# Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
snp_final <- snp_final[use,]

ld_matrix_2 <- read.table(paste0(path.folder_sim,"/",files_LD_pop1[1]),sep = ",",skip = variables_number,row.names = 1)
colnames(ld_matrix_2) <- rownames(ld_matrix_2)
rownames(ld_matrix_2) <- 1:nrow(ld_matrix_2)
colnames(ld_matrix_2) <- 1:ncol(ld_matrix_2)
ld_columns_2 <- as.data.frame(as.table(as.matrix(ld_matrix_2)))
ld_columns_2 <- ld_columns_2[-ld_columns_2$Freq < 0,] #remove cases where LD was not calculated
ld_columns_2$Var1 <- as.numeric(as.character(ld_columns_2$Var1))
ld_columns_2$Var2 <- as.numeric(as.character(ld_columns_2$Var2))
ld_columns_2 <- ld_columns_2[complete.cases(ld_columns_2),]
raster_haplo <- rasterFromXYZ(ld_columns_2)
polygon_haplo <- rasterToPolygons(raster_haplo, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=T)
polygon_haplo <- elide(polygon_haplo,rotate=45)
polygon_haplo$id <- rownames(as.data.frame(polygon_haplo))
polygon_haplo.pts <- fortify(polygon_haplo, polygon_haplo="id") #this only has the coordinates
polygon_haplo.df <- merge(polygon_haplo.pts, polygon_haplo, by="id", type='left') # add the attributes back 
width_poly <-  round(max(polygon_haplo.df$long),0)
height_poly <-  max(polygon_haplo.df$lat)

reduce_factor <- nrow(snp_final)/width_poly
enlarge_factor <- width_poly/nrow(snp_final)
#this is to correct the cell size unit when the polygon figure is more or less than 1000 units 
correction_factor <- (width_poly/1000)

ld_matrix_3 <- read.table(paste0(path.folder_sim,"/",files_LD_pop1[1]),sep = ",",skip = variables_number,row.names = 1)
ld_matrix_3 <- as.matrix(ld_matrix_3)
dimnames(ld_matrix_3)<-NULL 
matrix_rotate <- rotate.matrix_2(x=ld_matrix_3, angle=-45, method="simple")
matrix_rotate_2 <- matrix_rotate
matrix_rotate_2[matrix_rotate_2==0] <- NA
matrix_rotate_3 <- t(na.locf(t(matrix_rotate_2), fromLast = F, na.rm = T))

means_col <- apply(matrix_rotate_3,2,var,na.rm = T)
means_col[means_col==0] <- NA
means_col <- na.locf(means_col,fromLast =T)
means_col <- means_col * 2
df_col <- as.data.frame(matrix(ncol = 2 ,nrow = length(means_col))) 
df_col[,1] <- 1:nrow(df_col)
df_col[,2] <- means_col
df_col[,2] <-  rescale(df_col[,2],to = c(0,height_poly))
means_col_2 <- means_col
# means_col_2[means_col_order[head(means_col_order,2)]] <- mean(means_col,na.rm=T)
# means_col_2[means_col_order[tail(means_col_order,2)]] <- mean(means_col,na.rm=T)
mean_column <- as.numeric(summary(means_col_2,na.rm=T)[1:6])
mean_column <- rescale(mean_column,to = c(0,height_poly))

# as the matrix was rotated the position of the snps is not correct anymore. the following code reassigns the snp position based on the second row from bottom to top of the rotated matrix. the first two and the last snps are removed to take in account the snps that are not present in the second row
 second_row <- which(!is.na(matrix_rotate_3[,1])) -1
 second_row_2 <- matrix_rotate_2[second_row,]
 second_row_3 <- second_row_2
 reassign_loc <- snp_final$loc_bp[3:nrow(snp_final)]
 reassign_loc <- reassign_loc[1:length(reassign_loc)-1]
 
 element_reassign_loc <- 1
 for(i in 1:length(second_row_2)){
   if(is.na(second_row_2[i])){
     next
   }else{
     second_row_3[i] <- reassign_loc[element_reassign_loc]
      element_reassign_loc <-  element_reassign_loc + 1
   }
 }
 
 # putting back the first two snps and the last that were removed
 second_row_3[1:2] <- snp_final$loc_bp[1:2]
 second_row_3[length(second_row_3)] <- snp_final$loc_bp[length(snp_final$loc_bp)]
 #filling the NAs
 second_row_4 <- na.locf(second_row_3,fromLast=T)
 
 ld_threshold <- 0.5
  second_row_ver_2 <- na.locf(second_row_2,fromLast=T)
# 
 first_row <- which(!is.na(matrix_rotate_3[,1])) 
  first_row_2 <- matrix_rotate_3[first_row,]
 haplo_loc_test <- first_row_2 >= ld_threshold
 haplo_loc_test <- c(haplo_loc_test,F)
start_haplo <- NULL
end_haplo <- NULL
 for(i in 1:(length(haplo_loc_test)-1)){
   if( haplo_loc_test[i] == haplo_loc_test[i+1] ){
     next()
   }
   if( haplo_loc_test[i] == T ){
     end_haplo_temp <- i
     end_haplo <- c(end_haplo,end_haplo_temp )
   }
    if( haplo_loc_test[i] == F ){
     start_haplo_temp <- i
     start_haplo <- c(start_haplo,start_haplo_temp )
   }
 }

start_haplo <- c(1,start_haplo)
start_haplo_2 <-  unique(second_row_4[start_haplo])
# end_haplo <- c(end_haplo,length(first_row_2))
end_haplo_2 <-  unique(second_row_4[end_haplo])
haplo_1_ver_2 <- as.data.frame(cbind(start_haplo_2,end_haplo_2))


# haplo_loc_2_ver_2 <- unique(second_row_4[list_haplo])
# haplo_1_ver_2 <- as.data.frame(cbind(c(1,haplo_loc_2_ver_2),c(haplo_loc_2_ver_2,chromosome_length)))
# haplo_1_ver_2 <- haplo_1_ver_2[-1,]
# haplo_1_ver_2 <-  haplo_1_ver_2[seq(1,nrow(haplo_1_ver_2),2),]

haplo_1_ver_2$size <- (haplo_1_ver_2[,2]-haplo_1_ver_2[,1])
# haplo_1_ver_2 <- haplo_1_ver_2[which(haplo_1_ver_2$size>=100000),]
# mean_dis_2_ver_2 <- mean(diff(snp_final$loc_bp))
# haplo_1_ver_2$large <- haplo_1_ver_2$size > mean_dis_2_ver_2
# haplo_1_ver_2[1,4] <- TRUE
# haplo_1_ver_2 <- haplo_1_ver_2[which(haplo_1_ver_2$large==TRUE),]
# haplo_1 <- haplo_1_ver_2

# var_threshold <-  as.numeric(summary(means_col_2,na.rm=T)[2])
# haplo_loc <- which(means_col<var_threshold)
# haplo_loc_2 <- unique(second_row_4[haplo_loc])
# haplo_1 <- as.data.frame(cbind(c(1,haplo_loc_2),c(haplo_loc_2,chromosome_length)))
# haplo_1$size <- (haplo_1[,2]-haplo_1[,1])
# mean_dis_2 <- mean(diff(snp_final$loc_bp))
# haplo_1$large <- haplo_1$size > mean_dis_2
# haplo_1[1,4] <- TRUE

# counter <- 0
# new_haplo <- as.data.frame(matrix(ncol = 4))
# for(i in 1:nrow(haplo_1)){
#   if(haplo_1[i,4]==T){
#     counter <- counter + 1
#    new_haplo[counter,] <- haplo_1[i,]
#   }else{
#       new_haplo[counter,2] <- haplo_1[i,2]  
#     }
# }
# new_haplo <-  haplo_1 
# haplo_2 <- new_haplo[,1:3]
# haplo_2[,3] <- haplo_2[,2]- haplo_2[,1]
# # adding one to each number of one of the columns so they are different to be used as breaks
# haplo_2[,1] <- haplo_2[,1]+1
n_snps <- as.matrix(haplo_1_ver_2[,1:2])
n_snps[,1] <- n_snps[,1]+1
n_snps <- n_snps[which(n_snps[,1]!=n_snps[,2]),]
n_snps <- n_snps[1:(nrow(n_snps)-1),]
  
df.4.cut <- as.data.frame(table(cut(snp_final$loc_bp, breaks=n_snps)),stringsAsFactors=F)
# df.4.cut <-  df.4.cut[seq(1,nrow(df.4.cut),2),]
df.4.cut <- df.4.cut[which(df.4.cut$Freq>4),]
df.4.cut_3 <- gsub("[][()]", "", df.4.cut$Var1,",")
df.4.cut_3 <- strsplit(df.4.cut_3,",")
df.4.cut_4 <- lapply(df.4.cut_3, as.numeric)
df.4.cut_4 <- as.data.frame(laply(df.4.cut_4,rbind))
df.4.cut_4[,3] <- (df.4.cut_4[,2]-df.4.cut_4[,1])

# this is to calculate the real distance in bp of the polygon figure of LD
# NAs are replaced by 1 to not affect the real distance in bp
real_distance <- c(0,second_row_4)
real_distance_2 <- diff(real_distance)
real_distance_3 <- cumsum(real_distance_2)
real_distance_4 <- as.data.frame(cbind(1:length(real_distance_3),real_distance_3))

test_var <- unname(unlist(df.4.cut_4[,1:2]))

location_test <- lapply(test_var, findInterval,vec = as.numeric(paste(unlist(real_distance_4$real_distance_3))))
location_test_2 <-  unlist(location_test)
test_var_2 <- as.data.frame(cbind(1:length(location_test_2),location_test_2))

hap_blocks <- as.data.frame(cbind(paste0("CHR_",1:nrow(df.4.cut_4)),
                                  rep(1,times=nrow(df.4.cut_4)),
                                  df.4.cut_4[,1:3],
                                  df.4.cut$Freq))
colnames(hap_blocks) <- c("BLOCK", "CHR" ,  "BP1" ,  "BP2"   ,"SIZE",  "NSNP" )

offspring_ld_pop1 <- pop1
# offspring_ld_pop2 <- pop2

 offspring_ld_pop1$V1[offspring_ld_pop1$V1=="Male"]   <- 1
  offspring_ld_pop1$V1[offspring_ld_pop1$V1=="Female"] <- 2
  offspring_ld_pop1$V2 <- paste0("offspring_","pop1")
  offspring_ld_pop1$id <- paste0("pop1_",row.names(offspring_ld_pop1))
  
  # offspring_ld_pop2$V1[offspring_ld_pop2$V1=="Male"]   <- 1
  # offspring_ld_pop2$V1[offspring_ld_pop2$V1=="Female"] <- 2
  # offspring_ld_pop2$V2 <- paste0("offspring_","pop2")
  # offspring_ld_pop2$id <- paste0("pop2_",row.names(offspring_ld_pop2))
 # samples_haplotypes <-  rbind(offspring_ld_pop1,offspring_ld_pop2)
samples_haplotypes <-  offspring_ld_pop1
samples_haplo <- samples_haplotypes[,c("V2","id")]
chromosomes_pop <- c(samples_haplotypes$V3,samples_haplotypes$V4)

chromosomes_pop <- gsub("[-^1-9]", "o", chromosomes_pop)
chromosomes_pop <- gsub("a", "0", chromosomes_pop) # converting allele names to numbers
chromosomes_pop <- gsub("A", "1", chromosomes_pop)
find_haplo <- strsplit(chromosomes_pop,split = "")
find_haplo_2 <- data.frame((sapply(find_haplo,c)), stringsAsFactors = F)
find_haplo_2 <- find_haplo_2[snp_final$row_name,]
ref <- "A"
alt <- "a"
info_haplo <- as.data.frame(cbind(1,1:(nrow(recombination_map)-1),recombination_map$loc_bp[1:(nrow(recombination_map)-1)],ref,alt), stringsAsFactors = F)
info_haplo <- info_haplo[snp_final$row_name,]
colnames(info_haplo) <- c("chr","marker","pos","ref","alt")
find_haplo_3 <- cbind(info_haplo,find_haplo_2) 
snp_final_2 <- read.table(paste0(path.folder_sim,"/",files_snps_pop1[1]),header=T,sep = ",",skip = variables_number)

find_haplo_3$HW <- snp_final_2$z.HWE_pop1
find_haplo_3 <- find_haplo_3[complete.cases(find_haplo_3),]
find_haplo_3 <- find_haplo_3[,-ncol(find_haplo_3)]
markers_haplo <-  find_haplo_3[,1:5]
markers_haplo$pos <- as.numeric(as.character(markers_haplo$pos))

phase_haplo <- find_haplo_3[,6:ncol(find_haplo_3)]

phase_haplo_2 <- ghap.loadphase(
  samples.file=samples_haplo,
  markers.file=markers_haplo,
  phase.file=phase_haplo,
  verbose = TRUE
)

# blocks.mkr <- ghap.blockgen(phase_haplo_2, windowsize = 10, slide = 1, unit = "marker",nsnp = 4)

ghap.haplotyping(phase=phase_haplo_2, blocks=hap_blocks, outfile="res_haplo", freq = c(0,1), drop.minor = F, batchsize = 500, ncores = 1, verbose = TRUE)

haplo <- ghap.loadhaplo("res_haplo.hapsamples", "res_haplo.hapalleles", "res_haplo.hapgenotypes")
files_to_delete <- c("res_haplo.hapsamples", "res_haplo.hapalleles", "res_haplo.hapgenotypes")
files_to_delete <- paste0(getwd(),"/",files_to_delete)
unlink(files_to_delete)
hapstats <- ghap.hapstats(haplo, ncores = 2)
blockstats <- ghap.blockstats(hapstats, ncores = 2)

location_hap <- lapply(hapstats$BP1, findInterval, 
                   vec = as.numeric(paste(unlist(real_distance_4$real_distance_3))))
hapstats_2 <- hapstats[,c("BLOCK","BP1","ALLELE","FREQ")]
hapstats_2$BP1 <- c(location_hap)
hapstats_2$BLOCK <- as.factor(hapstats_2$BLOCK)

split_haplo <- split(hapstats_2,hapstats_2$BLOCK)
split_haplo <- lapply(split_haplo,function(x){x[order(x$FREQ,decreasing = T),]})
split_haplo2 <- lapply(split_haplo,"[",,c("BP1","ALLELE"))
position_hap <- lapply(split_haplo2,"[[",1,1)
position_hap <- unlist(position_hap)
position_hap <- position_hap 
# + correction_factor
# position_hap <-  floor(rescale(position_hap,to=c(1,floor(width_poly*0.98))))

split_haplo3 <- lapply(split_haplo2,function(x){strsplit(x$ALLELE,split = "") } )
split_haplo4 <- lapply(split_haplo3,sapply,rbind)
split_haplo5 <- lapply(split_haplo4,t)
for(i in 1:length(position_hap)){
  col_num <- ncol(split_haplo5[[i]]) -1
  row_num <- nrow(split_haplo5[[i]])
  col_names_haplo <- seq( position_hap[i] , position_hap[i]+(enlarge_factor*col_num) ,enlarge_factor)
    # position_hap[i] : (position_hap[i] + col_num)
  # col_names_haplo <- col_names_haplo * correction_factor
  colnames(split_haplo5[[i]]) <- col_names_haplo
  rownames(split_haplo5[[i]]) <-  seq( 1 , (enlarge_factor*row_num) ,enlarge_factor)
}
split_haplo6 <- lapply(split_haplo5,function(x){as.data.frame(as.table(x))})

plot_haplo <- rbindlist(split_haplo6)
plot_haplo$Var1 <- as.numeric(as.character(plot_haplo$Var1))
plot_haplo$Var1 <- 1-(plot_haplo$Var1)
plot_haplo$Var1 <- (plot_haplo$Var1)-2
plot_haplo$Var2 <- as.numeric(as.character(plot_haplo$Var2))
plot_haplo$Freq <- as.factor(plot_haplo$Freq)
# plot_haplo$Freq <- plot_haplo$Freq=="A"
# plot_haplo$Freq <- as.numeric(plot_haplo$Freq)


#to test selection coefficients and dominance per locus
HW_snp <- as.data.frame(snp_final$z.HWE_pop1)
HW_snp <- HW_snp + (-1*min(HW_snp))
HW_snp$loc_bp <- 1:nrow(HW_snp) 
HW_snp[,1] <- rescale(HW_snp[,1],to = c(0,width_poly))
HW_snp[,2] <-  rescale(HW_snp[,2],to = c(0,height_poly))

#to test selection coefficients and dominance per locus
test_sel <- snp_final
test_sel$hs <- test_sel$h * test_sel$s
test_sel <- test_sel[,c("loc_bp","hs")]
test_sel$loc_bp <- 1:nrow(test_sel)
test_sel[,1] <- rescale(test_sel[,1],to = c(0,width_poly))
test_sel[,2] <-  rescale(test_sel[,2],to = c(0,height_poly))

#this is Fst rescaled to 0 to 1 for the mixed statistic for SNP's
snp_fst <- snp_final[,c("loc_bp","Fst")]
snp_fst$loc_bp <- 1:nrow(snp_fst)
snp_fst[,1] <- rescale(snp_fst[,1],to = c(0,width_poly))
snp_fst[,2] <- snp_fst[,2] * height_poly
  # rescale(snp_fst[,2],to = c(0,height_poly))

#this is SNP's heterozygosity to be calculated alone
snp_het_alone <- snp_final[,c("loc_bp","P.AB_pop1")]
snp_het_alone$loc_bp <- 1:nrow(snp_het_alone)
snp_het_alone[,1] <- rescale(snp_het_alone[,1],to = c(0,width_poly))
snp_het_alone[,2] <- snp_het_alone[,2] * height_poly
  # rescale(snp_het_alone[,2],to = c(0,height_poly))

fst_msats <- as.data.frame(cbind(location_neutral_loci_analysis,unlist(mean_Fst[number_of_generations,])))
colnames(fst_msats) <- c("V1","V2")
location_msats_ld <- unlist(lapply(location_neutral_loci_analysis, findInterval,vec = as.numeric(paste(unlist(real_distance_4$real_distance_3)))))
fst_msats$V1 <- location_msats_ld
fst_msats <- fst_msats[complete.cases(fst_msats$V2),]
fst_msats <- unique(fst_msats)
fst_msats$V2 <- fst_msats$V2 * height_poly
  # rescale(fst_msats$V2,to=c(0,height_poly))

#this is msats heterozygosity to be calculated alone
het_msats_alone <- as.data.frame(cbind(location_neutral_loci_analysis,unlist(mean_het_pop1[number_of_generations,])))
colnames(het_msats_alone) <- c("V1","V2")
location_msats_ld <- unlist(lapply(location_neutral_loci_analysis, findInterval,vec = as.numeric(paste(unlist(real_distance_4$real_distance_3)))))
het_msats_alone$V1 <- location_msats_ld
het_msats_alone <- unique(het_msats_alone)
het_msats_alone <- het_msats_alone[complete.cases(het_msats_alone$V2),]
het_msats_alone$V2 <- het_msats_alone$V2 * height_poly
  # rescale(het_msats_alone$V2,to=c(0,height_poly))

exp_het <- as.data.frame(cbind(location_neutral_loci_analysis,unlist(mean_expected_het_pop1[number_of_generations,])))
colnames(exp_het) <- c("V1","V2")
location_msats_ld <- unlist(lapply(location_neutral_loci_analysis, findInterval,vec = as.numeric(paste(unlist(real_distance_4$real_distance_3)))))
exp_het$V1 <- location_msats_ld
exp_het <- unique(exp_het)
exp_het <- exp_het[complete.cases(exp_het$V2),]
exp_het$V2 <- exp_het$V2 * height_poly

recom_haplo <- as.data.frame(cbind((1:nrow(map)*map_resolution),map$cM))
colnames(recom_haplo) <- c("V1","V2")
location_recom <- unlist(lapply(recom_haplo$V1, findInterval,vec = as.numeric(paste(unlist(real_distance_4$real_distance_3)))))
recom_haplo$V1 <- location_recom
recom_haplo$V2 <- rescale(recom_haplo$V2,to=c(0,height_poly))
recom_haplo <- unique(recom_haplo)

ticks_df <- real_distance_4
ticks_df$real_distance_3 <- round(ticks_df$real_distance_3/1000000,1) 
ticks_df_2 <- ticks_df[which(ticks_df$real_distance_3 %% 0.5 == 0),]
colnames(ticks_df_2) <- c("V1","V2")
unique_val <- c(TRUE,diff(ticks_df_2$V2)!=0)
ticks_df_3 <- ticks_df_2[unique_val,]
ticks_lab <- as.character(ticks_df_3$V2)
ticks_breaks <- ticks_df_3$V1

colnames(plot_haplo) <- c("Var1", "Var2" ,"Variant")

cols <- c("SNP's"="lightgoldenrod","Msats observed"="deeppink","Msats expected"="cyan")

p1 <- ggplot() +
 geom_polygon(data=polygon_haplo.df, aes(long,lat,group=group, fill=Freq)) +
  scale_fill_viridis(name="R square",option= "viridis") +
  new_scale("fill") +   
  geom_vline(xintercept = test_var_2$location_test_2,color="red",size=1/2)+
 # geom_rect(data = plot_haplo, aes(xmin=Var2, xmax=Var2+enlarge_factor, ymin=Var1, ymax=Var1+enlarge_factor,fill=Variant),color="black",size=1/10,inherit.aes = F) +
  scale_fill_manual(values = c("red", "blue","black"), labels = c("Deleterious", "Alternative"))+
# geom_line(data=df_col,aes(x=df_col$V1,y=df_col$V2),color="orange",inherit.aes = F,size=1/2) +
# geom_hline(yintercept = mean_column[2],color="red",size=1/2)+

 # geom_step(data=recom_haplo,aes(x=recom_haplo$V1,y=recom_haplo$V2),color="white",inherit.aes = F,size=1/2) +
labs(x=NULL, y="Haplotypes", title="Simulation with a chromosome of 500 cM",subtitle = paste(nrow(snp_final),"SNP's ","mean Fst = ",round(final_res[1,18],2)," mean He = ",round(final_res[1,15],2),"s = ", s_gral,"h = ", h_gral,"q = ",q_gral, " with ",loci_number, " initial loci" ))+ 
# scale_x_continuous(breaks=ticks_breaks, labels=ticks_lab)+
theme_tufte(base_family="Helvetica") +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_text(angle = 0,face = "bold",size = 12,hjust=0,vjust=0.5),
        axis.ticks = element_blank(),
        axis.ticks.y=element_blank() ,
         axis.ticks.x=element_blank() ,
         axis.text.y = element_blank(),
        plot.margin=unit(c(0,1,-1,1), "in"),
        # legend.justification =c(0,0),
         # legend.margin=margin(0,0,0,0, unit = "in"),
        # legend.box.margin=margin(-1,0,0,0, unit = "in"),
         # legend.box.spacing= unit(0, "in"),
        plot.title = element_text(face = "bold",size = 14,hjust=0,vjust=0.5),
        axis.text.x = element_blank()) +
  coord_fixed(ratio = 1/1)

p2 <- ggplot() +
 geom_polygon(data=polygon_haplo.df, aes(long,lat,group=group, fill=Freq)) +
 scale_fill_viridis(name="R square",option= "viridis") +
 new_scale("fill") +   
 # geom_rect(data = plot_haplo, aes(xmin=Var2, xmax=Var2+enlarge_factor,ymin=Var1, ymax=Var1+enlarge_factor,fill=Variant),color="black",size=1/50,inherit.aes = F) +
  scale_fill_manual(values = c("red", "blue","black"), labels = c("Deleterious", "Alternative"))+
  # new_scale("color") +  
  geom_col(data=het_msats_alone,aes(x=V1,y=V2,color="Msats observed"),fill="deeppink",inherit.aes = F,position = position_dodge2(10),size=0.25) +
  geom_col(data = exp_het,aes(x=V1,y=V2,color="Msats expected"),fill="cyan",inherit.aes = F,position = position_dodge2(2),size=0.25)+
  geom_line(data=snp_het_alone,aes(x=loc_bp,y=P.AB_pop1),color= "lightgoldenrodyellow",inherit.aes = F,size=1/4)+
  geom_point(data=snp_het_alone,aes(x=loc_bp,y=P.AB_pop1,color="SNP's"),inherit.aes = F,size =1/4,stroke=0,shape=16)+
    scale_colour_manual(name="Markers",values=cols) +
  labs( x= "Chromosome location (Mbp)", y="He",title="")+ 
    scale_x_continuous(breaks=ticks_breaks, labels=ticks_lab)+

  # scale_x_continuous(breaks=ticks_breaks, labels=ticks_lab)+
  theme_tufte(base_family="Helvetica") +
   theme(legend.position = "left",
        # axis.title.x=element_blank(),
        axis.title.y=element_text(angle = 0,face = "bold",size = 12,hjust=0,vjust=0.5),
        # axis.ticks = element_blank(),
        axis.ticks.y=element_blank() ,
         axis.title.x= element_text(hjust =0.5),
        axis.text.x = element_text(face = "bold",size = 12),
         # axis.ticks.x=element_blank() ,
         axis.text.y = element_blank(),
        plot.margin=unit(c(-1,1,-1,1), "in"),
        legend.justification =c(0,0),
         # legend.margin=margin(0,0,0,0, unit = "in"),
        # legend.box.margin=margin(-1,0,0,0, unit = "in"),
         legend.box.spacing= unit(0, "in"))+
        # plot.title = element_text(face = "bold",size = 14,hjust=0,vjust=0.5),
        # axis.text.x = element_blank()) +
          # element_text(face = "bold",size = 4)) +
guides(fill = guide_legend(override.aes = list(colour = "blue","red","pink")))+
  coord_fixed(ratio = 1/1)

p3 <- ggplot() +
 geom_polygon(data=polygon_haplo.df, aes(long,lat,group=group, fill=Freq)) +
  scale_fill_viridis(name="R square",option= "viridis") +
  new_scale("fill") +   
 # geom_rect(data = plot_haplo, aes(xmin=Var2, xmax=Var2+enlarge_factor, ymin=Var1, ymax=Var1+enlarge_factor,fill=Variant),color="black",size=1/10,inherit.aes = F) +
   scale_fill_manual(values = c("red", "blue","black"), labels = c("Deleterious", "Alternative"))+
  # new_scale("color") +  
    geom_col(data=fst_msats,aes(x=V1,y=V2),color="deeppink",inherit.aes = F,position = position_dodge2(10),size=0.75) +
  geom_line(data=snp_fst, aes(x=loc_bp,y=Fst),color= "lightgoldenrodyellow",inherit.aes = F,size=1/2)+
  geom_point(data=snp_fst,aes(x=loc_bp,y=Fst,color="SNP's"),inherit.aes = F,size =1/2,stroke=0,shape=16)+
  geom_hline(yintercept = 0.6 * height_poly,color="cyan",size=1/2) +
    scale_colour_manual(name="Markers",values=cols) +
  labs(x= "Chromosome location (Mbp)", y="Fst", title="")+ 
  scale_x_continuous(breaks=ticks_breaks, labels=ticks_lab)+
  theme_tufte(base_family="Helvetica") +
  theme(legend.position = "none",
         # axis.title.x=element_blank(),
         axis.title.y=element_text(angle = 0,face = "bold",size = 12,hjust=0.5,vjust=0.5),
        axis.text.y=element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.y=element_blank() ,
       plot.margin=unit(c(-1,1,1,1), "in"),
      axis.title.x= element_text(hjust =0.5),
        axis.text.x = element_text(face = "bold",size = 12)) +
  coord_fixed(ratio = 1/1)

# p4 <- ggplot() +
# geom_polygon(data=polygon_haplo.df, aes(long,lat,group=group, fill=Freq)) +
# scale_fill_viridis(name="R square",option= "viridis") +
# new_scale("fill") +   
# geom_rect(data = plot_haplo, aes(xmin=Var2, xmax=Var2+enlarge_factor, ymin=Var1, ymax=Var1+enlarge_factor,fill=Variant),color="black",size=1/30,inherit.aes = F) +
# scale_fill_manual(values = c("red", "blue"))+
# geom_line(data=HW_snp ,aes(x=HW_snp[,1],y=HW_snp[,2]),color="deeppink",inherit.aes = F,size=1/2) +
# geom_point(data=HW_snp ,aes(x=HW_snp[,1],y=HW_snp[,2]),color="deeppink4",inherit.aes = F,size = 1/2, stroke = 0, shape = 16) +
# labs(x="Chromosome location (Mbp)", y="h*s", title=NULL)+ 
# scale_x_continuous(breaks=ticks_breaks, labels=ticks_lab)+
# theme_tufte(base_family="Helvetica") +
# theme(legend.position = "none",
#          # axis.title.x=element_blank(),
#          axis.title.y=element_text(angle = 0,face = "bold",size = 12,hjust=0,vjust=0.5),
#         axis.text.y=element_blank(),
#         axis.ticks = element_line(colour = "black", size = 0.25),
#         axis.ticks.y=element_blank() ,
#        plot.margin=unit(c(-1,1,1,1), "in"),
#       axis.title.x= element_text(hjust =0.5),
#         axis.text.x = element_text(face = "bold",size = 12)) +
#   coord_fixed(ratio = 1/1)
#  
 print(
     (p1 + p2 + p3 + plot_layout(ncol =1 , heights = c(1,1,1)))
     /
       (p4 + p5 + p6 ) 
 )

 ggsave("figure_5.pdf",  width = 30, height = 30, units = "in", dpi="retina" )
