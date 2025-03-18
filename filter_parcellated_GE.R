library("dplyr")
library("tidyr")
library("readr")
library('reticulate')

#unfiltered parcellated data pulled from Yohan
#can redo everything with sagittal ge data using find/replace feature
#some genes are only analyzed in coronal/sagittal slice (ex: Uba5 sagittal)
coronal_exp_unfiltered <- as.data.frame(read.csv('/data/chamal/projects/natvik/sir_extended/preprocessed/yohan_parcellated_GE/unfiltered/coronal_expression_data.csv'))
sagittal_exp_unfiltered <- as.data.frame(read.csv('/data/chamal/projects/natvik/sir_extended/preprocessed/yohan_parcellated_GE/unfiltered/sagittal_expression_data.csv'))
#Allen region IDs used by Shady/Steph for parcellated SIR code
steph_atrophy_rgn_names <- as.data.frame(read.csv('/data/chamal/projects/natvik/sir_extended/preprocessed/steph_inputs/merged_label_region_names.csv'))
steph_allen_rgn_names <- as.data.frame(read.csv('/data/chamal/projects/natvik/sir_extended/preprocessed/steph_inputs/ABAnewatlas_final_withacros.csv'))

###ISSUE - do not have data for 426 regions!###
#coronal_exp_filtered <- coronal_exp_unfiltered[
#  which(coronal_exp_unfiltered$region_id %in% 
#          c(steph_allen_rgn_names$rightLabel, steph_allen_rgn_names$leftLabel)),]

#length(unique(coronal_exp_filtered$region_id)) #should be 426 - is 132
yohan_allen_rgn_sagittal_unfilt <- gsub(",","",sagittal_exp_unfiltered$region_name)

yohan_allen_rgn_unfilt <- gsub(",","",coronal_exp_unfiltered$region_name)
coronal_exp_filtered <- coronal_exp_unfiltered[
  which(yohan_allen_rgn_unfilt %in% 
         steph_allen_rgn_names$Structure),]


setdiff(unique(yohan_allen_rgn_sagittal_unfilt), unique(yohan_allen_rgn_unfilt))
#same regions in both datasets. Can just use one of these vectors for future reference

length(unique(coronal_exp_filtered$region_id)) #should be 426 - is 208

#See which regions are not being included
yohan_allen_rgn_filt <- coronal_exp_filtered$region_name
setdiff(steph_allen_rgn_names$Structure, yohan_allen_rgn_filt)

#[1] "Median preoptic nucleus" "Nucleus of reunions"    
#[3] "Nucleus raphe magnus"    "Subiculum dorsal part"  
#[5] "Subiculum ventral part"


#see if these regions correspond to other names in Yohan's Allen data
yohan_allen_rgn_full <- sort(unique(yohan_allen_rgn_unfilt))
#MISSPELLINGS: Median --> Medial; reunions --> reuniens; 
#get rid of subiculum. Steph/Shady did this too

grep("raphe", yohan_allen_rgn_full, value = TRUE)

##ISSUE - "Medulla behavioral state related" and its children are not in Yohan's data
grep("pallidus", yohan_allen_rgn_full, value = TRUE)
grep("magnus", yohan_allen_rgn_full, value = TRUE) #in CCFv3 ontology but not in this dataset...
grep("obscurus", yohan_allen_rgn_full, value = TRUE)
grep("behavior", yohan_allen_rgn_full, value = TRUE)

#Change spelling discrepancies to line up with Steph's data 
sagittal_exp_unfiltered$region_name <- gsub("Nucleus of reuniens", "Nucleus of reunions", sagittal_exp_unfiltered$region_name)
sagittal_exp_unfiltered$region_name <- gsub(",","",sagittal_exp_unfiltered$region_name)

#NEW: COMBINE CORONAL AND SAGITTAL FOR OVERALL DATASET TO BE USED
#note: only ~1500 genes out of the 
combined_exp_unfiltered <- rbind(coronal_exp_unfiltered, sagittal_exp_unfiltered)
##################################################################################################

#Redo filtering and check that the only regions that don't match up are 
#replace combined with coronal or sagittal
yohan_allen_rgn_unfilt <- combined_exp_unfiltered$region_name
combined_exp_filtered <- combined_exp_unfiltered[
  which(yohan_allen_rgn_unfilt %in% 
          steph_allen_rgn_names$Structure),]


yohan_allen_rgn_filt <- sort(unique(combined_exp_filtered$region_name))
setdiff(steph_allen_rgn_names$Structure, yohan_allen_rgn_filt)

output_dir <- '/data/chamal/projects/natvik/sir_extended/derivatives/'
write.csv(combined_exp_filtered, file=paste0(output_dir,'combined_expression_rgn_min2.csv'), row.names=FALSE)


###PART 2: Create a parcellated dataset for each gene, making sure to exclude###
###genes with a low proportion of valid voxels across multiple regions##########

#sort by gene
output_dir <- '/data/chamal/projects/natvik/sir_extended/derivatives/'
combined_exp_filtered <- read.csv(paste0(output_dir,'combined_expression_rgn_min2.csv'))
yohan_allen_rgn_filt <- sort(unique(combined_exp_filtered$region_name))

gene_list_combined <- sort(unique(combined_exp_filtered$gene))
valid_voxel_thresholds <- c(0, 0.2, 0.5) 

for(valid_voxel_threshold in valid_voxel_thresholds)
{
  hist(combined_exp_filtered$voxel_valid_mean)
  mean(combined_exp_filtered$voxel_valid_mean < valid_voxel_threshold) #get rid of only 2% of data
  
  #Create an output data frame with regions x genes (209 x 4082)
  rgn_ge_df <- data.frame(matrix(ncol = length(gene_list_combined), nrow = length(yohan_allen_rgn_filt)))
  colnames(rgn_ge_df) <- gene_list_combined
  rownames(rgn_ge_df) <- yohan_allen_rgn_filt
  
  for(gene in gene_list_combined)
  {
    ge_rows <- combined_exp_filtered[which(combined_exp_filtered$gene == gene),]
    #DROP NAs
    ge_rows <- ge_rows[which(!is.na(ge_rows$expression_mean_normalized_brain)),]
    
    #Case 1: GE is invalid for rgn (< 0.5 voxel_valid_mean)
    #Case 2: GE is valid for experiment but is a repeat experiment: add to a LIST
    
    #Filter gene expression rows
    #Case 1
    ge_rows_f0 <- ge_rows[which(ge_rows[,"voxel_valid_mean"] >= valid_voxel_threshold),]
    
    #Case 2: average duplicate regions (for a single gene)
    rgn_mean_ge <- tapply(ge_rows_f0$expression_mean_normalized_brain, ge_rows_f0$region_name, mean)
    rgn_ge_df[names(rgn_mean_ge),gene] <- rgn_mean_ge
    
    
  }
  no_na_rows <- rownames(rgn_ge_df)[which(complete.cases(rgn_ge_df))]
  print(no_na_rows)
  
  write.csv(rgn_ge_df,paste0(output_dir,'combined_rgn_ge_df_full_vv',valid_voxel_threshold,'.csv'))
  
  #Find the number of valid regions per gene
  invalid_region_per_gene <- sapply(rgn_ge_df, function(x) sum(is.na(x)))
  print(invalid_region_per_gene["Snca"]) #only 2 regions
  rownames(rgn_ge_df[is.na(rgn_ge_df[,"Snca"]),]) #issue: olfactory bulb
  #RIP Irem
  
  invalid_gene_per_region <- apply(rgn_ge_df, 1, function(x) sum(is.na(x)))
  
  #discard genes with over 10 invalid regions from further analysis
  missing_rgn_threshold <- 20
  rgn_ge_df_filt <- rgn_ge_df[,which(invalid_region_per_gene <= missing_rgn_threshold)]
  #got rid of 68 genes
  
  write.csv(rgn_ge_df_filt,paste0(output_dir,'combined_rgn_ge_df_filt_mr',missing_rgn_threshold,
                                  'vv',valid_voxel_threshold,'.csv'))
  
}

#Write colmeans###################################################################
yohan_combined_ge <- as.data.frame(read_csv(paste0(output_dir, 'combined_rgn_ge_df_filt_mr20vv0.2.csv'), col_names=TRUE))
rownames(yohan_combined_ge) <- yohan_combined_ge[,1]
yohan_combined_ge <- yohan_combined_ge[,-1]

data <- yohan_combined_ge
mean_val <- colMeans(data,na.rm = TRUE)

# replacing NA with mean value of each column
for(i in colnames(data))
  data[,i][is.na(data[,i])] <- mean_val[i]
data

write.csv(data, "/data/chamal/projects/natvik/sir_extended/analysis/combined_rgn_ge_df_filt_mr20vv0.2_colmean.csv" , row.names=TRUE)


#################################################################################################
#postprocessing - see how Shady's SNCA correlates with Yohan's

#load in coronal, sagittal, and combined (need to process combined)
yohan_coronal_ge <- as.data.frame(read_csv(paste0(output_dir,'rgn_ge_df_filt_mr10vv0.2.csv'), col_names=TRUE))
rownames(yohan_coronal_ge) <- yohan_coronal_ge[,1]
yohan_coronal_ge <- yohan_coronal_ge[,-1]

yohan_sagittal_ge <- as.data.frame(read_csv(paste0(output_dir, 'sagittal_rgn_ge_df_filt_mr10vv0.2.csv'), col_names=TRUE))
rownames(yohan_sagittal_ge) <- yohan_sagittal_ge[,1]
yohan_sagittal_ge <- yohan_sagittal_ge[,-1]

yohan_combined_ge <- as.data.frame(read_csv(paste0(output_dir, 'combined_rgn_ge_df_filt_mr10vv0.2.csv'), col_names=TRUE))
rownames(yohan_combined_ge) <- yohan_combined_ge[,1]
yohan_combined_ge <- yohan_combined_ge[,-1]


yohan_coronal_snca <- yohan_coronal_ge["Snca"]
yohan_combined_snca <- yohan_combined_ge["Snca"]
yohan_sagittal_snca <- yohan_sagittal_ge["Snca"]

plot(yohan_coronal_snca[!is.na(yohan_sagittal_snca)], yohan_sagittal_snca[!is.na(yohan_sagittal_snca)],
     xlab ="Coronal GE (z score)",
     ylab="Sagittal GE (z score)",
     main="Coronal vs. Sagittal SNCA Expression, regions with > 20% valid voxels")

plot(pnorm(scale(yohan_coronal_snca[!is.na(yohan_sagittal_snca)])),
     pnorm(scale(yohan_sagittal_snca[!is.na(yohan_sagittal_snca)])),
     xlab ="Coronal GE (CDF)",
     ylab="Sagittal GE (CDF)",
     main="Coronal vs. Sagittal SNCA Expression, regions with > 20% valid voxels")

#Load Shady's data
pd <- import("pandas")
shady_snca <- pd$read_pickle("/data/chamal/projects/natvik/sir_extended/analysis/snca_norm.pickle")
shady_snca <- shady_snca[1:213] #equals shady_snca[214:426]
shady_snca <- data.frame(Snca=shady_snca)
rownames(shady_snca) <- steph_allen_rgn_names$Structure

plot(shady_snca[rownames(yohan_sagittal_snca)[which(!is.na(yohan_sagittal_snca))],], 
     pnorm(scale(yohan_sagittal_snca[!is.na(yohan_sagittal_snca)])),
     xlab ="Shady GE (z score)",
     ylab="Yohan Sagittal GE (z score)",
     main="Yohan Sagittal vs. Shady SNCA Expression")

plot(shady_snca[rownames(yohan_combined_snca)[which(!is.na(yohan_combined_snca))],], 
     pnorm(scale(yohan_combined_snca[!is.na(yohan_combined_snca)])),
     xlab ="Shady GE (z score)",
     ylab="Yohan combined GE (z score)",
     main="Yohan combined vs. Shady SNCA Expression")


total_genes_analyzed <- union(colnames(yohan_coronal_ge), colnames(yohan_sagittal_ge))

write.csv(yohan_combined_snca,"/data/chamal/projects/natvik/sir_extended/analysis/combined_snca.csv",row.names=TRUE)


