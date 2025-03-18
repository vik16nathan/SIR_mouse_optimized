library("dplyr")
library("tidyr")
library("readr")
###PART 2: Create a parcellated dataset for each gene, making sure to exclude###
###genes with a low proportion of valid voxels across multiple regions##########

#sort by gene
output_dir <- '/data/chamal/projects/natvik/sir_extended/derivatives/'
combined_exp_filtered <- read.csv(paste0(output_dir,'combined_expression_rgn_min2.csv')) ##takes a few mins to load
yohan_allen_rgn_filt <- sort(unique(combined_exp_filtered$region_name))

gene_list_combined <- sort(unique(combined_exp_filtered$gene))

#command-line arguments
args <- commandArgs(trailingOnly = TRUE)
missing_region_threshold <- args[1]
valid_voxel_threshold <- args[2]

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
  


#Write colmeans###################################################################
mr <- missing_region_threshold
vv <- valid_voxel_threshold
yohan_combined_ge <- as.data.frame(read_csv(paste0(output_dir, 
                               'combined_rgn_ge_df_filt_mr',mr,'vv',vv,'.csv'), col_names=TRUE))
rownames(yohan_combined_ge) <- yohan_combined_ge[,1]
yohan_combined_ge <- yohan_combined_ge[,-1]

data <- yohan_combined_ge
mean_val <- colMeans(data,na.rm = TRUE)

# replacing NA with mean value of each column
for(i in colnames(data))
  data[,i][is.na(data[,i])] <- mean_val[i]
data

write.csv(data, 
        paste0("/data/chamal/projects/natvik/sir_extended/analysis/combined_rgn_ge_df_filt_mr",mr,"vv",vv,"_colmean.csv") , row.names=TRUE)

