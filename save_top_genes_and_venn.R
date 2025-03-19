# Install ggVennDiagram if you haven't already
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  install.packages("ggVennDiagram")
}

###can always change CA1 to hipp to regenerate older plots
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggVennDiagram)

# Set the working directory
setwd("C:/Users/vik16/Documents/Masters/Figures 20241102")
baselines <- as.data.frame(read.csv("baselines.csv"))

# Load and format gene list data
list1_data <- read.csv("./rankings/cl_hipp_ret_results_eps1e-5.txt", header = FALSE)
list1 <- list1_data[, c(1, 2)]
colnames(list1) <- c("gene", "spearman_r")
list1["baseline"] <- baselines[baselines["pathology"] == "hipp_ret_redo40", "baseline"]

list2_data <- read.csv("./rankings/ca1_ret_ihc_clearance_24mpi_t1000.csv", header = FALSE)
list2 <- list2_data[, c(1, 2)]
colnames(list2) <- c("gene", "spearman_r")
list2["baseline"] <- baselines[baselines$pathology == "ca1_ret_ihc", "baseline"]

# Filter to keep only the genes that outperform the baseline
list1 <- list1 %>%
  filter(spearman_r > baseline)

list2 <- list2 %>%
  filter(spearman_r > baseline)

hipp_atrophy <- list1
hipp_ihc <- list2

# Identify genes in each set
genes <- list(
  "Atrophy" = list1$gene,
  "IHC" = list2$gene
)

# Create the Venn diagram with ggVennDiagram
venn_plot <- ggVennDiagram(genes, 
                           label = "count", 
                           label_alpha = 0, 
                           label_size = 5) + # Adjust label size here
  scale_fill_gradient(low = "lightblue", high = "steelblue") +
  
  # Title and theme adjustments for a cleaner look
  labs(title = "Hipp Ret. Atrophy vs. IHC Pathology") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(t = 15)),
    plot.margin = margin(30, 30, 20, 20), # Increased margin for better fit
    legend.position = "none",
    plot.background = element_blank()
  ) +  scale_x_continuous(expand = expansion(mult = .2))

# Save the plot with a wider, horizontal orientation
ggsave("hipp_ret_venn_diagram.png", plot = venn_plot, width = 10, height = 6, dpi = 300)

# Export lists for unique and shared genes
unique_list1_genes <- setdiff(list1$gene, list2$gene)
unique_list2_genes <- setdiff(list2$gene, list1$gene)
shared_genes <- intersect(list1$gene, list2$gene)

list1_only <- list1 %>%
  filter(gene %in% unique_list1_genes) %>%
  arrange(desc(spearman_r))
write.csv(list1_only, "hipp_ret_atrophy_unique_genes.csv", row.names = FALSE)

list2_only <- list2 %>%
  filter(gene %in% unique_list2_genes) %>%
  arrange(desc(spearman_r))
write.csv(list2_only, "ca1_ret_ihc_unique_genes.csv", row.names = FALSE)

shared_df <- list1 %>%
  filter(gene %in% shared_genes) %>%
  rename(spearman_r_list1 = spearman_r) %>%
  left_join(list2 %>%
              filter(gene %in% shared_genes) %>%
              rename(spearman_r_list2 = spearman_r),
            by = "gene")
write.csv(shared_df, "hipp_ret_shared_genes.csv", row.names = FALSE)


# Load and format gene list data for CP
cp_atrophy_data <- read.csv("./rankings/cl_cp_ant_results_eps1e-5.txt", header = FALSE)
cp_atrophy <- cp_atrophy_data[, c(1, 2)]
colnames(cp_atrophy) <- c("gene", "spearman_r")
cp_atrophy["baseline"] <- baselines[baselines["pathology"] == "cp_ant_redo40", "baseline"]

cp_ihc_data <- read.csv("./rankings/cp_ant_ihc_clearance_24mpi_t1000.csv", header = FALSE)
cp_ihc <- cp_ihc_data[, c(1, 2)]
colnames(cp_ihc) <- c("gene", "spearman_r")
cp_ihc["baseline"] <- baselines[baselines$pathology == "cp_ant_ihc", "baseline"]

# Filter to keep only the genes that outperform the baseline
cp_atrophy <- cp_atrophy %>%
  filter(spearman_r > baseline)

cp_ihc <- cp_ihc %>%
  filter(spearman_r > baseline)

# Define the gene sets
genes <- list(
  "CP Ant. Atrophy" = cp_atrophy$gene,
  "CP Ant. IHC" = cp_ihc$gene,
  "Hipp Ret. Atrophy" = hipp_atrophy$gene,
  "Hipp Ret. IHC" = hipp_ihc$gene
)
# Ensure all sets in `genes` are character vectors


genes <- lapply(genes, as.character)

# Recreate the gene sets for the Venn diagram
venn_genes <- list(
  "CP Ant. Atrophy" = genes[["CP Ant. Atrophy"]],
  "CP Ant. IHC" = genes[["CP Ant. IHC"]],
  "Hipp Ret. Atrophy" = genes[["Hipp Ret. Atrophy"]],
  "Hipp Ret. IHC" = genes[["Hipp Ret. IHC"]]
)

# Generate the Venn diagram
venn_plot <- ggVennDiagram(venn_genes, 
                           label = "count", 
                           label_alpha = 0, 
                           label_size = 4,  # Reduce count label size
                           set_label_size = 3) +  # Reduce category label size
  scale_fill_gradient(low = "lightblue", high = "steelblue") +
  labs(title = "CP & Hipp Ret. Atrophy vs. IHC Pathology") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(t = 10)),  # Increase title margin
    plot.margin = margin(150, 150, 150, 150),  # Add extra space around the diagram
    legend.position = "none",
    plot.background = element_blank()
  )

# Save the Venn diagram
ggsave("filtered_cp_hipp_ret_venn_diagram.png", plot = venn_plot, width = 18, height = 12, dpi = 300)  # Increased dimensions


# Extract intersections manually for saving
# Ensure all sets in `genes` are character vectors
library(dplyr)
cp_atrophy <- data.frame(gene = as.character(cp_atrophy$gene), spearman_r = cp_atrophy$spearman_r)
cp_ihc <- data.frame(gene = as.character(cp_ihc$gene), spearman_r = cp_ihc$spearman_r)
hipp_atrophy <- data.frame(gene = as.character(hipp_atrophy$gene), spearman_r = hipp_atrophy$spearman_r)
hipp_ihc <- data.frame(gene = as.character(hipp_ihc$gene), spearman_r = hipp_ihc$spearman_r)

# Define a function to calculate intersections with correlations
calculate_intersection <- function(df1, df2) {
  intersected_genes <- intersect(df1$gene, df2$gene)
  df1_intersect <- df1[df1$gene %in% intersected_genes, ]
  df2_intersect <- df2[df2$gene %in% intersected_genes, ]
  merge(df1_intersect, df2_intersect, by = "gene", suffixes = c("_df1", "_df2"))
}

# Define a function to calculate unique genes with correlations
calculate_unique <- function(df1, others) {
  unique_genes <- setdiff(df1$gene, unlist(lapply(others, function(df) df$gene)))
  df1[df1$gene %in% unique_genes, ]
}

# Calculate intersections
cp_atrophy_and_hipp_atrophy <- calculate_intersection(cp_atrophy, hipp_atrophy)
cp_ihc_and_hipp_ihc <- calculate_intersection(cp_ihc, hipp_ihc)
cp_atrophy_and_cp_ihc <- calculate_intersection(cp_atrophy, cp_ihc)
hipp_atrophy_and_hipp_ihc <- calculate_intersection(hipp_atrophy, hipp_ihc)
four_way_genes <- Reduce(intersect, list(cp_atrophy$gene, cp_ihc$gene, hipp_atrophy$gene, hipp_ihc$gene))
four_way <- cp_atrophy[cp_atrophy$gene %in% four_way_genes, ]

# Calculate unique sets
unique_cp_atrophy <- calculate_unique(cp_atrophy, list(cp_ihc))
unique_cp_ihc <- calculate_unique(cp_ihc, list(cp_atrophy))
unique_hipp_atrophy <- calculate_unique(hipp_atrophy, list(hipp_ihc))
unique_hipp_ihc <- calculate_unique(hipp_ihc, list(hipp_atrophy))

# Save results to CSV files
write.csv(unique_cp_atrophy, "unique_cp_atrophy_genes_with_correlation.csv", row.names = FALSE)
write.csv(unique_cp_ihc, "unique_cp_ihc_genes_with_correlation.csv", row.names = FALSE)
write.csv(unique_hipp_atrophy, "unique_hipp_atrophy_genes_with_correlation.csv", row.names = FALSE)
write.csv(unique_hipp_ihc, "unique_ca1_ihc_genes_with_correlation.csv", row.names = FALSE)
write.csv(cp_ihc_and_hipp_ihc, "cp_ihc_and_ca1_ihc_genes_with_correlation.csv", row.names = FALSE)
write.csv(cp_atrophy_and_hipp_atrophy, "cp_atrophy_and_hipp_atrophy_genes_with_correlation.csv", row.names = FALSE)
write.csv(cp_atrophy_and_cp_ihc, "cp_atrophy_and_cp_ihc_genes_with_correlation.csv", row.names = FALSE)
write.csv(hipp_atrophy_and_hipp_ihc, "hipp_atrophy_and_ca1_ihc_genes_with_correlation.csv", row.names = FALSE)
write.csv(four_way, "four_way_ca1_intersection_genes_with_correlation.csv", row.names = FALSE)
