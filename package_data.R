library(tidyverse)
library(glue)
library(RMINC)
library(data.tree)
library(data.table)

# Load all gene data
df <- read_csv("/projects/yyee/tools/Allen_brain/data/gene_data.csv") %>%
  mutate(fc=paste(msg.genes.acronym, msg.id, sep="_")) %>%
  select(filetag=fc, id=msg.id, 
         delegate=msg.delegate, failed=msg.failed, expression=msg.expression,
         gene_id=msg.genes.id, gene_acronym=msg.genes.acronym, gene_name=msg.genes.name,
         gene_original_symbol=msg.genes.original_symbol, gene_original_name=msg.genes.original_name,
         gene_aliases=msg.genes.alias_tags,
         gene_ensembl_id=msg.genes.ensembl_id, gene_legacy_ensembl_id=msg.genes.legacy_ensembl_gene_id, 
         gene_entrez_id=msg.genes.entrez_id, gene_homologene_id=msg.genes.homologene_id,
         gene_chromosome_id=msg.genes.chromosome_id,
         gene_version_status=msg.genes.version_status,
         imaging_plane_of_section=msg.plane_of_section_id, imaging_section_thickness=msg.section_thickness,
         imaging_reference_space=msg.reference_space_id,
         specimen_organism_id=msg.genes.organism_id, specimen_id=msg.specimen_id, specimen_weight=msg.weight,
         qc_date=msg.qc_date, 
  )

# Coronal files
coronal_files <- list.files("resources/expression/P56/coronal", full.names = T, pattern = "*.mnc")
df_coronal <- df %>%
  filter(!failed) %>%
  filter(filetag %in% gsub('.mnc', '', basename(coronal_files)))

# Sagittal files
sagittal_files <- list.files("resources/expression/P56/sagittal", full.names = T, pattern = "*.mnc")
df_sagittal <- df %>%
  filter(!failed) %>%
  filter(filetag %in% gsub('.mnc', '', basename(sagittal_files)))

# Load atlas
atlas <- "/projects/yyee/ip/2018-11-structural-covariance-genetics/data/outputs/mouse_atlas/atlas/AMBA_relabeled_backsampled_200um.mnc"
avol <- as.integer(round(mincGetVolume(atlas)))
mvol <- avol > 0.5
avec <- avol[mvol > 0.5]
all_labels <- sort(unique(avol))[-1]

##########
# Load tree
load("/projects/yyee/ip/2018-11-structural-covariance-genetics/data/outputs/mouse_atlas/atlas/AMBA_relabeled_backsampled_25um_split_hanat.RData")

# Drop split labels
Prune(tree_original, pruneFun = function(node) {node$is_unsplit})

# Given values detected in atlas, annotate which nodes are present and remove child nodes without any subtrees that have an atlas label
tree_original$Do(function(node) {ifelse(node$id %in% all_labels, node$total_child_labels <- 1, node$total_child_labels <- 0)})
tree_original$Do(function(node) {if(node$id %in% all_labels) {node$child_labels <- c(node$id)}})
tree_original$Do(function(node) {
  if (!is.null(node$parent)) {
    node$parent$total_child_labels <- node$parent$total_child_labels + node$total_child_labels
    node$parent$child_labels <- c(node$parent$child_labels, node$child_labels)
  } 
}, traversal = "post-order")
Prune(tree_original, pruneFun = function(node) {node$total_child_labels > 0})

# Make sagittal and coronal copies
tree_coronal <- Clone(tree_original)
tree_sagittal <- Clone(tree_original)

# Set filetag ordering
tree_coronal$Do(function(node) {
  node$filetag <- df_coronal$filetag
})
tree_sagittal$Do(function(node) {
  node$filetag <- df_sagittal$filetag
})

##########
# Compute expression voxel data

label_voxel_counts <- table(avec)
df_total_voxel <- tibble(label=as.integer(names(label_voxel_counts)), voxel_total_sum=label_voxel_counts)

get_gridded_gene_label_data <- function(filetag, dataset) {
  path <- glue("resources/expression/P56/{dataset}/{filetag}.mnc")
  gvec <- mincGetVolume(path)[mvol > 0.5]
  
  df_exp <- tibble(label=avec, expression=gvec) %>%
    mutate(valid=(gvec >= -0.5))
  
  df_exp_voxel <- df_exp %>%
    group_by(label) %>%
    summarize(voxel_valid_sum=sum(valid))
  
  df_exp_expression <- df_exp %>%
    filter(valid) %>%
    group_by(label) %>%
    summarize(expression_sum=sum(expression))
  
  df_out <- df_total_voxel %>%
    left_join(df_exp_voxel, by="label") %>%
    left_join(df_exp_expression, by="label") %>%
    mutate(filetag=filetag, dataset=dataset) %>%
    select(filetag, dataset, everything())
  
  return(df_out)
}

# Get data
if (file.exists("workspace/atlas_expression.RData")) {
  load("workspace/atlas_expression.RData")
} else {
  
  # Get data - coronal
  # ~5-10 minutes
  start_time <- Sys.time()
  atlas_expression_coronal_list <- list()
  prog <- txtProgressBar(max=nrow(df_coronal), style=3)
  for (i in 1:nrow(df_coronal)) {
    atlas_expression_coronal_list[[i]] <- get_gridded_gene_label_data(df_coronal$filetag[[i]], "coronal")
    setTxtProgressBar(prog, i)
  }
  close(prog)
  atlas_expression_coronal <- rbindlist(atlas_expression_coronal_list)
  rm(atlas_expression_coronal_list)
  end_time <- Sys.time()
  coronal_tdiff <- end_time - start_time
  print(glue("Start: {start_time} | End: {end_time}"))
  
  # Get data - sagittal
  # ~30-60 minutes
  start_time <- Sys.time()
  atlas_expression_sagittal_list <- list()
  prog <- txtProgressBar(max=nrow(df_sagittal), style=3)
  for (i in 1:nrow(df_sagittal)) {
    atlas_expression_sagittal_list[[i]] <- get_gridded_gene_label_data(df_sagittal$filetag[[i]], "sagittal")
    setTxtProgressBar(prog, i)
  }
  close(prog)
  atlas_expression_sagittal <- rbindlist(atlas_expression_sagittal_list)
  rm(atlas_expression_sagittal_list)
  end_time <- Sys.time()
  sagittal_tdiff <- end_time - start_time
  print(glue("Start: {start_time} | End: {end_time}"))
  
  # save
  save(list = c("atlas_expression_coronal", "atlas_expression_sagittal", "coronal_tdiff", "sagittal_tdiff"),
       file = "workspace/atlas_expression.RData")
  
}

# Add data to child nodes containing atlas labels
tree_coronal$Do(function(node) {
  print(glue("Working on node: {node$name}"))
  df_node <- tibble(filetag=node$filetag) %>%
    left_join(atlas_expression_coronal %>% filter(label==node$id), by="filetag") 
  
  node$voxel_total_sum <- df_node$voxel_total_sum
  node$voxel_valid_sum <- df_node$voxel_valid_sum
  node$expression_sum <- df_node$expression_sum
}, filterFun=function(node) {
  node$id %in% all_labels
})

tree_sagittal$Do(function(node) {
  print(glue("Working on node: {node$name}"))
  df_node <- tibble(filetag=node$filetag) %>%
    left_join(atlas_expression_sagittal %>% filter(label==node$id), by="filetag") 
  
  node$voxel_total_sum <- df_node$voxel_total_sum
  node$voxel_valid_sum <- df_node$voxel_valid_sum
  node$expression_sum <- df_node$expression_sum
}, filterFun=function(node) {
  node$id %in% all_labels
})

##########
# Aggregate: voxel_valid_sum, voxel_total_sum, expression_sum

add_to_parent <- function(node, attribute) {
  
  if (is.null(node$parent[[attribute]])) {
    node$parent[[attribute]] <- node[[attribute]]
  } else {
    node$parent[[attribute]] <- rowSums(
      data.frame(
        parent=as.numeric(node$parent[[attribute]]),
        child=as.numeric(node[[attribute]])
      ), na.rm=T)
  }
  
}

tree_coronal$Do(function(node) {
  if (!is.null(node$parent)) {
    print(glue("Working on node: {node$name}"))
    
    add_to_parent(node, "voxel_valid_sum")
    add_to_parent(node, "voxel_total_sum")
    add_to_parent(node, "expression_sum")
    
  } 
}, traversal = "post-order")

tree_sagittal$Do(function(node) {
  if (!is.null(node$parent)) {
    print(glue("Working on node: {node$name}"))
    
    add_to_parent(node, "voxel_valid_sum")
    add_to_parent(node, "voxel_total_sum")
    add_to_parent(node, "expression_sum")
    
  } 
}, traversal = "post-order")


##########
# At each node, compute: expression_mean (expression_sum/voxel_valid_sum), voxel_valid_mean (voxel_valid_sum/voxel_total_sum)

tree_coronal$Do(function(node) {
  node$expression_mean <- node$expression_sum / node$voxel_valid_sum
  node$voxel_valid_mean <- node$voxel_valid_sum / node$voxel_total_sum
})

tree_sagittal$Do(function(node) {
  node$expression_mean <- node$expression_sum / node$voxel_valid_sum
  node$voxel_valid_mean <- node$voxel_valid_sum / node$voxel_total_sum
})

# At each node: Normalize expression by parent, normalize by brain, divide by parent, divide by brain
tree_coronal$Do(function(node) {
  
  print(glue("Working on node: {node$name}"))
  
  root_node <- tree_coronal$root
  dataset_name <- "Coronal"
  
  if (!is.null(node$parent)) {
    # Has parent
    mod_parent <- lm(node$expression_mean ~ node$parent$expression_mean, na.action=na.exclude)
    mod_root <- lm(node$expression_mean ~ root_node$expression_mean, na.action=na.exclude)
    
    # Normalized expression
    node$expression_mean_normalized_parent <- residuals(mod_parent)
    node$expression_mean_normalized_brain <- residuals(mod_root)
    node$expression_mean_relative_parent <- node$expression_mean/node$parent$expression_mean
    node$expression_mean_relative_brain <- node$expression_mean/root_node$expression_mean
    
    # Normalization metrics
    mod <- mod_parent
    model_name <- "Parent"
    node$normalization_parent_model_coefficients <- 
      summary(mod)$coefficients %>%
      as.data.frame() %>%
      as.tibble() %>%
      mutate(term=rownames(summary(mod)$coefficients)) %>%
      mutate(term=case_when(
        term=="(Intercept)" ~ "Intercept",
        TRUE ~ "Slope")) %>%
      gather(metric, value, -term) %>%
      arrange(term, metric) %>%
      rbind(tibble(term="Model", metric="Adjusted R2", value=summary(mod)$adj.r.squared)) %>%
      mutate(model=model_name, dataset=dataset_name, 
             region_id=node$id, region_name=node$name, region_level=node$level, region_colour=node$color_hex_triplet, 
             region_parent_id=node$parent$id, region_children_ids=paste(sapply(node$children, "[[", "id"), collapse = ";"), region_in_atlas=(node$id %in% all_labels)) %>%
      select(model, dataset, term, metric, value, everything())
    
    mod <- mod_root
    model_name <- "Brain"
    node$normalization_brain_model_coefficients <- 
      summary(mod)$coefficients %>%
      as.data.frame() %>%
      as.tibble() %>%
      mutate(term=rownames(summary(mod)$coefficients)) %>%
      mutate(term=case_when(
        term=="(Intercept)" ~ "Intercept",
        TRUE ~ "Slope")) %>%
      gather(metric, value, -term) %>%
      arrange(term, metric) %>%
      rbind(tibble(term="Model", metric="Adjusted R2", value=summary(mod)$adj.r.squared)) %>%
      mutate(model=model_name, dataset=dataset_name, 
             region_id=node$id, region_name=node$name, region_level=node$level, region_colour=node$color_hex_triplet, 
             region_parent_id=node$parent$id, region_children_ids=paste(sapply(node$children, "[[", "id"), collapse = ";"), region_in_atlas=(node$id %in% all_labels)) %>%
      select(model, dataset, term, metric, value, everything())
    
  } else {
    # Is top of tree
    mod_parent <- lm(node$expression_mean ~ root_node$expression_mean, na.action=na.exclude)
    mod_root <- lm(node$expression_mean ~ root_node$expression_mean, na.action=na.exclude)
    
    # Normalized expression
    node$expression_mean_normalized_parent <- rep(0, length(root_node$expression_mean))
    node$expression_mean_normalized_brain <- rep(0, length(root_node$expression_mean))
    node$expression_mean_relative_parent <- rep(1, length(root_node$expression_mean))
    node$expression_mean_relative_brain <- rep(1, length(root_node$expression_mean))
    
    # Normalization metrics
    mod <- mod_parent
    model_name <- "Parent"
    node$normalization_parent_model_coefficients <- 
      summary(mod)$coefficients %>%
      as.data.frame() %>%
      as.tibble() %>%
      mutate(term=rownames(summary(mod)$coefficients)) %>%
      mutate(term=case_when(
        term=="(Intercept)" ~ "Intercept",
        TRUE ~ "Slope")) %>%
      gather(metric, value, -term) %>%
      arrange(term, metric) %>%
      rbind(tibble(term="Model", metric="Adjusted R2", value=summary(mod)$adj.r.squared)) %>%
      mutate(model=model_name, dataset=dataset_name, 
             region_id=node$id, region_name=node$name, region_level=node$level, region_colour=node$color_hex_triplet, 
             region_parent_id=NA, region_children_ids=paste(sapply(node$children, "[[", "id"), collapse = ";"), region_in_atlas=(node$id %in% all_labels)) %>%
      select(model, dataset, term, metric, value, everything())
    
    mod <- mod_root
    model_name <- "Brain"
    node$normalization_brain_model_coefficients <- 
      summary(mod)$coefficients %>%
      as.data.frame() %>%
      as.tibble() %>%
      mutate(term=rownames(summary(mod)$coefficients)) %>%
      mutate(term=case_when(
        term=="(Intercept)" ~ "Intercept",
        TRUE ~ "Slope")) %>%
      gather(metric, value, -term) %>%
      arrange(term, metric) %>%
      rbind(tibble(term="Model", metric="Adjusted R2", value=summary(mod)$adj.r.squared)) %>%
      mutate(model=model_name, dataset=dataset_name, 
             region_id=node$id, region_name=node$name, region_level=node$level, region_colour=node$color_hex_triplet, 
             region_parent_id=NA, region_children_ids=paste(sapply(node$children, "[[", "id"), collapse = ";"), region_in_atlas=(node$id %in% all_labels)) %>%
      select(model, dataset, term, metric, value, everything())
    
  }
})

tree_sagittal$Do(function(node) {
  
  print(glue("Working on node: {node$name}"))
  
  root_node <- tree_sagittal$root
  dataset_name <- "Sagittal"
  
  if (!is.null(node$parent)) {
    # Has parent
    mod_parent <- lm(node$expression_mean ~ node$parent$expression_mean, na.action=na.exclude)
    mod_root <- lm(node$expression_mean ~ root_node$expression_mean, na.action=na.exclude)
    
    # Normalized expression
    node$expression_mean_normalized_parent <- residuals(mod_parent)
    node$expression_mean_normalized_brain <- residuals(mod_root)
    node$expression_mean_relative_parent <- node$expression_mean/node$parent$expression_mean
    node$expression_mean_relative_brain <- node$expression_mean/root_node$expression_mean
    
    # Normalization metrics
    mod <- mod_parent
    model_name <- "Parent"
    node$normalization_parent_model_coefficients <- 
      summary(mod)$coefficients %>%
      as.data.frame() %>%
      as.tibble() %>%
      mutate(term=rownames(summary(mod)$coefficients)) %>%
      mutate(term=case_when(
        term=="(Intercept)" ~ "Intercept",
        TRUE ~ "Slope")) %>%
      gather(metric, value, -term) %>%
      arrange(term, metric) %>%
      rbind(tibble(term="Model", metric="Adjusted R2", value=summary(mod)$adj.r.squared)) %>%
      mutate(model=model_name, dataset=dataset_name, 
             region_id=node$id, region_name=node$name, region_level=node$level, region_colour=node$color_hex_triplet, 
             region_parent_id=node$parent$id, region_children_ids=paste(sapply(node$children, "[[", "id"), collapse = ";"), region_in_atlas=(node$id %in% all_labels)) %>%
      select(model, dataset, term, metric, value, everything())
    
    mod <- mod_root
    model_name <- "Brain"
    node$normalization_brain_model_coefficients <- 
      summary(mod)$coefficients %>%
      as.data.frame() %>%
      as.tibble() %>%
      mutate(term=rownames(summary(mod)$coefficients)) %>%
      mutate(term=case_when(
        term=="(Intercept)" ~ "Intercept",
        TRUE ~ "Slope")) %>%
      gather(metric, value, -term) %>%
      arrange(term, metric) %>%
      rbind(tibble(term="Model", metric="Adjusted R2", value=summary(mod)$adj.r.squared)) %>%
      mutate(model=model_name, dataset=dataset_name, 
             region_id=node$id, region_name=node$name, region_level=node$level, region_colour=node$color_hex_triplet, 
             region_parent_id=node$parent$id, region_children_ids=paste(sapply(node$children, "[[", "id"), collapse = ";"), region_in_atlas=(node$id %in% all_labels)) %>%
      select(model, dataset, term, metric, value, everything())
    
  } else {
    # Is top of tree
    mod_parent <- lm(node$expression_mean ~ root_node$expression_mean, na.action=na.exclude)
    mod_root <- lm(node$expression_mean ~ root_node$expression_mean, na.action=na.exclude)
    
    # Normalized expression
    node$expression_mean_normalized_parent <- rep(0, length(root_node$expression_mean))
    node$expression_mean_normalized_brain <- rep(0, length(root_node$expression_mean))
    node$expression_mean_relative_parent <- rep(1, length(root_node$expression_mean))
    node$expression_mean_relative_brain <- rep(1, length(root_node$expression_mean))
    
    # Normalization metrics
    mod <- mod_parent
    model_name <- "Parent"
    node$normalization_parent_model_coefficients <- 
      summary(mod)$coefficients %>%
      as.data.frame() %>%
      as.tibble() %>%
      mutate(term=rownames(summary(mod)$coefficients)) %>%
      mutate(term=case_when(
        term=="(Intercept)" ~ "Intercept",
        TRUE ~ "Slope")) %>%
      gather(metric, value, -term) %>%
      arrange(term, metric) %>%
      rbind(tibble(term="Model", metric="Adjusted R2", value=summary(mod)$adj.r.squared)) %>%
      mutate(model=model_name, dataset=dataset_name, 
             region_id=node$id, region_name=node$name, region_level=node$level, region_colour=node$color_hex_triplet, 
             region_parent_id=NA, region_children_ids=paste(sapply(node$children, "[[", "id"), collapse = ";"), region_in_atlas=(node$id %in% all_labels)) %>%
      select(model, dataset, term, metric, value, everything())
    
    mod <- mod_root
    model_name <- "Brain"
    node$normalization_brain_model_coefficients <- 
      summary(mod)$coefficients %>%
      as.data.frame() %>%
      as.tibble() %>%
      mutate(term=rownames(summary(mod)$coefficients)) %>%
      mutate(term=case_when(
        term=="(Intercept)" ~ "Intercept",
        TRUE ~ "Slope")) %>%
      gather(metric, value, -term) %>%
      arrange(term, metric) %>%
      rbind(tibble(term="Model", metric="Adjusted R2", value=summary(mod)$adj.r.squared)) %>%
      mutate(model=model_name, dataset=dataset_name, 
             region_id=node$id, region_name=node$name, region_level=node$level, region_colour=node$color_hex_triplet, 
             region_parent_id=NA, region_children_ids=paste(sapply(node$children, "[[", "id"), collapse = ";"), region_in_atlas=(node$id %in% all_labels)) %>%
      select(model, dataset, term, metric, value, everything())
    
  }
})

save(list=c("tree_coronal", "tree_sagittal"), file = "workspace/expression_trees.RData")

##########
# For Tomas:

tomas_coronal_list <- tree_coronal$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$Get(function(node) {
  tibble(
    dataset="Coronal",
    filetag=node$filetag,
    gene=sapply(strsplit(node$filetag, "_"), "[[", 1),
    id=sapply(strsplit(node$filetag, "_"), "[[", 2),
    region_id=node$id, 
    region_name=node$name, 
    region_level=node$level, 
    region_colour=node$color_hex_triplet, 
    region_parent_id=ifelse(is.null(node$parent), NA, node$parent$id), 
    region_children_ids=paste(sapply(node$children, "[[", "id"), collapse = ","), 
    region_in_atlas=(node$id %in% all_labels),
    voxel_total_sum=as.numeric(node$voxel_total_sum), 
    voxel_valid_sum=as.numeric(node$voxel_valid_sum), 
    voxel_valid_mean=as.numeric(node$voxel_valid_mean), 
    expression_sum=as.numeric(node$expression_sum), 
    expression_mean=as.numeric(node$expression_mean),
    expression_mean_relative_parent=as.numeric(node$expression_mean_relative_parent),
    expression_mean_relative_brain=as.numeric(node$expression_mean_relative_brain),
    expression_mean_normalized_parent=as.numeric(node$expression_mean_normalized_parent),
    expression_mean_normalized_brain=as.numeric(node$expression_mean_normalized_brain)
  )
}, simplify = F)

tomas_sagittal_list <- tree_sagittal$`Basic cell groups and regions`$Cerebrum$`Cerebral cortex`$`Cortical plate`$Isocortex$Get(function(node) {
  tibble(
    dataset="Sagittal",
    filetag=node$filetag,
    gene=sapply(strsplit(node$filetag, "_"), "[[", 1),
    id=sapply(strsplit(node$filetag, "_"), "[[", 2),
    region_id=node$id, 
    region_name=node$name, 
    region_level=node$level, 
    region_colour=node$color_hex_triplet, 
    region_parent_id=ifelse(is.null(node$parent), NA, node$parent$id), 
    region_children_ids=paste(sapply(node$children, "[[", "id"), collapse = ","), 
    region_in_atlas=(node$id %in% all_labels),
    voxel_total_sum=as.numeric(node$voxel_total_sum), 
    voxel_valid_sum=as.numeric(node$voxel_valid_sum), 
    voxel_valid_mean=as.numeric(node$voxel_valid_mean), 
    expression_sum=as.numeric(node$expression_sum), 
    expression_mean=as.numeric(node$expression_mean),
    expression_mean_relative_parent=as.numeric(node$expression_mean_relative_parent),
    expression_mean_relative_brain=as.numeric(node$expression_mean_relative_brain),
    expression_mean_normalized_parent=as.numeric(node$expression_mean_normalized_parent),
    expression_mean_normalized_brain=as.numeric(node$expression_mean_normalized_brain)
  )
}, simplify = F)


tomas_coronal_data <- tomas_coronal_list %>%
  rbindlist()
rm(tomas_coronal_list)

tomas_sagittal_data <- tomas_sagittal_list %>%
  rbindlist()
rm(tomas_sagittal_list)

dir.create("workspace/for_tomas")
tomas_coronal_data %>% write_csv(file="workspace/for_tomas/coronal_expression_data.csv")
tomas_sagittal_data %>% write_csv(file="workspace/for_tomas/sagittal_expression_data.csv")
df_coronal %>% write_csv(file="workspace/for_tomas/coronal_gene_info.csv")
df_sagittal %>% write_csv(file="workspace/for_tomas/sagittal_gene_info.csv")

##########
# For Andrea:
andrea_coronal_list <- tree_coronal$root$Get(function(node) {
  tibble(
    dataset="Coronal",
    filetag=node$filetag,
    gene=sapply(strsplit(node$filetag, "_"), "[[", 1),
    id=sapply(strsplit(node$filetag, "_"), "[[", 2),
    region_id=node$id, 
    region_name=node$name, 
    region_level=node$level, 
    region_colour=node$color_hex_triplet, 
    region_parent_id=ifelse(is.null(node$parent), NA, node$parent$id), 
    region_children_ids=paste(sapply(node$children, "[[", "id"), collapse = ","), 
    region_in_atlas=(node$id %in% all_labels),
    voxel_total_sum=as.numeric(node$voxel_total_sum), 
    voxel_valid_sum=as.numeric(node$voxel_valid_sum), 
    voxel_valid_mean=as.numeric(node$voxel_valid_mean), 
    expression_sum=as.numeric(node$expression_sum), 
    expression_mean=as.numeric(node$expression_mean),
    expression_mean_relative_parent=as.numeric(node$expression_mean_relative_parent),
    expression_mean_relative_brain=as.numeric(node$expression_mean_relative_brain),
    expression_mean_normalized_parent=as.numeric(node$expression_mean_normalized_parent),
    expression_mean_normalized_brain=as.numeric(node$expression_mean_normalized_brain)
  )
}, simplify = F)

andrea_sagittal_list <- tree_sagittal$root$Get(function(node) {
  tibble(
    dataset="Sagittal",
    filetag=node$filetag,
    gene=sapply(strsplit(node$filetag, "_"), "[[", 1),
    id=sapply(strsplit(node$filetag, "_"), "[[", 2),
    region_id=node$id, 
    region_name=node$name, 
    region_level=node$level, 
    region_colour=node$color_hex_triplet, 
    region_parent_id=ifelse(is.null(node$parent), NA, node$parent$id), 
    region_children_ids=paste(sapply(node$children, "[[", "id"), collapse = ","), 
    region_in_atlas=(node$id %in% all_labels),
    voxel_total_sum=as.numeric(node$voxel_total_sum), 
    voxel_valid_sum=as.numeric(node$voxel_valid_sum), 
    voxel_valid_mean=as.numeric(node$voxel_valid_mean), 
    expression_sum=as.numeric(node$expression_sum), 
    expression_mean=as.numeric(node$expression_mean),
    expression_mean_relative_parent=as.numeric(node$expression_mean_relative_parent),
    expression_mean_relative_brain=as.numeric(node$expression_mean_relative_brain),
    expression_mean_normalized_parent=as.numeric(node$expression_mean_normalized_parent),
    expression_mean_normalized_brain=as.numeric(node$expression_mean_normalized_brain)
  )
}, simplify = F)


andrea_coronal_data <- andrea_coronal_list %>%
  rbindlist()
rm(andrea_coronal_list)

andrea_sagittal_data <- andrea_sagittal_list %>%
  rbindlist()
rm(andrea_sagittal_list)

dir.create("workspace/for_andrea")
andrea_coronal_data %>% write_csv(file="workspace/for_andrea/coronal_expression_data.csv")
andrea_sagittal_data %>% write_csv(file="workspace/for_andrea/sagittal_expression_data.csv")
df_coronal %>% write_csv(file="workspace/for_andrea/coronal_gene_info.csv")
df_sagittal %>% write_csv(file="workspace/for_andrea/sagittal_gene_info.csv")

# Package nifti files for Andrea

# Package
# - provide link to template and atlas
# - describe expression format of nifti files + missing data
# - general overview of process to download and convert to RAW
# - note some repeats, and other caveats