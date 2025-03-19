# Read the CSV file (assuming no header)

####sample code for atrophies############
#setwd("/data/chamal/projects/natvik/sir_extended/analysis")
setwd("C:/Users/vik16/Documents/Masters/Figures 20241102")
data <- read.csv("./repeat_results/hipp_ant_redo40_t3000.csv", header = FALSE)

# Remove the first column
data <- data[ , -1]

# Rename the columns
colnames(data) <- c("Spearman_r", "v", "spread_rate", "injection_amount", "k1_atrophy", "k2_atrophy", "peak_timestep")
# Define a function to calculate mean and 95% confidence interval
mean_ci <- function(x, confidence = 0.95) {
  n <- length(x)
  mean_x <- mean(x)
  std_err <- sd(x) / sqrt(n)  # Standard Error
  error_margin <- qt((1 + confidence) / 2, df = n-1) * std_err
  lower <- mean_x - error_margin
  upper <- mean_x + error_margin
  return(c(mean = mean_x, lower_CI = lower, upper_CI = upper))
}

# Apply the function to each column
results <- data.frame(t(apply(data, 2, mean_ci)))

# Display the results
print(round(results,3))



####repeat for IHC########################################
library(ggplot2)
library(dplyr)

# Set working directory
# Function to calculate mean and 95% confidence interval (pre-calculated outside ggplot)
mean_ci <- function(x, confidence = 0.95) {
  n <- length(x)
  mean_x <- mean(x)
  std_err <- sd(x) / sqrt(n)
  error_margin <- qt((1 + confidence) / 2, df = n-1) * std_err
  lower <- mean_x - error_margin
  upper <- mean_x + error_margin
  return(data.frame(mean = mean_x, ymin = lower, ymax = upper))
}

# Function to read and process data, including calculating mean and CI
read_and_process_data <- function(file_prefix, times, mpi_ext) {
  combined_data <- data.frame()
  
  for (time in times) {
    file_path <- paste0("./repeat_results/", file_prefix, "_", mpi_ext, "_t", time, ".csv")
    data <- read.csv(file_path, header = FALSE)
    
    # Remove the first column
    data <- data[ , -1]
    
    # Rename columns
    colnames(data) <- c("Spearman_r", "v", "spread_rate", "injection_amount", "peak_timestep")
    
    # Loop over each column to calculate stats
    for (col in colnames(data)) {
      stats <- mean_ci(data[[col]])
      
      # Create a new dataframe with the raw values and summary statistics
      temp_data <- data.frame(
        values = data[[col]],
        time = time,
        mpi_ext = mpi_ext,
        column = col,
        mean = stats$mean,
        ymin = stats$ymin,
        ymax = stats$ymax
      )
      
      # Combine with the existing data
      combined_data <- rbind(combined_data, temp_data)
    }
  }
  
  return(combined_data)
}

# Function to generate the final plot using facet_grid for efficient multi-panel plotting
generate_multi_panel_plot <- function(file_prefix, times, mpi_exts) {
  all_data <- data.frame()
  
  # Process all mpi extensions and times
  for (mpi in mpi_exts) {
    data <- read_and_process_data(file_prefix, times, mpi)
    all_data <- rbind(all_data, data)
  }
  
  # Set factor levels to place "Spearman_r" at the top and order mpi_exts
  all_data$column <- factor(all_data$column, levels = c("Spearman_r", "v", "spread_rate", "injection_amount", "peak_timestep"))
  all_data$mpi_ext <- factor(all_data$mpi_ext, levels = mpi_exts)  # Order mpi_ext
  
  # Create the plot with facet_grid (5 rows for each column, 7 columns for mpi_ext)
  p <- ggplot(all_data, aes(x = factor(time), y = values, color = factor(time))) +
    geom_point(alpha = 0.5) +
    geom_pointrange(aes(y = mean, ymin = ymin, ymax = ymax), color = "black") +
    facet_grid(column ~ mpi_ext, scales = "free_y", switch = "y") +  # Switch labels to left
    labs(x = "Time") +  # Add the x-axis label
    theme_minimal() +
    theme(
      strip.text = element_text(size = 14),           # Increase size of facet labels
      strip.placement = "outside",                    # Move column labels to outside
      axis.title.y = element_blank(),                 # Remove y-axis title
      axis.title.x = element_text(size = 16, hjust = 0.05),  # Left-align x-axis title
      axis.text = element_text(size = 12),            # Increase size of axis text
      legend.title = element_text(size = 14),         # Increase size of legend title
      legend.text = element_text(size = 12),          # Increase size of legend text
      plot.margin = margin(t = 10, r = 10, b = 30, l = 10)  # Increase bottom margin
    )
  
  # Save the multi-panel plot as a PNG file
  ggsave(filename = paste0(file_prefix, "_multi_panel_plot.jpg"), plot = p, width = 14, height = 10)
}


# Define the file name prefixes, timepoints, and mpi extensions
file_prefixes <- c("hipp_ant_ihc", "hipp_ret_ihc", "cp_ant_ihc", "cp_ret_ihc", "ca1_ret_ihc")
times <- c(1000, 3000, 5000)
mpi_exts <- c("0.5mpi", "1mpi", "3mpi", "6mpi", "12mpi", "18mpi", "24mpi")  # Specify the desired order here

# Process the files and generate the multi-panel plots
for(file_prefix in file_prefixes) {
  generate_multi_panel_plot(file_prefix, times, mpi_exts)
  
}

###save all the means across 40 sets of 200 Optuna trials as baselines before 
###clearance genes
baselines <- data.frame(matrix(nrow=8, ncol=2))
colnames(baselines) <- c("pathology","baseline")
all_data <- data.frame()
i <- 1
for (file_prefix in file_prefixes) {
  for (mpi in c("24mpi")) {
    data <- read_and_process_data(file_prefix, times, mpi)
  }
  baselines[i, "pathology"] <- file_prefix
  ##row 1 is t=1000
  baselines[i, "baseline"] <- data[1,"mean"] ####mean across 40 repeats of 200 optuna trials
  i <- i + 1
}

file_prefixes <- c("hipp_ant_redo40", "hipp_ret_redo40", "cp_ant_redo40", "cp_ret_redo40")
for(file_prefix in file_prefixes) {
  data <- as.data.frame(read.csv(paste0("./repeat_results/",file_prefix,"_t3000.csv"), header=FALSE))
  baselines[i, "pathology"] <- file_prefix
  baselines[i, "baseline"] <- mean(data[,2])
  i <- i+1
  
}

write.csv(baselines, "baselines.csv")