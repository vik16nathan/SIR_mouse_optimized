# Load required libraries
library(ggplot2)  # For plotting
library(dplyr)    # For data manipulation
library(readr)    # For reading CSV files
library(scales)

# Set the directory containing the CSV files
directory <- "C:/Users/vik16/Documents/Masters/Figures 20241102/optuna_csvs/35/ihc/"  # Replace with your directory path

# List all CSV files in the directory
csv_files <- list.files(directory, pattern = "^retFalse.*\\.csv$", full.names = TRUE)

# Initialize an empty list to store data frames
data_list <- list()

# Loop through each file, read the data, and store it in the list
for (file in csv_files) {
  data <- read_csv(file)
  data['simulation'] <- 'cp, ant, top 50 clearance'
  if(file == paste0(directory, "retFalse_Zmat4_ihc_24mpi.csv")) {
    data['simulation'] <- "cp, Zmat4, ant"
  }
  data_list[[file]] <- data
}

##repeat for ret
#csv_files <- list.files(directory, pattern = "^retTrue.*\\.csv$", full.names = TRUE)

# Loop through each file, read the data, and store it in the list
#for (file in csv_files) {
#  data <- read_csv(file)
#  data['simulation'] <- 'cp, ret, top 50 clearance'
#  if(file == paste0(directory, "retTrue_Zmat4_ihc_24mpi.csv")) {
#    data['simulation'] <- "cp, Zmat4, ret"
#  }
#  data_list[[file]] <- data
#}


##and for baseline
directory <- "C:/Users/vik16/Documents/Masters/Figures 20241102/optuna_csvs/35/baseline/"  # Replace with your directory path
csv_files <- list.files(directory, pattern = "^retFalse.*\\.csv$", full.names = TRUE)

for (file in csv_files) {
  data <- read_csv(file)
  data['simulation'] <- 'cp, ant, baseline'
  data_list[[file]] <- data
}

#csv_files <- list.files(directory, pattern = "^retTrue.*\\.csv$", full.names = TRUE)

#for (file in csv_files) {
#  data <- read_csv(file)
#  data['simulation'] <- 'cp, ret, baseline'
#  data_list[[file]] <- data
#}

#######repeat everything above for hipp. 
# Set the directory containing the CSV files
directory <- "C:/Users/vik16/Documents/Masters/Figures 20241102/optuna_csvs/24/ihc/"  # Replace with your directory path

# List all CSV files in the directory
#csv_files <- list.files(directory, pattern = "^retFalse.*\\.csv$", full.names = TRUE)

# Initialize an empty list to store data frames

# Loop through each file, read the data, and store it in the list
#for (file in csv_files) {
#  data <- read_csv(file)
#  data['simulation'] <- 'hipp, ant, top 50 clearance'
#  if(file == paste0(directory, "retFalse_Ctnnbip1_ihc_24mpi.csv")) {
#    data['simulation'] <- "hipp, Ctnnbip1, ant"
#  }
#  data_list[[file]] <- data
#}

##repeat for ret
csv_files <- list.files(directory, pattern = "^retTrue.*\\.csv$", full.names = TRUE)

# Loop through each file, read the data, and store it in the list
for (file in csv_files) {
  data <- read_csv(file)
  data['simulation'] <- 'hipp, ret, top 80 clearance'
  if(file == paste0(directory, "retTrue_Spp1_ihc_24mpi.csv")) {
    data['simulation'] <- "hipp, Spp1, ret"
  }
  if(file %in% c(paste0(directory, "retTrue_Acot7_ihc_24mpi.csv"),
                 paste0(directory, "retTrue_Acsl6_ihc_24mpi.csv"),
                 paste0(directory, "retTrue_Elovl5_ihc_24mpi.csv"),
                 paste0(directory, "retTrue_Elovl6_ihc_24mpi.csv"))
     ){
    
    data['simulation'] <- "hipp, long-chain fatty-acyl-CoA\n metabolic process"
  }
  if(file %in% c(paste0(directory, "retTrue_Arl2_ihc_24mpi.csv"),
                 paste0(directory, "retTrue_Capn2_ihc_24mpi.csv"),
                 paste0(directory, "retTrue_Itgb1_ihc_24mpi.csv"),
                 paste0(directory, "retTrue_Limk1_ihc_24mpi.csv"),
                 paste0(directory, "retTrue_Tns1_ihc_24mpi.csv"))
  ){
    data['simulation'] <- "hipp, focal adhesion"
  }
  
  if(file %in% c(paste0(directory, "retTrue_Klk6_ihc_24mpi.csv"),
                 paste0(directory, "retTrue_Pde6g_ihc_24mpi.csv"),
                 paste0(directory, "retTrue_Slc39a14_ihc_24mpi.csv"))
  ){
    data['simulation'] <- "hipp, pos. reg. GPCR signaling"
  }

  data_list[[file]] <- data
}


##and for baseline
directory <- "C:/Users/vik16/Documents/Masters/Figures 20241102/optuna_csvs/24/baseline/"  # Replace with your directory path
#csv_files <- list.files(directory, pattern = "^retFalse.*\\.csv$", full.names = TRUE)

#for (file in csv_files) {
#  data <- read_csv(file)
#  data['simulation'] <- 'hipp, ant, baseline'
#  data_list[[file]] <- data
#}

csv_files <- list.files(directory, pattern = "^retTrue.*\\.csv$", full.names = TRUE)

for (file in csv_files) {
  data <- read_csv(file)
  data['simulation'] <- 'hipp, ret, baseline'
  data_list[[file]] <- data
}

# Combine all data frames into a single data frame
combined_data <- bind_rows(data_list, .id = "source")

# Create a scatterplot
#Load 'no clearance' csv files
#setwd("C:/Users/vik26/Documents/Masters/optuna_csvs/")
#hipp_baseline <- read_csv("hipp_retro_redo_eps1e-5.csv") ##just one baseline file for now
#hipp_baseline['simulation'] <- 'hipp, no clearance'

#cp_baseline <- read_csv("cp_antero_redo_eps1e-5.csv")
#cp_baseline['simulation'] <- 'cp, no clearance'
#combined <- rbind(hipp_baseline, cp_baseline)

#custom_colors <- c("hipp, no clearance" = "#ff969c", 
#                   "hipp, top 50 clearance" = "#65de82",
#                   "hipp, chaperone" = "#d359dd",
#                   "hipp, Pfkfb2" ="#bc2f00",
#                   
#                   "cp, no clearance" = "#0374ec",
#                   "cp, top 50 clearance" =  "#b5c910",
#                   "cp, Zfhx2os"= "#019272")


custom_colors <- c("cp, ant, baseline" = "#b88500",
                   "cp, ant, top 50 clearance" = "#882987",
                   "cp, Zmat4, ant" = "#01ca80",
                   "hipp, focal adhesion" = "#b00046",
                   "hipp, long-chain fatty-acyl-CoA\n metabolic process" = "#595419",
                   "hipp, pos. reg. GPCR signaling" = "#008ECC",
                   "hipp, ret, baseline" = "#e8a7ff",
                   "hipp, ret, top 50 clearance" =   "#ff5255",
                   "hipp, Spp1, ret" = "#ff9b83")
# Plotting
ggplot() +
  #geom_point(data=combined_data,aes(x = params_spread_rate, y = value, colour = simulation), alpha=0.2) + 
  # Smoothed curve with lighter shade of the corresponding line color for combined
  #geom_smooth(data = combined, aes(x = `Param v`, y = Value, colour = simulation, fill = simulation), 
  #            method = "loess", span = 0.3, se = TRUE, alpha = 0.2) +  # alpha controls the transparency of the fill
  
  # Smoothed curve with lighter shade of the corresponding line color for combined_data
  geom_smooth(data = combined_data, aes(x = params_v, y = value, colour = simulation, fill = simulation), 
              method = "loess", span = 0.3, se = TRUE, linetype = "dashed", alpha = 0.2) +  # same transparency for the second dataset +
  # Custom color scheme for lines and fill
  scale_color_manual(values = custom_colors) + 
  scale_fill_manual(values = custom_colors) +  # use the same custom colors for the fill
  ylim(0.4, 0.75) + 
  # Logarithmic scale for x-axis
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  
  # Labels and theme
  labs(title = "24 MPI pSyn Correlation vs aSyn Velocity",
       x = "aSyn Velocity",
       y = "Spearman r") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  ) + 
  guides(colour = guide_legend(override.aes = list(size = 7))) 


# Plotting
ggplot() +
  # Smoothed curve with lighter shade of the corresponding line color for `combined_data` subset to include only "hipp" simulations
  geom_smooth(data = subset(combined_data, grepl("^hipp", simulation)), 
              aes(x = params_spread_rate, y = value, colour = simulation, fill = simulation), 
              method = "loess", span = 0.1, se = TRUE, linetype = "dashed", alpha = 0.2) +  # same transparency for the second dataset

  # Logarithmic scale for x-axis
  #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #              labels = trans_format("log10", math_format(10^.x))) +
  
  # Labels and theme
  labs(title = "Hipp., 24 MPI pSyn Correlation vs aSyn Spreading Rate",
       x = "aSyn spread rate",
       y = "Spearman r (Simulated I fraction vs. IHC pathology)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  ) + 
  guides(colour = guide_legend(override.aes = list(size = 7)))
