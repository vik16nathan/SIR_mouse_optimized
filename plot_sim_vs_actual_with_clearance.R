library(dplyr)
library(ggplot2)

#setwd('C:/Users/vik16/Documents/Masters/Figures 20241102/') #for  = IHC
setwd('C:/Users/vik16/Documents/Masters/ADPD 2025/') #for atrophy
ABA_region_names <- read.csv('./yohan_source_full.csv', header=TRUE)
ABA_region_names<-rbind(ABA_region_names,ABA_region_names)
ABA_region_names <- ABA_region_names[,-1]

#line 1 for IHC; line 2 for atrophy
#sim_atrophy_data <- read.csv("~/Masters/Figures 20241102/simulated_ihc_12012024/abm_spread_v.231.71887821552662.spread_rate.0.015447945715065374.dt.0.1.seed.35.injection_amount.63.92107951571679.clearance_gene.Pfkfb2.csv", header=FALSE)
sim_atrophy_data <- read.csv("~/Masters/ADPD 2025/SIR_clearance_ADPD_20240909/retro_v.0.008943021004069162.spread_rate.0.27372277130315303.dt.0.1.seed.40.injection_amount.84.3914540615607.clearance_gene.Pfkfb2.k1.0.3778796660334787.k2.0.5304646596391381.csv", header = FALSE)

#results <- read.csv('./simulated_ihc_12012024/CP injection.csv')
#rownames(results) <- results[,1]
#results <- results['Total.Pathology_24.MPI']
results <- read.csv('./hipp_rgn_t_stats_hemiMsPff_hemiPBS.csv')
result_rownames <- results[,1]
results <- as.data.frame(results[,"Time.3"])
rownames(results) <- result_rownames
# targets_data <- read.csv('/Users/stephanietullo/Documents/VSC/SIR_mouse-main/model/model/targets.csv', header = FALSE)
ABA_region_names <- read.csv('./yohan_source_full.csv', header = TRUE)
ABA_region_names<-rbind(ABA_region_names,ABA_region_names)
ABA_region_names <- ABA_region_names[,-1]

sim_atrophy_data <- cbind(ABA_region_names, sim_atrophy_data) 

###add hemisphere for the simulated atrophy####
nreg <- 209
alternating_values <- rep(c("right", "left"), each = nreg)
sim_atrophy_data$hemisphere <- alternating_values

# Create unique row names based on first and last columns 
rownames(sim_atrophy_data) <- paste(sim_atrophy_data[, ncol(sim_atrophy_data)], sim_atrophy_data[, 1], sep = " ")

# delete first columns and last column
sim_atrophy_data <- sim_atrophy_data[, -1]
sim_atrophy_data <- sim_atrophy_data[, -ncol(sim_atrophy_data)]

#PLOT
library(ggplot2)
time_steps <- 1:ncol(sim_atrophy_data)  # Assuming the time steps are represented by column indices

# plot(time_steps, sim_atrophy_data["right Caudoputamen",], type = "l", xlab = "Time Steps", ylab = "Atrophy")
# plot(1:4, results["right Caudoputamen",], type = "l", xlab = "Time Steps", ylab = "Atrophy")


# Get the common row names
common_row_names <- intersect(rownames(results), rownames(sim_atrophy_data))

# Sort matrices based on common row names
empirical_common <- results[common_row_names, ]
simulated_common <- sim_atrophy_data[common_row_names, ]

# # Sort matrices based on row names
# empirical_common <- empirical_common[order(rownames(empirical_common)), ]
# simulated_common <- simulated_common[order(rownames(simulated_common)), ]

# Find rows not common between the two dataframes
# not_common_rows <- setdiff(rownames(results), rownames(sim_atrophy_data))
not_common_rows <- setdiff(rownames(sim_atrophy_data), rownames(results))


####find correlations #####

# Initialize a vector to store the correlations
correlations <- vector("numeric", length = ncol(simulated_common))

# Calculate Spearman correlation with each column in sim_atrophy_data
for (i in 1:ncol(simulated_common)) {
  correlations[i] <- cor(empirical_common, simulated_common[, i], method = "spearman")
}

peak_timestep <- which.max(correlations)
peak_corr <- max(correlations)
simulated_peak <- simulated_common[, peak_timestep]
sim_v_emp <- data.frame(simulated = simulated_peak, empirical = empirical_common)
rownames(sim_v_emp) <- rownames(sim_atrophy_data)
###########based on asyn gene expression########
#replace aSyn with clearance gene for visualization!! 

#aSyn <- read.csv("./simulated_ihc_12012024/Pfkfb2.csv", header = TRUE)
aSyn <- read.csv("./Pfkfb2.csv")
aSyn$Pfkfb2 <- scale(aSyn$Pfkfb2)
aSyn <- rbind(aSyn, aSyn) 

###add hemisphere for the simulated atrophy####
aSyn$hemisphere <- alternating_values

# Create unique row names based on first and last columns
#double chek
rownames(aSyn) <- paste(aSyn[, ncol(aSyn)], ABA_region_names, sep = " ")

# delete first columns and last column
aSyn <- aSyn[, -1]

# Get the common row names and not common rows -- double check
common_row_names <- intersect(rownames(sim_v_emp), rownames(aSyn))
not_common_rows <- setdiff(rownames(aSyn), rownames(sim_v_emp))

# Sort matrices based on common row names
sim_v_emp_aSyn <- sim_v_emp[common_row_names, ]
aSyn_common <- aSyn[common_row_names, ]

color_palette <- c("blue", "white", "red")

# Assuming Pfkfb2 has matching rows in both datasets
sim_v_emp_aSyn$Pfkfb2 <- aSyn_common$Pfkfb2

# Pre-scale the simulated and empirical variables to use throughout the plot
sim_v_emp_aSyn$scaled_simulated <- scale(sim_v_emp_aSyn$simulated)
sim_v_emp_aSyn$scaled_empirical <- scale(sim_v_emp_aSyn$empirical)
sim_v_emp_aSyn$scaled_Pfkfb2 <- scale(sim_v_emp_aSyn$Pfkfb2)

# Create a new column to flag the points that meet the criteria based on the pre-scaled columns
sim_v_emp_aSyn$highlight <- with(sim_v_emp_aSyn, ifelse(scaled_empirical > 5.5 | scaled_simulated > 5.5, rownames(sim_v_emp_aSyn), NA))
# Use pre-scaled variables for consistent plotting
library(ggrepel)

ggplot(data = sim_v_emp_aSyn, aes(y = scaled_simulated, x = scaled_empirical, fill = scaled_Pfkfb2)) +
  geom_point(size = 3, shape = 21, stroke = 1, color = "black") +  # Black border, fill based on scaled Pfkfb2
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  
  # Highlight points with different border color based on pre-scaled empirical and simulated columns
  geom_point(data = subset(sim_v_emp_aSyn, scaled_empirical > 5.5 | scaled_simulated > 5.5), 
             aes(y = scaled_simulated, x = scaled_empirical, fill = scaled_Pfkfb2), 
             size = 4, shape = 21, stroke = 2, color = "black") +
  
  # Add labels to the highlighted points with ggrepel to avoid overlap
  geom_text_repel(data = subset(sim_v_emp_aSyn, scaled_empirical > 5.5 | scaled_simulated > 5.5),
                  aes(x = scaled_empirical, y = scaled_simulated, label = highlight), 
                  size = 6, color = "black", max.overlaps = Inf) +
  
  # Continue with the rest of your ggplot setup
  labs(y = "Simulated Atrophy", x = "Empirical Atrophy") +
  scale_fill_gradient2(
    low = color_palette[1],   # Color for low values
    mid = color_palette[2],   # Color for middle value (zero)
    high = color_palette[3],  # Color for high values
    midpoint = 0,             # The central point of the scale
    limits = range(sim_v_emp_aSyn$Pfkfb2),
    breaks = c(min(sim_v_emp_aSyn$Pfkfb2), max(sim_v_emp_aSyn$Pfkfb2))
  ) +
  ggtitle("DG Simulated vs. Empirical Atrophy") + 
  labs(fill = 'Pfkfb2') +  # Adjusted legend label
  theme(
    plot.title = element_text(hjust = 0.5, size = 28, face = "bold"),  # Center and enlarge the title
    axis.title.x = element_text(size = 24, face = "bold"),  # Enlarge x-axis label
    axis.title.y = element_text(size = 24, face = "bold"),  # Enlarge y-axis label
    legend.title = element_text(size = 20, face = "bold"),  # Enlarge legend title
    legend.text = element_text(size = 20),  # Enlarge legend text
    legend.key.size = unit(2, "cm"),  # Enlarge legend key (symbol) size
    axis.text.x = element_text(size = 20),  # Enlarge x-axis tick labels
    axis.text.y = element_text(size = 20)   # Enlarge y-axis tick labels
  )
###################check correlation between asyn gene expression and simulated atrophy##########

###sim atrophy vs Pfkfb2 expression####
correlations_sim_Pfkfb2 <- cor(aSyn_common$Pfkfb2, sim_v_emp_aSyn[,1], method = "spearman")
 print(correlations_sim_Pfkfb2)

# Create data frame for plotting
plot_data <- data.frame(
  x = aSyn_common$Pfkfb2,
  y = sim_v_emp_aSyn[,1],
  region = rownames(sim_v_emp_aSyn)
)

# Scatter plot with ggplot
ggplot(data = plot_data, aes(x = x, y = y)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, color = "#0072B2") +
  geom_text(x = min(plot_data$x), y = max(plot_data$y), 
            label = paste0("Spearman correlation = ", round(correlations_sim_Pfkfb2, 3)), hjust = 0, vjust = 1,color="black",size=6) +
  labs(x = "Pfkfb2 expression", y = "Simulated Atrophy", title = "") +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 18),
    plot.title = element_text(size = 13)
  )
