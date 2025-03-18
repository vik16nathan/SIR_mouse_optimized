library("readxl")
library("writexl")

# Define the path to the .xlsx file
input_dir <- "/data/chamal/projects/natvik/sir_extended/preprocessed/shady_ihc/"
output_dir <- "/data/chamal/projects/natvik/sir_extended/preprocessed/shady_ihc/"
ABA_region_names <- as.data.frame(read.csv('/data/chamal/projects/natvik/sir_extended/preprocessed/steph_inputs/ABAnewatlas_final_withacros.csv'))
xlsx_file <- paste0(input_dir,'Supp Material - dataset.xlsx')

# Read the names of the sheets in the Excel file
sheet_names <- excel_sheets(xlsx_file)

rep_right_left <- c(rep("right ", 213), rep("left ", 213))
# Loop through each sheet and write it to a CSV file
for (sheet in sheet_names) {
  # Read the current sheet
  data <- read_excel(xlsx_file, sheet = sheet)
  
  
  ####COMMENT OUT STARTING HERE#################################
  # Get the special names from the first row (superheader denoting
  #total/somatic/neuronal aSyn burden in Shady's xlsx file)
  special_name_e <- as.character(colnames(data[1, 5]))  # Column E
  special_name_l <- as.character(colnames(data[1, 12]))
  special_name_s <- as.character(colnames(data[1, 19]))
  bound_1 <- 11
  bound_2 <- 18
  
  ###edge case
  if(sheet == "ACB injection") {
    print("edge case")
    special_name_e <- as.character(colnames(data[1, 5]))  # Column E
    special_name_l <- as.character(colnames(data[1, 11]))
    special_name_s <- as.character(colnames(data[1, 17]))
    bound_1 <- 10
    bound_2 <- 16
  }
  
  # Prepend the special names to the entries in columns E-K
  for (i in 5:bound_1) {
    if (!is.na(data[1, i])) {
      data[1, i] <- paste(special_name_e, data[1, i], sep = "_")
    }
  }
  
  # Prepend the special names to the entries in columns L-end
  for (i in (bound_1+1):bound_2) {
    if (!is.na(data[1, i])) {
      data[1, i] <- paste(special_name_l, data[1, i], sep = "_")
    }
  }
  
  for (i in (bound_2+1):ncol(data)) {
    if (!is.na(data[1, i])) {
      data[1, i] <- paste(special_name_s, data[1, i], sep = "_")
    }
  }
  
  
  # Set the new column names
  colnames(data) <- as.character(data[1, ])
  
  # Remove the first two rows (original headers and new header row)
  ##remove the first two columns and replace with a column with full names
  
  data <- data[-c(1), ]
  #convert rownames to full region names instead of abbreviated region names
  #use key from Steph/Yohan's data to avoid discarding regions due to spelling issues, etc.
  old_names <- data['Abbreviation']
  print(length(old_names$Abbreviation))
  new_names <- ABA_region_names$Structure[match(old_names$Abbreviation, ABA_region_names$acronym)]
  new_names <- paste0(rep_right_left, new_names)
  data <- data[,-c(1,2)]
  
  data <- cbind(new_names, data)
  

  ##############COMMENT OUT ENDING HERE##########################
  # Create a CSV file name based on the sheet name
  csv_file <- paste0(sheet, ".csv")
  
  # Write the data to a CSV file
  write.csv(data, file = paste0(output_dir,csv_file), row.names = FALSE)
  
  # Print a message indicating that the file has been created
  cat("Created:", csv_file, "\n")
}

