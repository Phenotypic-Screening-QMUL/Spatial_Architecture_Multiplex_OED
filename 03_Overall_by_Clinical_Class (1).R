# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Convert Percentage_Data to long format
percentage_data_long <- convert_to_long(Percentage_Data, c("Slide", "Region", "Analysis Region", "Transformation", "Clinical_Class"))

# Filter and shorten the variable strings
percentage_data_long <- percentage_data_long %>% 
  filter((`Analysis Region` == "Epi" & grepl("^Epi: ", variable)) |
           (`Analysis Region` != "Epi" & grepl("^Stroma: ", variable))) %>%
  mutate(variable = gsub("^(Epi: |Stroma: |% )", "", variable)) %>%  # Shorten strings and remove percentage sign
  rename(Cell_Type = variable) %>%  # Rename variable to Cell_Type
  select(Slide, Region, `Analysis Region`, Transformation, Clinical_Class, Cell_Type, value)


# Calculate the mean percentage of positive cells for each Cell_Type within each Clinical_Class
mean_percentage_data <- percentage_data_long %>%
  group_by(Slide, Cell_Type, Clinical_Class) %>%
  summarise(Mean_Value = mean(value, na.rm = TRUE), .groups = 'drop') %>%  # Ensure groups are dropped after summarizing
  ungroup()

mean_percentage_data <- mean_percentage_data %>%
  mutate(Clinical_Class = factor(Clinical_Class, 
                                 levels = c("Non-dysplastic (Control)", 
                                            "Hyperkeratosis", "Mild OED", "Moderate OED", "Severe OED", 
                                            "HPV OED", 
                                            "Veruccous OED","OSCC")))

# Remove the "%" from the values in the "Cell_Type" column without adding spaces
mean_percentage_data$Cell_Type <- gsub("%", "", mean_percentage_data$Cell_Type, fixed = TRUE)

# Trim whitespace from the beginning and end of the values
mean_percentage_data$Cell_Type <- trimws(mean_percentage_data$Cell_Type)

# Define distinct colors for each Clinical_Class, including hyperkeratosis
clinical_class_colors <- c(
  "Non-dysplastic (Control)" = "#FF9999",            # Light red
  "Hyperkeratosis" = "#FF8080",                      # Soft red
  "Mild OED" = "#FF4D4D",                            # Medium red
  "Moderate OED" = "#FF1A1A",                       # Darker red
  "Severe OED" = "#CC0000",                         # Strong red
  "OSCC" = "#7F0000",                               # Very dark red
  "HPV OED" = "#A30000",                            # Dark red
  "Veruccous OED" = "#990000"                       # Deepest red
)

# Define categories to exclude
excluded_categories <- c("Hyperkeratosis", "HPV OED", "Veruccous OED")

# Filter out the excluded categories from the data
mean_percentage_data <- mean_percentage_data %>%
  filter(!Clinical_Class %in% excluded_categories)

# Define categories to exclude
excluded_regions <- c("Sub002", "Sub003", "SubIN")


# Create the directory if it doesn't exist for Clinical Class plots
output_dir <- "Clinical_Class_Percentage_Positive_Overall"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Get unique Cell_Types
unique_analysis_regions <- unique(mean_percentage_data$`Cell_Type`)

# Loop through each unique Cell_Type
for (analysis_region in unique_analysis_regions) {
  # Filter the data for the current Cell_Type
  filtered_data <- mean_percentage_data[mean_percentage_data$`Cell_Type` == analysis_region, ]
  
  # Loop through each unique Cell_Type for the current Cell_Type
  unique_cell_types <- unique(filtered_data$Cell_Type)
  
  for (cell_type in unique_cell_types) {
    # Filter the data for the current Cell_Type
    cell_type_data <- filtered_data[filtered_data$Cell_Type == cell_type, ]
    
    # Create the bar plot for the current Cell_Type and Cell_Type
    bar_plot <- ggplot(cell_type_data, aes(x = Clinical_Class, y = Mean_Value, fill = Clinical_Class)) +
      geom_bar(stat = "summary", fun = "mean", position = "dodge") +  # Use identity for direct values
      geom_point(aes(color = Clinical_Class), position = position_dodge(width = 0.9), size = 1, shape = 21, fill = "black") +  # Add points for individual reps
      labs(title = paste("Mean", cell_type, "in", analysis_region), 
           x = "Clinical Class", 
           y = "Mean Percentage of Total Cells (%)") +
      scale_fill_manual(values = clinical_class_colors) +  # Use custom red gradient colors
      scale_color_manual(values = clinical_class_colors) +  # Use the same colors for points
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  # Rotate x-axis labels 90 degrees and make them smaller
            axis.text.y = element_text(size = 10),              # Adjust y-axis label size
            panel.background = element_rect(fill = "white", color = NA),  # Set background color to white without border
            plot.background = element_rect(fill = "white", color = NA))   # Set plot background color to white without border
    
    # Save the plot in the specified directory
    ggsave(filename = file.path(output_dir, paste0(cell_type, "_", analysis_region, "_bar_plot.png")), plot = bar_plot)
  }
}



# Create the directory for saving Clinical Class CSV files (single folder)
clinical_class_dir <- "Clinical_Class_Percentage_Positive_CSV_Overall"
if (!dir.exists(clinical_class_dir)) {
  dir.create(clinical_class_dir)
}

#Save the mean percentage data to CSVs for each Clinical Class in the single folder
for (Cell_Type in unique(mean_percentage_data$Cell_Type)) {
  # Filter data for the current Clinical Class
  clinical_class_data <- mean_percentage_data[mean_percentage_data$Cell_Type == Cell_Type, ]
  
  # Sort the data by Analysis Region and then Transformation
  sorted_by_transform <- clinical_class_data %>%
    arrange(Clinical_Class)
  
  # Save the filtered data as a CSV file
  write.csv(sorted_by_transform, file = file.path(clinical_class_dir, paste0(Cell_Type, "_mean_percentage_overall.csv")), row.names = FALSE)
}
