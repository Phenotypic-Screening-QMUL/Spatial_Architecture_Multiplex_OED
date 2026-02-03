# Load necessary libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(viridis)
library(scales)
library(tidyr)

# Load the data
Data <- fread("Neighbourhood Analysis.csv")

# Select relevant columns for analysis
data_for_plot <- Data[, c(25:34)]



# Reshape data for plotting
data_for_plot_long <- data_for_plot %>%
  pivot_longer(
    cols = c("Vasculature","CD4", "CD8", "Epithelial", "Negative","Stroma", "B-Cell", "TREG"), 
    names_to = "Cell_Type", 
    values_to = "Average_Count"
  )

# Aggregate data by ClusterNN to get the average count per cell type
aggregated_data <- data_for_plot_long %>%
  group_by(ClusterNN, Cell_Type) %>%
  summarise(Average_Count = mean(Average_Count, na.rm = TRUE), .groups = 'drop')

# Define a custom color palette with specified colors for each cell type
# Define custom colors
cell_type_colors <- c(
  "Epithelial" = "#27b913",
  "Vasculature" = "#b91313",
  "CD4" = "#ff51fc",
  "CD8" = "#e0dd28",
  "Negative" = "grey",
  "Unknown" = "beige",
  "B-Cell" = "blue",
  "TREG" = "purple",
  "Stroma" = "black"
)



# Plotting
plot_avg_counts <- ggplot(aggregated_data, aes(x = factor(ClusterNN), y = Average_Count, fill = Cell_Type)) +
  geom_col(position = "stack", width = 0.7) +
  labs(x = "ClusterNN", y = "Average Count", fill = "Cell Type") +
  scale_fill_manual(values = cell_type_colors) +  # Use custom colors
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if necessary

# Print the plot
print(plot_avg_counts)

# Create the directory if it doesn't exist
if (!dir.exists("Neighbourhoods")) {
  dir.create("Neighbourhoods")
}

# Save the cluster plot
ggsave("Neighbourhoods/Cluster_Compositions_Aggregated.png", plot_avg_counts, width = 12, height = 8, dpi = 600)

# Create the directory if it doesn't exist
if (!dir.exists("Neighbourhoods")) {
  dir.create("Neighbourhoods")
}

# Save the cluster plot
ggsave("Neighbourhoods/Cluster_Compositions_Aggregated.png", plot_avg_counts, width = 12, height = 8, dpi = 600)

# Reshape data to wide format for heatmap
heatmap_data <- aggregated_data %>%
  pivot_wider(names_from = Cell_Type, values_from = Average_Count, values_fill = 0)

# Convert ClusterNN to a factor for the heatmap
heatmap_data$ClusterNN <- factor(heatmap_data$ClusterNN) 

# Convert the data into a matrix (row = ClusterNN, column = Cell Type)
heatmap_matrix <- as.matrix(heatmap_data[,-1])  # Exclude the first column (ClusterNN)

heatmap_plot <- ggplot(aggregated_data, aes(x = Cell_Type, y = factor(ClusterNN), fill = Average_Count)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C") +  # Use a perceptually uniform color scale with log transformation
  labs(x = "Cell Type", y = "ClusterNN", fill = "Average Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Increase x-axis text size
    axis.text.y = element_text(size = 14),  # Increase y-axis text size
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 12)    # Increase legend text size
  )

# Print the heatmap plot
print(heatmap_plot)

# Save the heatmap plot
ggsave("Neighbourhoods/Cluster_Compositions_Heatmap.png", heatmap_plot, width = 12, height = 8, dpi = 600)

