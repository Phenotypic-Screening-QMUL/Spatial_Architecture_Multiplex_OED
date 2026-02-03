# Useful Packages
library("gplots")       # Graphical plots
library("RColorBrewer") # Color palettes
library("matrixStats")  # Efficient matrix computations
library("plyr")        # Data manipulation
library("dplyr")       # Data manipulation
library("data.table")  # High-performance data manipulation 
library("stringr")     # String manipulation
library("ggplot2")     # Elegant plots
library(dplyr)
library(ggplot2)
library(data.table)
library(viridis)
library(scales)
library(FNN)
library(tidyr)

Data_Intensity <- fread("Data_Intensity_Scaled.csv", encoding = "Latin-1")

Data_Intensity <- Data_Intensity[,-c(1,2)]


#Data_Intensity <- Data_Intensity %>% filter(`Analysis Region` == "Sub001")

  
Meta_Data <- Data_Intensity[,c(1,11:22)]


Data_Intensity <- Data_Intensity[,-c(1,11:22)]

library(umap)

# Set up the main clustering and dimensionality reduction for the data
Data_For_Clustering <- Data_Intensity

# Check for missing values
if(any(is.na(Data_For_Clustering))) {
  Data_For_Clustering <- na.omit(Data_For_Clustering)
}

# Select only numeric columns
numeric_data <- Data_For_Clustering  %>% 
  select_if(is.numeric)

#Scale the data 
scaled_data <- scale(numeric_data)

# Elbow method for optimal clusters 
wss <- function(k) {
  kmeans(scaled_data, centers = k, nstart = 25 ,iter.max=30)$tot.withinss
}

# Compute and plot WSS for different values of k
k_values <- 1:15
wss_values <- sapply(k_values, wss)

# Create the elbow plot
elbow_plot <- data.frame(k = k_values, wss = wss_values)

ggplot(elbow_plot, aes(x = k, y = wss)) +
  geom_line() +
  geom_point() +
  ggtitle("Elbow Method for Determining Optimal Number of Clusters") +
  xlab("Number of Clusters (k)") +
  ylab("Total Within-Cluster Sum of Squares (WSS)")

ggsave("Elbow_Plot_Cell_Types.png", plot = last_plot())

# Perform K-means clustering with 12 clusters
set.seed(123)
kmeans_result <- kmeans(scaled_data, centers = 20, nstart = 25, iter.max=30)

# Add cluster labels to the original data
Data_Intensity$Cluster <- kmeans_result$cluster
Data_Intensity <- cbind(Data_Intensity,Meta_Data)

write.csv(Data_Intensity, "Clustered_Data_Tidy.csv")

Data_Intensity <- fread("Clustered_Data_Tidy.csv", encoding = "Latin-1")

# UMAP dimensionality reduction
umap_result <- umap(scaled_data, n_neighbors = 15, min_dist = 0.1, n_components = 2)
umap_data <- data.frame(umap_result$layout)

# Merge UMAP results with the original data
umap_data_merged <- cbind(Data_Intensity, umap_data)
umap_data_merged$Cluster <- as.factor(umap_data_merged$Cluster)


# Plot UMAP results colored by cluster
cluster_plot <- ggplot(umap_data_merged, aes(x = X1, y = X2, color = Cluster)) +
  geom_point(shape = ".") +
  ggtitle("UMAP Visualization of Clusters") +
  xlab("UMAP Dimension 1") +
  ylab("UMAP Dimension 2") +
  theme_bw() + xlim (-25,25) + ylim(-20,20)

# Display the cluster plot
print(cluster_plot)

# Intensity-based UMAP Plots
create_and_save_umap_plots <- function(data_merged, variables, folder_path) {
  # Create folder if it doesn't exist
  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
  }
  
  for (variable in variables) {
    sanitized_variable <- gsub("[^[:alnum:]_]", "_", variable) # Replace problematic characters with '_'
    
    plot <- ggplot(data_merged, aes(x = X1, y = X2)) +
      geom_point(shape = ".", aes(color = !!sym(variable))) +
      scale_color_viridis(trans = "log10") +  # Use viridis color scale with log transformation
      ggtitle(paste("UMAP Visualization Colored by", variable)) +
      xlab("UMAP Dimension 1") +
      ylab("UMAP Dimension 2") +
      theme_bw() + xlim (-25,25) + ylim(-20,20)
    
    # Save the plot
    ggsave(paste0(folder_path, "/", sanitized_variable, "_UMAP_plot.png"), plot)
  }
}

# Create intensity plots
intensity_variables <- colnames(numeric_data)
intensity_folder_path <- "UMAP_plots_intensity"
create_and_save_umap_plots(umap_data_merged, intensity_variables, intensity_folder_path)

# Save the intensity cluster plot
ggsave("UMAP_plots_intensity/Cluster_UMAP_plot.png", cluster_plot)

##

library(ggplot2)
library(dplyr)

# Update cell types based on clusters
Update_Cell_Types <- umap_data_merged %>%  
  mutate(Cell_Type = case_when(
    Cluster == 1  ~ "B-Cell",
    Cluster == 2  ~ "Negative",
    Cluster == 3  ~ "Epithelial",
    Cluster == 4  ~ "Stroma",
    Cluster == 5  ~ "Epithelial",
    Cluster == 6  ~ "Vasculature",
    Cluster == 7  ~ "Stroma",
    Cluster == 8  ~ "PDL1 High Epithelial",
    Cluster == 9  ~ "Artefact",
    Cluster == 10 ~ "Vasculature",
    Cluster == 11 ~ "CD8",
    Cluster == 12 ~ "CD4",
    Cluster == 13 ~ "CD8",
    Cluster == 14 ~ "TREG",
    Cluster == 15 ~ "CD4",
    Cluster == 16 ~ "CD8",
    Cluster == 17 ~ "TREG",
    Cluster == 18 ~ "Epithelial",
    Cluster == 19 ~ "CD4",
    Cluster == 20 ~ "Epithelial",
    TRUE ~ "Unknown"
  ))

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

# Create UMAP plot
plot <- ggplot(Update_Cell_Types, aes(x = X1, y = X2)) +
  geom_point(shape = ".", aes(color = Cell_Type)) +
  scale_color_manual(values = cell_type_colors) +  # Use custom colors
  ggtitle("UMAP Visualization Colored by New Cell Type") +
  xlab("UMAP Dimension 1") +
  ylab("UMAP Dimension 2") +
  theme_bw() +
  xlim(-25, 25) + 
  ylim(-20, 20)

# Save plot to file
ggsave("Merged_Cell_Types_UMAP_plot.png", plot, width = 6, height = 4, dpi = 300)


