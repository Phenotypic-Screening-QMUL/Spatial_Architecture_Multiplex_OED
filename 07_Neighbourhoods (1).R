# Load necessary libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(viridis)
library(scales)
library(umap)

# Load the dataset
Nearest_Neighbours <- fread("Nearest_Neighbours.csv")
Nearest_Neighbours <- Nearest_Neighbours [,-c(1,2)] 

# Select relevant columns 
selected_columns <- Nearest_Neighbours[, c(24:32)]  

selected_columns <- selected_columns %>% na.omit()

# Function to compute total within-cluster sum of squares
wss <- function(k) {
  kmeans(selected_columns, centers = k, nstart = 25, iter.max = 30)$tot.withinss
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

# Save the elbow plot
ggsave("Elbow_Plot_NN.png", plot = last_plot())

# Determine optimal number of clusters based on elbow method (e.g., from plot inspection)
optimal_k <- 8  # Adjust based on the elbow plot

# Perform K-means clustering with optimal number of clusters
set.seed(123)
kmeans_result <- kmeans(selected_columns, centers = optimal_k, nstart = 25, iter.max = 30)

# Add cluster labels to the original data
Nearest_Neighbours$ClusterNN <- kmeans_result$cluster

write.csv(Nearest_Neighbours,"Neighbourhood Analysis.csv")

# Perform UMAP for dimensionality reduction and visualization
umap_result <- umap(selected_columns, n_neighbors = 15, min_dist = 0.1, n_components = 2)
umap_data <- data.frame(umap_result$layout)

# Merge UMAP results with the original data
Nearest_Neighbours <- cbind(Nearest_Neighbours, umap_data)

Nearest_Neighbours <- Nearest_Neighbours[,-c(1)]
Nearest_Neighbours$ClusterNN <- as.factor(Nearest_Neighbours$ClusterNN)

# Plot UMAP results colored by Cluster and save it
umap_cluster_plot <- ggplot(Nearest_Neighbours, aes(x = X1, y = X2, colour = as.factor(ClusterNN))) +
  geom_point(size = 0.1) +
  scale_color_manual(values = custom_colors) +
  ggtitle("UMAP Visualization of Clusters") +
  xlab("UMAP Dimension 1") +
  ylab("UMAP Dimension 2") +
  theme_bw()

# Display the UMAP cluster plot
print(umap_cluster_plot)

# Save the UMAP cluster plot
ggsave("Neighbourhoods/Cluster_UMAP_Neighbours_plot.png", umap_cluster_plot, width = 8, height = 6, dpi = 300)


##################

#Skip UMAP due to time it takes to run

#Reload Data if Required
Nearest_Neighbours <- fread("Neighbourhood Analysis.csv")
Nearest_Neighbours <- Nearest_Neighbours[,-c(1)]


#  Define a custom color palette



custom_colors <- c("#e0dd28", "blue", "red", "grey", "#00BA38", "#ff51fc", "darkgreen","black")

# Define custom colors for Cell_Type
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

# Example Slides - Adjusted to filter by the correct Slide names
Slide_Examples <- Nearest_Neighbours %>% filter(Slide %in% c("92865_R0" , "92882_R1", "95470_R0","92879_R0"))


all_clusters <- sort(as.numeric(unique(Nearest_Neighbours$ClusterNN)))


##################################################################################################################



process_and_plot <- function(sample_name, all_clusters, data) {
  
  # Filter data for the specific sample
  plot_data_NN <- data %>%
    filter(Slide == sample_name)
  
  # Define the number of bins for higher resolution
  num_bins <- 20
  
  # Bin the XCentroid and YCentroid variables with higher resolution
  plot_data_NN <- plot_data_NN %>%
    mutate(
      bin_x = cut(XCentroid, breaks = num_bins, labels = FALSE),
      bin_y = cut(YCentroid, breaks = num_bins, labels = FALSE)
    )
  
  # Convert ClusterNN column to a factor with all clusters to maintain color consistency
  plot_data_NN$ClusterNN <- factor(plot_data_NN$ClusterNN, levels = all_clusters)
  
  # Plot the cluster data using geom_tile() with higher resolution
  cluster_plot_Tile <- ggplot(plot_data_NN, aes(x = bin_x, y = bin_y, fill = ClusterNN)) +
    geom_tile() +
    scale_fill_manual(values = custom_colors, drop = FALSE) +  # Use custom color palette and prevent dropping of unused levels
    labs(title = paste("Cluster Visualization - Neighbourhoods Tiled for", sample_name), 
         x = "X Centroid", y = "Y Centroid") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white", color = NA),
          axis.text = element_blank(),
          axis.title = element_text(size = 12),
          axis.ticks = element_blank(),
          legend.position = "right") 
  
  # Display and save the tile plot
  print(cluster_plot_Tile)
  ggsave(paste0("Neighbourhoods/Cluster_Visualization_Tiled_", sample_name, ".png"), 
         plot = cluster_plot_Tile, width = 10, height = 8, dpi = 300)
  
  # Plot the centroids using geom_point() (ClusterNN)
  cluster_plot <- ggplot(plot_data_NN, aes(x = XCentroid, y = YCentroid, colour = ClusterNN)) +
    geom_point(shape = ".") +
    scale_color_manual(values = custom_colors, drop = FALSE) +  # Use custom color palette and prevent dropping of unused levels
    labs(title = paste("Cell Map - Neighbourhoods for", sample_name), 
         x = "X Centroid", y = "Y Centroid") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white", color = NA),
          axis.text = element_blank(),
          axis.title = element_text(size = 12),
          axis.ticks = element_blank()) +
    scale_x_continuous(labels = NULL) +
    scale_y_continuous(labels = NULL) 
  
  # Display and save the ClusterNN point plot
  print(cluster_plot)
  ggsave(paste0("Neighbourhoods/Cluster_plot_Neighbours_", sample_name, ".png"), 
         plot = cluster_plot, width = 8, height = 6, dpi = 600)
  
  # Plot the centroids using geom_point() (Cell_Type)
  cell_type_plot <- ggplot(plot_data_NN, aes(x = XCentroid, y = YCentroid, colour = Cell_Type)) +
    geom_point(shape = ".") +
    scale_color_manual(values = cell_type_colors, drop = FALSE) +  # Use custom colors for Cell_Type
    labs(title = paste("Cell Type Visualization for", sample_name), 
         x = "X Centroid", y = "Y Centroid") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white", color = NA),
          axis.text = element_blank(),
          axis.title = element_text(size = 12),
          axis.ticks = element_blank()) +
    scale_x_continuous(labels = NULL) +
    scale_y_continuous(labels = NULL) 
  
  # Display and save the Cell_Type plot
  print(cell_type_plot)
  ggsave(paste0("Neighbourhoods/Cell_Type_plot_", sample_name, ".png"), 
         plot = cell_type_plot, width = 8, height = 6, dpi = 600)
}


# Process and plot for each of the sample names in Slide_Examples
for (sample in unique(Slide_Examples$Slide)) {
  process_and_plot(sample_name = sample, all_clusters = all_clusters, data = Slide_Examples)
}



##############

plot_cluster_nn_for_selected_slides <- function(data, output_folder = "Neighbourhoods") {
  
  # Ensure output folder exists
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Filter dataset for Analysis Region == "Sub001"
  filtered_data <- data %>% filter(`Analysis Region` == "Sub001")
  
  # Get unique slide names from the filtered dataset
  selected_slides <- unique(filtered_data$Slide)
  
  # Loop through selected slides
  for (slide in selected_slides) {
    
    slide_data <- filtered_data %>% filter(Slide == slide)
    
    # Convert ClusterNN to factor to maintain color consistency
    slide_data$ClusterNN <- factor(slide_data$ClusterNN, levels = all_clusters)
    
    # Generate plot
    cluster_plot <- ggplot(slide_data, aes(x = XCentroid, y = YCentroid, colour = ClusterNN)) +
      geom_point(shape = ".") +  # Smallest possible point size
      scale_color_manual(values = custom_colors, drop = FALSE) +  # Use custom color palette
      labs(title = paste0("Cell Map - Neighbourhoods for ", slide, " (Sub001)"), 
           x = "X Centroid", y = "Y Centroid") +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white"),
            panel.background = element_rect(fill = "white", color = NA),
            legend.background = element_rect(fill = "white", color = NA),
            legend.key = element_rect(fill = "white", color = NA),
            axis.text = element_blank(),
            axis.title = element_text(size = 12),
            axis.ticks = element_blank()) +
      scale_x_continuous(labels = NULL) +
      scale_y_continuous(labels = NULL)
    
    # Save plot
    plot_filename <- file.path(output_folder, paste0("Cluster_plot_Neighbours_", slide, "_Sub001.png"))
    ggsave(plot_filename, plot = cluster_plot, width = 8, height = 6, dpi = 600)
    
    # Print confirmation message
    message("Saved: ", plot_filename)
  }
}

# Run the function on the already-filtered dataset
plot_cluster_nn_for_selected_slides(Slide_Examples)

