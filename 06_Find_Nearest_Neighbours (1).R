# Load necessary libraries - this script takes approx 2hrs to run
library(dplyr)
library(ggplot2)
library(data.table)
library(viridis)
library(scales)
library(FNN)
library(tidyr)

# Load the dataset - Note Slide = Slide + Region so FNN is within tissue sections 
Data <- fread("Clustered_Data_Tidy.csv")

# Label
Data <- Data %>%
  mutate(
    Cell_Type = case_when(
      Cluster == 1  ~ "CD8",
      Cluster == 2  ~ "CD4",
      Cluster == 3  ~ "TREG",
      Cluster == 4  ~ "PDL1 High",
      Cluster == 5  ~ "B-Cell",
      Cluster == 6  ~ "Epithelial",
      Cluster == 7  ~ "Stroma",
      Cluster == 8  ~ "CD4",
      Cluster == 9  ~ "Epithelial",
      Cluster == 10 ~ "Vasculature",
      Cluster == 11 ~ "CD8",
      Cluster == 12 ~ "Epithelial",
      Cluster == 13 ~ "CD4",
      Cluster == 14 ~ "Stroma",
      Cluster == 15 ~ "B-Cell",
      Cluster == 16 ~ "Vasculature",
      Cluster == 17 ~ "PDL1 High",
      Cluster == 18 ~ "Artefact",
      Cluster == 19 ~ "Epithelial",
      Cluster == 20 ~ "Negative",
      TRUE ~ "Unknown"
    )
  )
# Define all possible cell types in the default cluster order for ggplot colors
all_cell_types <- c("Vasculature", "CD4", "CD8", "Epithelial", "Negative", "Stroma", "B-Cell", "TREG", "PDL1 High")


# Create a function to count cell types among nearest neighbors
count_cell_types <- function(indices, Data) {
  cell_types <- Data$Cell_Type[indices]
  counts <- table(factor(cell_types, levels = all_cell_types))
  counts_df <- as.data.frame(t(as.matrix(counts)))
  colnames(counts_df) <- all_cell_types
  return(counts_df)
}

# Initialize an empty dataframe to store results
final_results <- data.frame(matrix(ncol = length(all_cell_types), nrow = 0))
colnames(final_results) <- all_cell_types


# Loop tsample()# Loop through each Slide and calculate nearest neighbors within each Slide
for(slide in unique(Data$Slide)) {
  # Subset the data for the current Slide
  Slide_data <- Data %>% filter(Slide == slide)
  
  # Extract the coordinates for the current Slide
  Slide_coords <- Slide_data %>% select(XCentroid, YCentroid)
  
  # Calculate the 20 nearest neighbors within the Slide
  Slide_nearest_neighbors <- get.knn(Slide_coords, k = 10)
  
  # Initialize a temporary dataframe to store results for this Slide
  Slide_results <- data.frame(matrix(ncol = length(all_cell_types), nrow = nrow(Slide_data)))
  colnames(Slide_results) <- all_cell_types
  
  # Populate the results for this Slide
  for (i in 1:nrow(Slide_data)) {
    indices <- Slide_nearest_neighbors$nn.index[i,]
    counts_df <- count_cell_types(indices, Slide_data)
    Slide_results[i, ] <- counts_df
  }
  
  # Combine the results with the original Slide data using cbind
  Slide_combined <- cbind(Slide_data, Slide_results)
  
  # Append the results to the final dataframe
  final_results <- rbind(final_results, Slide_combined)
}

# View the final combined dataframe
print(final_results)

# Save the final dataframe
write.csv(final_results, "Nearest_Neighbours_with_Valid_Slides.csv")

# Calculate average counts of each cell type among nearest neighbors
avg_counts <- final_results %>%
  group_by(Cell_Type) %>%
  summarise(across(all_of(all_cell_types), mean, na.rm = TRUE))

# Rename Cell_Type column temporarily for reshaping
names(avg_counts)[1] <- "Reference_Cell"

# Reshape the data for plotting
avg_counts_long <- avg_counts %>%
  pivot_longer(cols = -Reference_Cell, names_to = "Cell_Type", values_to = "Average_Count") %>%
  mutate(Cell_Type = factor(Cell_Type, levels = unique(Cell_Type)))  # Ensure correct factor levels

# Define custom colors
cell_type_colors <- c(
  "Epithelial" = "#27b913",
  "Vasculature" = "#b91313",
  "CD4" = "#ff51fc",
  "CD8" = "#e0dd28",
  "Negative" = "grey",
  "PDL1 High" = "beige",
  "B-Cell" = "blue",
  "TREG" = "purple",
  "Stroma" = "black"
)


# Use the custom color palette in the plot
Neighbours_plot <- ggplot(avg_counts_long, aes(x = Cell_Type, y = Average_Count, fill = Cell_Type)) +
  geom_col() +
  facet_wrap(~ Reference_Cell, scales = "free_y", ncol = 2) +
  labs(x = "Cell Type", y = "Average Count among Nearest Neighbors") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = cell_type_colors)  # Apply custom colors

# Save the plot
ggsave("NN/Neighbours_plot.png", Neighbours_plot, width = 8, height = 6, dpi = 600)

# Save the final dataframe to a CSV file
write.csv(final_results, "Nearest_Neighbours.csv")

