# Load libraries
library(data.table)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)  # optional

# Create folder "heatmaps" if it doesn't exist
if (!dir.exists("heatmaps")) {
  dir.create("heatmaps")
}

# --- Load and prepare data ---
Heat_Map_Data <- fread("Clustered_Data_Tidy.csv", encoding = "Latin-1")
Heat_Map_Data$Cluster <- as.factor(Heat_Map_Data$Cluster)

# Select relevant intensity columns, remove one unwanted column (6th)
intensity_columns <- grep("Intensity$", colnames(Heat_Map_Data), value = TRUE)
#intensity_columns <- intensity_columns[-6]

# Calculate mean intensity per cluster
cluster_means <- Heat_Map_Data %>%
  group_by(Cluster) %>%
  summarise(across(all_of(intensity_columns), mean, na.rm = TRUE))

# Convert to matrix and label rows
heatmap_matrix <- as.matrix(cluster_means[, -1])
rownames(heatmap_matrix) <- paste0("Cluster_", cluster_means$Cluster)

# Reorder columns
new_order <- c("PANCK Cell Intensity", "SMA Cell Intensity", "CD31 Cell Intensity", 
               "CD3d Cell Intensity", "CD4 Cell Intensity", "CD8a Cell Intensity", 
               "CD20 Cell Intensity", "FOXP3 Nucleus Intensity", "PD-L1 Cell Intensity")
heatmap_matrix <- heatmap_matrix[, new_order]

# Z-score scaling by column
heatmap_matrix <- scale(heatmap_matrix)

# --- Define marker group annotations (Epithelial / Stromal / Immune) ---
marker_groups <- data.frame(
  Marker_Group = c("Epithelial", "Stromal", "Stromal",
                   "Immune", "Immune", "Immune",
                   "Immune", "Immune", "Epithelial")
)
rownames(marker_groups) <- new_order  # Row names must match markers (rows after transpose)

# Define custom colors for marker groups
annotation_colors <- list(
  Marker_Group = c(
    Epithelial = "#33a02c",
    Stromal    = "lightgoldenrod",
    Immune     = "#e31a1c"
  )
)

# --- Define heatmap color palette and breaks ---
breaks <- unique(c(seq(-2, -0.01, length = 100),
                   seq(-0.01, 0.01, length = 100),
                   seq(0.01, 2, length = 100)))
my_palette <- colorRampPalette(c("#2166AC", "white", "white", "#B2182B"))(length(breaks) - 1)

# --- Plot and save cluster heatmap PNG ---
png("heatmaps/pheatmap_cluster_biomarkers.png", width = 2000, height = 1000, res = 200)
pheatmap(
  t(heatmap_matrix),  # Transpose so markers are rows
  color             = my_palette,
  breaks            = breaks,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  border_color      = "black",
  cellwidth         = 10,
  cellheight        = 10,
  fontsize          = 8,
  treeheight_row    = 50,
  treeheight_col    = 50,
  legend            = TRUE,
  angle_col         = 90,
  main              = "Cluster Heatmap (Mean Z-Score)",
  annotation_row    = marker_groups,
  annotation_colors = annotation_colors
)
dev.off()

# --- Bar plot: number of cells per cluster ---
cluster_counts <- Heat_Map_Data %>%
  group_by(Cluster) %>%
  summarise(Cell_Count = n())
cluster_counts$Cluster <- factor(cluster_counts$Cluster)

# Save bar plot PNG in heatmaps folder
png("heatmaps/cluster_counts_barplot.png", width = 1600, height = 1200, res = 150)
ggplot(cluster_counts, aes(x = Cluster, y = Cell_Count, fill = Cluster)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  theme_minimal(base_size = 22) +
  labs(title = "Cell Counts per Cluster", x = "Cluster", y = "Cell Count") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 22, hjust = 0.5)
  )
dev.off()

# --- Cell Type heatmap ---

# Add Cell_Type column based on Cluster
Heat_Map_Data <- Heat_Map_Data %>%
  mutate(Cell_Type = case_when(
    Cluster == 1  ~ "CD8",
    Cluster == 2  ~ "CD4",
    Cluster == 3  ~ "TREG",
    Cluster == 4  ~ "PDL1 High",
    Cluster == 5  ~ "B-Cell",
    Cluster == 6  ~ "Epithelium",
    Cluster == 7  ~ "Stroma",
    Cluster == 8  ~ "CD4",
    Cluster == 9  ~ "Epithelium",
    Cluster == 10 ~ "Vasculature",
    Cluster == 11 ~ "CD8",
    Cluster == 12 ~ "Epithelium",
    Cluster == 13 ~ "CD4",
    Cluster == 14 ~ "Stroma",
    Cluster == 15 ~ "B-Cell",
    Cluster == 16 ~ "Vasculature",
    Cluster == 17 ~ "PDL1 High",
    Cluster == 18 ~ "Artefact",
    Cluster == 19 ~ "Epithelium",
    Cluster == 20 ~ "Negative",
    TRUE ~ "Unknown"
  ))

Heat_Map_Data$Cell_Type <- as.factor(Heat_Map_Data$Cell_Type)

# Calculate mean intensities per Cell_Type
celltype_means <- Heat_Map_Data %>%
  group_by(Cell_Type) %>%
  summarise(across(all_of(intensity_columns), mean, na.rm = TRUE))

# Convert to matrix, rownames = Cell_Type
heatmap_matrix <- as.matrix(celltype_means[, -1])
rownames(heatmap_matrix) <- celltype_means$Cell_Type

# Reorder columns
heatmap_matrix <- heatmap_matrix[, new_order]

# Z-score scale by column
heatmap_matrix <- scale(heatmap_matrix)

# Marker group annotation for markers (rows after transpose)
marker_groups <- data.frame(
  Marker_Group = c("Epithelial", "Stromal", "Stromal",
                   "Immune", "Immune", "Immune",
                   "Immune", "Immune", "Epithelial")
)
rownames(marker_groups) <- new_order

# Plot and save Cell Type heatmap PNG
png("heatmaps/pheatmap_celltype_biomarkers.png", width = 2000, height = 1200, res = 200)
pheatmap(
  t(heatmap_matrix),
  color             = my_palette,
  breaks            = breaks,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  border_color      = "black",
  cellwidth         = 10,
  cellheight        = 10,
  fontsize          = 8,
  treeheight_row    = 50,
  treeheight_col    = 50,
  legend            = TRUE,
  angle_col         = 90,
  main              = "Cell Type Heatmap (Mean Z-Score)",
  annotation_row    = marker_groups,
  annotation_colors = annotation_colors
)
dev.off()


# Group by Cell_Type and Transformation, count cells
celltype_counts <- Heat_Map_Data %>%
  group_by(Cell_Type, Transformation) %>%
  summarise(Cell_Count = n(), .groups = "drop")

# Optional: order Cell_Type by total counts
celltype_counts$Cell_Type <- factor(celltype_counts$Cell_Type,
                                    levels = celltype_counts %>%
                                      group_by(Cell_Type) %>%
                                      summarise(Total = sum(Cell_Count)) %>%
                                      arrange(desc(Total)) %>%
                                      pull(Cell_Type))

# Save plot to PNG in heatmaps folder
png("heatmaps/celltype_counts_stacked_barplot.png", width = 1800, height = 1200, res = 150)

ggplot(celltype_counts, aes(x = Cell_Type, y = Cell_Count, fill = Transformation)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  theme_minimal(base_size = 20) +
  labs(
    title = "Cell Counts per Cell Type (by Transformation)",
    x = "Cell Type",
    y = "Cell Count",
    fill = "Transformation"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  scale_fill_brewer(palette = "Set1")  # Or choose another palette

dev.off()
