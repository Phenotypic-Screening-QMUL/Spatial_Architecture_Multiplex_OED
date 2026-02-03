# =========================
# Load libraries
# =========================
library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
library(tidyr)
library(pheatmap)

# =========================
# 1. Load data
# =========================
Data <- fread("Neighbourhood Analysis.csv", encoding = "Latin-1") %>%
  mutate(
    Cluster        = as.factor(Cluster),
    ClusterNN      = as.integer(ClusterNN),
    Slide          = sub("_.*", "", Slide),
    Neighbourhood  = case_when(
      ClusterNN == 1  ~ "N1_CD8_Mixed",
      ClusterNN == 2  ~ "N2_Immune_Bndry",
      ClusterNN == 3  ~ "N3_Vasculature",
      ClusterNN == 4  ~ "N4_Negative",
      ClusterNN == 5  ~ "N5_Epi",
      ClusterNN == 6  ~ "N6_CD4",
      ClusterNN == 7  ~ "N7_Epi_Bndry",
      ClusterNN == 8  ~ "N8_Stroma",
      TRUE            ~ "Unknown"
    )
  ) %>%
  filter(Cell_Type != "Artefact")   # <-- remove Artefact cells

# =========================
# 2. Simplified metadata
# =========================
Simplified_Data <- Data %>%
  select(Region, Clinical_Class, Slide, Transformation, Cluster, Neighbourhood, Cell_Type)

# =========================
# 3. Count per Slide/Region/Neighbourhood
# =========================
Neighbourhood_Counts <- Simplified_Data %>%
  group_by(Slide, Region, Neighbourhood) %>%
  summarise(Cell_Count = n(), .groups = "drop")

# =========================
# 4. Pivot wider so each Neighbourhood is a column
# =========================
Neighbourhood_Wide <- Neighbourhood_Counts %>%
  pivot_wider(
    names_from  = Neighbourhood,
    values_from = Cell_Count,
    values_fill = 0
  )

# =========================
# 5. Compute total cells & percentages
# =========================
Neighbourhood_Wide <- Neighbourhood_Wide %>%
  mutate(Total_Cells = rowSums(select(., -Slide, -Region))) %>%
  mutate(across(
    -c(Slide, Region, Total_Cells),
    ~ .x / Total_Cells * 100,
    .names = "Pct_{.col}"
  ))

# =========================
# 6. Average percentages per slide
# =========================
Slide_Level_Percent <- Neighbourhood_Wide %>%
  group_by(Slide) %>%
  summarise(across(starts_with("Pct_"), mean, na.rm = TRUE)) %>%
  ungroup()

# =========================
# 7. Add metadata
# =========================
meta_per_slide <- Simplified_Data %>%
  select(Slide, Clinical_Class, Transformation) %>%
  distinct()

Slide_Level_Percent <- Slide_Level_Percent %>%
  left_join(meta_per_slide, by = "Slide")

# =========================
# 8. Long format for plotting
# =========================
long_df <- Slide_Level_Percent %>%
  pivot_longer(
    cols = starts_with("Pct_"),
    names_to = "Neighbourhood",
    values_to = "Percent"
  ) %>%
  mutate(
    Neighbourhood = str_remove(Neighbourhood, "^Pct_"),
    Transformation = factor(Transformation, levels = c("Normal Control", "Cancer Control", "0", "1"))
  ) %>%
  mutate(Neighbourhood = factor(Neighbourhood, levels = sort(unique(Neighbourhood))))

# =========================
# 9. ANOVA + TukeyHSD per Neighbourhood
# =========================
tukey_results <- long_df %>%
  group_by(Neighbourhood) %>%
  do({
    fit   <- aov(Percent ~ Transformation, data = .)
    tukey <- TukeyHSD(fit, "Transformation")
    df    <- as.data.frame(tukey$Transformation)
    df$Comparison    <- rownames(df)
    df$Neighbourhood <- unique(.$Neighbourhood)
    df
  }) %>%
  ungroup()

sig_results <- tukey_results %>%
  filter(`p adj` <= 0.05) %>%
  mutate(
    annotation = ifelse(`p adj` < 0.001, "<0.001",
                        paste0("p = ", formatC(`p adj`, format = "f", digits = 3)))
  )

# =========================
# 10. Build annotation df for geom_signif
# =========================
annotation_df <- sig_results %>%
  rowwise() %>%
  mutate(
    grp1 = str_split(Comparison, "-", simplify = TRUE)[1],
    grp2 = str_split(Comparison, "-", simplify = TRUE)[2]
  ) %>%
  ungroup() %>%
  mutate(
    start = as.numeric(factor(grp1, levels = levels(long_df$Transformation))),
    end   = as.numeric(factor(grp2, levels = levels(long_df$Transformation)))
  ) %>%
  left_join(
    long_df %>%
      group_by(Neighbourhood) %>%
      summarise(max_y = max(Percent, na.rm = TRUE), .groups = "drop"),
    by = "Neighbourhood"
  ) %>%
  group_by(Neighbourhood) %>%
  arrange(Comparison) %>%
  mutate(
    gap        = max_y * 0.2,
    y_position = max_y + row_number() * gap
  ) %>%
  ungroup() %>%
  select(Neighbourhood, start, end, y_position, annotation)

# =========================
# 11. Define colors and make boxplot
# =========================
trans_colors <- c(
  "Normal Control" = "grey",
  "Cancer Control" = "#838383",
  "0"               = "#89CFF0",
  "1"               = "red"
)

p <- ggplot(long_df, aes(Transformation, Percent, fill = Transformation)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, colour = "black") +
  facet_wrap(~ Neighbourhood, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = trans_colors) +
  labs(
    title = "Average % Neighbourhood Composition by Transformation (with Tukey p‑values)",
    x     = "Transformation Status",
    y     = "Average %"
  ) +
  geom_signif(
    data        = annotation_df,
    aes(xmin = start, xmax = end, annotations = annotation, y_position = y_position),
    manual      = TRUE,
    inherit.aes = FALSE,
    textsize    = 4,
    tip_length  = 0.01
  )+
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y  = element_text(color = "black"),
    axis.title   = element_text(color = "black"),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black")
  )

# =========================
# 12. Save boxplot (safe path)
# =========================
if (!dir.exists("Cluster_Percentage_Summaries")) dir.create("Cluster_Percentage_Summaries")
ggsave("Cluster_Percentage_Summaries/Percent_CN_by_Transformation.png",
       plot = p, width = 12, height = 8, dpi = 300)

# =========================
# 13. Create CellType_by_Neighbourhood for heatmap (filtered)
# =========================
CellType_by_Neighbourhood <- Data %>%
  filter(Cell_Type != "Artefact") %>%      # <-- remove Artefact cells
  group_by(Neighbourhood, Cell_Type) %>%
  summarise(Cell_Count = n(), .groups = "drop")

heatmap_df <- CellType_by_Neighbourhood %>%
  pivot_wider(names_from = Cell_Type, values_from = Cell_Count, values_fill = 0)

heatmap_matrix <- as.matrix(heatmap_df[, -1])
rownames(heatmap_matrix) <- heatmap_df$Neighbourhood

# =========================
# 14. Scale (z-score) and define colors
# =========================
scaled_matrix <- scale(heatmap_matrix, center = TRUE, scale = TRUE)
breaks <- unique(c(seq(-3, -0.01, length = 100),
                   seq(-0.01, 0.01, length = 100),
                   seq(0.01, 3, length = 100)))
my_palette <- colorRampPalette(c("#2166AC", "white", "white", "#B2182B"))(length(breaks) - 1)

# =========================
# 15. Save heatmap (safe path)
# =========================
if(!dir.exists("Neighbourhoods_Heatmap")) dir.create("Neighbourhoods_Heatmap")
output_file <- "Neighbourhoods_Heatmap/CellType_by_Neighbourhood_Heatmap.png"

pheatmap(
  t(scaled_matrix),
  color        = my_palette,
  breaks       = breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = "black",
  cellwidth    = 10,
  cellheight   = 10,
  fontsize     = 8,
  treeheight_row = 50,
  treeheight_col = 50,
  legend       = TRUE,
  angle_col    = 90,
  main         = "Raw Cell Counts \n Per Neighbourhood (Scaled)",
  filename     = output_file
)

# Save slide-level percentages to CSV
write.csv(Slide_Level_Percent, 
          "Slide_Level_Percentages_NoArtefact.csv", 
          row.names = FALSE)


# =========================
# 16. Raw cell counts per Neighbourhood, stacked by Cell_Type
# =========================
cluster_celltype_counts <- Data %>%
  filter(Cell_Type != "Artefact") %>%
  group_by(Neighbourhood, Cell_Type) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  mutate(Neighbourhood = factor(Neighbourhood, levels = sort(unique(Neighbourhood))))

stacked_plot <- ggplot(cluster_celltype_counts, aes(x = Neighbourhood, y = Cell_Count, fill = Cell_Type)) +
  geom_col() +
  labs(
    title = "Raw Cell Counts per Neighbourhood (stacked by Cell Type)",
    x = "Neighbourhood",
    y = "Raw Cell Count",
    fill = "Cell Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y  = element_text(color = "black"),
    axis.title   = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save the stacked bar plot
if (!dir.exists("Cluster_Percentage_Summaries")) dir.create("Cluster_Percentage_Summaries")
ggsave("Cluster_Percentage_Summaries/Raw_CellCounts_by_Neighbourhood_CellType.png",
       plot = stacked_plot, width = 12, height = 8, dpi = 300)


Neighbourhood_Region_Data <- Neighbourhood_Wide[,c(12:19)]

write.csv(Neighbourhood_Region_Data, "Neighbourhood_Region_Data.csv")
