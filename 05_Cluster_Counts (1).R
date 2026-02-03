# --- Script A: Load, clean, count, percent, average per slide  ----

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(stringr)

# 1. Load and clean data
Heat_Map_Data <- fread("Clustered_Data_Tidy.csv", encoding = "Latin-1") %>%
  mutate(
    Cluster       = as.factor(Cluster),
    Slide         = sub("_.*", "", Slide),       # strip "_R0" etc.
    Cell_Type     = case_when(
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

# 2. Simplified metadata-only dataframe
Simplified_Data <- Heat_Map_Data %>%
  select(Region, Clinical_Class, Slide, Transformation, Cluster, Cell_Type)

# 3. Count per Slide/Region/Cell_Type
CellType_Counts <- Simplified_Data %>%
  group_by(Slide, Region, Cell_Type) %>%
  summarise(Cell_Count = n(), .groups = "drop")

# 4. Pivot wider so each Cell_Type is a column
CellType_Wide <- CellType_Counts %>%
  pivot_wider(
    names_from  = Cell_Type,
    values_from = Cell_Count,
    values_fill = 0
  )

# 5. Compute total cells & percentages per Slide/Region
CellType_Wide <- CellType_Wide %>%
  mutate(Total_Cells = rowSums(select(., -Slide, -Region))) %>%
  mutate(across(
    -c(Slide, Region, Total_Cells),
    ~ .x / Total_Cells * 100,
    .names = "Pct_{.col}"
  ))

# 6. Average those percentages per Slide
Slide_Level_Percent <- CellType_Wide %>%
  group_by(Slide) %>%
  summarise(across(starts_with("Pct_"), mean, na.rm = TRUE)) %>%
  ungroup()

# 7. Bring back Clinical_Class & Transformation
meta_per_slide <- Simplified_Data %>%
  select(Slide, Clinical_Class, Transformation) %>%
  distinct()

Slide_Level_Percent <- Slide_Level_Percent %>%
  left_join(meta_per_slide, by = "Slide")

# --- Script B–style: reshape, stats, and plot ----

# 8. Pivot and prepare long format for plotting
long_df <- Slide_Level_Percent %>%
  pivot_longer(
    cols = starts_with("Pct_"),
    names_to = "Cell_Type",
    values_to = "Percent"
  ) %>%
  mutate(
    Cell_Type = str_remove(Cell_Type, "^Pct_"),
    Transformation = factor(Transformation)
  ) %>%
  filter(Cell_Type != "Artefact") %>%
  mutate(
    Transformation = factor(Transformation, levels = c("Normal Control", "Cancer Control", "0", "1")),
    Cell_Type = factor(Cell_Type, levels = c(
      "Epithelial", "Stroma", "Vasculature",
      "CD4", "CD8", "TREG", "B-Cell", "PDL1 High", "Negative"
    ))
  )

# 9. ANOVA + TukeyHSD per Cell_Type
tukey_results <- long_df %>%
  group_by(Cell_Type) %>%
  do({
    fit   <- aov(Percent ~ Transformation, data = .)
    tukey <- TukeyHSD(fit, "Transformation")
    df    <- as.data.frame(tukey$Transformation)
    df$Comparison <- rownames(df)
    df$Cell_Type  <- unique(.$Cell_Type)
    df
  }) %>%
  ungroup()

sig_results <- tukey_results %>%
  filter(`p adj` <= 0.05) %>%
  mutate(
    annotation = ifelse(`p adj` < 0.001, "<0.001",
                        paste0("p = ", formatC(`p adj`, format = "f", digits = 3)))
  )

# 10. Build annotation df for geom_signif
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
      group_by(Cell_Type) %>%
      summarise(max_y = max(Percent, na.rm = TRUE), .groups = "drop"),
    by = "Cell_Type"
  ) %>%
  group_by(Cell_Type) %>%
  arrange(Comparison) %>%
  mutate(
    gap        = max_y * 0.2,
    y_position = max_y + row_number() * gap
  ) %>%
  ungroup() %>%
  select(Cell_Type, start, end, y_position, annotation)

# 11. Define colors and make plot
trans_colors <- c(
  "Normal Control" = "grey",
  "Cancer Control" = "#838383",
  "0"               = "#89CFF0",
  "1"               = "red"
)

p <- ggplot(long_df, aes(Transformation, Percent, fill = Transformation)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, colour = "black") +
  facet_wrap(~ Cell_Type, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = trans_colors) +
  labs(
    title = "Average % Cluster-Based Cell Type by Transformation (with Tukey p‑values)",
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
  ) +
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

# 12. Save output
if (!dir.exists("Cluster_Percentage_Summaries")) {
  dir.create("Cluster_Percentage_Summaries")
}
ggsave("Cluster_Percentage_Summaries/Percent_CellType_by_Transformation.png",
       plot = p, width = 12, height = 8, dpi = 300)


Cluster_Cells_Region <- CellType_Wide[,c(13:21)]

write.csv(Cluster_Cells_Region, "Cluster_Cells_Region.csv")
