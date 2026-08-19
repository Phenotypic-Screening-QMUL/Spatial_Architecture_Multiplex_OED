
# Exploratory PCA Only 


# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble) 


# 1. Load and Prepare Slide Summary

Data <- fread("Neighbourhood Analysis.csv", encoding = "Latin-1") 

Slide_Summary <- Data %>%
  group_by(Slide, Region, Transformation, Clinical_Class) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  select(-V1, -ClusterNN)

# Remove unwanted columns by index
Slide_Summary <- Slide_Summary[, -c(14, 15:29)]

# Separate meta and intensity
Meta_Data      <- Slide_Summary[, c(1:4)]
Intensity_Data <- Slide_Summary[, c(5:13)]


# 2. Load Other Tables

HALO_Data     <- fread("Percentage_Positive_Counts.csv", encoding = "Latin-1") %>% select(-V1)
Cluster_Data  <- fread("Cluster_Cells_Region.csv", encoding = "Latin-1") %>% select(-V1)
CNN_Data      <- fread("Neighbourhood_Region_Data.csv", encoding = "Latin-1") %>% select(-V1)
Distance_Data <- fread("Distance_DF.csv", encoding = "Latin-1") %>% select(-Slide)


# 3. Combine All Data

Final_Data <- cbind(
  Meta_Data,
  Intensity_Data,
  HALO_Data,
  Cluster_Data,
  CNN_Data,
  Distance_Data
)

# Filter Transformation to 0 and 1
Final_Data <- Final_Data %>% filter(Transformation %in% c(0, 1))


# 4. Aggregate to Slide-Level Means

All_Means <- Final_Data %>%
  mutate(Slide = sub("_.*", "", Slide)) %>%
  group_by(Slide) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    across(where(~ !is.numeric(.)), first),
    .groups = "drop"
  ) %>%
  mutate(Clinical_Class = as.factor(Clinical_Class))


# 5. Prepare predictors (cleaned)

meta <- All_Means %>% select(Slide, Transformation)
preds <- All_Means %>% select(-Slide, -Region, -Clinical_Class, -Transformation)

# Drop zero-variance columns
vars_var  <- sapply(preds, function(x) var(x, na.rm = TRUE))
keep_cols <- which(!is.na(vars_var) & vars_var > 0)
preds_clean <- preds[, keep_cols]

# Drop rows with any NA
keep_rows  <- complete.cases(preds_clean)
preds_clean <- preds_clean[keep_rows, , drop = FALSE]
meta_clean  <- meta[keep_rows, ]


# 6. PCA on full dataset (for visualization & loadings only)

pca_res_viz <- prcomp(preds_clean, scale. = TRUE)
pc_scores_viz <- as.data.frame(pca_res_viz$x[, c(1:4)])

viz_df <- bind_cols(
  meta_clean %>% mutate(
    Transformation = factor(Transformation, levels = c(0, 1), labels = c("Class0","Class1"))
  ),
  pc_scores_viz
)


# 7. PCA Scatter Plot: PC1 vs PC3

pc_plot <- ggplot(viz_df, aes(x = PC1, y = PC3, color = Transformation)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, aes(fill = Transformation), geom = "polygon", alpha = 0.2, color = NA) +
  scale_color_manual(values = c("Class1" = "red", "Class0" = "#89CFF0")) +
  scale_fill_manual(values = c("Class1" = "red", "Class0" = "#89CFF0")) +
  labs(
    title = "PCA: PC1 vs PC3",
    x = "PC1",
    y = "PC3",
    color = "Class",
    fill = "Class"
  ) +
  theme_minimal(base_size = 14)

print(pc_plot)
ggsave("Plots/PCA_PC1_PC3.png", plot = pc_plot, width = 8, height = 6, dpi = 300)


# 8. PC1 Loadings (Top 10 Contributing Features)

loadings_df <- pca_res_viz$rotation %>%
  as.data.frame() %>%
  rownames_to_column("Feature")

pc1_loadings <- loadings_df %>%
  transmute(
    Feature    = Feature,
    Loading    = PC1,
    AbsLoading = abs(PC1)
  )

Top10_PC1_Features <- pc1_loadings %>%
  arrange(desc(AbsLoading)) %>%
  slice(1:10)

write.csv(
  Top10_PC1_Features,
  "Top10_PC1_Features.csv",
  row.names = FALSE
)

print(Top10_PC1_Features)

cat("\nExploratory PCA complete. Plots and loadings saved.\n")


# 9. PC1 Loadings Bar Plot (Top 20 Features)

top20_loadings <- pc1_loadings %>%
  arrange(desc(AbsLoading)) %>%
  slice(1:10) %>%
  mutate(
    Feature = factor(Feature, levels = rev(Feature))
  )

pc1_bar_plot <- ggplot(top20_loadings, aes(x = Feature, y = Loading, fill = Loading > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "#D73027", "FALSE" = "#4575B4"),
    guide = "none"
  ) +
  labs(
    title = "Top 10 Features Contributing to PC1",
    x = "Feature",
    y = "PC1 Loading"
  ) +
  theme_minimal(base_size = 14)

print(pc1_bar_plot)

ggsave("Plots/PCA_PC1_Top20_Features.png", plot = pc1_bar_plot, width = 10, height = 8, dpi = 300)

cat("\nExploratory PCA complete. Plots and loadings saved.\n")
