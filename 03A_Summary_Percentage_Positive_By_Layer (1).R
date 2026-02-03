# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(stringr)
library(ggsignif)

# 1. Load data
Positive_Cell_Data <- fread("Positive_Cell_Data.csv")

# 2. Identify marker columns
positive_cols <- c(
  "PANCK Positive Classification",
  "SMA Positive Classification",
  "CD31 Positive Classification",
  "CD3d Positive Classification",
  "CD4 Positive Classification",
  "CD8a Positive Classification",
  "FOXP3 Positive Classification",
  "CD20 Positive Classification",
  "PD-L1 Positive Classification"
)

# 3. Summarise per Slide and AnalysisRegion: count positives and calculate percentages
Positive_Summary <- Positive_Cell_Data %>%
  group_by(Slide, `Analysis Region`) %>%
  summarise(
    Total_Count = n(),
    across(all_of(positive_cols), ~ sum(. == 1, na.rm = TRUE), .names = "Pos_{col}"),
    .groups = 'drop'
  ) %>%
  mutate(across(starts_with("Pos_"), ~ . / Total_Count * 100, .names = "Pct_{col}"))

# 4. Add slide-level metadata
meta <- Positive_Cell_Data %>%
  select(Slide, Transformation, Clinical_Class) %>%
  distinct()

df <- Positive_Summary %>%
  left_join(meta, by = "Slide")

# 5. Pivot longer and set explicit Marker order
long_df <- df %>%
  pivot_longer(
    cols = starts_with("Pct_Pos_"),
    names_to = "Marker",
    values_to = "PercentPositive"
  ) %>%
  mutate(
    Marker = str_remove(Marker, "^Pct_Pos_"),
    Marker = factor(Marker, levels = c(
      "PANCK Positive Classification",
      "SMA Positive Classification",
      "CD31 Positive Classification",
      "CD3d Positive Classification",
      "CD4 Positive Classification",
      "CD8a Positive Classification",
      "FOXP3 Positive Classification",
      "CD20 Positive Classification",
      "PD-L1 Positive Classification"
    )),
    AnalysisRegion = factor(`Analysis Region`, levels = c("Epi", "Stroma", "Sub001"))
  )

# 6. ANOVA + TukeyHSD
tukey_results <- long_df %>%
  group_by(Marker) %>%
  do({
    fit <- aov(PercentPositive ~ AnalysisRegion, data = .)
    tuk <- TukeyHSD(fit, "AnalysisRegion")
    res <- as.data.frame(tuk$AnalysisRegion)
    res$Comparison <- rownames(res)
    res$Marker <- unique(.$Marker)
    res
  }) %>%
  ungroup()

# 7. Filter and annotate p-values
sig_results <- tukey_results %>%
  filter(`p adj` <= 0.05) %>%
  mutate(
    annotation = ifelse(`p adj` < 0.001, "<0.001",
                        paste0("p = ", formatC(`p adj`, format = "f", digits = 3)))
  )

# 8. Build annotation dataframe
annotation_df <- sig_results %>%
  rowwise() %>%
  mutate(
    grp1 = str_split(Comparison, "-", simplify = TRUE)[1],
    grp2 = str_split(Comparison, "-", simplify = TRUE)[2]
  ) %>%
  ungroup() %>%
  mutate(
    start = as.numeric(factor(grp1, levels = levels(long_df$AnalysisRegion))),
    end   = as.numeric(factor(grp2, levels = levels(long_df$AnalysisRegion)))
  ) %>%
  left_join(
    long_df %>% group_by(Marker) %>% summarise(max_y = max(PercentPositive, na.rm = TRUE), .groups = 'drop'),
    by = "Marker"
  ) %>%
  group_by(Marker) %>%
  arrange(Comparison) %>%
  mutate(
    gap        = max_y * 0.2,
    y_position = max_y + row_number() * gap
  ) %>%
  ungroup() %>%
  mutate(row_id = row_number()) %>%
  select(Marker, start, end, y_position, annotation, row_id)

# 9. Define colors for regions
region_colors <- c(
  "Epi"    = "darkgreen",
  "Stroma" = "darkgoldenrod1",
  "Sub001" = "darkred"
)

# 10. Plot
p <- ggplot(long_df, aes(x = AnalysisRegion, y = PercentPositive, fill = AnalysisRegion)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, colour = "black") +
  facet_wrap(~ Marker, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = region_colors) +
  labs(
    title = "Percent Positive Markers by Tissue Layer (with Tukey p-values)",
    x = "Tissue Layer",
    y = "Percent Positive"
  ) +
  geom_signif(
    data = annotation_df,
    aes(
      xmin = start,
      xmax = end,
      annotations = annotation,
      y_position = y_position,
      group = row_id
    ),
    manual = TRUE,
    inherit.aes = FALSE,
    textsize = 4,
    tip_length = 0.01
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    legend.position = "none",
    plot.margin = margin(5, 5, 20, 5),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Print and save
print(p)
ggsave("HALO_Percentage_Summarys/Percent_Positive_by_AnalysisRegion.png", plot = p, width = 12, height = 8, dpi = 300)
