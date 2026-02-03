# -------------------------------
# Updated Script: Slide-level ROI mean
# -------------------------------

# Load libraries
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(ggsignif)
library(ggplot2)
library(grid)

# 1. Load data
Positive_Cell_Data <- fread("Positive_Cell_Data.csv")

# 2. Define marker columns
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

# 3. Summarise per Slide × Region × Analysis Region × Transformation (ROI-level)
Positive_Summary <- Positive_Cell_Data %>%
  group_by(Slide, Region, `Analysis Region`, Transformation) %>%
  summarise(
    Total_Count = n(),
    across(all_of(positive_cols), ~ sum(. == 1, na.rm = TRUE), .names = "Pos_{col}"),
    .groups = "drop"
  ) %>%
  mutate(
    across(starts_with("Pos_"), ~ . / Total_Count * 100, .names = "Pct_{.col}")
  )

# 4. Slide-level metadata
meta_per_slide <- Positive_Cell_Data %>%
  select(Slide, Clinical_Class) %>%
  distinct()

# -------------------------------
# NEW STEP: Average across ROIs within each Slide × Analysis Region × Transformation
# -------------------------------
Slide_ROI_Mean <- Positive_Summary %>%
  group_by(Slide, `Analysis Region`, Transformation) %>%
  summarise(
    across(starts_with("Pct_Pos_"), mean, na.rm = TRUE),
    .groups = "drop"
  )

# 5. Pivot longer for stats & plotting
long_ar <- Slide_ROI_Mean %>%
  pivot_longer(
    cols = starts_with("Pct_Pos_"),
    names_to = "Marker",
    values_to = "PercentPositive"
  ) %>%
  mutate(
    Marker         = str_remove(Marker, "^Pct_Pos_"),
    Marker         = factor(Marker, levels = positive_cols),
    AnalysisRegion = factor(`Analysis Region`, levels = c("Epi", "Stroma", "Sub001")),
    Transformation = factor(Transformation, levels = c("Normal Control", "Cancer Control", "0", "1"))
  )

# 6. ANOVA + TukeyHSD per Marker × Analysis Region
tukey_ar <- long_ar %>%
  group_by(Marker, AnalysisRegion) %>%
  do({
    fit <- aov(PercentPositive ~ Transformation, data = .)
    tk  <- TukeyHSD(fit, "Transformation")
    df  <- as.data.frame(tk$Transformation)
    df$Comparison     <- rownames(df)
    df$Marker         <- unique(.$Marker)
    df$AnalysisRegion <- unique(.$AnalysisRegion)
    df
  }) %>%
  ungroup()

# 7. Filter significant comparisons
sig_ar <- tukey_ar %>%
  filter(`p adj` <= 0.05) %>%
  mutate(
    annotation = ifelse(`p adj` < 0.001, "<0.001",
                        paste0("p = ", formatC(`p adj`, format = "f", digits = 3)))
  )

# 8. Prepare annotation dataframe for geom_signif
annotation_ar <- sig_ar %>%
  rowwise() %>%
  mutate(
    grp1 = str_split(Comparison, "-", simplify = TRUE)[1],
    grp2 = str_split(Comparison, "-", simplify = TRUE)[2]
  ) %>%
  ungroup() %>%
  mutate(
    start = as.numeric(factor(grp1, levels = levels(long_ar$Transformation))),
    end   = as.numeric(factor(grp2, levels = levels(long_ar$Transformation)))
  ) %>%
  left_join(
    long_ar %>%
      group_by(Marker, AnalysisRegion) %>%
      summarise(max_y = max(PercentPositive, na.rm = TRUE), .groups = "drop"),
    by = c("Marker", "AnalysisRegion")
  ) %>%
  group_by(Marker, AnalysisRegion) %>%
  arrange(Comparison) %>%
  mutate(
    gap        = max_y * 0.1,
    y_position = max_y + row_number() * gap
  ) %>%
  ungroup() %>%
  mutate(row_id = row_number()) %>%
  select(Marker, AnalysisRegion, start, end, y_position, annotation, row_id)

# 9. Plot: Transformation within each Analysis Region
trans_colors <- c(
  "Normal Control" = "grey70",
  "Cancer Control" = "grey40",
  "0"               = "#89CFF0",
  "1"               = "red"
)

p_ar <- ggplot(long_ar,
               aes(x = Transformation, y = PercentPositive, fill = Transformation)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, colour = "black") +
  facet_grid(AnalysisRegion ~ Marker, scales = "free_y") +
  scale_fill_manual(values = trans_colors) +
  labs(
    title = "Percent Positive Markers by Transformation\n(ROI-averaged per Slide)",
    x     = "Transformation Status",
    y     = "Percent Positive"
  ) +
  geom_signif(
    data       = annotation_ar,
    aes(xmin        = start,
        xmax        = end,
        annotations = annotation,
        y_position  = y_position,
        group       = row_id),
    manual      = TRUE,
    inherit.aes = FALSE,
    textsize    = 3.5,
    tip_length  = 0.01
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, color = "black"),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    plot.margin      = ggplot2::margin(5, 5, 20, 5),
    panel.grid       = element_blank(),
    axis.line        = element_line(color = "black"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
    panel.spacing    = unit(0.5, "lines")
  )

# 10. Print & save
print(p_ar)
if (!dir.exists("HALO_Percentage_Summarys")) dir.create("HALO_Percentage_Summarys")
ggsave("HALO_Percentage_Summarys/Percent_Positive_Trans_by_AR_ROImean.png",
       plot = p_ar, width = 16, height = 10, dpi = 300)
