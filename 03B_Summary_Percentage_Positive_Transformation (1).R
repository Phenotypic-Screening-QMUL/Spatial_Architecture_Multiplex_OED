library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(stringr)
library(ggsignif)
library(grid)


# Load CSV
Positive_Cell_Data <- fread("Positive_Cell_Data.csv")

# Extract metadata
meta_per_slide <- Positive_Cell_Data %>%
  select(Slide, Clinical_Class, Transformation) %>%
  distinct()

# Identify marker columns (those containing "Positive Classification")
positive_cols <- grep("Positive Classification", colnames(Positive_Cell_Data), value = TRUE)

# Summarise per Slide and Region: count positives and calculate percentages
Positive_Summary <- Positive_Cell_Data %>%
  group_by(Slide, Region) %>%
  summarise(
    Total_Count = n(),
    across(all_of(positive_cols), ~ sum(. == 1, na.rm = TRUE), .names = "Pos_{col}")
  ) %>%
  ungroup() %>%
  mutate(across(starts_with("Pos_"), ~ . / Total_Count * 100, .names = "Pct_{.col}"))

# Average percentages per Slide
Positive_Slide_Summary <- Positive_Summary %>%
  group_by(Slide) %>%
  summarise(across(starts_with("Pct_Pos_"), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(meta_per_slide, by = "Slide")

# Pivot longer and set Marker factor with your custom order right here!
long_df <- Positive_Slide_Summary %>%
  pivot_longer(
    cols = starts_with("Pct_Pos_"),
    names_to = "Marker",
    values_to = "PercentPositive"
  ) %>%
  mutate(
    Marker = gsub("^Pct_Pos_", "", Marker)
  ) %>%
  mutate(
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
    Transformation = factor(Transformation, levels = c("Normal Control", "Cancer Control", "0", "1"))
  )

# Verify factor levels
print(levels(long_df$Marker))
print(unique(long_df$Marker))

# Run ANOVA + TukeyHSD, collect p-values
tukey_results <- long_df %>%
  group_by(Marker) %>%
  do({
    fit <- aov(PercentPositive ~ Transformation, data = .)
    tukey <- TukeyHSD(fit, "Transformation")
    tukey_df <- as.data.frame(tukey$Transformation)
    tukey_df$Comparison <- rownames(tukey_df)
    tukey_df$Marker <- unique(.$Marker)
    tukey_df
  }) %>%
  ungroup()

# Filter significant results with p adj <= 0.05 and prepare annotations
sig_results <- tukey_results %>%
  filter(`p adj` <= 0.05) %>%
  mutate(
    annotation = ifelse(`p adj` < 0.001, "<0.001",
                        paste0("p = ", formatC(`p adj`, format = "f", digits = 3)))
  )

# Prepare annotation dataframe for ggsignif manual mode
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
      group_by(Marker) %>%
      summarise(max_y = max(PercentPositive, na.rm = TRUE), .groups = "drop"),
    by = "Marker"
  ) %>%
  group_by(Marker) %>%
  arrange(Comparison) %>%
  mutate(
    gap = max_y * 0.2,
    y_position = max_y + row_number() * gap
  ) %>%
  ungroup() %>%
  select(Marker, start, end, y_position, annotation)

# Colors for transformation groups
trans_colors <- c(
  "Normal Control" = "grey",
  "Cancer Control" = "#838383",
  "0" = "#89CFF0",
  "1" = "red"
)

# Plot with custom facet order and Tukey p-values
p <- ggplot(long_df, aes(Transformation, PercentPositive, fill = Transformation)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, colour = "black") +
  facet_wrap(~ Marker, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = trans_colors) +
  labs(
    title = "Percent Positive Markers by Transformation Status (with Tukey p-values)",
    x = "Transformation Status",
    y = "Percent Positive"
  ) +
  geom_signif(
    data = annotation_df,
    aes(xmin = start, xmax = end, annotations = annotation, y_position = y_position),
    manual = TRUE,
    inherit.aes = FALSE,
    textsize = 4,
    tip_length = 0.01
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none",
    plot.margin = unit(c(5, 5, 20, 5), "pt"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )


# Create folder if needed and saveå
if (!dir.exists("HALO_Percentage_Summarys")) {
  dir.create("HALO_Percentage_Summarys")
}
ggsave("HALO_Percentage_Summarys/Percent_Positive_by_Transformation.png", plot = p, width = 12, height = 8, dpi = 300)


Percentage_Positive_Counts <- Positive_Summary [,c(13:21)]
write.csv(Percentage_Positive_Counts, "Percentage_Positive_Counts.csv")
