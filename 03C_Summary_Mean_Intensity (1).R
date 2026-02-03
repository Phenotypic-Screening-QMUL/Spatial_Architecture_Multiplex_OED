library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggsignif)

# --------------------------
# STEP 1: Pivot to long format
# --------------------------
long_df <- Intensity_Slide_Summary %>%
  pivot_longer(
    cols = starts_with("Mean_"),
    names_to = "Marker",
    values_to = "MeanIntensity"
  ) %>%
  mutate(
    Marker = gsub("^Mean_", "", Marker),
    Transformation = factor(Transformation, levels = c("Normal Control", "Cancer Control", "0", "1"))
  )

# --------------------------
# STEP 2: Manually set marker display order
# --------------------------
# Actual marker column names after gsub
# e.g., "CD20 Cell Intensity", etc.
desired_marker_order <- c(
  "PANCK Cell Intensity",
  "SMA Cell Intensity",
  "CD31 Cell Intensity",
  "CD3d Cell Intensity",
  "CD4 Cell Intensity",
  "CD8a Cell Intensity",
  "FOXP3 Nucleus Intensity",
  "CD20 Cell Intensity",
  "PD-L1 Cell Intensity"
)

# Set Marker as a factor with the desired order
long_df <- long_df %>%
  mutate(
    Marker = factor(Marker, levels = desired_marker_order)
  )

# --------------------------
# STEP 3: Tukey HSD per marker
# --------------------------
tukey_results <- long_df %>%
  group_by(Marker) %>%
  do({
    fit <- aov(MeanIntensity ~ Transformation, data = .)
    tukey <- TukeyHSD(fit, "Transformation")
    tukey_df <- as.data.frame(tukey$Transformation)
    tukey_df$Comparison <- rownames(tukey_df)
    tukey_df$Marker <- unique(.$Marker)
    tukey_df
  }) %>%
  ungroup()

# --------------------------
# STEP 4: Filter significant results
# --------------------------
sig_results <- tukey_results %>%
  filter(`p adj` <= 0.05) %>%
  mutate(
    annotation = ifelse(`p adj` < 0.001, "<0.001",
                        paste0("p = ", formatC(`p adj`, format = "f", digits = 3)))
  )

# --------------------------
# STEP 5: Create annotation dataframe
# --------------------------
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
      summarise(max_y = max(MeanIntensity, na.rm = TRUE), .groups = "drop"),
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

# --------------------------
# STEP 6: Define fill colors
# --------------------------
trans_colors <- c(
  "Normal Control" = "grey",
  "Cancer Control" = "#838383",
  "0" = "#89CFF0",
  "1" = "red"
)

# --------------------------
# STEP 7: Plot
# --------------------------
p <- ggplot(long_df, aes(Transformation, MeanIntensity, fill = Transformation)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1, colour = "black") +
  facet_wrap(~ Marker, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = trans_colors) +
  labs(
    title = "Mean Marker Intensity by Transformation Status (with Tukey p-values)",
    x = "Transformation Status",
    y = "Mean Intensity"
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
    plot.margin = ggplot2::margin(5, 5, 20, 5),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )

# --------------------------
# STEP 8: Save the plot
# --------------------------
if (!dir.exists("HALO_Percentage_Summarys")) {
  dir.create("HALO_Percentage_Summarys")
}
ggsave("HALO_Percentage_Summarys/Mean_Intensity_by_Transformation.png", plot = p, width = 12, height = 8, dpi = 300)