# =========================
# Load libraries
# =========================
library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
library(tidyr)

# =========================
# 1. Load and preprocess data
# =========================
Data <- fread("Neighbourhood Analysis.csv", encoding = "Latin-1") %>%
  mutate(
    Cluster       = as.factor(Cluster),
    ClusterNN     = as.integer(ClusterNN),
    Slide         = sub("_.*", "", Slide),
    AnalysisRegion= `Analysis Region`,
    Neighbourhood = case_when(
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
  filter(Cell_Type != "Artefact")   # Remove Artefact

# =========================
# 2. Prepare simplified dataframe
# =========================
Simplified_Data <- Data %>%
  select(Slide, Region, AnalysisRegion, Transformation, Neighbourhood, Cell_Type)

# =========================
# 3. Count per ROI (Slide × Region × AnalysisRegion × Neighbourhood)
# =========================
ROI_Counts <- Simplified_Data %>%
  group_by(Slide, Region, AnalysisRegion, Neighbourhood) %>%
  summarise(Cell_Count = n(), .groups = "drop")

# =========================
# 4. Pivot wider and compute percentages per ROI
# =========================
ROI_Wide <- ROI_Counts %>%
  pivot_wider(
    names_from  = Neighbourhood,
    values_from = Cell_Count,
    values_fill = 0
  ) %>%
  mutate(Total_Cells = rowSums(select(., -Slide, -Region, -AnalysisRegion))) %>%
  mutate(across(
    -c(Slide, Region, AnalysisRegion, Total_Cells),
    ~ .x / Total_Cells * 100,
    .names = "Pct_{.col}"
  ))

# =========================
# 5. Average across ROIs to get Slide-level percentages per AnalysisRegion
# =========================
Slide_Level_Percent <- ROI_Wide %>%
  group_by(Slide, AnalysisRegion) %>%
  summarise(across(starts_with("Pct_"), mean, na.rm = TRUE), .groups = "drop")

# =========================
# 6. Add Transformation metadata
# =========================
meta_per_slide <- Simplified_Data %>%
  select(Slide, Transformation) %>%
  distinct()

Slide_Level_Percent <- Slide_Level_Percent %>%
  left_join(meta_per_slide, by = "Slide")

# =========================
# 7. Long format for plotting
# =========================
long_df <- Slide_Level_Percent %>%
  pivot_longer(
    cols = starts_with("Pct_"),
    names_to = "Neighbourhood",
    values_to = "Percent"
  ) %>%
  mutate(
    Neighbourhood   = str_remove(Neighbourhood, "^Pct_"),
    Transformation  = factor(Transformation, levels = c("0", "1")),
    Neighbourhood   = factor(Neighbourhood, levels = sort(unique(Neighbourhood))),
    AnalysisRegion  = factor(AnalysisRegion, levels = c("Epi", "Stroma", "Sub001"))
  )

# =========================
# 8. Compute mean difference per Neighbourhood × AnalysisRegion
# =========================
diff_df <- long_df %>%
  group_by(AnalysisRegion, Neighbourhood) %>%
  summarise(
    mean_0    = mean(Percent[Transformation=="0"], na.rm = TRUE),
    mean_1    = mean(Percent[Transformation=="1"], na.rm = TRUE),
    diff      = mean_1 - mean_0,
    size      = abs(diff),
    direction = ifelse(diff > 0, "Up", "Down"),
    .groups = "drop"
  )

# =========================
# 9. T-test for significance (slide-level)
# =========================
stat_summary <- long_df %>%
  group_by(AnalysisRegion, Neighbourhood) %>%
  summarise(
    p.value = tryCatch(t.test(Percent ~ Transformation)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  )

diff_df <- diff_df %>%
  left_join(stat_summary, by = c("AnalysisRegion","Neighbourhood")) %>%
  mutate(sig = !is.na(p.value) & p.value < 0.05)

# =========================
# 10. Dotplot: size = abs difference, color = direction, outline = significance
# =========================
p_dot <- ggplot(diff_df, aes(x = Neighbourhood, y = AnalysisRegion)) +
  # Non-significant points
  geom_point(aes(size = size, fill = direction),
             shape = 21, color = "grey70", stroke = 0.4) +
  # Significant points with thicker outlines
  geom_point(data = filter(diff_df, sig),
             aes(size = size, fill = direction),
             shape = 21, color = "black", stroke = 1.2) +
  scale_fill_manual(values = c(Up = "red", Down = "#89CFF0")) +
  scale_size_continuous(range = c(4, 16), name = "Δ in % Positive\n(|1 − 0|)") +
  labs(
    title = "Change in Neighbourhood Composition with Transformation",
    x     = "Neighbourhood",
    y     = "Analysis Region",
    fill  = "Direction"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y      = element_text(size = 14),
    axis.title       = element_text(size = 16, face = "bold"),
    plot.title       = element_text(size = 18, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# =========================
# 11. Save plot
# =========================
if(!dir.exists("Neighbourhoods")) dir.create("Neighbourhoods")
ggsave("Neighbourhoods/Dotplot_Transformation_0_vs_1_ROImean.png", p_dot, width = 12, height = 8, dpi = 300)
