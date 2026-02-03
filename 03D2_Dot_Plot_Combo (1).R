# Load libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

# 1. Compute mean % positive for Transformation = 0 and 1,
#    then calculate the difference (1 − 0)
dot_df <- long_ar %>%
  filter(Transformation %in% c("0", "1")) %>%  # keep only 0 and 1
  group_by(AnalysisRegion, Marker, Transformation) %>%
  summarise(meanPct = mean(PercentPositive, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Transformation, values_from = meanPct, names_prefix = "T") %>%
  mutate(
    diffPct = T1 - T0   # difference = mean(1) − mean(0)
  )

# 2. Extract the 0–1 Tukey comparison p-values
p_0_1 <- tukey_ar %>%
  rowwise() %>%
  mutate(
    grp1 = str_split(Comparison, "-", simplify = TRUE)[1],
    grp2 = str_split(Comparison, "-", simplify = TRUE)[2]
  ) %>%
  ungroup() %>%
  filter((grp1 == "0" & grp2 == "1") | (grp1 == "1" & grp2 == "0")) %>%
  select(Marker, AnalysisRegion, p.adj = `p adj`)

# 3. Join the p-values to the difference data
dot_df2 <- dot_df %>%
  left_join(p_0_1, by = c("AnalysisRegion", "Marker")) %>%
  mutate(
    p.adj       = ifelse(is.na(p.adj), 1, p.adj),
    neg_log_p   = -log10(p.adj),
    significant = p.adj <= 0.05,
    Marker      = str_remove(Marker, " Positive Classification$"),
    Marker      = factor(Marker, levels = c("PANCK", "SMA", "CD31", "CD3d",
                                            "CD4", "CD8a", "FOXP3", "CD20", "PD-L1"))
  )

# 4. Dot plot:
#    - Bubble size = difference (abs value of 1 − 0)
#    - Bubble colour = significance level (-log10 p)
#    - Black outline = significant
p_dot <- ggplot(dot_df2, aes(x = Marker, y = AnalysisRegion)) +
  geom_point(aes(size = abs(diffPct), colour = neg_log_p), alpha = 0.9) +
  
  # Add black outline for significant comparisons
  geom_point(
    data = subset(dot_df2, significant),
    aes(size = abs(diffPct)),
    shape = 21, stroke = 2, colour = "black"
  ) +
  
  scale_size_continuous(
    range = c(4, 16),
    name  = "Δ in % Positive\n(|1 − 0|)"
  ) +
  
  scale_colour_gradient(
    low = "grey85", high = "red",
    name = expression(-log[10](italic(p))),
    limits = c(0, max(dot_df2$neg_log_p, na.rm = TRUE)),
    guide = guide_colourbar(barwidth = 1.2, barheight = 10)
  ) +
  
  labs(
    title = "Cell Type Composition by Analysis Region and Transformation",
    x     = "Marker",
    y     = "Analysis Region"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                    size = 14, face = "bold", color = "black"),
    axis.text.y      = element_text(size = 14, face = "bold", color = "black"),
    axis.title       = element_text(size = 16, face = "bold", color = "black"),
    legend.position  = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line        = element_line(color = "black", size = 0.5)
  )

# 5. Display and save
print(p_dot)
ggsave("Mean_Positive_DotPlot_Diff.png", p_dot, width = 10, height = 6, dpi = 300)
