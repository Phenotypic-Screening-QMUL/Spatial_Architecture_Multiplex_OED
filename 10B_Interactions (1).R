# --- Load Libraries ---
library(dplyr)
library(tidyr)
library(purrr)
library(data.table)  # for fread if CSV is large

# --- Load the Data ---
Data <- fread("Neighbourhood Analysis.csv", encoding = "Latin-1") %>%
  mutate(
    Transformation = factor(Transformation, levels = c("Normal Control", "Cancer Control", "0", "1")),
    ClusterNN = factor(ClusterNN),
    XMicron = XCentroid / 3,
    YMicron = YCentroid / 3
  )

# --- Define neighbor count columns ---
cell_types <- c("Vasculature","CD4","CD8","Epithelial","Negative","Stroma","B-Cell","TREG")

# --- Prepare slide-level data ---
slide_data <- Data %>%
  select(Slide, Cell_Type, all_of(cell_types))

# --- Function to compute observed counts per slide ---
compute_observed_counts <- function(df, cell_types) {
  df %>%
    group_by(Cell_Type) %>%
    summarise(across(all_of(cell_types), sum), .groups = "drop") %>%
    pivot_longer(cols = all_of(cell_types), names_to = "Neighbor_Type", values_to = "Observed_Count")
}

# --- Permutation function ---
permute_counts <- function(df, cell_types, n_perm = 1000) {
  # Store results
  perm_results <- vector("list", length = n_perm)
  
  for (i in 1:n_perm) {
    # Shuffle cell types
    df_shuffled <- df %>% mutate(Cell_Type = sample(Cell_Type))
    
    # Compute counts for shuffled labels
    counts <- compute_observed_counts(df_shuffled, cell_types)
    counts$Perm <- i
    perm_results[[i]] <- counts
  }
  
  bind_rows(perm_results)
}

# --- Function to compute enrichment/depletion P-values ---
compute_p_values <- function(observed_df, perm_df) {
  observed_df %>%
    rowwise() %>%
    mutate(
      p_enrich  = mean(perm_df$Observed_Count[perm_df$Cell_Type == Cell_Type &
                                                perm_df$Neighbor_Type == Neighbor_Type] >= Observed_Count),
      p_deplete = mean(perm_df$Observed_Count[perm_df$Cell_Type == Cell_Type &
                                                perm_df$Neighbor_Type == Neighbor_Type] <= Observed_Count)
    ) %>%
    ungroup()
}

# --- Run per slide ---
slide_list <- unique(slide_data$Slide)

all_p_values <- list()

for (slide in slide_list) {
  cat("Processing slide:", slide, "\n")
  df_slide <- slide_data %>% filter(Slide == slide)
  
  # Observed counts
  obs <- compute_observed_counts(df_slide, cell_types)
  
  # Permutations
  perms <- permute_counts(df_slide, cell_types, n_perm = 1000)  # can increase to 1000 if memory allows
  
  # Compute p-values
  pvals <- compute_p_values(obs, perms)
  pvals$Slide <- slide
  all_p_values[[slide]] <- pvals
}

# Combine all slides
p_values_df <- bind_rows(all_p_values)

# --- View results ---
head(p_values_df)

# Extract slide -> transformation mapping
slide_trans <- Data %>%
  select(Slide, Transformation) %>%
  distinct()

# Add Transformation to the p-values dataframe
p_values_with_trans <- p_values_df %>%
  left_join(slide_trans, by = "Slide")

# View
head(p_values_with_trans)

# --- Set significance threshold ---
sig_thresh <- 0.05

# --- Assign interaction scores and keep Transformation ---
p_values_scored <- p_values_df %>%
  left_join(
    Data %>% select(Slide, Transformation) %>% distinct(),
    by = "Slide"
  ) %>%
  mutate(
    Interaction_Score = case_when(
      p_enrich < sig_thresh  ~ 1,   # Interaction
      p_deplete < sig_thresh ~ -1,  # Avoidance
      TRUE                    ~ 0   # Neutral
    )
  )

# --- Remove artefact comparisons ---
p_values_scored <- p_values_scored %>%
  filter(Cell_Type != "artefact", Neighbor_Type != "artefact")


# --- Optional: view results ---
head(p_values_scored)

# --- Aggregate interaction scores per transformation ---
interaction_by_trans <- p_values_scored %>%
  group_by(Cell_Type, Neighbor_Type, Transformation) %>%
  summarise(
    Mean_Interaction = mean(Interaction_Score, na.rm = TRUE),
    .groups = "drop"
  )

# --- Prepare summary for plotting ---
interaction_summary <- p_values_scored %>%
  group_by(Cell_Type, Neighbor_Type, Transformation) %>%
  summarise(
    Mean_Score = mean(Interaction_Score, na.rm = TRUE),
    .groups = "drop"
  )

interaction_plot_df <- interaction_summary %>%
  filter(
    Transformation %in% c("0", "1"),
    Cell_Type != Neighbor_Type,
    Cell_Type != "Artefact",
    Neighbor_Type != "Artefact"
  ) %>%
  mutate(Pair = paste(Cell_Type, "→", Neighbor_Type))


# --- Prepare wide table with signed differences ---
interaction_wide_signed <- interaction_plot_df %>%
  pivot_wider(
    names_from = Transformation,
    values_from = Mean_Score,
    names_prefix = "Trans_"
  ) %>%
  mutate(
    Score_Diff = Trans_1 - Trans_0,   # signed difference
    Pair = paste(Cell_Type, "→", Neighbor_Type)
  )

# --- Top 10 up in Trans_1 ---
top_up <- interaction_wide_signed %>%
  arrange(desc(Score_Diff)) %>%
  slice_head(n = 10) %>%
  pivot_longer(
    cols = starts_with("Trans_"),
    names_to = "Transformation",
    values_to = "Mean_Score"
  ) %>%
  mutate(Transformation = ifelse(Transformation == "Trans_0", "0", "1"))

# --- Reorder top_up pairs by Trans_1 score (unique only) ---
top_up <- top_up %>%
  mutate(Pair = factor(Pair, levels = top_up %>% 
                         filter(Transformation == "1") %>% 
                         arrange(Mean_Score) %>% 
                         pull(Pair) %>% unique()))

# --- Top 10 down in Trans_1 ---
top_down <- interaction_wide_signed %>%
  arrange(Score_Diff) %>%
  slice_head(n = 10) %>%
  pivot_longer(
    cols = starts_with("Trans_"),
    names_to = "Transformation",
    values_to = "Mean_Score"
  ) %>%
  mutate(Transformation = ifelse(Transformation == "Trans_0", "0", "1"))

# --- Reorder top_down pairs by Trans_1 score (unique only) ---
top_down <- top_down %>%
  mutate(Pair = factor(Pair, levels = top_down %>% 
                         filter(Transformation == "1") %>% 
                         arrange(Mean_Score) %>% 
                         pull(Pair) %>% unique()))

# --- Plotting function ---
# --- Plotting function with vertical dashed line at 0 ---
plot_top <- function(df, title_text) {
  ggplot(df, aes(x = Mean_Score, y = Pair, color = Transformation, group = Pair)) +
    geom_point(size = 5, shape = 21, stroke = 0.5, fill = NA) +  # bigger circles
    geom_line(color = "black", linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +  # vertical dashed line at 0
    scale_color_manual(values = c("0" = "#89CFF0", "1" = "red")) +
    scale_x_continuous(limits = c(-1, 1)) +   # x-axis from -1 to 1
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),                  # remove grid
      axis.line = element_line(color = "black"),     # show axis lines
      axis.text = element_text(color = "black", size = 14),
      axis.title = element_text(color = "black", size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold")
    ) +
    labs(
      title = title_text,
      x = "Mean Interaction Score",
      y = "Cell → Neighbor",
      color = "Transformation"
    )
}


# --- Generate plots ---
plot_up <- plot_top(top_up, "Top 10 Interactions Up in Trans 1")
plot_down <- plot_top(top_down, "Top 10 Interactions Down in Trans 1")

# --- Display plots ---
plot_up
plot_down

# --- Save plots ---
ggsave("Top_10_Interactions_Up_Trans1.png", plot = plot_up, 
       width = 10, height = 6, dpi = 300)

ggsave("Top_10_Interactions_Down_Trans1.png", plot = plot_down, 
       width = 10, height = 6, dpi = 300)

