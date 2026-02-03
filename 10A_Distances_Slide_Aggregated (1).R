# --- Load Required Libraries ---
library(dplyr)
library(data.table)
library(FNN)
library(tidyr)
library(ggplot2)

#Note Slide Var = Slide+Region - then split into Slide_Base for stats

# --- Load the Data ---
Data <- fread("Neighbourhood Analysis.csv", encoding = "Latin-1") %>%
  filter(!Cell_Type %in% c("Artefact", "PDL1 High", "Negative")) %>%
  mutate(
    Transformation = factor(Transformation, 
                            levels = c("Normal Control", "Cancer Control", "0", "1")),
    XMicron = XCentroid / 3,
    YMicron = YCentroid / 3
  )

# --- Define NN Distance Function (Within-Slide) ---
find_nn_distances_wide_by_Slide <- function(df) {
  results_list <- list()
  
  for (slide_id in unique(df$Slide)) {
    df_slide <- df %>% filter(Slide == slide_id)
    unique_types <- unique(df_slide$Cell_Type)
    
    # Initialize distance columns
    for (other_type in unique_types) {
      df_slide[[paste0("Dist_to_", other_type)]] <- NA
    }
    
    # Loop through each target cell type
    for (other_type in unique_types) {
      other_cells <- df_slide %>% filter(Cell_Type == other_type)
      
      if (nrow(other_cells) > 0) {
        knn_result <- get.knnx(
          data  = other_cells[, c("XMicron", "YMicron")],
          query = df_slide[, c("XMicron", "YMicron")],
          k     = 1
        )
        # Assign distances
        df_slide[[paste0("Dist_to_", other_type)]] <- as.vector(knn_result$nn.dist)
      }
    }
    
    results_list[[slide_id]] <- df_slide
  }
  
  return(bind_rows(results_list))
}

# --- Run NN Calculation on the FULL Dataset (Within-Slide) ---
nn_df <- find_nn_distances_wide_by_Slide(Data)

# --- Group by Slide and Cell_Type, then average distances ---
summary_distances <- nn_df %>%
  group_by(Slide, Cell_Type, Transformation) %>%
  summarise(across(starts_with("Dist_to_"), mean, na.rm = TRUE), .groups = "drop")

# --- Extract base Slide name (remove everything after first underscore) ---
summary_distances <- summary_distances %>%
  mutate(Slide_Base = sub("_.*", "", Slide))

# --- Group by cleaned Slide, Cell_Type, and Transformation, then average distances ---
summary_distances_aggregated <- summary_distances %>%
  group_by(Slide_Base, Cell_Type, Transformation) %>%
  summarise(across(starts_with("Dist_to_"), mean, na.rm = TRUE), .groups = "drop")

# --- Filter to Transformation values "0" and "1" ---
filtered_data <- summary_distances_aggregated %>%
  filter(Transformation %in% c("0", "1"))

# --- Get distance columns ---
dist_cols <- names(filtered_data)[grepl("^Dist_to_", names(filtered_data))]

# --- Initialize results list for t-tests ---
results_list <- list()

# --- Loop over Cell_Types and Distance columns for t-tests ---
for (cell_type in unique(filtered_data$Cell_Type)) {
  df_sub <- filtered_data %>% filter(Cell_Type == cell_type)
  
  for (dist_col in dist_cols) {
    test_data <- df_sub %>% select(Transformation, !!sym(dist_col))
    test_data <- test_data %>% filter(!is.na(.data[[dist_col]]))
    
    if (n_distinct(test_data$Transformation) == 2) {
      n0 <- sum(test_data$Transformation == "0")
      n1 <- sum(test_data$Transformation == "1")
      
      if (n0 >= 3 & n1 >= 3) {
        test_data$Transformation <- factor(test_data$Transformation)
        t_result <- t.test(as.formula(paste0("`", dist_col, "` ~ Transformation")), data = test_data)
        
        results_list[[length(results_list) + 1]] <- data.frame(
          Cell_Type   = cell_type,
          Distance_To = gsub("^Dist_to_", "", dist_col),
          Mean        = mean(test_data[[dist_col]], na.rm = TRUE),
          p_value     = t_result$p.value,
          Status      = "OK"
        )
      } else {
        results_list[[length(results_list) + 1]] <- data.frame(
          Cell_Type   = cell_type,
          Distance_To = gsub("^Dist_to_", "", dist_col),
          Mean        = mean(test_data[[dist_col]], na.rm = TRUE),
          p_value     = NA,
          Status      = "Too few obs"
        )
      }
    }
  }
}

# --- Combine t-test results ---
t_summary <- bind_rows(results_list) %>%
  mutate(
    Significant = case_when(
      Status == "Too few obs" ~ "Too few obs",
      p_value < 0.05           ~ "Yes",
      TRUE                    ~ "No"
    )
  )

# --- Calculate mean distances by Cell_Type, Distance_To, and Transformation ---
mean_values <- filtered_data %>%
  pivot_longer(cols = starts_with("Dist_to_"), names_to = "Distance_To", values_to = "Distance") %>%
  mutate(Distance_To = sub("^Dist_to_", "", Distance_To)) %>%
  group_by(Cell_Type, Distance_To, Transformation) %>%
  summarise(Mean_Distance = mean(Distance, na.rm = TRUE), .groups = "drop")

# --- Reshape means to wide format ---
mean_wide <- mean_values %>%
  pivot_wider(names_from = Transformation, values_from = Mean_Distance, names_prefix = "Trans_") %>%
  mutate(mean_diff = Trans_1 - Trans_0)

# --- Join mean_diff back to t_summary ---
t_summary <- t_summary %>%
  left_join(mean_wide %>% select(Cell_Type, Distance_To, mean_diff), by = c("Cell_Type", "Distance_To"))

# --- Clip mean_diff for bubble size ---
t_summary <- t_summary %>%
  mutate(
    bubble_size = pmin(abs(mean_diff), 100)  # Clip at 100
  )

# --- Plot with clipped bubble sizes and legend ---
ggplot(t_summary, aes(y = Distance_To, x = Cell_Type)) +
  geom_point(
    aes(
      size = bubble_size,
      fill = case_when(
        mean_diff > 0 ~ "Up in Transformed",
        mean_diff < 0 ~ "Down in Transformed",
        TRUE          ~ "No Change"
      )
    ),
    shape = 21, color = "black", stroke = 0.8
  ) +
  geom_point(
    data = filter(t_summary, Significant == "Yes"),
    aes(x = Cell_Type, y = Distance_To),
    shape = 19, color = "black", size = 4
  ) +
  scale_fill_manual(
    name   = "Direction",
    values = c(
      "Up in Transformed"   = "red",
      "Down in Transformed" = "turquoise",
      "No Change"           = "grey70"
    )
  ) +
  scale_size_continuous(
    name = "Absolute\nMean Diff",
    range = c(6, 18),
    breaks = c(25, 50, 75, 100)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 8)),
    size = guide_legend()
  ) +
  theme_minimal(base_size = 20) +
  labs(
    title = "Distance Changes (Nearest Neighbour, Within Slide)",
    y     = "Distance To Cell Type",
    x     = "Cell Type"
  ) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y     = element_text(size = 18),
    axis.title.x    = element_text(size = 20, face = "bold"),
    axis.title.y    = element_text(size = 20, face = "bold"),
    plot.title      = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.key.size = unit(2, "lines"),
    legend.text     = element_text(size = 18),
    legend.title    = element_text(size = 20, face = "bold")
  )

# --- Export per-slide distance matrix (no Region) ---
per_slide_dist <- nn_df %>%
  group_by(Slide, Cell_Type) %>%
  summarise(across(starts_with("Dist_to_"), mean, na.rm = TRUE), .groups = "drop")

wide_df <- per_slide_dist %>%
  pivot_longer(
    cols = starts_with("Dist_to_"),
    names_to = "Target_Type",
    values_to = "Distance"
  ) %>%
  mutate(
    Target_Type = sub("^Dist_to_", "", Target_Type),
    VarName = paste0(Cell_Type, "_DistTo_", Target_Type)
  ) %>%
  select(Slide, VarName, Distance) %>%
  pivot_wider(
    names_from = VarName,
    values_from = Distance
  )

write.csv(wide_df, "Distance_DF.csv", row.names = FALSE)

# Save plot
if (!dir.exists("Neighbourhoods/Plots")) dir.create("Neighbourhoods/Plots", recursive = TRUE)
ggsave("Neighbourhoods/Plots/NN_Distance_Changes.png", plot = last_plot(), width = 12, height = 8, dpi = 300)
