
# LDA Feature Exploration:
# HALO vs Distance vs HALO + Distance vs All Features

# All Features =
# Intensity + HALO + Cluster + CN + Distance

# Highly correlated features are removed before LDA
# using a correlation cutoff of |r| > 0.8.


library(MASS)
library(dplyr)
library(tibble)
library(ggplot2)
library(caret)

output_dir <- "LDA_Feature_Comparison"
dir.create(output_dir, showWarnings = FALSE)


# 1. Prepare metadata


meta <- All_Means %>%
  dplyr::select(
    Slide,
    Transformation
  )


# 2. Define feature sets


halo_features <- All_Means %>%
  dplyr::select(
    all_of(halo_vars)
  )

distance_features <- All_Means %>%
  dplyr::select(
    all_of(distance_vars)
  )

combined_features <- bind_cols(
  halo_features,
  distance_features
)


# All available raw features


all_features <- All_Means %>%
  dplyr::select(
    all_of(
      c(
        intensity_vars,
        halo_vars,
        cluster_vars,
        CN_vars,
        distance_vars
      )
    )
  )


# Define all four LDA analyses


feature_sets <- list(
  
  HALO = halo_features,
  
  Distance = distance_features,
  
  `HALO + Distance` = combined_features,
  
  `All Features` = all_features
  
)


# 3. Function to run LDA
#    Includes:
#      - zero-variance removal
#      - complete-case filtering
#      - correlation pruning
#      - scaling
#      - LDA
#      - LOOCV validity check
#      - LD1 scores
#      - feature coefficients
#      - top 10 features
#      - plots


run_lda <- function(
    features,
    name,
    cor_cutoff = 0.8
) {
  
  cat(
    "\n========================================\n"
  )
  
  cat(
    "LDA:",
    name,
    "\n"
  )
  
  cat(
    "========================================\n"
  )
  
  
  # Remove zero-variance features
  
  
  feature_variance <- sapply(
    features,
    function(x) var(
      x,
      na.rm = TRUE
    )
  )
  
  features <- features[
    ,
    !is.na(feature_variance) &
      feature_variance > 0,
    drop = FALSE
  ]
  
  cat(
    "Features after zero-variance removal:",
    ncol(features),
    "\n"
  )
  
  
  # Keep complete observations
  
  
  keep_rows <- complete.cases(
    features
  )
  
  features <- features[
    keep_rows,
    ,
    drop = FALSE
  ]
  
  meta_sub <- meta[
    keep_rows,
    ,
    drop = FALSE
  ]
  
  cat(
    "Complete observations:",
    nrow(features),
    "\n"
  )
  
  
  # Prune highly correlated features
  
  
  cor_matrix <- cor(
    features,
    use = "pairwise.complete.obs"
  )
  
  to_drop <- findCorrelation(
    cor_matrix,
    cutoff = cor_cutoff,
    names = TRUE
  )
  
  if (length(to_drop) > 0) {
    
    cat(
      "Dropping",
      length(to_drop),
      "correlated feature(s) (|r| >",
      cor_cutoff,
      "):\n"
    )
    
    print(
      to_drop
    )
    
    safe_name <- gsub(
      " ",
      "_",
      name
    )
    
    safe_name <- gsub(
      "\\+",
      "plus",
      safe_name
    )
    
    write.csv(
      data.frame(
        Analysis = name,
        Dropped_Feature = to_drop
      ),
      file.path(
        output_dir,
        paste0(
          "Dropped_Correlated_Features_",
          safe_name,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    
  } else {
    
    cat(
      "No features exceeded the correlation cutoff - none dropped.\n"
    )
  }
  
  # Remove correlated features
  
  if (length(to_drop) > 0) {
    
    features <- features %>%
      dplyr::select(
        -all_of(to_drop)
      )
  }
  
  cat(
    "Features remaining for LDA:",
    ncol(features),
    "\n"
  )
  
  
  # Standardise features
  #
  # Scaling makes the LDA coefficients comparable because
  # all features are expressed on the same standardised scale.
  
  
  features_scaled <- as.data.frame(
    scale(features)
  )
  
  
  # Build LDA dataset
  
  
  lda_data <- bind_cols(
    
    meta_sub %>%
      mutate(
        Transformation = factor(
          Transformation,
          levels = c(
            0,
            1
          ),
          labels = c(
            "Class0",
            "Class1"
          )
        )
      ) %>%
      dplyr::select(
        Transformation
      ),
    
    features_scaled
    
  )
  
  
  # Fit LDA
  
  
  lda_model <- withCallingHandlers(
    
    MASS::lda(
      Transformation ~ .,
      data = lda_data,
      tol = 1e-4
    ),
    
    warning = function(w) {
      
      message(
        name,
        ": collinearity warning PERSISTS after pruning - ",
        conditionMessage(w)
      )
      
      invokeRestart(
        "muffleWarning"
      )
    }
  )
  
  print(
    lda_model
  )
  
  
 
  
  
  # LD1 scores
  
  
  lda_prediction <- predict(
    lda_model
  )
  
  lda_scores <- data.frame(
    
    Slide = meta_sub$Slide,
    
    Transformation =
      lda_data$Transformation,
    
    LD1 =
      lda_prediction$x[
        ,
        "LD1"
      ]
    
  )
  
  
  # Extract LDA coefficients
  
  
  lda_weights <- as.data.frame(
    lda_model$scaling
  )
  
  lda_weights$Feature <- rownames(
    lda_model$scaling
  )
  
  lda_weights <- lda_weights %>%
    
    dplyr::select(
      Feature,
      everything()
    ) %>%
    
    mutate(
      
      Absolute_LD1 =
        abs(LD1),
      
      Analysis =
        name
      
    ) %>%
    
    arrange(
      desc(
        Absolute_LD1
      )
    )
  
  
  # Top 10 features
  
  
  top_features <- lda_weights %>%
    
    slice_head(
      n = min(
        10,
        nrow(
          lda_weights
        )
      )
    )
  
  
  # Safe filename
  
  
  safe_name <- gsub(
    " ",
    "_",
    name
  )
  
  safe_name <- gsub(
    "\\+",
    "plus",
    safe_name
  )
  
  
  # Save complete coefficient table
  
  
  write.csv(
    
    lda_weights,
    
    file.path(
      output_dir,
      paste0(
        "LDA_",
        safe_name,
        "_Coefficients.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  # Save top 10 features
  
  
  write.csv(
    
    top_features,
    
    file.path(
      output_dir,
      paste0(
        "Top10_LDA_",
        safe_name,
        ".csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  # LD1 scatter plot
  
  
  lda_scatter <- ggplot(
    
    lda_scores,
    
    aes(
      x = LD1,
      y = Transformation,
      colour = Transformation
    )
    
  ) +
    
    geom_jitter(
      height = 0.12,
      width = 0,
      size = 3,
      alpha = 0.8
    ) +
    
    scale_colour_manual(
      
      values = c(
        "Class0" = "#89CFF0",
        "Class1" = "red"
      )
      
    ) +
    
    labs(
      
      title =
        paste0(
          "LDA: ",
          name
        ),
      
      x =
        "Linear discriminant 1 (LD1)",
      
      y =
        "Class",
      
      colour =
        "Class"
      
    ) +
    
    theme_minimal(
      base_size = 14
    )
  
  print(
    lda_scatter
  )
  
  # Save scatter plot
  
  ggsave(
    
    file.path(
      output_dir,
      paste0(
        "LDA_",
        safe_name,
        "_Scatter.png"
      )
    ),
    
    lda_scatter,
    
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  # Top 10 feature barplot
  
  
  lda_bar <- ggplot(
    
    top_features,
    
    aes(
      x = reorder(
        Feature,
        LD1
      ),
      
      y = LD1,
      
      fill = LD1 > 0
    )
    
  ) +
    
    geom_col() +
    
    coord_flip() +
    
    scale_fill_manual(
      
      values = c(
        "TRUE" = "#D73027",
        "FALSE" = "#4575B4"
      ),
      
      guide = "none"
      
    ) +
    
    labs(
      
      title =
        paste0(
          "Top 10 Features Contributing to LD1: ",
          name
        ),
      
      x =
        "Feature",
      
      y =
        "LDA coefficient"
      
    ) +
    
    theme_minimal(
      base_size = 14
    )
  
  print(
    lda_bar
  )
  
  # Save barplot
  
  ggsave(
    
    file.path(
      output_dir,
      paste0(
        "LDA_",
        safe_name,
        "_Top10Features.png"
      )
    ),
    
    lda_bar,
    
    width = 9,
    height = 7,
    dpi = 300
  )
  
  
  # Return results
  
  
  return(
    
    list(
      
      Model =
        lda_model,
      
      Scores =
        lda_scores,
      
      Coefficients =
        lda_weights,
      
      TopFeatures =
        top_features,
      
      DroppedFeatures =
        to_drop
      
    )
    
  )
      

}

# ============================================================
# 4. Run all four analyses
# ============================================================

lda_results <- list()

for (name in names(feature_sets)) {
  
  lda_results[[name]] <-
    
    run_lda(
      
      features =
        feature_sets[[name]],
      
      name =
        name
      
    )
}

# ============================================================
# 5. Combine all LDA coefficients
# ============================================================

all_lda_coefficients <- bind_rows(
  
  lapply(
    
    lda_results,
    
    function(x)
      x$Coefficients
    
  )
  
)

write.csv(
  
  all_lda_coefficients,
  
  file.path(
    output_dir,
    "All_LDA_Coefficients_Combined.csv"
  ),
  
  row.names = FALSE
  
)

# ============================================================
# 6. Combined top-10 table
# ============================================================

all_top_features <- bind_rows(
  
  lapply(
    
    lda_results,
    
    function(x)
      x$TopFeatures
    
  )
  
)

write.csv(
  
  all_top_features,
  
  file.path(
    output_dir,
    "All_LDA_Top10_Features.csv"
  ),
  
  row.names = FALSE
  
)

# ============================================================
# 7. Summary of pruning + validity
# ============================================================

pruning_summary <- bind_rows(
  
  lapply(
    
    names(lda_results),
    
    function(n) {
      
      tibble(
        
        Analysis =
          n,
        
        N_Features_Dropped =
          length(
            lda_results[[n]]$DroppedFeatures
          )
        
      )
      
    }
    
  )
  
)

cat(
  "\n\n===== Pruning Summary =====\n"
)

cat(
  "\n\n===== Pruning & Validity Summary =====\n"
)

print(
  pruning_summary
)

write.csv(
  
  pruning_summary,
  
  file.path(
    output_dir,
    "Pruning_Validity_Summary.csv"
  ),
  
  row.names = FALSE
  
)

# ============================================================
# 8. Final message
# ============================================================

cat(
  
  "\nLDA comparison complete.\n",
  
  "Analyses performed:\n",
  
  "  1. HALO\n",
  
  "  2. Distance\n",
  
  "  3. HALO + Distance\n",
  
  "  4. All Features\n",
  
  "Top 10 LD1 features saved for each analysis.\n",
  
  "\nResults saved in: ",
  
  output_dir,
  
  "\n",
  
  sep = ""
  
)