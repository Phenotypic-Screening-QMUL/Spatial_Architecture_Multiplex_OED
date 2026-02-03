# ------------------------------
# 0. Load Libraries
# ------------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(caret)
library(pROC)
library(randomForest)
library(tibble)

set.seed(123)

# ------------------------------
# 1. Load and Prepare Slide Summary
# ------------------------------
Data <- fread("Neighbourhood Analysis.csv", encoding = "Latin-1") 

Slide_Summary <- Data %>%
  group_by(Slide, Region, Transformation, Clinical_Class) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  select(-V1, -ClusterNN)

Slide_Summary <- Slide_Summary[, -c(14, 15:29)]

Meta_Data      <- Slide_Summary[, c(1:4)]
Intensity_Data <- Slide_Summary[, c(5:13)] %>% scale() %>% as.data.frame()

# ------------------------------
# 2. Load Other Tables
# ------------------------------
HALO_Data     <- fread("Percentage_Positive_Counts.csv", encoding = "Latin-1") %>% select(-V1)
Cluster_Data  <- fread("Cluster_Cells_Region.csv", encoding = "Latin-1") %>% select(-V1)
CNN_Data      <- fread("Neighbourhood_Region_Data.csv", encoding = "Latin-1") %>% select(-V1)
Distance_Data <- fread("Distance_DF.csv", encoding = "Latin-1") %>% select(-Slide)

# Drop unwanted cols from cluster data
Cluster_Data <- Cluster_Data[,-c(1,2)]

# ------------------------------
# 3. Combine All Data
# ------------------------------
Final_Data <- cbind(
  Meta_Data,
  Intensity_Data,
  HALO_Data,
  Cluster_Data,
  CNN_Data,
  Distance_Data
)

Final_Data <- Final_Data %>% filter(Transformation %in% c(0,1))

# ------------------------------
# 4. Aggregate to Slide-Level Means
# ------------------------------
All_Means <- Final_Data %>%
  mutate(Slide = sub("_.*", "", Slide)) %>%
  group_by(Slide) %>%
  summarise(
    across(where(is.numeric), mean, na.rm = TRUE),
    across(where(~ !is.numeric(.)), first),
    .groups = "drop"
  ) %>%
  mutate(Clinical_Class = as.factor(Clinical_Class))

# ------------------------------
# 5. Define Subsets
# ------------------------------
intensity_vars <- colnames(Intensity_Data)
halo_vars      <- colnames(HALO_Data)
cluster_vars   <- colnames(Cluster_Data)
cnn_vars       <- colnames(CNN_Data)
distance_vars  <- colnames(Distance_Data)

subset_list <- list(
  Intensity = All_Means %>% select(Slide, Transformation, Clinical_Class, all_of(intensity_vars)),
  HALO      = All_Means %>% select(Slide, Transformation, Clinical_Class, all_of(halo_vars)),
  Cluster   = All_Means %>% select(Slide, Transformation, Clinical_Class, all_of(cluster_vars)),
  CNN       = All_Means %>% select(Slide, Transformation, Clinical_Class, all_of(cnn_vars)),
  Distance  = All_Means %>% select(Slide, Transformation, Clinical_Class, all_of(distance_vars))
)

# ------------------------------
# 6. k-Fold CV Random Forest with CV Probabilities
# ------------------------------
k <- 10
all_results_rf <- list()

for (subset_name in names(subset_list)) {
  
  cat("\n========== Subset:", subset_name, "==========\n")
  
  df <- subset_list[[subset_name]] %>%
    filter(Transformation %in% c(0,1)) %>%
    mutate(Transformation = factor(Transformation, levels = c(0,1)))
  
  predictors <- df %>% select(-Slide, -Transformation, -Clinical_Class) %>% select(where(is.numeric))
  
  # Drop zero variance predictors
  zero_var <- sapply(predictors, function(x) var(x, na.rm = TRUE) == 0)
  predictors <- predictors[, !zero_var]
  
  # Drop rows with NA
  df_clean <- bind_cols(df %>% select(Slide, Transformation), predictors) %>% drop_na()
  n <- nrow(df_clean)
  
  if (n < k) {
    cat("Not enough samples after cleaning. Skipping.\n")
    next
  }
  
  # Create folds
  folds <- caret::createFolds(df_clean$Transformation, k = k, list = TRUE, returnTrain = FALSE)
  
  all_pred_class <- rep(NA, n)
  all_pred_prob  <- rep(NA, n)
  feat_imp_accum <- setNames(rep(0, ncol(predictors)), colnames(predictors))
  
  # --- k-Fold CV Loop ---
  for (i in 1:k) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(1:n, test_idx)
    
    train_data <- df_clean[train_idx, ]
    test_data  <- df_clean[test_idx, ]
    
    rf_model <- randomForest(
      x = train_data %>% select(-Slide, -Transformation),
      y = train_data$Transformation,
      importance = TRUE,
      ntree = 500
    )
    
    # Predict class and probability for "1"
    all_pred_class[test_idx] <- as.character(predict(rf_model, newdata = test_data %>% select(-Slide, -Transformation)))
    all_pred_prob[test_idx]  <- predict(rf_model, newdata = test_data %>% select(-Slide, -Transformation), type = "prob")[, "1"]
    
    # Accumulate feature importance
    feat_imp_accum <- feat_imp_accum + rf_model$importance[, "MeanDecreaseGini"]
  }
  
  predicted_classes <- factor(all_pred_class, levels = levels(df_clean$Transformation))
  actual_classes    <- df_clean$Transformation
  
  # --- Confusion Matrix ---
  cm <- caret::confusionMatrix(predicted_classes, actual_classes)
  TP <- sum(predicted_classes == "1" & actual_classes == "1")
  TN <- sum(predicted_classes == "0" & actual_classes == "0")
  FP <- sum(predicted_classes == "1" & actual_classes == "0")
  FN <- sum(predicted_classes == "0" & actual_classes == "1")
  
  # --- Metrics ---
  roc_obj <- roc(response = as.numeric(actual_classes)-1, predictor = all_pred_prob, levels = c(0,1))
  
  metrics <- tibble(
    Subset      = subset_name,
    Accuracy    = cm$overall["Accuracy"],
    Kappa       = cm$overall["Kappa"],
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    Precision   = cm$byClass["Pos Pred Value"],
    Recall      = cm$byClass["Sensitivity"],
    F1          = 2*(cm$byClass["Sensitivity"]*cm$byClass["Pos Pred Value"])/(cm$byClass["Sensitivity"]+cm$byClass["Pos Pred Value"]),
    AUC         = as.numeric(auc(roc_obj)),
    TP = TP, TN = TN, FP = FP, FN = FN
  )
  
  # --- Feature Importance (average over folds) ---
  feat_imp_df <- data.frame(
    Feature = names(feat_imp_accum),
    MeanDecreaseGini = feat_imp_accum / k,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(MeanDecreaseGini)) %>%
    mutate(Subset = subset_name)
  
  # --- Save Results ---
  all_results_rf[[subset_name]] <- list(
    ConfusionMatrix   = as.data.frame(cm$table),
    Metrics           = metrics,
    FeatureImportance = feat_imp_df,
    CV_Probabilities  = data.frame(Slide = df_clean$Slide, Actual = actual_classes, Pred_Prob = all_pred_prob)
  )
  
  # --- Plot ROC Curve (CV-based) ---
  roc_df <- data.frame(
    specificity = rev(roc_obj$specificities),
    sensitivity = rev(roc_obj$sensitivities)
  )
  
  roc_plot <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(color = "#1c61b6", size = 1.2) +
    geom_abline(linetype = "dashed", color = "gray") +
    labs(title = paste0("ROC Curve (CV): ", subset_name),
         x = "1 - Specificity", y = "Sensitivity") +
    annotate("text", x = 0.7, y = 0.1, label = paste0("AUC = ", round(auc(roc_obj),3)), size = 5) +
    theme_minimal(base_size = 14) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          axis.line        = element_line(color = "black", size = 0.8))
  
  ggsave(paste0("Plots/ROC_RF_CV_", subset_name, ".png"), plot = roc_plot, width = 8, height = 6, dpi = 300)
}

# ------------------------------
# 7. Save All Results to CSV
# ------------------------------
flat_cm <- bind_rows(lapply(all_results_rf, function(x) x$ConfusionMatrix))
fwrite(flat_cm, "All_RF_Subsets_ConfusionMatrix.csv")

flat_metrics <- bind_rows(lapply(all_results_rf, function(x) x$Metrics))
fwrite(flat_metrics, "All_RF_Subsets_Metrics.csv")

flat_feat_imp <- bind_rows(lapply(all_results_rf, function(x) x$FeatureImportance))
fwrite(flat_feat_imp, "All_RF_Subsets_FeatureImportance.csv")

flat_probs <- bind_rows(lapply(all_results_rf, function(x) x$CV_Probabilities))
fwrite(flat_probs, "All_RF_Subsets_CV_Probabilities.csv")

cat("Random Forest k-Fold CV complete. Metrics, feature importance, and CV-based ROC curves saved.\n")

# ------------------------------
# 8. Plot Top 5 Feature Importances
# ------------------------------
for (subset_name in names(all_results_rf)) {
  
  feat_imp_df <- all_results_rf[[subset_name]]$FeatureImportance
  
  top5 <- feat_imp_df %>%
    top_n(5, wt = MeanDecreaseGini) %>%
    arrange(MeanDecreaseGini)
  
  p <- ggplot(top5, aes(x = reorder(Feature, MeanDecreaseGini), y = MeanDecreaseGini)) +
    geom_col(fill = "#1c61b6") +
    coord_flip() +
    labs(
      title = paste0("Top 5 Feature Importance: ", subset_name),
      x = "Feature",
      y = "Mean Decrease Gini"
    ) +
    theme_minimal(base_size = 14)
  
  print(p)
  
  ggsave(paste0("Plots/Top5_FeatureImportance_", subset_name, ".png"),
         plot = p, width = 8, height = 6, dpi = 300)
}
