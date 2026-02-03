# Load libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------
# 1. Load and Prepare Slide Summary
# ------------------------------
Data <- fread("Neighbourhood Analysis.csv", encoding = "Latin-1") 

Slide_Summary <- Data %>%
  group_by(Slide, Region, Transformation, Clinical_Class) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
  select(-V1, -ClusterNN)

# Remove unwanted columns by index
Slide_Summary <- Slide_Summary[, -c( 14, 15:29)]

# Separate meta and scaled intensity
Meta_Data      <- Slide_Summary[, c(1:4)]
Intensity_Data <- Slide_Summary[, c(5:13)] %>% scale() %>% as.data.frame()

# ------------------------------
# 2. Load Other Tables
# ------------------------------
HALO_Data     <- fread("Percentage_Positive_Counts.csv", encoding = "Latin-1") %>% select(-V1)
Cluster_Data  <- fread("Cluster_Cells_Region.csv", encoding = "Latin-1") %>% select(-V1)
CNN_Data      <- fread("Neighbourhood_Region_Data.csv", encoding = "Latin-1") %>% select(-V1)
Distance_Data <- fread("Distance_DF.csv", encoding = "Latin-1") %>% select( -Slide )

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

# Filter Transformation to 0 and 1
Final_Data <- Final_Data %>% filter(Transformation %in% c(0, 1))

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



#### Double Check - Something Has changed

# ─── Libraries ───────────────────────────────────────────────────────────────
library(dplyr)
library(caret)
library(randomForest)
library(pROC)
library(ggplot2)

set.seed(42) #Note this can run very slightly differently on different hardward / R versions

# ─── 1. Prepare predictors + metadata from your existing All_Means ─────────
# (no assumptions—just use All_Means as is)
meta <- All_Means %>% select(Slide, Transformation)
preds <- All_Means %>% select(-Slide, -Region, -Clinical_Class, -Transformation)

# ─── 2. Drop zero-variance columns safely ─────────────────────────────────
vars_var  <- sapply(preds, function(x) var(x, na.rm = TRUE))
keep_cols <- which(!is.na(vars_var) & vars_var > 0)
preds_clean <- preds[, keep_cols]

# ─── 3. Drop rows with any NA ──────────────────────────────────────────────
keep_rows  <- complete.cases(preds_clean)
preds_clean <- preds_clean[keep_rows, , drop = FALSE]
meta_clean  <- meta[keep_rows, ]

# ─── 4. PCA on cleaned predictors ─────────────────────────────────────────
pca_res    <- prcomp(preds_clean, scale. = TRUE)
pc_scores  <- as.data.frame(pca_res$x[, c(1:4)])  # select PCs 1–4

# ─── 5. Build modeling data.frame ─────────────────────────────────────────
model_df <- bind_cols(
  meta_clean %>% mutate(
    Transformation = factor(Transformation, levels = c(0, 1), labels = c("Class0","Class1"))
  ),
  pc_scores
)

# ─── 6. 10‑fold CV setup ───────────────────────────────────────────────────
ctrl <- trainControl(
  method           = "cv",
  number           = 10,
  classProbs       = TRUE,
  summaryFunction  = twoClassSummary,
  savePredictions  = TRUE
)

# ─── 7. Train RF with 10‑fold CV, optimizing AUC ───────────────────────────
rf_cv <- train(
  Transformation ~ .,
  data    = model_df %>% select(-Slide),
  method  = "rf",
  metric  = "ROC",
  trControl = ctrl,
  tuneLength = 5
)

print(rf_cv)  # shows CV ROC and Accuracy

# ─── 8. Confusion matrix (pooled CV predictions) ─────────────────────────
best_pred <- rf_cv$pred %>% filter(mtry == rf_cv$bestTune$mtry)

# Confusion matrix
cm <- confusionMatrix(best_pred$pred, best_pred$obs)
print(cm)

# Extract Accuracy
accuracy <- cm$overall["Accuracy"]
cat("Pooled CV Accuracy:", round(accuracy, 4), "\n")

# Also view other metrics
rf_cv$results       # ROC, Sens, Spec, etc. for all tuning parameters
rf_cv$bestTune      # best mtry

library(pROC)
library(ggplot2)

# ─── 1. Get predicted probabilities for the positive class (Class1)
best_pred <- rf_cv$pred %>% 
  filter(mtry == rf_cv$bestTune$mtry)

roc_obj <- roc(
  response = best_pred$obs, 
  predictor = best_pred$Class1,   # probabilities for "Class1"
  levels = c("Class0", "Class1"), # ensure correct order
  direction = "<"
)

library(ggplot2)

roc_df <- data.frame(
  specificity = rev(roc_obj$specificities),
  sensitivity = rev(roc_obj$sensitivities)
)

roc_plot <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#1c61b6", size = 1.2) +
  geom_abline(linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Model 6: PCA",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  annotate(
    "text", x = 0.7, y = 0.1, 
    label = paste0("AUC = ", round(auc(roc_obj), 3)),
    size = 5
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90"),  # optional light grid
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black", size = 0.8)
  )

# Save
ggsave("Plots/ROC_RF_CV_ggplot.png", plot = roc_plot, width = 8, height = 6)


# Get variable importance
rf_importance <- varImp(rf_cv, scale = T)$importance

# Add PC names
rf_importance$PC <- rownames(rf_importance)

# Arrange by importance
rf_importance <- rf_importance %>% arrange(desc(Overall))

# Plot all PCs (or top n if you have >5)
ggplot(rf_importance, aes(x = reorder(PC, Overall), y = Overall)) +
  geom_col(fill = "#1c61b6") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Random Forest: PC Importance",
    x = "Principal Component",
    y = "Importance (scaled)"
  )

# Save the RF importance plot
ggsave("Plots/RF_PC_Importance.png", plot = last_plot(), width = 8, height = 6, dpi = 300)


# ------------------------------
# PCA Scatter Plot: PC1 vs PC3
# ------------------------------
pc_plot <- ggplot(model_df, aes(x = PC1, y = PC3, color = Transformation)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, aes(fill = Transformation), geom = "polygon", alpha = 0.2, color = NA) +
  scale_color_manual(values = c("Class1" = "red", "Class0" = "#89CFF0")) +
  scale_fill_manual(values = c("Class1" = "red", "Class0" = "#89CFF0")) +
  labs(
    title = "PCA: PC1 vs PC3",
    x = "PC1",
    y = "PC3",
    color = "Class",
    fill = "Class"
  ) +
  theme_minimal(base_size = 14) 

# Display plot
print(pc_plot)

# Save plot
ggsave("Plots/PCA_PC1_PC3.png", plot = pc_plot, width = 8, height = 6, dpi = 300)


# ------------------------------
# 10. Compute metrics for PCA RF with TP/FP/TN/FN
# ------------------------------
library(dplyr)
library(pROC)
library(tibble)

# Use pooled CV predictions from caret
best_pred <- rf_cv$pred %>% filter(mtry == rf_cv$bestTune$mtry)

predicted_classes <- best_pred$pred
actual_classes    <- best_pred$obs

# Confusion matrix
cm <- caret::confusionMatrix(predicted_classes, actual_classes)

# Extract TP, TN, FP, FN
TP <- sum(predicted_classes == "Class1" & actual_classes == "Class1")
TN <- sum(predicted_classes == "Class0" & actual_classes == "Class0")
FP <- sum(predicted_classes == "Class1" & actual_classes == "Class0")
FN <- sum(predicted_classes == "Class0" & actual_classes == "Class1")

# Metrics tibble
metrics_pca_rf <- tibble(
  Accuracy    = cm$overall["Accuracy"],
  Kappa       = cm$overall["Kappa"],
  Sensitivity = cm$byClass["Sensitivity"],
  Specificity = cm$byClass["Specificity"],
  Precision   = cm$byClass["Pos Pred Value"],
  Recall      = cm$byClass["Sensitivity"],
  F1          = 2 * (cm$byClass["Sensitivity"] * cm$byClass["Pos Pred Value"]) /
    (cm$byClass["Sensitivity"] + cm$byClass["Pos Pred Value"]),
  AUC         = tryCatch(
    as.numeric(
      roc(
        response  = as.numeric(actual_classes)-1,
        predictor = as.numeric(predicted_classes)-1
      )$auc
    ),
    error = function(e) NA
  ),
  TP = TP,
  TN = TN,
  FP = FP,
  FN = FN
)

print(metrics_pca_rf)


#### Double Check - Something Has changed

# ─── Libraries ───────────────────────────────────────────────────────────────
library(dplyr)
library(caret)
library(randomForest)
library(pROC)
library(ggplot2)

set.seed(42)

# ─── 1. Prepare predictors + metadata from your existing All_Means ─────────
# (no assumptions—just use All_Means as is)
meta <- All_Means %>% select(Slide, Transformation)
preds <- All_Means %>% select(-Slide, -Region, -Clinical_Class, -Transformation)

# ─── 2. Drop zero-variance columns safely ─────────────────────────────────
vars_var  <- sapply(preds, function(x) var(x, na.rm = TRUE))
keep_cols <- which(!is.na(vars_var) & vars_var > 0)
preds_clean <- preds[, keep_cols]

# ─── 3. Drop rows with any NA ──────────────────────────────────────────────
keep_rows  <- complete.cases(preds_clean)
preds_clean <- preds_clean[keep_rows, , drop = FALSE]
meta_clean  <- meta[keep_rows, ]

# ─── 4. PCA on cleaned predictors ─────────────────────────────────────────
pca_res    <- prcomp(preds_clean, scale. = TRUE)
pc_scores  <- as.data.frame(pca_res$x[, c(1:4)])  # select PCs 1–4

# ─── 5. Build modeling data.frame ─────────────────────────────────────────
model_df <- bind_cols(
  meta_clean %>% mutate(
    Transformation = factor(Transformation, levels = c(0, 1), labels = c("Class0","Class1"))
  ),
  pc_scores
)

# ─── 6. 10‑fold CV setup ───────────────────────────────────────────────────
ctrl <- trainControl(
  method           = "cv",
  number           = 10,
  classProbs       = TRUE,
  summaryFunction  = twoClassSummary,
  savePredictions  = TRUE
)

# ─── 7. Train RF with 10‑fold CV, optimizing AUC ───────────────────────────
rf_cv <- train(
  Transformation ~ .,
  data    = model_df %>% select(-Slide),
  method  = "rf",
  metric  = "ROC",
  trControl = ctrl,
  tuneLength = 5
)

print(rf_cv)  # shows CV ROC and Accuracy

# ─── 8. Confusion matrix (pooled CV predictions) ─────────────────────────
best_pred <- rf_cv$pred %>% filter(mtry == rf_cv$bestTune$mtry)

# Confusion matrix
cm <- confusionMatrix(best_pred$pred, best_pred$obs)
print(cm)

# Extract Accuracy
accuracy <- cm$overall["Accuracy"]
cat("Pooled CV Accuracy:", round(accuracy, 4), "\n")

# Also view other metrics
rf_cv$results       # ROC, Sens, Spec, etc. for all tuning parameters
rf_cv$bestTune      # best mtry

library(pROC)
library(ggplot2)

# ─── 1. Get predicted probabilities for the positive class (Class1)
best_pred <- rf_cv$pred %>% 
  filter(mtry == rf_cv$bestTune$mtry)

roc_obj <- roc(
  response = best_pred$obs, 
  predictor = best_pred$Class1,   # probabilities for "Class1"
  levels = c("Class0", "Class1"), # ensure correct order
  direction = "<"
)

library(ggplot2)

roc_df <- data.frame(
  specificity = rev(roc_obj$specificities),
  sensitivity = rev(roc_obj$sensitivities)
)

roc_plot <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#1c61b6", size = 1.2) +
  geom_abline(linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Model 6: PCA",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  annotate(
    "text", x = 0.7, y = 0.1, 
    label = paste0("AUC = ", round(auc(roc_obj), 3)),
    size = 5
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90"),  # optional light grid
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black", size = 0.8)
  )

# Save
ggsave("Plots/ROC_RF_CV_ggplot.png", plot = roc_plot, width = 8, height = 6)


# Get variable importance
rf_importance <- varImp(rf_cv, scale = T)$importance

# Add PC names
rf_importance$PC <- rownames(rf_importance)

# Arrange by importance
rf_importance <- rf_importance %>% arrange(desc(Overall))

# Plot all PCs (or top n if you have >5)
ggplot(rf_importance, aes(x = reorder(PC, Overall), y = Overall)) +
  geom_col(fill = "#1c61b6") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Random Forest: PC Importance",
    x = "Principal Component",
    y = "Importance (scaled)"
  )

# Save the RF importance plot
ggsave("Plots/RF_PC_Importance.png", plot = last_plot(), width = 8, height = 6, dpi = 300)


# ------------------------------
# PCA Scatter Plot: PC1 vs PC3
# ------------------------------
pc_plot <- ggplot(model_df, aes(x = PC1, y = PC3, color = Transformation)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, aes(fill = Transformation), geom = "polygon", alpha = 0.2, color = NA) +
  scale_color_manual(values = c("Class1" = "red", "Class0" = "#89CFF0")) +
  scale_fill_manual(values = c("Class1" = "red", "Class0" = "#89CFF0")) +
  labs(
    title = "PCA: PC1 vs PC3",
    x = "PC1",
    y = "PC3",
    color = "Class",
    fill = "Class"
  ) +
  theme_minimal(base_size = 14) 

# Display plot
print(pc_plot)

# Save plot
ggsave("Plots/PCA_PC1_PC3.png", plot = pc_plot, width = 8, height = 6, dpi = 300)


# ------------------------------
# 10. Compute metrics for PCA RF with TP/FP/TN/FN
# ------------------------------
library(dplyr)
library(pROC)
library(tibble)

# Use pooled CV predictions from caret
best_pred <- rf_cv$pred %>% filter(mtry == rf_cv$bestTune$mtry)

predicted_classes <- best_pred$pred
actual_classes    <- best_pred$obs

# Confusion matrix
cm <- caret::confusionMatrix(predicted_classes, actual_classes)

# Extract TP, TN, FP, FN
TP <- sum(predicted_classes == "Class1" & actual_classes == "Class1")
TN <- sum(predicted_classes == "Class0" & actual_classes == "Class0")
FP <- sum(predicted_classes == "Class1" & actual_classes == "Class0")
FN <- sum(predicted_classes == "Class0" & actual_classes == "Class1")

# Metrics tibble
metrics_pca_rf <- tibble(
  Accuracy    = cm$overall["Accuracy"],
  Kappa       = cm$overall["Kappa"],
  Sensitivity = cm$byClass["Sensitivity"],
  Specificity = cm$byClass["Specificity"],
  Precision   = cm$byClass["Pos Pred Value"],
  Recall      = cm$byClass["Sensitivity"],
  F1          = 2 * (cm$byClass["Sensitivity"] * cm$byClass["Pos Pred Value"]) /
    (cm$byClass["Sensitivity"] + cm$byClass["Pos Pred Value"]),
  AUC         = tryCatch(
    as.numeric(
      roc(
        response  = as.numeric(actual_classes)-1,
        predictor = as.numeric(predicted_classes)-1
      )$auc
    ),
    error = function(e) NA
  ),
  TP = TP,
  TN = TN,
  FP = FP,
  FN = FN
)

print(metrics_pca_rf)

library(dplyr)
library(tibble)

# ------------------------------
# 1. Extract PCA loadings
# ------------------------------
loadings_df <- pca_res$rotation %>%
  as.data.frame() %>%
  rownames_to_column("Feature")

# ------------------------------
# 2. Build PC1 loading table
# ------------------------------
pc1_loadings <- loadings_df %>%
  transmute(
    Feature    = Feature,
    Loading    = PC1,
    AbsLoading = abs(PC1)
  )

# ------------------------------
# 3. Select Top 10 by absolute loading
# ------------------------------
Top10_PC1_Features <- pc1_loadings %>%
  arrange(desc(AbsLoading)) %>%
  slice(1:10)

# ------------------------------
# 4. Save to CSV
# ------------------------------
write.csv(
  Top10_PC1_Features,
  "Top10_PC1_Features.csv",
  row.names = FALSE
)

# View
print(Top10_PC1_Features)


