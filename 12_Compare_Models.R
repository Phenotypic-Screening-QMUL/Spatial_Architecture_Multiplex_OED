
# DeLong AUC Comparisons: All Models vs HALO - include test with now redundant PCA model


library(data.table)
library(dplyr)
library(pROC)
library(tibble)


# 1. Load per-slide predictions for every model into one named list

subset_probs <- fread("All_RF_Subsets_CV_Probabilities.csv")
pca_probs    <- fread("PCA_CV_Probabilities.csv")

# Standardise labels to 0/1 numeric across both sources, since the subset
# script uses factor levels "0"/"1" and the PCA script uses "Class0"/"Class1"
subset_probs <- subset_probs %>%
  mutate(Actual = as.numeric(as.character(Actual)))

pca_probs <- pca_probs %>%
  mutate(Actual = ifelse(Actual == "Class1", 1, 0)) %>%
  mutate(Subset = "PCA") %>%
  select(Subset, Slide, Actual, Pred_Prob)

all_probs <- bind_rows(subset_probs, pca_probs)

model_names <- unique(all_probs$Subset)

# Build one roc() object per model, keeping the per-slide data alongside it
# so we can align on common slides before each pairwise DeLong test.
model_data <- list()
for (m in model_names) {
  d <- all_probs %>% filter(Subset == m) %>% arrange(Slide)
  model_data[[m]] <- d
}


# 2. Helper: run DeLong test for one baseline vs every other model

run_delong_vs_baseline <- function(baseline_name, model_data) {
  
  if (!(baseline_name %in% names(model_data))) {
    stop(paste0("Baseline '", baseline_name, "' not found. Available: ",
                paste(names(model_data), collapse = ", ")))
  }
  
  baseline_df <- model_data[[baseline_name]]
  results <- list()
  
  for (m in names(model_data)) {
    if (m == baseline_name) next
    
    comp_df <- model_data[[m]]
    
    common_slides <- intersect(baseline_df$Slide, comp_df$Slide)
    
    if (length(common_slides) < 10) {
      cat("\nSkipping", m, "vs", baseline_name, "-", length(common_slides),
          "common slides, too few for a meaningful DeLong test.\n")
      next
    }
    
    b <- baseline_df %>% filter(Slide %in% common_slides) %>% arrange(Slide)
    c <- comp_df      %>% filter(Slide %in% common_slides) %>% arrange(Slide)
    
    if (!all(b$Actual == c$Actual)) {
      cat("\nWARNING:", m, "vs", baseline_name,
          "- ground-truth labels disagree on common slides. Skipping.\n")
      next
    }
    
    roc_b <- roc(response = b$Actual, predictor = b$Pred_Prob, levels = c(0,1), quiet = TRUE)
    roc_c <- roc(response = c$Actual, predictor = c$Pred_Prob, levels = c(0,1), quiet = TRUE)
    
    dl <- roc.test(roc_b, roc_c, method = "delong", paired = TRUE)
    
    results[[m]] <- tibble(
      Baseline        = baseline_name,
      Comparison      = m,
      N_Common_Slides = length(common_slides),
      AUC_Baseline    = as.numeric(auc(roc_b)),
      AUC_Comparison  = as.numeric(auc(roc_c)),
      AUC_Diff        = as.numeric(auc(roc_c)) - as.numeric(auc(roc_b)),
      Z_Statistic     = as.numeric(dl$statistic),
      P_Value         = dl$p.value
    )
    
    cat("\n---", m, "vs", baseline_name, "---\n")
    print(dl)
  }
  
  if (length(results) == 0) return(NULL)
  
  bind_rows(results) %>%
    mutate(P_Adj_BH = p.adjust(P_Value, method = "BH"))
}


# 3. Run both sets of comparisons

cat("\n========== DeLong: All Models vs HALO ==========\n")
vs_halo <- run_delong_vs_baseline("HALO", model_data)



# 4. Save results

if (!is.null(vs_halo)) {
  print(vs_halo)
  fwrite(vs_halo, "DeLong_AUC_Comparisons_vs_HALO.csv")
  cat("\nSaved: DeLong_AUC_Comparisons_vs_HALO.csv\n")
}

