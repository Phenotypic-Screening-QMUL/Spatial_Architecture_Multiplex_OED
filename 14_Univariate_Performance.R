library(dplyr)
library(pROC)
library(tibble)

output_dir <- "Univariate_Feature_AUC"
dir.create(output_dir, showWarnings = FALSE)


# 1. Build the full feature list with subset labels

feature_map <- bind_rows(
  tibble(Feature = intensity_vars, Subset = "Intensity"),
  tibble(Feature = halo_vars,      Subset = "HALO"),
  tibble(Feature = cluster_vars,   Subset = "Cluster"),
  tibble(Feature = CN_vars,        Subset = "CN"),
  tibble(Feature = distance_vars,  Subset = "Distance")
)

outcome <- All_Means %>%
  filter(Transformation %in% c(0, 1)) %>%
  mutate(Transformation = factor(Transformation, levels = c(0, 1)))


# 2. Compute univariate AUC + CI + Wilcoxon p for each feature

run_univariate <- function(feature_name, subset_name) {
  
  vals <- outcome[[feature_name]]
  keep <- !is.na(vals) & is.finite(vals)
  
  y <- outcome$Transformation[keep]
  x <- vals[keep]
  n <- length(x)
  
  # Need both classes represented and some variance to fit
  if (n < 10 || length(unique(y)) < 2 || var(x) == 0) {
    return(tibble(
      Feature = feature_name, Subset = subset_name, N = n,
      AUC = NA, AUC_CI = NA, Wilcoxon_P = NA
    ))
  }
  
  roc_obj <- tryCatch(
    roc(response = y, predictor = x, levels = c(0, 1), direction = "auto", quiet = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(roc_obj)) {
    return(tibble(
      Feature = feature_name, Subset = subset_name, N = n,
      AUC = NA, AUC_CI = NA, Wilcoxon_P = NA
    ))
  }
  
  auc_val <- as.numeric(auc(roc_obj))
  
  auc_ci <- tryCatch(
    ci.auc(roc_obj, conf.level = 0.95, method = "bootstrap",
           boot.n = 1000, stratified = TRUE, progress = "none"),
    error = function(e) c(NA, NA, NA)
  )
  auc_ci_str <- if (!is.na(auc_ci[1])) paste0("+/-", round((auc_ci[3] - auc_ci[1]) / 2, 3)) else NA
  
  wilcox_p <- tryCatch(
    wilcox.test(x ~ y)$p.value,
    error = function(e) NA
  )
  
  tibble(
    Feature = feature_name, Subset = subset_name, N = n,
    AUC = auc_val, AUC_CI = auc_ci_str, Wilcoxon_P = wilcox_p
  )
}

cat("Computing univariate AUC for", nrow(feature_map), "features",
    "(this may take a minute due to bootstrap CIs)...\n")

univariate_results <- purrr::pmap_dfr(
  list(feature_map$Feature, feature_map$Subset),
  run_univariate
)


# 3. Multiple-comparison correction + ranking

univariate_results <- univariate_results %>%
  mutate(Wilcoxon_P_Adj_BH = p.adjust(Wilcoxon_P, method = "BH")) %>%
  arrange(desc(AUC))

write.csv(
  univariate_results,
  file.path(output_dir, "Univariate_Feature_AUC_Ranked.csv"),
  row.names = FALSE
)

cat("\nTop 15 individual features by AUC:\n")
print(univariate_results %>% slice_head(n = 15))

n_sig <- sum(univariate_results$Wilcoxon_P_Adj_BH < 0.05, na.rm = TRUE)
cat("\n", n_sig, "of", nrow(univariate_results),
    "features remain significant after BH correction (adj. p < 0.05).\n")


# 4. Plot: top 20 features by AUC, with CI error bars, colored by subset

library(ggplot2)

top20 <- univariate_results %>%
  filter(!is.na(AUC)) %>%
  slice_head(n = 20) %>%
  mutate(
    CI_HalfWidth = as.numeric(gsub("\\+/-", "", AUC_CI)),
    Feature = factor(Feature, levels = rev(Feature))
  )

univariate_plot <- ggplot(top20, aes(x = Feature, y = AUC, color = Subset)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = AUC - CI_HalfWidth, ymax = AUC + CI_HalfWidth), width = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  coord_flip() +
  labs(
    title = "Top 20 Individual Features by Univariate AUC",
    subtitle = "Error bars = 95% bootstrap CI. Dashed line = chance (AUC 0.5)",
    x = NULL, y = "AUC"
  ) +
  theme_minimal(base_size = 13)

print(univariate_plot)
ggsave(file.path(output_dir, "Top20_Univariate_Feature_AUC.png"),
       plot = univariate_plot, width = 9, height = 8, dpi = 300)

cat("\nUnivariate feature analysis complete. Results saved in:", output_dir, "\n")