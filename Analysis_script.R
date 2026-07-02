# ===========================================================================================
# PAPER: Body Adiposity Index Outperforms Central Obesity Measures in Predicting HbA1c Levels
# among Hausa Women across Menopausal Status
# ===========================================================================================
#
# Scope: descriptive statistics, correlation analysis, regression model
# comparison (Linear Regression, LASSO, Random Forest, SVR), feature
# importance (LASSO coefficients + RF impurity), and subgroup regression.
#
# HOW TO RUN IN RStudio:
# 1. Place 'clean_data.csv' in your working directory.
# 2. Run once:
#      install.packages(c("tidymodels", "tidyverse", "vip", "corrplot",
#                          "patchwork", "glmnet", "ranger", "kernlab"))
# 3. Source this script top to bottom. Tables print to console;
#    figures save to ./figures_paperA/
#
# ==============================================================================

library(tidymodels)
library(tidyverse)
library(vip)
library(corrplot)
library(patchwork)
library(glmnet)
library(ranger)
library(kernlab)

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------
DATA_PATH    <- "clean_data.csv"
FIGURES_DIR  <- "figures_paperA"
RANDOM_STATE <- 42
N_FOLDS      <- 5

dir.create(FIGURES_DIR, showWarnings = FALSE)
set.seed(RANDOM_STATE)
theme_set(theme_minimal(base_size = 12))

FEATURE_NAMES <- c("Age", "BMI", "Neck Circ.", "Hip Circ.", "Waist Circ.",
                    "WHR", "WHtR", "BAI", "Menopausal Status")
PREDICTOR_COLS <- c("Age", "BMI", "NeckCirc", "HipCirc", "WaistCirc", "WHR", "WHtR", "BAI")

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================
df_raw <- read_csv(DATA_PATH, show_col_types = FALSE)

df <- df_raw %>%
  mutate(MENO_STATUS = if_else(str_trim(`Menopausal status`) == "Postmenopausal", 1, 0)) %>%
  rename(Age = AGE, BMI = `BMI(kg/m²)`, NeckCirc = `NC (cm)`,
         HipCirc = `HC (cm)`, WaistCirc = `WC (cm)`, WHtR = WHTR, BAI = `BAI (%)`) %>%
  select(Age, BMI, NeckCirc, HipCirc, WaistCirc, WHR, WHtR, BAI, MENO_STATUS, HbA1c) %>%
  drop_na() %>%
  mutate(MenoGroup = factor(if_else(MENO_STATUS == 1, "Postmenopausal", "Premenopausal")))

cat(sprintf("\nDataset loaded: n = %d\n", nrow(df)))

# ==============================================================================
# 2. DESCRIPTIVE STATISTICS (Table 1, Figure 1)
# ==============================================================================
cat("\n", strrep("=", 70), "\n")
cat("TABLE 1. DESCRIPTIVE STATISTICS BY MENOPAUSAL STATUS\n")
cat(strrep("=", 70), "\n")

pre  <- df %>% filter(MenoGroup == "Premenopausal")
post <- df %>% filter(MenoGroup == "Postmenopausal")
cat(sprintf("Premenopausal  n = %d\n", nrow(pre)))
cat(sprintf("Postmenopausal n = %d\n", nrow(post)))

desc_vars  <- c(PREDICTOR_COLS, "HbA1c")
desc_names <- c(FEATURE_NAMES[1:8], "HbA1c")

desc_table <- map2_dfr(desc_vars, desc_names, function(v, vname) {
  t_result <- t.test(pre[[v]], post[[v]])
  p <- t_result$p.value
  sig <- case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
  tibble(Variable = vname,
         Premenopausal  = sprintf("%.2f ± %.2f", mean(pre[[v]]), sd(pre[[v]])),
         Postmenopausal = sprintf("%.2f ± %.2f", mean(post[[v]]), sd(post[[v]])),
         p_value = sprintf("%.4f %s", p, sig))
})
print(desc_table, n = Inf)
write_csv(desc_table, file.path(FIGURES_DIR, "table-1.csv"))

# --- Figure 1: Boxplot comparison ---
plot_vars <- c("BAI", "WHtR", "WHR", "NeckCirc", "HipCirc", "HbA1c")

box_plots <- map(plot_vars, function(v) {
  ggplot(df, aes(x = MenoGroup, y = .data[[v]], fill = MenoGroup)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
    scale_fill_manual(values = c("Premenopausal" = "#4C72B0", "Postmenopausal" = "#DD8452")) +
    labs(title = v, x = NULL, y = NULL) +
    theme(
      legend.position = "none", 
      plot.title = element_text(face = "bold"),
      # FIX 1: Rotate text 30 degrees and shift them slightly right to clear the overlap
      axis.text.x = element_text(angle = 30, hjust = 1) 
    )
})

fig1 <- wrap_plots(box_plots, ncol = 3) +
  plot_annotation(
    theme = theme(plot.title = element_text(face = "bold", size = 15))
  )

# FIX 2: Explicitly define a wider width and height in ggsave so the columns aren't squished
ggsave(
  filename = file.path(FIGURES_DIR, "Figure-1.png"), 
  plot = fig1,
  width = 10,       # Wider page dimension in inches
  height = 7,       # Balanced height for a 3x2 grid
  dpi = 300
)

cat(sprintf("Saved: %s/Figure-1.png\n", FIGURES_DIR))

# ==============================================================================
# 3. CORRELATION ANALYSIS (Table 2, Figures 2-3)
# ==============================================================================
cat("\n", strrep("=", 70), "\n")
cat("TABLE 2. PEARSON CORRELATIONS WITH HbA1c\n")
cat(strrep("=", 70), "\n")

corr_table <- map2_dfr(PREDICTOR_COLS, FEATURE_NAMES[1:8], function(v, vname) {
  ra <- cor.test(df[[v]],  df$HbA1c)
  rp <- cor.test(pre[[v]], pre$HbA1c)
  rq <- cor.test(post[[v]], post$HbA1c)
  tibble(Variable = vname,
         r_all = ra$estimate, p_all = ra$p.value,
         r_pre = rp$estimate, p_pre = rp$p.value,
         r_post = rq$estimate, p_post = rq$p.value)
})
print(corr_table, n = Inf)
write_csv(corr_table, file.path(FIGURES_DIR, "table-2.csv"))

# --- Figure 2: Correlation heatmap ---
corr_matrix <- df %>%
  select(all_of(PREDICTOR_COLS), HbA1c) %>%
  rename_with(~ c(FEATURE_NAMES[1:8], "HbA1c"), everything()) %>%
  cor()

png(file.path(FIGURES_DIR, "Figure-2.png"), width = 9, height = 7, units = "in", res = 300)
corrplot(corr_matrix, method = "color", type = "upper", addCoef.col = "black",
         tl.col = "black", tl.srt = 45, col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
         title = "Figure 2. Correlation Matrix: Adiposity Indices and HbA1c", mar = c(0, 0, 2, 0))
dev.off()
cat(sprintf("Saved: %s/Figure-2.png\n", FIGURES_DIR))

# --- Figure 3: BAI vs HbA1c scatter ---
fig3 <- ggplot(df, aes(x = BAI, y = HbA1c, color = MenoGroup)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
  scale_color_manual(values = c("Premenopausal" = "#4C72B0", "Postmenopausal" = "#DD8452")) +
  labs(
       x = "Body Adiposity Index (BAI, %)", y = "HbA1c (%)", color = NULL) +
  theme(plot.title = element_text(face = "bold", size = 14))
ggsave(file.path(FIGURES_DIR, "Figure-3.png"), fig3, width = 8, height = 6, dpi = 300)
cat(sprintf("Saved: %s/Figure-3.png\n", FIGURES_DIR))
# ==============================================================================
# 4. REGRESSION MODELS (Table 3 and Figure 4)
# ==============================================================================
cat("\n", strrep("=", 70), "\n")
cat("TABLE 3. REGRESSION MODEL PERFORMANCE (5-fold Cross-Validation)\n")
cat(strrep("=", 70), "\n")

reg_data <- df %>% select(all_of(PREDICTOR_COLS), MENO_STATUS, HbA1c)

set.seed(RANDOM_STATE)
reg_folds <- vfold_cv(reg_data, v = N_FOLDS)

reg_recipe <- recipe(HbA1c ~ ., data = reg_data) %>%
  step_normalize(all_numeric_predictors())

lm_spec    <- linear_reg() %>% set_engine("lm") %>% set_mode("regression")
lasso_spec <- linear_reg(penalty = tune(), mixture = 1) %>% set_engine("glmnet") %>% set_mode("regression")
rf_spec    <- rand_forest(trees = 300) %>%
  set_engine("ranger", importance = "impurity", seed = RANDOM_STATE) %>% set_mode("regression")
svr_spec   <- svm_rbf(cost = 1) %>% set_engine("kernlab") %>% set_mode("regression")

reg_metrics <- metric_set(rmse, mae, rsq)

# Linear Regression
lm_wf  <- workflow() %>% add_recipe(reg_recipe) %>% add_model(lm_spec)
lm_res <- fit_resamples(lm_wf, resamples = reg_folds, metrics = reg_metrics)

# LASSO — tune penalty
lasso_wf    <- workflow() %>% add_recipe(reg_recipe) %>% add_model(lasso_spec)
lasso_grid  <- grid_regular(penalty(range = c(-4, 0)), levels = 20)
lasso_tuned <- tune_grid(lasso_wf, resamples = reg_folds, grid = lasso_grid, metrics = reg_metrics)
best_lasso_penalty <- select_best(lasso_tuned, metric = "rmse")
lasso_final_wf <- finalize_workflow(lasso_wf, best_lasso_penalty)
lasso_res <- fit_resamples(lasso_final_wf, resamples = reg_folds, metrics = reg_metrics)

# Random Forest
rf_wf  <- workflow() %>% add_recipe(reg_recipe) %>% add_model(rf_spec)
rf_res <- fit_resamples(rf_wf, resamples = reg_folds, metrics = reg_metrics)

# SVR
svr_wf  <- workflow() %>% add_recipe(reg_recipe) %>% add_model(svr_spec)
svr_res <- fit_resamples(svr_wf, resamples = reg_folds, metrics = reg_metrics)

collect_reg_metrics <- function(res, name) {
  m <- collect_metrics(res)
  tibble(Model = name,
         RMSE = m$mean[m$.metric == "rmse"], RMSE_sd = m$std_err[m$.metric == "rmse"],
         MAE  = m$mean[m$.metric == "mae"],  MAE_sd  = m$std_err[m$.metric == "mae"],
         R2   = m$mean[m$.metric == "rsq"],  R2_sd   = m$std_err[m$.metric == "rsq"])
}

reg_results_df <- bind_rows(
  collect_reg_metrics(lm_res, "Linear Regression"),
  collect_reg_metrics(lasso_res, "LASSO"),
  collect_reg_metrics(rf_res, "Random Forest"),
  collect_reg_metrics(svr_res, "SVR (RBF)")
)

for (i in seq_len(nrow(reg_results_df))) {
  r <- reg_results_df[i, ]
  cat(sprintf("%-20s RMSE=%.3f(±%.3f)  MAE=%.3f(±%.3f)  R²=%.3f(±%.3f)\n",
              r$Model, r$RMSE, r$RMSE_sd, r$MAE, r$MAE_sd, r$R2, r$R2_sd))
}
write_csv(reg_results_df, file.path(FIGURES_DIR, "table-3.csv"))

# --- Figure 4: Regression comparison ---
fig4a <- ggplot(reg_results_df, aes(x = Model, y = RMSE)) +
  geom_col(fill = "#4C72B0") +
  geom_errorbar(aes(ymin = RMSE - RMSE_sd, ymax = RMSE + RMSE_sd), width = 0.2) +
  labs(title = "Regression RMSE (lower is better)", x = NULL) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), plot.title = element_text(face = "bold"))

fig4b <- ggplot(reg_results_df, aes(x = Model, y = R2)) +
  geom_col(fill = "#55A868") +
  geom_errorbar(aes(ymin = R2 - R2_sd, ymax = R2 + R2_sd), width = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  labs(title = "Regression R² (higher is better)", x = NULL) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), plot.title = element_text(face = "bold"))

fig4 <- (fig4a + fig4b) +
  plot_annotation(
                   theme = theme(plot.title = element_text(face = "bold", size = 14)))
ggsave(file.path(FIGURES_DIR, "Figure-4.png"), fig4,
       width = 13, height = 5, dpi = 300)
cat(sprintf("Saved: %s/Figure-4.png\n", FIGURES_DIR))

# ==============================================================================
# 5. FEATURE IMPORTANCE: LASSO COEFFICIENTS + RF IMPURITY (Table 4, Figure 5)
# ==============================================================================
cat("\n", strrep("=", 70), "\n")
cat("TABLE 4. LASSO COEFFICIENTS AND RANDOM FOREST IMPORTANCE (Regression)\n")
cat(strrep("=", 70), "\n")

rf_final  <- rf_wf %>% fit(data = reg_data)
rf_imp    <- rf_final %>% extract_fit_parsnip() %>% vip::vi()

lasso_final_fit <- lasso_final_wf %>% fit(data = reg_data)
lasso_coefs <- lasso_final_fit %>% extract_fit_parsnip() %>% tidy() %>%
  filter(term != "(Intercept)")

name_map <- tibble(
  term = c("Age", "BMI", "NeckCirc", "HipCirc", "WaistCirc", "WHR", "WHtR", "BAI", "MENO_STATUS"),
  Feature = FEATURE_NAMES
)

importance_df <- name_map %>%
  left_join(rf_imp %>% rename(term = Variable, RF_MDI_Importance = Importance), by = "term") %>%
  left_join(lasso_coefs %>% select(term, estimate) %>% rename(LASSO_Coefficient = estimate), by = "term") %>%
  mutate(LASSO_Coefficient = replace_na(LASSO_Coefficient, 0)) %>%
  arrange(desc(RF_MDI_Importance)) %>%
  select(Feature, LASSO_Coefficient, RF_MDI_Importance)

cat(sprintf("\nLASSO penalty (lambda): %.4f\n", best_lasso_penalty$penalty))
print(importance_df, n = Inf)
write_csv(importance_df, file.path(FIGURES_DIR, "table-4.csv"))

# --- Figure 5b: RF feature importance ---
fig5b <- importance_df %>%
  mutate(Feature = fct_reorder(Feature, RF_MDI_Importance)) %>%
  ggplot(aes(x = Feature, y = RF_MDI_Importance)) +
  geom_col(fill = "#4C72B0") + coord_flip() +
  labs(title = "(B) RF Feature Importance",
       x = NULL, y = "Importance (impurity)") +
  theme(plot.title = element_text(face = "bold"))

# --- Figure 5a: LASSO coefficients ---
fig5a <- importance_df %>%
  mutate(Feature = fct_reorder(Feature, LASSO_Coefficient),
         Retained = abs(LASSO_Coefficient) > 1e-6) %>%
  ggplot(aes(x = Feature, y = LASSO_Coefficient, fill = Retained)) +
  geom_col() + coord_flip() +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#CCCCCC"), guide = "none") +
  labs(title = "(A) LASSO Regression Coefficients",
       subtitle = "Grey = eliminated by regularization", x = NULL, y = "Standardized Coefficient") +
  theme(plot.title = element_text(face = "bold", size = 13))

fig5 = fig5a + fig5b
ggsave(file.path(FIGURES_DIR, "Figure-5.png"), fig5,
       width = 13, height = 5, dpi = 300)
cat(sprintf("Saved: %s/figure-5.png\n", FIGURES_DIR))

# ==============================================================================
# 6. SUBGROUP REGRESSION ANALYSIS (Table 5)
# ==============================================================================
cat("\n", strrep("=", 70), "\n")
cat("TABLE 5. SUBGROUP REGRESSION ANALYSIS BY MENOPAUSAL STATUS\n")
cat(strrep("=", 70), "\n")

run_subgroup_reg <- function(subdf, label) {
  sub_reg_data <- subdf %>% select(all_of(PREDICTOR_COLS), HbA1c)
  set.seed(RANDOM_STATE)
  sub_folds <- vfold_cv(sub_reg_data, v = N_FOLDS)

  sub_recipe <- recipe(HbA1c ~ ., data = sub_reg_data) %>% step_normalize(all_numeric_predictors())
  sub_rf_wf <- workflow() %>% add_recipe(sub_recipe) %>%
    add_model(rand_forest(trees = 200) %>% set_engine("ranger", seed = RANDOM_STATE) %>% set_mode("regression"))

  sub_res <- fit_resamples(sub_rf_wf, resamples = sub_folds, metrics = metric_set(rmse, rsq))
  m <- collect_metrics(sub_res)
  rmse_mean <- m$mean[m$.metric == "rmse"]; rmse_sd <- m$std_err[m$.metric == "rmse"]
  r2_mean   <- m$mean[m$.metric == "rsq"];  r2_sd   <- m$std_err[m$.metric == "rsq"]

  sub_rf_full <- rand_forest(trees = 500) %>%
    set_engine("ranger", importance = "impurity", seed = RANDOM_STATE) %>%
    set_mode("regression") %>% fit(HbA1c ~ ., data = sub_reg_data)
  top3 <- sub_rf_full %>% vip::vi() %>% arrange(desc(Importance)) %>% slice_head(n = 3) %>%
    left_join(name_map %>% rename(Variable = term), by = "Variable") %>% pull(Feature)

  cat(sprintf("\n%s (n=%d):\n", label, nrow(subdf)))
  cat(sprintf("  RMSE=%.3f±%.3f, R²=%.3f±%.3f\n", rmse_mean, rmse_sd, r2_mean, r2_sd))
  cat(sprintf("  Top 3 predictors: %s\n", paste(top3, collapse = ", ")))

  tibble(Group = label, n = nrow(subdf), RMSE = rmse_mean, RMSE_sd = rmse_sd,
         R2 = r2_mean, R2_sd = r2_sd, Top3 = paste(top3, collapse = ", "))
}

subgroup_df <- bind_rows(
  run_subgroup_reg(pre, "Premenopausal"),
  run_subgroup_reg(post, "Postmenopausal")
)
write_csv(subgroup_df, file.path(FIGURES_DIR, "table-5.csv"))

cat("\n", strrep("=", 70), "\n")
cat("PAPER A ANALYSIS COMPLETE\n")
cat(sprintf("All figures and tables saved to: ./%s/\n", FIGURES_DIR))
cat(strrep("=", 70), "\n")
