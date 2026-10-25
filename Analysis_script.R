#============================================================
# Analysis_script.R
# BAI vs waist-based indices for dysglycemia in Hausa women
# Sourced by Article.Rmd; creates every table, figure and
# inline number used in the manuscript.
#============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(splines)
  library(pROC)
  library(glmnet)
  library(ranger)
  library(kernlab)
  library(patchwork)
})

set.seed(2026)
N_BOOT     <- 1000   # bootstrap CIs for delta R2
N_OPT      <- 300    # bootstrap optimism correction for AUC
N_LASSO    <- 500    # bootstrap LASSO selection stability
CV_FOLDS   <- 5
CV_REPEATS <- 10

#------------------------------------------------------------
# 1. Data
#    clean_data.csv holds the measured variables and the derived
#    indices (BMI, WHR, WHtR, BAI) calculated from them:
#    BAI = HC (cm) / height (m)^1.5 - 18 (Bergman et al., 2011)
#------------------------------------------------------------
d <- read.csv("clean_data.csv") %>%
  transmute(
    id,
    Age      = age_years,
    Status   = factor(menopausal_status, levels = c("Premenopausal", "Postmenopausal")),
    Menarche = age_at_menarche_years,
    Ht       = height_cm,
    Wt       = weight_kg,
    WC       = waist_circumference_cm,
    HC       = hip_circumference_cm,
    NC       = neck_circumference_cm,
    BMI      = bmi_kg_m2,
    WHR      = whr,
    WHtR     = whtr,
    BAI      = bai_pct,
    HbA1c    = hba1c_pct,
    post     = as.integer(Status == "Postmenopausal"),
    dysgly   = as.integer(HbA1c >= 5.7),   # prediabetes or diabetes range (ADA)
    dmr      = as.integer(HbA1c >= 6.5)    # diabetes range
  )

IDX <- c("BAI", "WHtR", "WC", "HC", "NC", "BMI", "WHR")
IDX_LAB <- c(BAI = "BAI", WHtR = "WHtR", WC = "WC",
             HC = "HC", NC = "NC", BMI = "BMI", WHR = "WHR")
ABBREV <- "BAI, body adiposity index; WHtR, waist-to-height ratio; WC, waist circumference; HC, hip circumference; NC, neck circumference; BMI, body mass index; WHR, waist-to-hip ratio."
for (v in IDX) d[[paste0(v, "_z")]] <- as.numeric(scale(d[[v]]))

pre  <- filter(d, post == 0)
postm <- filter(d, post == 1)


#------------------------------------------------------------
# helpers
#------------------------------------------------------------
fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
fmt   <- function(x, k = 2) formatC(x, format = "f", digits = k)
msd   <- function(x, k = 2) paste0(fmt(mean(x), k), " ± ", fmt(sd(x), k))
ci    <- function(est, lo, hi, k = 2) paste0(fmt(est, k), " (", fmt(lo, k), " to ", fmt(hi, k), ")")

# Meng, Rosenthal & Rubin (1992) test for two dependent correlations sharing a variable
mrr_test <- function(r12, r13, r23, n) {
  rm2 <- (r12^2 + r13^2) / 2
  f   <- min((1 - r23) / (2 * (1 - rm2)), 1)
  h   <- (1 - f * rm2) / (1 - rm2)
  Z   <- (atanh(r12) - atanh(r13)) * sqrt((n - 3) / (2 * (1 - r23) * h))
  c(Z = Z, p = 2 * pnorm(-abs(Z)))
}

nri_idi <- function(y, p0, p1) {
  up <- p1 > p0; dn <- p1 < p0
  nri <- (mean(up[y == 1]) - mean(dn[y == 1])) + (mean(dn[y == 0]) - mean(up[y == 0]))
  idi <- (mean(p1[y == 1]) - mean(p0[y == 1])) - (mean(p1[y == 0]) - mean(p0[y == 0]))
  c(NRI = nri, IDI = idi)
}

quiet_roc <- function(y, p) roc(y, p, quiet = TRUE, direction = "<", levels = c(0, 1))

#------------------------------------------------------------
# 2. Table 1: descriptives
#------------------------------------------------------------
desc_vars <- c(Age = "Age (years)", Ht = "Height (cm)", Wt = "Weight (kg)",
               BMI = "BMI (kg/m²)", WC = "Waist circumference (cm)",
               HC = "Hip circumference (cm)", NC = "Neck circumference (cm)",
               WHR = "Waist-to-hip ratio", WHtR = "Waist-to-height ratio",
               BAI = "BAI (%)", HbA1c = "HbA1c (%)")
desc_digits <- c(Age = 1, Ht = 1, Wt = 1, BMI = 1, WC = 1, HC = 1, NC = 1,
                 WHR = 2, WHtR = 2, BAI = 1, HbA1c = 2)

table1 <- map_dfr(names(desc_vars), function(v) {
  k <- desc_digits[[v]]
  tibble(Variable = desc_vars[[v]],
         `All (n = 249)` = msd(d[[v]], k),
         Premenopausal   = msd(pre[[v]], k),
         Postmenopausal  = msd(postm[[v]], k),
         `P value`       = fmt_p(t.test(pre[[v]], postm[[v]])$p.value))
})
cat_row <- function(lab, var) {
  pct <- function(x) sprintf("%d (%.1f)", sum(x), 100 * mean(x))
  tibble(Variable = lab, `All (n = 249)` = pct(d[[var]]),
         Premenopausal = pct(pre[[var]]), Postmenopausal = pct(postm[[var]]),
         `P value` = fmt_p(chisq.test(table(d$post, d[[var]]))$p.value))
}
table1 <- bind_rows(table1,
                    cat_row("HbA1c ≥5.7% (dysglycemia), n (%)", "dysgly"),
                    cat_row("HbA1c ≥6.5% (diabetes range), n (%)", "dmr"))
names(table1)[3:4] <- c(sprintf("Premenopausal (n = %d)", nrow(pre)),
                        sprintf("Postmenopausal (n = %d)", nrow(postm)))

#------------------------------------------------------------
# 3. Correlations, dependent-correlation tests, Fisher z
#------------------------------------------------------------
corr_res <- map_dfr(list(All = d, Premenopausal = pre, Postmenopausal = postm), .id = "group", function(s) {
  r_bai <- cor(s$BAI, s$HbA1c)
  map_dfr(IDX, function(v) {
    ct <- cor.test(s[[v]], s$HbA1c); ct$estimate <- unname(ct$estimate)
    t  <- if (v == "BAI") c(Z = NA, p = NA) else mrr_test(r_bai, ct$estimate, cor(s$BAI, s[[v]]), nrow(s))
    tibble(index = v, r = unname(ct$estimate), lo = ct$conf.int[1], hi = ct$conf.int[2],
           p = ct$p.value, p_vs_BAI = unname(t["p"]))
  })
})

fisher_res <- map_dfr(IDX, function(v) {
  a <- filter(corr_res, group == "Premenopausal", index == v)
  b <- filter(corr_res, group == "Postmenopausal", index == v)
  Z <- (atanh(a$r) - atanh(b$r)) / sqrt(1 / (nrow(pre) - 3) + 1 / (nrow(postm) - 3))
  tibble(index = v, Z = Z, p = 2 * pnorm(-abs(Z)))
})

table2 <- corr_res %>%
  mutate(cell = paste0(fmt(r, 2), " (", fmt(lo, 2), ", ", fmt(hi, 2), ")"),
         pv   = ifelse(index == "BAI", "Reference", fmt_p(p_vs_BAI))) %>%
  select(group, index, cell, pv) %>%
  pivot_wider(names_from = group, values_from = c(cell, pv)) %>%
  left_join(fisher_res %>% transmute(index, `Pre vs post, P` = fmt_p(p)), by = "index") %>%
  transmute(Index = IDX_LAB[index],
            `All: r (95% CI)` = cell_All, `All: P vs BAI` = pv_All,
            `Pre: r (95% CI)` = cell_Premenopausal, `Pre: P vs BAI` = pv_Premenopausal,
            `Post: r (95% CI)` = cell_Postmenopausal, `Post: P vs BAI` = pv_Postmenopausal,
            `Pre vs post, P`)

#------------------------------------------------------------
# 4. Incremental value over age + menopausal status (linear)
#------------------------------------------------------------
m_base <- lm(HbA1c ~ Age + post, d)

boot_dr2 <- function(v) {
  f1 <- as.formula(paste0("HbA1c ~ Age + post + ", v, "_z"))
  replicate(N_BOOT, {
    s <- d[sample(nrow(d), replace = TRUE), ]
    summary(lm(f1, s))$r.squared - summary(lm(HbA1c ~ Age + post, s))$r.squared
  })
}

inc_res <- map_dfr(IDX, function(v) {
  term <- paste0(v, "_z")
  m  <- lm(as.formula(paste0("HbA1c ~ Age + post + ", term)), d)
  mi <- lm(as.formula(paste0("HbA1c ~ Age + post * ", term)), d)
  cf <- confint(m)[term, ]
  b  <- boot_dr2(v)
  tibble(index = v, beta = coef(m)[term], lo = cf[1], hi = cf[2],
         p = summary(m)$coefficients[term, 4],
         R2 = summary(m)$r.squared, dR2 = R2 - summary(m_base)$r.squared,
         dR2_lo = quantile(b, .025), dR2_hi = quantile(b, .975),
         AIC = AIC(m),
         p_int = summary(mi)$coefficients[paste0("post:", term), 4])
})

m_bai <- lm(HbA1c ~ Age + post + BAI_z, d)
beyond_bai <- map_dfr(setdiff(IDX, "BAI"), function(v) {
  term <- paste0(v, "_z")
  m <- lm(as.formula(paste0("HbA1c ~ Age + post + BAI_z + ", term)), d)
  m0 <- lm(as.formula(paste0("HbA1c ~ Age + post + ", term)), d)
  tibble(index = v,
         p_index_given_BAI = summary(m)$coefficients[term, 4],
         dR2_index_given_BAI = summary(m)$r.squared - summary(m_bai)$r.squared,
         p_BAI_given_index = summary(m)$coefficients["BAI_z", 4],
         dR2_BAI_given_index = summary(m)$r.squared - summary(m0)$r.squared)
})

table3 <- inc_res %>%
  left_join(beyond_bai, by = "index") %>%
  transmute(Index = IDX_LAB[index],
            `β per SD (95% CI)` = ci(beta, lo, hi),
            P = fmt_p(p),
            `ΔR² vs age + status (95% CI)` = ci(dR2, dR2_lo, dR2_hi, 3),
            AIC = fmt(AIC, 1),
            `ΔR² when added to BAI model` = ifelse(is.na(dR2_index_given_BAI), "—",
                                                  paste0(fmt(dR2_index_given_BAI, 3), " (P = ", fmt_p(p_index_given_BAI), ")")),
            `Index × status, P` = fmt_p(p_int))

#------------------------------------------------------------
# 5. Discrimination of dysglycemia (HbA1c >=5.7%) and diabetes range
#------------------------------------------------------------
# 5a. adjusted models in the whole sample
auc_adj <- function(outcome) {
  y  <- d[[outcome]]
  g0 <- glm(as.formula(paste(outcome, "~ Age + post")), binomial, d)
  p0 <- fitted(g0); r0 <- quiet_roc(y, p0)
  gb <- glm(as.formula(paste(outcome, "~ Age + post + BAI_z")), binomial, d)
  rb <- quiet_roc(y, fitted(gb))
  map_dfr(IDX, function(v) {
    term <- paste0(v, "_z")
    g  <- glm(as.formula(paste(outcome, "~ Age + post +", term)), binomial, d)
    p1 <- fitted(g); r1 <- quiet_roc(y, p1)
    aci <- ci.auc(r1, method = "delong")
    or  <- exp(coef(g)[term]); orci <- exp(confint.default(g)[term, ])
    ni  <- nri_idi(y, p0, p1)
    tibble(outcome = outcome, index = v, events = sum(y),
           OR = or, OR_lo = orci[1], OR_hi = orci[2],
           AUC_base = as.numeric(auc(r0)), AUC = as.numeric(aci[2]),
           AUC_lo = aci[1], AUC_hi = aci[3],
           p_vs_base = roc.test(r1, r0, method = "delong", paired = TRUE)$p.value,
           dAUC_BAI_minus = as.numeric(auc(rb)) - as.numeric(auc(r1)),
           p_vs_BAI = if (v == "BAI") NA_real_ else roc.test(rb, r1, method = "delong", paired = TRUE)$p.value,
           NRI = ni["NRI"], IDI = ni["IDI"])
  })
}
auc_dys <- auc_adj("dysgly")
auc_dm  <- auc_adj("dmr")

# 5b. index alone within strata (DeLong vs BAI, Youden cut-off)
auc_alone <- function(s, outcome) {
  y <- s[[outcome]]; rb <- quiet_roc(y, s$BAI)
  map_dfr(IDX, function(v) {
    r <- quiet_roc(y, s[[v]]); a <- ci.auc(r, method = "delong")
    yo <- coords(r, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))[1, ]
    tibble(index = v, events = sum(y), n = length(y),
           AUC = as.numeric(a[2]), lo = a[1], hi = a[3],
           p_vs_BAI = if (v == "BAI") NA_real_ else roc.test(rb, r, method = "delong", paired = TRUE)$p.value,
           cut = yo$threshold, sens = yo$sensitivity, spec = yo$specificity)
  })
}
auc_pre_dys  <- auc_alone(pre, "dysgly")
auc_post_dm  <- auc_alone(postm, "dmr")
roc_pre <- map(set_names(IDX), ~ quiet_roc(pre$dysgly, pre[[.x]]))

# 5c. bootstrap optimism-corrected AUC (Harrell)
opt_auc <- function(formula, outcome) {
  fit <- glm(formula, binomial, d)
  app <- as.numeric(auc(quiet_roc(d[[outcome]], fitted(fit))))
  opt <- replicate(N_OPT, {
    s <- d[sample(nrow(d), replace = TRUE), ]
    g <- suppressWarnings(glm(formula, binomial, s))
    as.numeric(auc(quiet_roc(s[[outcome]], fitted(g)))) -
      as.numeric(auc(quiet_roc(d[[outcome]], predict(g, d, type = "response"))))
  })
  c(apparent = app, corrected = app - mean(opt))
}
opt_models <- c("Age + status" = "~ Age + post", "Age + status + BAI" = "~ Age + post + BAI_z",
                "Age + status + WHtR" = "~ Age + post + WHtR_z", "Age + status + WC" = "~ Age + post + WC_z")
opt_res <- map_dfr(c("dysgly", "dmr"), function(o) map_dfr(names(opt_models), function(nm) {
  r <- opt_auc(as.formula(paste(o, opt_models[[nm]])), o)
  tibble(outcome = o, model = nm, apparent = r["apparent"], corrected = r["corrected"])
}))

table4 <- auc_dys %>%
  left_join(auc_dm %>% select(index, AUC_dm = AUC, lo_dm = AUC_lo, hi_dm = AUC_hi, p_dm = p_vs_BAI), by = "index") %>%
  left_join(auc_pre_dys %>% select(index, AUC_pre = AUC, lo_pre = lo, hi_pre = hi, p_pre = p_vs_BAI), by = "index") %>%
  transmute(Index = IDX_LAB[index],
            `Dysglycemia: OR per SD (95% CI)` = ci(OR, OR_lo, OR_hi),
            `Dysglycemia: AUC (95% CI)` = ci(AUC, AUC_lo, AUC_hi, 3),
            `P vs BAI` = ifelse(index == "BAI", "Reference", fmt_p(p_vs_BAI)),
            `Diabetes range: AUC (95% CI)` = ci(AUC_dm, lo_dm, hi_dm, 3),
            `P vs BAI ` = ifelse(index == "BAI", "Reference", fmt_p(p_dm)),
            `Premenopausal, index alone: AUC (95% CI)` = ci(AUC_pre, lo_pre, hi_pre, 3),
            `P vs BAI  ` = ifelse(index == "BAI", "Reference", fmt_p(p_pre)))

#------------------------------------------------------------
# 6. Menopausal status, age and adiposity
#------------------------------------------------------------
status_fit <- function(label, formula, data) {
  m <- lm(formula, data); cf <- confint(m)["post", ]
  tibble(Model = label, n = nrow(data), beta = coef(m)["post"], lo = cf[1], hi = cf[2],
         p = summary(m)$coefficients["post", 4])
}
ov <- filter(d, Age >= 45, Age <= 55)
status_res <- bind_rows(
  status_fit("Unadjusted", HbA1c ~ post, d),
  status_fit("Adjusted for age (linear)", HbA1c ~ post + Age, d),
  status_fit("Adjusted for age (natural spline, 3 df)", HbA1c ~ post + ns(Age, 3), d),
  status_fit("Adjusted for age and BAI", HbA1c ~ post + Age + BAI_z, d),
  status_fit("Adjusted for age, BAI and neck circumference", HbA1c ~ post + Age + BAI_z + NC_z, d),
  status_fit("Age-overlap subsample (45–55 y), adjusted for age and BAI", HbA1c ~ post + Age + BAI_z, ov)
)
table5 <- status_res %>% transmute(Model, n, `Difference in HbA1c, % (95% CI)` = ci(beta, lo, hi), P = fmt_p(p))

age_slopes <- map_dfr(list(Premenopausal = pre, Postmenopausal = postm), .id = "group", function(s) {
  m <- lm(HbA1c ~ Age, s); cf <- confint(m)["Age", ] * 10
  tibble(slope10 = coef(m)["Age"] * 10, lo = cf[1], hi = cf[2], p = summary(m)$coefficients["Age", 4])
})
d$Ht_z <- as.numeric(scale(d$Ht))
hcht <- summary(lm(HbA1c ~ Age + post + HC_z + Ht_z, d))$coefficients
vif_post <- 1 / (1 - summary(lm(post ~ Age + BAI_z, d))$r.squared)
ov_n <- table(ov$Status)

#------------------------------------------------------------
# 7. Supplementary: repeated CV model comparison
#------------------------------------------------------------
FEAT <- c("BAI", "WHtR", "WC", "HC", "NC", "BMI", "WHR", "Age", "post")
folds <- map(seq_len(CV_REPEATS), ~ sample(rep(1:CV_FOLDS, length.out = nrow(d))))
metric <- function(y, p) c(RMSE = sqrt(mean((y - p)^2)), MAE = mean(abs(y - p)), R2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2))
fit_predict <- list(
  `Linear regression (all predictors)` = function(tr, te) predict(lm(reformulate(FEAT, "HbA1c"), tr), te),
  `LASSO` = function(tr, te) {
    X <- as.matrix(tr[FEAT]); cvf <- cv.glmnet(X, tr$HbA1c, alpha = 1, nfolds = 5, standardize = TRUE)
    as.numeric(predict(cvf, as.matrix(te[FEAT]), s = "lambda.min"))
  },
  `Random forest` = function(tr, te) predict(ranger(reformulate(FEAT, "HbA1c"), tr, num.trees = 500, min.node.size = 5), te)$predictions,
  `SVR (RBF kernel)` = function(tr, te) as.numeric(predict(ksvm(as.matrix(scale(tr[FEAT])), tr$HbA1c, type = "eps-svr", kernel = "rbfdot"),
                                                         scale(as.matrix(te[FEAT]), center = colMeans(tr[FEAT]), scale = apply(tr[FEAT], 2, sd)))),
  `Age + status (linear)` = function(tr, te) predict(lm(HbA1c ~ Age + post, tr), te),
  `Age + status + BAI (linear)` = function(tr, te) predict(lm(HbA1c ~ Age + post + BAI, tr), te)
)
cv_raw <- map_dfr(seq_len(CV_REPEATS), function(r) map_dfr(1:CV_FOLDS, function(k) {
  tr <- d[folds[[r]] != k, ]; te <- d[folds[[r]] == k, ]
  map_dfr(names(fit_predict), function(nm) {
    m <- metric(te$HbA1c, fit_predict[[nm]](tr, te))
    tibble(rep = r, fold = k, model = nm, RMSE = m["RMSE"], MAE = m["MAE"], R2 = m["R2"])
  })
}))
ref <- filter(cv_raw, model == "Age + status + BAI (linear)")
n_te <- nrow(d) / CV_FOLDS; n_tr <- nrow(d) - n_te
nb_p <- function(nm) {   # Nadeau–Bengio corrected resampled t-test vs reference
  if (nm == "Age + status + BAI (linear)") return(NA_real_)
  dl <- filter(cv_raw, model == nm)$R2 - ref$R2; J <- length(dl)
  t <- mean(dl) / sqrt(var(dl) * (1 / J + n_te / n_tr)); 2 * pt(-abs(t), J - 1)
}
tableS1 <- cv_raw %>% group_by(model) %>%
  summarise(across(c(RMSE, MAE, R2), list(m = mean, s = sd)), .groups = "drop") %>%
  mutate(order = match(model, names(fit_predict))) %>% arrange(order) %>%
  rowwise() %>% mutate(p = nb_p(model)) %>% ungroup() %>%
  transmute(Model = model, RMSE = paste0(fmt(RMSE_m, 3), " ± ", fmt(RMSE_s, 3)),
            MAE = paste0(fmt(MAE_m, 3), " ± ", fmt(MAE_s, 3)),
            `R²` = paste0(fmt(R2_m, 3), " ± ", fmt(R2_s, 3)),
            `P vs age + status + BAI` = ifelse(is.na(p), "Reference", fmt_p(p)))
cv_sum <- cv_raw %>% group_by(model) %>% summarise(R2 = mean(R2), .groups = "drop")

#------------------------------------------------------------
# 8. Supplementary: LASSO stability and RF permutation importance
#------------------------------------------------------------
X_all <- as.matrix(d[FEAT])
lasso_full <- cv.glmnet(X_all, d$HbA1c, alpha = 1, nfolds = 10)
lasso_boot <- replicate(N_LASSO, {
  i <- sample(nrow(d), replace = TRUE)
  f <- cv.glmnet(X_all[i, ], d$HbA1c[i], alpha = 1, nfolds = 5)
  as.numeric(coef(f, s = "lambda.min"))[-1] != 0
})
rf_full <- ranger(reformulate(FEAT, "HbA1c"), d, num.trees = 2000, min.node.size = 5,
                  importance = "permutation", seed = 2026)
FEAT_LAB <- c(BAI = "BAI", WHtR = "WHtR", WC = "Waist circumference", HC = "Hip circumference",
              NC = "Neck circumference", BMI = "BMI", WHR = "WHR", Age = "Age", post = "Menopausal status")
tableS2 <- tibble(Predictor = FEAT_LAB[FEAT],
                  `LASSO selection frequency (%)` = fmt(100 * rowMeans(lasso_boot), 1),
                  `RF permutation importance` = fmt(rf_full$variable.importance[FEAT], 3)) %>%
  arrange(desc(as.numeric(`LASSO selection frequency (%)`)))
lasso_freq <- setNames(100 * rowMeans(lasso_boot), FEAT)

#------------------------------------------------------------
# 9. Supplementary: decision curve (dysglycemia)
#------------------------------------------------------------
dca <- function(outcome) {
  y <- d[[outcome]]; prev <- mean(y)
  mods <- list(`Age + status` = "~ Age + post", `Age + status + BAI` = "~ Age + post + BAI_z",
               `Age + status + WHtR` = "~ Age + post + WHtR_z")
  P <- map(mods, ~ fitted(glm(as.formula(paste(outcome, .x)), binomial, d)))
  map_dfr(seq(0.05, 0.80, 0.01), function(t) {
    nb <- map_dbl(P, function(p) mean(p >= t & y == 1) - mean(p >= t & y == 0) * t / (1 - t))
    tibble(threshold = t, `Treat all` = prev - (1 - prev) * t / (1 - t), `Treat none` = 0, !!!as.list(nb))
  })
}
dca_dys <- dca("dysgly")

#------------------------------------------------------------
# 10. Figures
#------------------------------------------------------------
theme_set(theme_minimal(base_size = 10) +
            theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10)))
COL <- c(Premenopausal = "#4C72B0", Postmenopausal = "#DD8452")
IDX_COL <- c(BAI = "#C0392B", WHtR = "#2E86C1", WC = "#17A589",
             HC = "#8E44AD", NC = "#D68910", BMI = "#7F8C8D", WHR = "#34495E")

# Figure 1: HbA1c against age by status (step at menopause, flat within groups)
fig1 <- ggplot(d, aes(Age, HbA1c, colour = Status)) +
  annotate("rect", xmin = 45, xmax = 55, ymin = -Inf, ymax = Inf, alpha = .08) +
  geom_hline(yintercept = c(5.7, 6.5), linetype = "dashed", colour = "grey55", linewidth = .3) +
  geom_point(alpha = .6, size = 1.4) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, linewidth = .8) +
  scale_colour_manual(values = COL, name = NULL) +
  labs(x = "Age (years)", y = "HbA1c (%)") +
  theme(legend.position = "top")

# Figure 2: delta R2 forest plot
fig2 <- inc_res %>% mutate(lab = factor(IDX_LAB[index], levels = rev(IDX_LAB[index[order(-dR2)]]))) %>%
  ggplot(aes(dR2, lab)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_errorbarh(aes(xmin = dR2_lo, xmax = dR2_hi), height = .2) +
  geom_point(aes(colour = index == "BAI"), size = 2.5) +
  scale_colour_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "black"), guide = "none") +
  labs(x = "Increase in R² over age + menopausal status (95% bootstrap CI)", y = NULL)

# Figure 3: ROC curves
roc_df <- function(rlist, labs) map_dfr(names(rlist), function(v) {
  r <- rlist[[v]]
  tibble(index = paste0(labs[[v]], " (AUC ", fmt(as.numeric(auc(r)), 2), ")"),
         key = labs[[v]], fpr = 1 - r$specificities, tpr = r$sensitivities)
}) %>% arrange(index, fpr, tpr)
show <- c("BAI", "WHtR", "WC", "HC", "NC")
roc_adj <- map(set_names(show), ~ quiet_roc(d$dysgly, fitted(glm(as.formula(paste0("dysgly ~ Age + post + ", .x, "_z")), binomial, d))))
roc_adj[["Base"]] <- quiet_roc(d$dysgly, fitted(glm(dysgly ~ Age + post, binomial, d)))
labs_adj <- c(as.list(IDX_LAB[show]), Base = "Age + status only")
roc_plot <- function(df, title) {
  cols <- c(IDX_COL, `Age + status only` = "grey50")
  ggplot(df, aes(fpr, tpr, colour = key, group = index)) +
    geom_abline(linetype = "dotted", colour = "grey70") +
    geom_step(direction = "vh", linewidth = .6) +
    scale_colour_manual(values = cols, breaks = unique(df$key),
                        labels = unique(df$index[match(unique(df$key), df$key)]), name = NULL) +
    coord_equal() + labs(x = "1 – specificity", y = "Sensitivity", title = title) +
    theme(legend.position = c(.68, .25), legend.text = element_text(size = 7),
          legend.key.height = unit(.35, "cm"), legend.background = element_rect(fill = "white", colour = NA))
}
fig3 <- roc_plot(roc_df(roc_adj, labs_adj), "(A) All women, models with age + status") +
  roc_plot(roc_df(roc_pre[show], as.list(IDX_LAB[show])), "(B) Premenopausal women, index alone")

# Supplementary Figure S1: decision curve
figS1 <- dca_dys %>% pivot_longer(-threshold, names_to = "Strategy", values_to = "nb") %>%
  ggplot(aes(threshold, nb, colour = Strategy, linetype = Strategy)) +
  geom_line(linewidth = .6) +
  coord_cartesian(ylim = c(-0.05, max(dca_dys$`Treat all`) + .05)) +
  scale_colour_manual(values = c(`Treat all` = "grey60", `Treat none` = "black", `Age + status` = "grey30",
                                 `Age + status + BAI` = "#C0392B", `Age + status + WHtR` = "#2E86C1")) +
  scale_linetype_manual(values = c(`Treat all` = "dashed", `Treat none` = "solid", `Age + status` = "dotdash",
                                   `Age + status + BAI` = "solid", `Age + status + WHtR` = "solid")) +
  labs(x = "Threshold probability", y = "Net benefit")

#------------------------------------------------------------
# 11. Inline values used in the text
#------------------------------------------------------------
g <- function(df, ...) filter(df, ...)
v_corr  <- function(grp, v) g(corr_res, group == grp, index == v)
v_inc   <- function(v) g(inc_res, index == v)
v_dys   <- function(v) g(auc_dys, index == v)
v_dm    <- function(v) g(auc_dm, index == v)
v_pre   <- function(v) g(auc_pre_dys, index == v)
v_stat  <- function(i) status_res[i, ]
v_opt   <- function(o, m) g(opt_res, outcome == o, model == m)
nb_at   <- function(t, col) dca_dys[[col]][which.min(abs(dca_dys$threshold - t))]

#------------------------------------------------------------
# 12. Export figures as separate files for submission
#------------------------------------------------------------
dir.create("figures", showWarnings = FALSE)
save_fig <- function(plot, name, w, h) {
  ggsave(file.path("figures", paste0(name, ".pdf")), plot, width = w, height = h, device = cairo_pdf)
  ggsave(file.path("figures", paste0(name, ".tiff")), plot, width = w, height = h, dpi = 300, compression = "lzw")
}
save_fig(fig3,  "Figure_1", 8, 4.3)            # ROC curves
save_fig(fig1,  "Figure_2", 6, 3.8)            # HbA1c by age and status
save_fig(fig2,  "Supplementary_Figure_S1", 6, 3)
save_fig(figS1, "Supplementary_Figure_S2", 6, 3.5)
