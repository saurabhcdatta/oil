# ============================================================
#  LASSO Logistic Regression with Parallel Processing
#  Binary Y | X1-X10 | Post-LASSO significance filtering
#  + Final AUC and ROC curve
# ============================================================

# ── 0. Install / load packages ───────────────────────────────
pkgs <- c("glmnet", "doParallel", "foreach", "dplyr", "broom", "pROC", "ggplot2")
install.packages(setdiff(pkgs, rownames(installed.packages())))

library(glmnet)
library(doParallel)
library(foreach)
library(dplyr)
library(broom)
library(pROC)
library(ggplot2)

set.seed(42)


# ── 1. Simulate data ──────────────────────────────────────────
# True signal: X1, X3, X5, X7 are associated with Y
# X2, X4, X6, X8, X9, X10 are noise
n <- 500

X <- matrix(rnorm(n * 10), nrow = n,
            dimnames = list(NULL, paste0("X", 1:10)))

log_odds <- -0.5 +
  1.2 * X[, "X1"] +
  -0.8 * X[, "X3"] +
  0.6 * X[, "X5"] +
  1.0 * X[, "X7"]

Y <- rbinom(n, size = 1, prob = plogis(log_odds))

cat("Class balance — Y=1:", round(mean(Y), 3),
    "| Y=0:", round(1 - mean(Y), 3), "\n")


# ── 2. Set up parallel backend ────────────────────────────────
n_cores <- max(1, parallel::detectCores() - 1)
cl      <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Using", n_cores, "cores\n")


# ── 3. Cross-validated LASSO (parallel fold evaluation) ───────
# Three CV replicates run in parallel with different fold seeds
# for a more stable lambda estimate
cv_fits <- foreach(
  seed      = c(42, 7, 123),
  .packages = c("glmnet"),
  .combine  = list,
  .multicombine = TRUE
) %dopar% {
  set.seed(seed)
  cv.glmnet(
    x            = X,
    y            = Y,
    family       = "binomial",   # logistic LASSO
    alpha        = 1,            # 1 = LASSO, 0 = Ridge, (0,1) = Elastic Net
    nfolds       = 10,
    type.measure = "auc",        # optimise on AUC for binary outcome
    standardize  = TRUE          # standardise predictors (recommended)
  )
}

stopCluster(cl)


# ── 4. Choose best lambda ─────────────────────────────────────
# Average lambda.1se across replicates for stability
# lambda.1se = most regularised model within 1 SE of minimum CV error
lambda_1se_vals <- sapply(cv_fits, `[[`, "lambda.1se")
best_lambda     <- mean(lambda_1se_vals)

cat("\nLambda range across replicates:",
    round(range(lambda_1se_vals), 4), "\n")
cat("Averaged lambda.1se used:", round(best_lambda, 4), "\n")

# Final LASSO model on full data with averaged lambda
final_lasso <- glmnet(
  x           = X,
  y           = Y,
  family      = "binomial",
  alpha       = 1,
  lambda      = best_lambda,
  standardize = TRUE
)


# ── 5. Extract non-zero LASSO coefficients ─────────────────────
lasso_coef <- coef(final_lasso)
lasso_df   <- data.frame(
  variable    = rownames(lasso_coef),
  coefficient = as.numeric(lasso_coef)
) %>%
  filter(variable != "(Intercept)",
         coefficient != 0)

cat("\n── Variables selected by LASSO ──────────────────────────\n")
print(lasso_df)


# ── 6. Post-LASSO refit: logistic regression for p-values ──────
# LASSO shrinkage biases standard errors → refit standard logistic
# GLM on selected variables to obtain valid inference
if (nrow(lasso_df) == 0) {
  stop("LASSO selected no variables. Try a smaller lambda.")
}

selected_vars <- lasso_df$variable

refit_formula <- as.formula(
  paste("Y ~", paste(selected_vars, collapse = " + "))
)

refit_data <- as.data.frame(X) %>%
  mutate(Y = Y)

refit_model <- glm(
  formula = refit_formula,
  data    = refit_data,
  family  = binomial(link = "logit")
)


# ── 7. Filter to p < 0.05 ──────────────────────────────────────
refit_summary <- tidy(refit_model, exponentiate = FALSE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    odds_ratio = exp(estimate),
    sig        = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

final_vars <- refit_summary %>%
  filter(p.value < 0.05) %>%
  arrange(p.value)

cat("\n── Full post-LASSO logistic model ───────────────────────\n")
print(refit_summary %>%
        select(term, estimate, std.error, statistic, p.value, odds_ratio, sig),
      digits = 4)

cat("\n── Significant variables (p < 0.05) ─────────────────────\n")
print(final_vars %>%
        select(term, estimate, odds_ratio, std.error, p.value, sig),
      digits = 4)


# ── 8. Rebuild final model using only significant variables ─────
if (nrow(final_vars) > 0) {
  
  sig_vars      <- final_vars$term
  final_formula <- as.formula(
    paste("Y ~", paste(sig_vars, collapse = " + "))
  )
  
  final_model <- glm(
    formula = final_formula,
    data    = refit_data,
    family  = binomial(link = "logit")
  )
  
  cat("\n── Final parsimonious model summary ─────────────────────\n")
  print(tidy(final_model, exponentiate = TRUE, conf.int = TRUE) %>%
          select(term, estimate, conf.low, conf.high, p.value),
        digits = 4)
  
} else {
  cat("No variables met p < 0.05 after refitting.\n")
}


# ── 9. AUC ────────────────────────────────────────────────────
# Use final_model (sig. vars only) if it exists, else refit_model
model_for_roc <- if (exists("final_model")) final_model else refit_model

model_label <- if (exists("final_model")) {
  paste0("Selected variables: ", paste(sig_vars, collapse = ", "))
} else {
  paste0("All LASSO-selected variables: ", paste(selected_vars, collapse = ", "))
}

pred_probs <- predict(model_for_roc, type = "response")

roc_obj <- roc(
  response  = refit_data$Y,
  predictor = pred_probs,
  quiet     = TRUE
)

auc_val <- auc(roc_obj)
ci_val  <- ci.auc(roc_obj, conf.level = 0.95)

cat("\n── Final Model AUC ───────────────────────────────────────\n")
cat(sprintf("AUC:      %.4f\n",         as.numeric(auc_val)))
cat(sprintf("95%% CI:  [%.4f, %.4f]\n", ci_val[1], ci_val[3]))


# ── 10. ROC plot ───────────────────────────────────────────────
roc_df <- data.frame(
  specificity = roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

# Optimal threshold via Youden's J (maximises sensitivity + specificity - 1)
best_idx   <- which.max(roc_obj$sensitivities + roc_obj$specificities - 1)
best_point <- data.frame(
  specificity = roc_obj$specificities[best_idx],
  sensitivity = roc_obj$sensitivities[best_idx],
  threshold   = roc_obj$thresholds[best_idx]
)

auc_label <- sprintf(
  "AUC = %.4f\n95%% CI [%.4f, %.4f]",
  as.numeric(auc_val), ci_val[1], ci_val[3]
)

roc_plot <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  
  # Chance diagonal
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60", linewidth = 0.5) +
  
  # Shaded area under curve
  geom_ribbon(aes(ymin = 0, ymax = sensitivity),
              fill = "#4E6EBF", alpha = 0.08) +
  
  # ROC curve
  geom_line(colour = "#4E6EBF", linewidth = 1) +
  
  # Optimal threshold point
  geom_point(
    data = best_point,
    aes(x = 1 - specificity, y = sensitivity),
    colour = "#C0392B", size = 3
  ) +
  geom_label(
    data = best_point,
    aes(
      x     = 1 - specificity + 0.04,
      y     = sensitivity - 0.07,
      label = sprintf("Optimal threshold\n%.3f", threshold)
    ),
    size = 3, colour = "#C0392B", fill = "white",
    label.size = 0.3, hjust = 0
  ) +
  
  # AUC annotation
  annotate("label",
           x = 0.58, y = 0.18,
           label      = auc_label,
           size       = 3.5,
           colour     = "#2C3E50",
           fill       = "white",
           label.size = 0.4,
           hjust      = 0
  ) +
  
  scale_x_continuous(
    name   = "1 − Specificity (False Positive Rate)",
    limits = c(0, 1), expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    name   = "Sensitivity (True Positive Rate)",
    limits = c(0, 1), expand = c(0.01, 0.01)
  ) +
  labs(
    title    = "ROC Curve — Post-LASSO Logistic Model",
    subtitle = model_label,
    caption  = "Optimal cut-point via Youden's J statistic"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(colour = "grey40", size = 11),
    plot.caption     = element_text(colour = "grey50", size = 9),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "grey85", fill = NA)
  )

print(roc_plot)

# ggsave("roc_curve.png", roc_plot, width = 7, height = 6, dpi = 300)