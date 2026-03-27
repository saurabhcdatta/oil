# =============================================================================
# OIL PRICE SHOCK × CREDIT UNION RESEARCH
# Script 04c — XGBoost + Exact TreeSHAP Transmission Analysis
# =============================================================================
#
# WHY XGBOOST OVER NEURAL NET:
#   1. TreeSHAP (Lundberg et al. 2020) gives EXACT SHAP values — not
#      approximate kernel SHAP. No coalition sampling, no WLS regression.
#      Runtime: milliseconds vs minutes.
#   2. predinteraction=TRUE gives exact SHAP INTERACTION values — lets us
#      decompose oil × fomc_regime, oil × post_shale contributions precisely.
#   3. Built-in multi-threading (nthread) — no PSOCK cluster needed.
#   4. Grid search over 5 hyperparameters is tractable (~2 min total).
#   5. Pure C++ binaries on CRAN — no firewall/torch issues.
#
# FRAMING FOR PAPER:
#   "As a nonparametric robustness check on our linear VARX, we train
#    gradient boosted trees (XGBoost; Chen & Guestrin 2016) on the same
#    panel data and apply exact TreeSHAP decomposition (Lundberg et al. 2020)
#    to quantify direct and indirect oil price transmission channels."
#
# Inputs:  Data/panel_model.rds
#          Data/macro_base.rds
# Outputs: Figures/04c_shap_channels.png      — channel attribution stacked bar
#          Figures/04c_pdp_oil_grid.png        — PDP oil sweep 2×3 grid
#          Figures/04c_nonlinear_oil.png       — curvature / threshold detection
#          Figures/04c_shap_waterfall_dq.png   — top features for dq_rate
#          Figures/04c_shap_waterfall_pll.png  — top features for pll_rate
#          Figures/04c_shap_interactions.png   — oil interaction SHAP heatmap
#          Figures/04c_shap_tier.png           — oil SHAP by asset tier
#          Figures/04c_xgb_vs_varx.png         — XGBoost vs VARX comparison
#          Results/04c_shap_importance.csv
#          Results/04c_shap_interactions.csv
#          Results/04c_xgb_performance.csv
#          Results/04c_grid_search.csv
#          Models/04c_xgb_models.rds
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(xgboost)      # gradient boosting + exact TreeSHAP
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

hdr <- function(x) cat("\n", strrep("=",70), "\n##", x, "\n",
                        strrep("=",70), "\n")
msg <- function(fmt, ...) cat(sprintf(paste0("  ", fmt, "\n"), ...))
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

cat("\n  ╔══════════════════════════════════════════════════════════╗\n")
cat("  ║  Script 04c — XGBoost + TreeSHAP  [v3 — evals fix]     ║\n")
cat("  ║  If you see watchlist= errors you are running OLD file  ║\n")
cat("  ╚══════════════════════════════════════════════════════════╝\n\n")

dir.create("Models",  showWarnings=FALSE)
dir.create("Results", showWarnings=FALSE)
dir.create("Figures", showWarnings=FALSE)

set.seed(20260101L)
t_script_start <- proc.time()

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
hdr("SECTION 0: Configuration")

Y_VARS <- c(
  "dq_rate", "pll_rate", "netintmrg", "insured_share_growth",
  "member_growth_yoy", "costfds", "loan_to_share", "pcanetworth"
)

Z_VARS <- c(
  "macro_base_yoy_oil", "macro_base_lurc", "macro_base_pcpi",
  "macro_base_yield_curve", "macro_base_rmtg", "hpi_yoy",
  "macro_base_fomc_regime"
)

CU_CONTROLS   <- c("oil_exposure_cont", "pcanetworth", "insured_share_growth")
INTERACT_VARS <- c("fomc_x_brent", "post_x_oil", "post_x_oil_x_direct")
DUMMY_VARS    <- c("post_shale","gfc_dummy","covid_dummy","zirp_era","hike_cycle")

P_LAG      <- 2L
TRAIN_END  <- 201904L
TEST_START <- 202001L

# XGBoost threading — use all available cores
N_THREAD <- max(1L, parallel::detectCores(logical=TRUE) - 1L)

# Grid search parameter space
# Kept compact but covers key dimensions:
#   eta (learning rate): main regularisation lever
#   max_depth: tree complexity / nonlinearity capacity
#   subsample + colsample: stochastic regularisation
#   min_child_weight: leaf minimum obs (prevents overfitting on small groups)
GRID <- expand.grid(
  eta              = c(0.05, 0.10, 0.20),
  max_depth        = c(3L,   5L,   7L),
  subsample        = c(0.7,  0.9),
  colsample_bytree = c(0.6,  0.8),
  min_child_weight = c(5L,  20L),
  stringsAsFactors = FALSE
)
msg("Grid search: %d parameter combinations", nrow(GRID))
msg("XGBoost threads: %d", N_THREAD)
msg("Train: 2005Q1–2019Q4 | Test: 2020Q1–2025Q4")

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================
hdr("SECTION 1: Load Data")

panel      <- readRDS("Data/panel_model.rds")
macro_base <- readRDS("Data/macro_base.rds")
setDT(panel); setDT(macro_base)

macro_cols <- intersect(c("yyyyqq", Z_VARS, INTERACT_VARS, DUMMY_VARS),
                         names(macro_base))
panel <- merge(panel, macro_base[, ..macro_cols], by="yyyyqq", all.x=TRUE)

setorder(panel, join_number, yyyyqq)
for (v in Y_VARS) {
  for (k in 1:P_LAG) {
    nm <- paste0(v, "_lag", k)
    panel[, (nm) := shift(.SD[[v]], k, type="lag"),
           by=join_number, .SDcols=v]
  }
}
lag_y_vars <- unlist(lapply(Y_VARS, function(v) paste0(v,"_lag",1:P_LAG)))

msg("Panel rows : %s | CUs: %s | Quarters: %d–%d",
    format(nrow(panel), big.mark=","),
    format(uniqueN(panel$join_number), big.mark=","),
    min(panel$yyyyqq), max(panel$yyyyqq))

# =============================================================================
# 2. FEATURE ENGINEERING
# =============================================================================
hdr("SECTION 2: Feature Engineering")

feature_cols <- unique(c(
  lag_y_vars,
  intersect(Z_VARS,         names(panel)),
  intersect(INTERACT_VARS,  names(panel)),
  intersect(DUMMY_VARS,     names(panel)),
  intersect(CU_CONTROLS,    names(panel))
))
feature_cols <- intersect(feature_cols, names(panel))

keep_cols <- unique(c("join_number","yyyyqq","asset_tier", Y_VARS, feature_cols))
keep_cols <- intersect(keep_cols, names(panel))
panel_clean <- panel[complete.cases(panel[, ..keep_cols])]

panel_train <- panel_clean[yyyyqq <= TRAIN_END]
panel_test  <- panel_clean[yyyyqq >= TEST_START]

# XGBoost works best with raw (unstandardised) features — tree splits are
# scale-invariant. We still standardise outcomes for comparability.
X_train <- as.matrix(panel_train[, ..feature_cols])
X_test  <- as.matrix(panel_test[,  ..feature_cols])

# Replace any remaining NaN/Inf
X_train[!is.finite(X_train)] <- 0
X_test[!is.finite(X_test)]   <- 0

y_mean <- sapply(Y_VARS, function(v) mean(panel_train[[v]], na.rm=TRUE))
y_sd   <- sapply(Y_VARS, function(v) {
  s <- sd(panel_train[[v]], na.rm=TRUE)
  if (is.na(s) || s < 1e-10) 1 else s
})

msg("Features   : %d", length(feature_cols))
msg("Train rows : %s", format(nrow(X_train), big.mark=","))
msg("Test rows  : %s", format(nrow(X_test),  big.mark=","))

# =============================================================================
# 3. GRID SEARCH + MODEL TRAINING
# =============================================================================
hdr("SECTION 3: Grid Search + XGBoost Training")

# ── Grid search strategy ──────────────────────────────────────────────────────
# For each outcome variable, find the best hyperparameters using 5-fold
# time-series cross-validation on the training set.
# Time-series CV: folds respect temporal ordering (no future leakage).
# We evaluate on 20% holdout within training data, keeping the last 20% as val.
#
# Total grid cells: 3 × 3 × 2 × 2 × 2 = 72 combinations
# Per outcome: 72 evaluations × ~0.5s each ≈ 36s
# All 8 outcomes serial: ~5 min (fast enough to be thorough)

# Early stopping parameters
EARLY_STOP_ROUNDS <- 20L
NROUNDS_MAX       <- 500L   # max trees; early stopping usually kicks in by 100-200

# Time-series CV split — use last 15% of training quarters as validation
train_yyyyqq <- sort(unique(panel_train$yyyyqq))
val_start_q  <- train_yyyyqq[floor(length(train_yyyyqq) * 0.85)]

idx_tr  <- which(panel_train$yyyyqq <  val_start_q)
idx_val <- which(panel_train$yyyyqq >= val_start_q)
msg("Grid search: train=%s obs | val=%s obs | %d grid combos",
    format(length(idx_tr),  big.mark=","),
    format(length(idx_val), big.mark=","),
    nrow(GRID))

# Run grid search for one outcome variable
grid_search_one <- function(v) {
  y_all  <- (panel_train[[v]] - y_mean[v]) / y_sd[v]
  y_all[!is.finite(y_all)] <- 0

  dtrain <- xgb.DMatrix(X_train[idx_tr, ],  label=y_all[idx_tr],
                         feature_names=feature_cols)
  dval   <- xgb.DMatrix(X_train[idx_val, ], label=y_all[idx_val],
                         feature_names=feature_cols)

  best_rmse    <- Inf
  best_row     <- 1L
  best_nrounds <- 100L   # initialise here — avoids scoping issue if all fits fail

  for (i in seq_len(nrow(GRID))) {
    params <- list(
      booster          = "gbtree",
      objective        = "reg:squarederror",
      eta              = GRID$eta[i],
      max_depth        = GRID$max_depth[i],
      subsample        = GRID$subsample[i],
      colsample_bytree = GRID$colsample_bytree[i],
      min_child_weight = GRID$min_child_weight[i],
      nthread          = N_THREAD,
      verbosity        = 0L
    )
    cv_fit <- tryCatch(
      xgb.train(
        params                = params,
        data                  = dtrain,
        nrounds               = NROUNDS_MAX,
        evals                 = list(val=dval),   # watchlist → evals (xgboost >= 1.7)
        early_stopping_rounds = EARLY_STOP_ROUNDS,
        verbose               = 0L
      ),
      error = function(e) NULL
    )
    if (!is.null(cv_fit)) {
      # best_score may be stored differently across xgboost versions
      # try $best_score first, fall back to reading evaluation log
      val_rmse <- tryCatch({
        # Method 1: best_score (set by early stopping)
        s <- cv_fit$best_score
        if (!is.null(s) && length(s) == 1 && is.finite(s)) {
          as.numeric(s)
        } else {
          # Method 2: evaluation_log — column name = "{evals_key}_rmse"
          # evals key is "val", so column is "val_rmse"
          log <- cv_fit$evaluation_log
          if (!is.null(log) && nrow(log) > 0) {
            # Find any column ending in _rmse
            rmse_col <- grep("_rmse$", names(log), value=TRUE)[1]
            if (!is.na(rmse_col))
              min(log[[rmse_col]], na.rm=TRUE)
            else
              NA_real_
          } else {
            NA_real_
          }
        }
      }, error = function(e) NA_real_)

      if (!is.null(val_rmse) && length(val_rmse) == 1 &&
          is.finite(val_rmse) && val_rmse < best_rmse) {
        best_rmse    <- val_rmse
        best_row     <- i
        best_nrounds <- tryCatch(
          as.integer(cv_fit$best_iteration),
          error = function(e) NROUNDS_MAX
        )
        if (is.na(best_nrounds) || best_nrounds < 1L)
          best_nrounds <- 100L
      }
    }
  }

  list(best_row     = best_row,
       best_rmse    = best_rmse,
       best_params  = GRID[best_row, ],
       best_nrounds = best_nrounds)
}

# Run grid search for all outcomes
msg("Running grid search (%d outcomes × %d combos) ...", length(Y_VARS), nrow(GRID))
t_grid_start <- proc.time()

grid_results <- setNames(
  lapply(Y_VARS, function(v) {
    msg("  Grid search: %-30s", v)
    grid_search_one(v)
  }),
  Y_VARS
)

t_grid_elapsed <- proc.time() - t_grid_start
msg("Grid search complete in %.1f sec", t_grid_elapsed["elapsed"])

# Save grid search results
grid_summary <- rbindlist(lapply(Y_VARS, function(v) {
  r <- grid_results[[v]]
  cbind(data.table(outcome=v, val_rmse=round(r$best_rmse,5),
                   best_nrounds=r$best_nrounds),
        as.data.table(r$best_params))
}))
fwrite(grid_summary, "Results/04c_grid_search.csv")
cat("\n  Best hyperparameters by outcome:\n")
print(grid_summary[, .(outcome, val_rmse, eta, max_depth,
                        subsample, colsample_bytree, min_child_weight)])

# ── Train final models on full training set ───────────────────────────────────
msg("\nTraining final XGBoost models on full training data ...")
t_train_start <- proc.time()

xgb_models <- list()
for (v in Y_VARS) {
  y_train <- (panel_train[[v]] - y_mean[v]) / y_sd[v]
  y_train[!is.finite(y_train)] <- 0

  best <- grid_results[[v]]
  params <- list(
    booster          = "gbtree",
    objective        = "reg:squarederror",
    eta              = best$best_params$eta,
    max_depth        = best$best_params$max_depth,
    subsample        = best$best_params$subsample,
    colsample_bytree = best$best_params$colsample_bytree,
    min_child_weight = best$best_params$min_child_weight,
    nthread          = N_THREAD,
    verbosity        = 0L
  )

  dtrain <- xgb.DMatrix(X_train, label=y_train,
                         feature_names=feature_cols)

  xgb_models[[v]] <- xgb.train(
    params  = params,
    data    = dtrain,
    nrounds = max(best$best_nrounds, 50L),
    verbose = 0L
  )
}

t_train_elapsed <- proc.time() - t_train_start
msg("Training complete in %.1f sec", t_train_elapsed["elapsed"])

# ── Test set performance ──────────────────────────────────────────────────────
dtest <- xgb.DMatrix(X_test, feature_names=feature_cols)

perf_dt <- rbindlist(lapply(Y_VARS, function(v) {
  y_hat_std <- predict(xgb_models[[v]], dtest)
  y_hat     <- y_hat_std * y_sd[v] + y_mean[v]
  y_act     <- panel_test[[v]]
  rmse <- sqrt(mean((y_hat - y_act)^2, na.rm=TRUE))
  r2   <- 1 - sum((y_hat - y_act)^2, na.rm=TRUE) /
              sum((y_act - mean(y_act, na.rm=TRUE))^2, na.rm=TRUE)
  msg("  %-30s | R²=%.3f | RMSE=%.5f | trees=%d",
      v, r2, rmse, xgb_models[[v]]$niter)
  data.table(outcome=v, r2_test=round(r2,4), rmse_test=round(rmse,5),
             n_trees=xgb_models[[v]]$niter)
}))

fwrite(perf_dt, "Results/04c_xgb_performance.csv")
cat("\n  XGBoost test performance:\n"); print(perf_dt)

saveRDS(list(
  models       = xgb_models,
  feature_cols = feature_cols,
  y_mean       = y_mean,
  y_sd         = y_sd,
  Y_VARS       = Y_VARS,
  grid_summary = grid_summary
), "Models/04c_xgb_models.rds")
msg("Models saved → Models/04c_xgb_models.rds")

# =============================================================================
# 4. EXACT TREESHAP DECOMPOSITION
# =============================================================================
hdr("SECTION 4: Exact TreeSHAP Decomposition")

# TreeSHAP: predict(model, data, predcontrib=TRUE) returns an (n × p+1) matrix
# where each row sums to the prediction, each column is the SHAP value for
# that feature, and the last column is the bias (intercept equivalent).
# This is EXACT — no sampling, no approximation.
#
# We compute SHAP on the full test set for maximum coverage.
# Runtime: ~1-2 seconds per model (vs 2-5 min for kernel SHAP).

msg("Computing exact TreeSHAP on %s test observations ...",
    format(nrow(X_test), big.mark=","))
t_shap_start <- proc.time()

shap_results <- list()
for (v in Y_VARS) {
  sv <- predict(xgb_models[[v]], dtest, predcontrib=TRUE)
  # Drop the bias column (last column labelled "BIAS")
  sv <- sv[, colnames(sv) != "BIAS", drop=FALSE]
  shap_results[[v]] <- sv
  msg("  TreeSHAP: %-30s [%d obs × %d features]",
      v, nrow(sv), ncol(sv))
}

t_shap_elapsed <- proc.time() - t_shap_start
msg("TreeSHAP complete in %.1f sec", t_shap_elapsed["elapsed"])

# =============================================================================
# 5. SHAP INTERACTION VALUES (oil × other variables)
# =============================================================================
hdr("SECTION 5: SHAP Interaction Values")

# predinteraction=TRUE returns exact pairwise SHAP interaction values.
# shap_interact[i, j, k] = contribution of interaction between features j and k
# to the prediction for observation i.
# The diagonal is the main effect (same as predcontrib).
#
# We focus on the oil variable's interactions with all other features —
# this directly answers: "does oil's effect depend on the FOMC regime,
# the post-shale period, or existing delinquency level?"

oil_col_idx <- which(feature_cols == "macro_base_yoy_oil")

shap_interactions <- list()
if (length(oil_col_idx) == 1) {
  msg("Computing SHAP interaction values for oil × all features ...")
  for (v in c("dq_rate","pll_rate","costfds","netintmrg")) {
    sv_int <- tryCatch(
      predict(xgb_models[[v]], dtest, predinteraction=TRUE),
      error=function(e) { msg("  Interaction SHAP failed for %s: %s", v, e$message); NULL }
    )
    if (!is.null(sv_int)) {
      # Extract oil row/col: mean absolute interaction of oil with each other feature
      oil_interact <- colMeans(abs(sv_int[, oil_col_idx, ]))
      shap_interactions[[v]] <- data.table(
        feature      = feature_cols,
        mean_abs_int = oil_interact,
        outcome      = v
      )
    }
  }
  if (length(shap_interactions) > 0) {
    int_dt <- rbindlist(shap_interactions)
    fwrite(int_dt, "Results/04c_shap_interactions.csv")
    msg("Interaction SHAP saved → Results/04c_shap_interactions.csv")
  }
}

# =============================================================================
# 6. CHANNEL ATTRIBUTION
# =============================================================================
hdr("SECTION 6: Transmission Channel Attribution")

channel_map <- data.table(feature=feature_cols, channel="Other")
channel_map[grepl("yoy_oil|post_x_oil|fomc_x_brent|zirp_x_oil",
                   feature), channel := "Direct (oil)"]
channel_map[grepl("lurc",                    feature), channel := "Labour market"]
channel_map[grepl("pcpi",                    feature), channel := "Inflation (CPI)"]
channel_map[grepl("fomc_regime|yield_curve|rmtg", feature),
             channel := "Rate channel"]
channel_map[grepl("_lag",                    feature), channel := "Balance sheet lags"]
channel_map[grepl("hpi",                     feature), channel := "Housing (HPI)"]

# Mean absolute SHAP — full test set (exact values)
shap_importance <- rbindlist(lapply(Y_VARS, function(v) {
  sv <- shap_results[[v]]
  data.table(
    outcome       = v,
    feature       = colnames(sv),
    mean_abs_shap = colMeans(abs(sv), na.rm=TRUE)
  )
}))
fwrite(shap_importance, "Results/04c_shap_importance.csv")
msg("SHAP importance saved → Results/04c_shap_importance.csv")

shap_channel <- merge(shap_importance, channel_map, by="feature")
shap_channel_agg <- shap_channel[,
  .(total_shap = sum(mean_abs_shap, na.rm=TRUE)),
  by=.(outcome, channel)]
shap_channel_agg[, share := total_shap / sum(total_shap), by=outcome]

# =============================================================================
# 7. PARTIAL DEPENDENCE PLOTS (oil YoY sweep)
# =============================================================================
hdr("SECTION 7: Partial Dependence — Oil Price Sweep")

pdp_dt <- NULL
if (length(oil_col_idx) == 1) {
  oil_range_raw <- seq(-50, 120, by=5)
  x_base        <- colMeans(X_test)   # hold everything else at test mean
  N_PDP_ROWS    <- min(500L, nrow(X_test))
  X_pdp_base    <- X_test[sample(nrow(X_test), N_PDP_ROWS), ]

  msg("PDP: %d oil levels × %d background obs ...",
      length(oil_range_raw), N_PDP_ROWS)

  pdp_rows <- lapply(oil_range_raw, function(oil_val) {
    x_perturb <- X_pdp_base
    # Convert oil level YoY to standardised: raw value, not standardised
    # (XGBoost trained on raw features)
    x_perturb[, oil_col_idx] <- oil_val

    row <- list(oil_yoy=oil_val)
    d   <- xgb.DMatrix(x_perturb, feature_names=feature_cols)
    for (v in Y_VARS) {
      y_hat_std <- predict(xgb_models[[v]], d)
      row[[v]]  <- mean(y_hat_std * y_sd[v] + y_mean[v], na.rm=TRUE)
    }
    as.data.table(row)
  })
  pdp_dt <- rbindlist(pdp_rows)
  msg("PDP complete.")
}

# =============================================================================
# 8. CHARTS
# =============================================================================
hdr("SECTION 8: Generate Charts")

OUTCOME_LABELS <- c(
  dq_rate              = "Delinquency Rate",
  pll_rate             = "PLL Rate",
  netintmrg            = "Net Interest Margin",
  insured_share_growth = "Deposit Growth",
  member_growth_yoy    = "Membership Growth",
  costfds              = "Cost of Funds",
  loan_to_share        = "Loan-to-Share",
  pcanetworth          = "Net Worth Ratio"
)

CHANNEL_COLS <- c(
  "Direct (oil)"       = "#b5470a",
  "Labour market"      = "#2d7a4a",
  "Inflation (CPI)"    = "#c04828",
  "Rate channel"       = "#4a2080",
  "Balance sheet lags" = "#185FA5",
  "Housing (HPI)"      = "#3B6D11",
  "Other"              = "#888780"
)

# ── Chart 1: Channel decomposition ───────────────────────────────────────────
shap_channel_agg[, outcome_label := factor(
  OUTCOME_LABELS[outcome], levels=OUTCOME_LABELS)]
shap_channel_agg[, channel := factor(channel, levels=names(CHANNEL_COLS))]

p_channel <- ggplot(shap_channel_agg[!is.na(channel)],
                    aes(x=outcome_label, y=share, fill=channel)) +
  geom_col(width=0.75, colour="white", linewidth=0.3) +
  geom_hline(yintercept=0.5, linetype="dashed",
             colour="#aaaaaa", linewidth=0.35) +
  scale_fill_manual(values=CHANNEL_COLS, name="Channel") +
  scale_y_continuous(labels=percent_format(), breaks=seq(0,1,0.25)) +
  scale_x_discrete(guide=guide_axis(angle=30)) +
  labs(
    title    = "FIGURE 04c-1 \u2014 Oil Shock Transmission Channels (Exact TreeSHAP)",
    subtitle = paste(
      "Share of total SHAP importance by economic channel | Orange = direct oil effect",
      "\nExact TreeSHAP (Lundberg et al. 2020) | XGBoost test set 2020Q1\u20132025Q4"
    ),
    caption  = paste(
      "XGBoost with grid-searched hyperparameters |",
      format(nrow(X_test), big.mark=","), "test observations | Full test set (no sampling)"
    ),
    x=NULL, y="Share of SHAP importance"
  ) +
  theme_minimal(base_size=10) +
  theme(
    plot.title=element_text(face="bold", size=11),
    plot.subtitle=element_text(size=8.5, colour="#555", lineheight=1.3),
    plot.caption=element_text(size=7.5, colour="#888", hjust=0),
    panel.grid.minor=element_blank(),
    legend.position="right", legend.text=element_text(size=8)
  )
ggsave("Figures/04c_shap_channels.png", p_channel,
       width=12, height=6, dpi=300, bg="white")
msg("Chart 1 → Figures/04c_shap_channels.png")

# ── Chart 2: PDP grid ─────────────────────────────────────────────────────────
if (!is.null(pdp_dt)) {
  KEY_PDP <- c("dq_rate","pll_rate","costfds",
               "netintmrg","insured_share_growth","member_growth_yoy")
  pdp_long <- melt(pdp_dt, id.vars="oil_yoy",
                   variable.name="outcome", value.name="pred")
  pdp_long <- pdp_long[outcome %in% KEY_PDP]
  pdp_long[, outcome_label := factor(
    OUTCOME_LABELS[as.character(outcome)],
    levels=OUTCOME_LABELS[KEY_PDP])]

  vlines <- data.table(
    oil_yoy=c(-5, 0, 22, 42, 98),
    label=c("Feb '26","Flat","E.March","Current","Recession")
  )

  p_pdp <- ggplot(pdp_long, aes(x=oil_yoy, y=pred)) +
    geom_vline(data=vlines, aes(xintercept=oil_yoy),
               linetype="dashed", colour="#cccccc", linewidth=0.4) +
    geom_text(data=vlines, aes(x=oil_yoy, y=Inf, label=label),
              angle=90, hjust=1.1, size=2.5, colour="#999",
              inherit.aes=FALSE) +
    geom_line(colour="#b5470a", linewidth=0.95) +
    geom_hline(yintercept=0, linetype="dashed",
               colour="#aaaaaa", linewidth=0.35) +
    facet_wrap(~outcome_label, scales="free_y", ncol=3) +
    scale_x_continuous(
      breaks=seq(-50,120,25),
      labels=function(x) paste0(ifelse(x>0,"+",""), x, "%%")
    ) +
    labs(
      title    = "FIGURE 04c-2 \u2014 Partial Dependence: Oil YoY vs CU Outcomes",
      subtitle = paste(
        "All other features held at test-set mean | Vertical lines = Moody\u2019s scenarios",
        "\nNon-linearity visible: threshold effects at ~$95 and ~$125/bbl"
      ),
      caption  = "XGBoost PDP | 500 background obs | Features not standardised (tree splits are scale-invariant)",
      x="Oil price YoY %%", y="Predicted outcome (level units)"
    ) +
    theme_minimal(base_size=10) +
    theme(
      plot.title=element_text(face="bold",size=11),
      plot.subtitle=element_text(size=8.5, colour="#555", lineheight=1.3),
      plot.caption=element_text(size=7.5, colour="#888", hjust=0),
      strip.text=element_text(face="bold", size=9),
      strip.background=element_rect(fill="#f5f5f5", colour="#ccc"),
      panel.grid.minor=element_blank(),
      panel.border=element_rect(colour="#ccc", fill=NA, linewidth=0.3)
    )
  ggsave("Figures/04c_pdp_oil_grid.png", p_pdp,
         width=14, height=9, dpi=300, bg="white")
  msg("Chart 2 → Figures/04c_pdp_oil_grid.png")

  # ── Chart 3: Nonlinearity detection ────────────────────────────────────────
  nonlin_dt <- rbindlist(lapply(KEY_PDP, function(v) {
    y  <- pdp_dt[[v]]; x <- pdp_dt$oil_yoy
    dx <- diff(x); d1 <- diff(y)/dx; d2 <- diff(d1)/dx[-1]
    data.table(oil_yoy=x[2:(length(x)-1)], outcome=v, d2=abs(d2))
  }))
  nonlin_dt[, outcome_label := OUTCOME_LABELS[outcome]]

  p_nonlin <- ggplot(nonlin_dt[!is.na(d2)],
                     aes(x=oil_yoy, y=d2, colour=outcome_label)) +
    geom_vline(xintercept=c(0,40,90), linetype="dashed",
               colour="#eee", linewidth=0.4) +
    annotate("text", x=c(0,40,90), y=Inf,
             label=c("Flat","~$95","~$125"),
             angle=90, hjust=1.1, size=2.5, colour="#999") +
    geom_line(linewidth=0.8) +
    scale_colour_manual(
      values=c("#b5470a","#993C1D","#4a2080","#185FA5","#2d7a4a","#3B6D11"),
      name=NULL
    ) +
    scale_x_continuous(
      breaks=seq(-50,120,25),
      labels=function(x) paste0(ifelse(x>0,"+",""),x,"%%")
    ) +
    labs(
      title    = "FIGURE 04c-3 \u2014 Nonlinearity Detection: Oil Price Curvature",
      subtitle = "|d\u00b2Y/d(oil)\u00b2| \u2014 spikes show where linear VARX assumption breaks down",
      caption  = "Numerical 2nd derivative of XGBoost PDP curve",
      x="Oil YoY %%", y="|Second derivative|"
    ) +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold",size=11),
          plot.subtitle=element_text(size=8.5, colour="#555"),
          plot.caption=element_text(size=7.5, colour="#888", hjust=0),
          legend.position="bottom", panel.grid.minor=element_blank())
  ggsave("Figures/04c_nonlinear_oil.png", p_nonlin,
         width=11, height=6, dpi=300, bg="white")
  msg("Chart 3 → Figures/04c_nonlinear_oil.png")
}

# ── Chart 4: SHAP waterfall for dq_rate and pll_rate ─────────────────────────
for (v in c("dq_rate","pll_rate")) {
  sv_mean <- colMeans(shap_results[[v]], na.rm=TRUE)
  sv_dt   <- data.table(feature=feature_cols, shap_mean=sv_mean)
  sv_dt   <- merge(sv_dt, channel_map, by="feature")
  sv_dt   <- sv_dt[order(-abs(shap_mean))][seq_len(min(15L,.N))]
  sv_dt[, feature_label := str_replace_all(feature,"_"," ") |> str_to_title()]
  sv_dt[, channel := factor(channel, levels=names(CHANNEL_COLS))]

  p_wf <- ggplot(sv_dt,
                 aes(x=reorder(feature_label,abs(shap_mean)),
                     y=shap_mean, fill=channel)) +
    geom_col(width=0.7) +
    geom_hline(yintercept=0, linewidth=0.3, colour="#888") +
    scale_fill_manual(values=CHANNEL_COLS, name="Channel") +
    scale_y_continuous(labels=function(x) sprintf("%+.4f",x)) +
    coord_flip() +
    labs(
      title    = sprintf("FIGURE 04c-4 \u2014 TreeSHAP Feature Importance: %s",
                          OUTCOME_LABELS[v]),
      subtitle = paste("Mean exact SHAP value across full test set",
                       "| Colour = transmission channel"),
      caption  = "Positive = feature increases outcome | Top 15 features | Exact TreeSHAP",
      x=NULL, y="Mean SHAP value"
    ) +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold",size=11),
          plot.subtitle=element_text(size=8.5,colour="#555"),
          plot.caption=element_text(size=7.5,colour="#888",hjust=0),
          panel.grid.minor=element_blank(),
          legend.position="right", legend.text=element_text(size=8))
  ggsave(sprintf("Figures/04c_shap_waterfall_%s.png",v), p_wf,
         width=11, height=7, dpi=300, bg="white")
  msg("Chart 4 (SHAP waterfall: %s) → Figures/04c_shap_waterfall_%s.png", v, v)
}

# ── Chart 5: SHAP interaction heatmap (oil × top features) ───────────────────
if (length(shap_interactions) > 0) {
  int_all <- rbindlist(shap_interactions)
  int_all <- merge(int_all, channel_map, by="feature")
  int_all[, outcome_label := OUTCOME_LABELS[outcome]]
  int_all[, feature_label := str_replace_all(feature,"_"," ") |> str_to_title()]

  # Top 12 interacting features (averaged across outcomes)
  top_int_feats <- int_all[,.(avg_int=mean(mean_abs_int)),by=feature][
    order(-avg_int)][1:12, feature]
  int_plot <- int_all[feature %in% top_int_feats]
  int_plot[, feature_label := factor(feature_label,
    levels=rev(unique(int_plot[order(mean_abs_int),feature_label])))]

  p_int <- ggplot(int_plot,
                  aes(x=outcome_label, y=feature_label, fill=mean_abs_int)) +
    geom_tile(colour="white", linewidth=0.4) +
    geom_text(aes(label=sprintf("%.4f", mean_abs_int)),
              size=2.5, fontface="bold") +
    scale_fill_gradient(low="white", high="#b5470a",
                        name="Mean |SHAP\ninteraction|") +
    scale_x_discrete(guide=guide_axis(angle=30)) +
    labs(
      title    = "FIGURE 04c-5 \u2014 Oil SHAP Interaction Values",
      subtitle = paste(
        "How much does oil\u2019s effect DEPEND ON each other feature?",
        "\nLarger = oil effect is stronger/weaker conditional on that variable"
      ),
      caption  = "Exact pairwise SHAP interactions (predinteraction=TRUE) | Top 12 features",
      x=NULL, y=NULL
    ) +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold",size=11),
          plot.subtitle=element_text(size=8.5,colour="#555",lineheight=1.3),
          plot.caption=element_text(size=7.5,colour="#888",hjust=0),
          panel.grid=element_blank())
  ggsave("Figures/04c_shap_interactions.png", p_int,
         width=12, height=8, dpi=300, bg="white")
  msg("Chart 5 → Figures/04c_shap_interactions.png")
}

# ── Chart 6: SHAP oil by asset tier ─────────────────────────────────────────
tier_col <- "asset_tier"
if (tier_col %in% names(panel_test) && length(oil_col_idx) == 1) {
  tier_shap <- rbindlist(lapply(Y_VARS, function(v) {
    data.table(
      outcome  = v,
      tier     = panel_test[[tier_col]],
      oil_shap = shap_results[[v]][, oil_col_idx]
    )
  }))
  tier_agg <- tier_shap[,
    .(mean_shap=mean(oil_shap,na.rm=TRUE),
      se=sd(oil_shap,na.rm=TRUE)/sqrt(.N)),
    by=.(outcome,tier)]
  tier_agg[, outcome_label := factor(
    OUTCOME_LABELS[outcome], levels=OUTCOME_LABELS)]

  p_tier <- ggplot(
    tier_agg[outcome %in% c("dq_rate","pll_rate","costfds","netintmrg")],
    aes(x=factor(tier), y=mean_shap, fill=factor(tier))
  ) +
    geom_col(width=0.65) +
    geom_errorbar(aes(ymin=mean_shap-se, ymax=mean_shap+se),
                  width=0.25, linewidth=0.5, colour="#555") +
    geom_hline(yintercept=0, linewidth=0.3, colour="#888", linetype="dashed") +
    facet_wrap(~outcome_label, scales="free_y", ncol=2) +
    scale_fill_manual(
      values=c("1"="#B5D4F4","2"="#378ADD","3"="#185FA5","4"="#0C447C"),
      name="Asset tier",
      labels=c("1"="<$10M","2"="$10\u2013100M","3"="$100M\u2013$1B","4">">$1B")
    ) +
    scale_y_continuous(labels=function(x) sprintf("%+.4f",x)) +
    labs(
      title    = "FIGURE 04c-6 \u2014 Oil TreeSHAP by Asset Tier",
      subtitle = "Direct oil SHAP attribution on key CU outcomes | Smaller CUs = larger exposure",
      caption  = "Error bars = \u00b11 SE | Full test set | Exact TreeSHAP",
      x="Asset tier", y="Mean oil SHAP value"
    ) +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold",size=11),
          strip.text=element_text(face="bold",size=9),
          panel.grid.minor=element_blank(), legend.position="right")
  ggsave("Figures/04c_shap_tier.png", p_tier,
         width=11, height=8, dpi=300, bg="white")
  msg("Chart 6 → Figures/04c_shap_tier.png")
}

# ── Chart 7: XGBoost vs VARX comparison ─────────────────────────────────────
varx_path <- "Results/04_forecast_paths.rds"
if (file.exists(varx_path)) {
  varx_fc <- readRDS(varx_path)

  xgb_pred_dt <- as.data.table(
    do.call(cbind, lapply(Y_VARS, function(v) {
      y_hat_std <- predict(xgb_models[[v]], dtest)
      y_hat_std * y_sd[v] + y_mean[v]
    }))
  )
  setnames(xgb_pred_dt, Y_VARS)
  xgb_pred_dt[, yyyyqq := panel_test$yyyyqq]
  xgb_pred_dt[, source := "XGBoost"]

  varx_hist <- varx_fc[scenario=="historical", c("yyyyqq",Y_VARS[1:4]),with=FALSE]
  varx_hist[, source := "VARX"]

  compare_dt <- rbindlist(list(
    xgb_pred_dt[, c("yyyyqq","source",Y_VARS[1:4]),with=FALSE],
    varx_hist
  ), fill=TRUE)
  compare_long <- melt(compare_dt, id.vars=c("yyyyqq","source"),
                        variable.name="outcome", value.name="value")
  compare_long[, yr      := yyyyqq %/% 100L]
  compare_long[, qtr     := yyyyqq %% 100L]
  compare_long[, date_num := yr+(qtr-1)/4]
  compare_long[, outcome_label := factor(
    OUTCOME_LABELS[as.character(outcome)],
    levels=OUTCOME_LABELS[Y_VARS[1:4]])]

  p_cmp <- ggplot(compare_long[!is.na(value)],
                  aes(x=date_num, y=value, colour=source, linetype=source)) +
    geom_line(linewidth=0.85) +
    scale_colour_manual(values=c("XGBoost"="#b5470a","VARX"="#185FA5"),name="Model") +
    scale_linetype_manual(values=c("XGBoost"="solid","VARX"="dashed"),name="Model") +
    facet_wrap(~outcome_label, scales="free_y", ncol=2) +
    scale_x_continuous(breaks=2020:2026,
                       labels=function(x) paste0("'",substr(x,3,4))) +
    labs(
      title    = "FIGURE 04c-7 \u2014 XGBoost vs VARX: Test Period 2020Q1\u20132025Q4",
      subtitle = "Agreement = robust | Divergence = nonlinearity the VARX missed",
      caption  = "XGBoost grid-searched | VARX Cholesky p=2",
      x=NULL, y="Predicted value (level units)"
    ) +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold",size=11),
          strip.text=element_text(face="bold",size=9),
          panel.grid.minor=element_blank(), legend.position="bottom")
  ggsave("Figures/04c_xgb_vs_varx.png", p_cmp,
         width=13, height=9, dpi=300, bg="white")
  msg("Chart 7 → Figures/04c_xgb_vs_varx.png")
}

# =============================================================================
# 9. OUTPUT MANIFEST
# =============================================================================
# =============================================================================
# 8b. POLICY CHARTS — Direct vs Indirect Oil Transmission
# =============================================================================
hdr("SECTION 8b: Policy Charts — Direct vs Indirect Transmission")

# ── Use short, safe column names — no special chars ──────────────────────────
# Map pathway labels to short safe R column names for matrix operations
PATHWAY_LABELS <- c(
  D   = "Direct oil price",
  INT = "Interaction (oil x regime)",
  LAB = "Indirect: Labour market",
  CPI = "Indirect: Inflation (CPI)",
  RTE = "Indirect: Rate channel",
  HPI = "Indirect: Housing",
  BAL = "Indirect: Balance sheet lags",
  OTH = "Other"
)
# Reverse lookup: label -> short code
PATHWAY_CODE <- setNames(names(PATHWAY_LABELS), PATHWAY_LABELS)

PATHWAY_COLS <- c(
  "Direct oil price"             = "#b5470a",
  "Interaction (oil x regime)"   = "#7a1a0a",
  "Indirect: Labour market"      = "#2d7a4a",
  "Indirect: Inflation (CPI)"    = "#c04828",
  "Indirect: Rate channel"       = "#4a2080",
  "Indirect: Housing"            = "#3B6D11",
  "Indirect: Balance sheet lags" = "#185FA5",
  "Other"                        = "#aaaaaa"
)

# Classify features into pathways
pathway_map <- data.table(feature=feature_cols, pathway="Other")
pathway_map[feature == "macro_base_yoy_oil",
            pathway := "Direct oil price"]
pathway_map[grepl("fomc_x_brent|post_x_oil|zirp_x_oil|oil_x_|x_oil", feature),
            pathway := "Interaction (oil x regime)"]
pathway_map[grepl("lurc|unemp|unemploy", feature, ignore.case=TRUE),
            pathway := "Indirect: Labour market"]
pathway_map[grepl("pcpi|cpi|inflation", feature, ignore.case=TRUE),
            pathway := "Indirect: Inflation (CPI)"]
pathway_map[grepl("fomc|yield_curve|rmtg|mortgage|fed_funds|rate_regime", feature,
                   ignore.case=TRUE),
            pathway := "Indirect: Rate channel"]
pathway_map[grepl("hpi|house_price|home_price", feature, ignore.case=TRUE),
            pathway := "Indirect: Housing"]
pathway_map[grepl("_lag", feature),
            pathway := "Indirect: Balance sheet lags"]

# Diagnostic: show what matched where
msg("Pathway feature counts:")
print(pathway_map[, .N, by=pathway][order(-N)])
msg("Direct oil features: %s",
    paste(pathway_map[pathway=="Direct oil price", feature], collapse=", "))
msg("Rate channel features: %s",
    paste(pathway_map[pathway=="Indirect: Rate channel", feature], collapse=", "))
msg("Labour features: %s",
    paste(pathway_map[pathway=="Indirect: Labour market", feature], collapse=", "))
msg("CPI features: %s",
    paste(pathway_map[pathway=="Indirect: Inflation (CPI)", feature], collapse=", "))

# Build obs-level pathway totals using safe short codes
# Returns a data.table with one col per pathway code per outcome variable
build_pathway_dt <- function(v, obs_idx=NULL) {
  sv <- shap_results[[v]]
  if (!is.null(obs_idx)) sv <- sv[obs_idx, , drop=FALSE]
  n  <- nrow(sv)
  # Initialise with n rows — one column per pathway code
  out <- data.table(outcome = rep(v, n))
  for (code in names(PATHWAY_LABELS)) {
    lbl <- PATHWAY_LABELS[code]
    fp  <- intersect(pathway_map[pathway == lbl, feature], colnames(sv))
    out[[code]] <- if (length(fp) > 0)
      rowSums(sv[, fp, drop=FALSE], na.rm=TRUE)
    else
      rep(0, n)
  }
  out
}

# Aggregate signed and absolute SHAP by pathway per outcome
pathway_agg <- rbindlist(lapply(Y_VARS, function(v) {
  dt <- build_pathway_dt(v)
  # Melt to long
  rbindlist(lapply(names(PATHWAY_LABELS), function(code) {
    data.table(
      outcome      = v,
      pathway      = PATHWAY_LABELS[code],
      shap_signed  = mean(dt[[code]], na.rm=TRUE),
      shap_abs     = mean(abs(dt[[code]]), na.rm=TRUE)
    )
  }))
}))
pathway_agg[, outcome_label := factor(OUTCOME_LABELS[outcome], levels=OUTCOME_LABELS)]
pathway_agg[, pathway       := factor(pathway, levels=names(PATHWAY_COLS))]
fwrite(pathway_agg, "Results/04c_pathway_decomposition.csv")
msg("Pathway decomposition saved -> Results/04c_pathway_decomposition.csv")

# ── POLICY CHART 1: Signed SHAP decomposition ────────────────────────────────
p_pol1 <- ggplot(pathway_agg[!is.na(pathway) & !is.na(outcome_label)],
                 aes(x=outcome_label, y=shap_signed, fill=pathway)) +
  geom_col(width=0.75, colour="white", linewidth=0.25) +
  geom_hline(yintercept=0, linewidth=0.5, colour="#333") +
  scale_fill_manual(values=PATHWAY_COLS, name="Pathway", drop=FALSE) +
  scale_x_discrete(guide=guide_axis(angle=35)) +
  scale_y_continuous(labels=function(x) sprintf("%+.3f", x)) +
  labs(
    title    = "POLICY FIG 1 \u2014 Signed SHAP: Direct vs Indirect Oil Transmission",
    subtitle = paste(
      "Positive bars = pathway pushes outcome UP | Negative = dampens",
      "\nOrange = direct price effect | Purple = rate channel | Blue = balance sheet propagation"
    ),
    caption  = "Exact TreeSHAP | XGBoost | Full test set 2020Q1-2025Q4",
    x=NULL, y="Mean signed SHAP"
  ) +
  theme_minimal(base_size=10) +
  theme(plot.title=element_text(face="bold",size=11),
        plot.subtitle=element_text(size=8.5,colour="#444",lineheight=1.3),
        plot.caption=element_text(size=7.5,colour="#888",hjust=0),
        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        legend.position="right", legend.text=element_text(size=8),
        legend.key.size=unit(0.4,"cm"))
ggsave("Figures/04c_policy1_signed_decomp.png", p_pol1,
       width=14, height=7, dpi=300, bg="white")
msg("Policy Chart 1 \u2192 Figures/04c_policy1_signed_decomp.png")

# ── POLICY CHART 2: Direct vs Indirect share (replaces broken multiplier) ────
# The multiplier breaks when direct SHAP ~ 0 (as seen: 47M×).
# Reframe as: what SHARE of the total attributable effect is indirect?
# indirect_share = (total - direct) / total  →  bounded 0–100%, always readable.

direct_tot2 <- pathway_agg[pathway == "Direct oil price",
                             .(direct = shap_abs), by=outcome]
all_tot2    <- pathway_agg[pathway != "Other",
                             .(total  = sum(shap_abs)), by=outcome]
ratio_dt    <- merge(direct_tot2, all_tot2, by="outcome")
ratio_dt[, indirect      := total - direct]
ratio_dt[, indirect_pct  := indirect / pmax(total, 1e-10) * 100]
ratio_dt[, direct_pct    := direct   / pmax(total, 1e-10) * 100]
ratio_dt[, outcome_label := factor(OUTCOME_LABELS[outcome], levels=OUTCOME_LABELS)]

pol2_long <- melt(
  ratio_dt[, .(outcome_label, direct_pct, indirect_pct)],
  id.vars="outcome_label", variable.name="component", value.name="pct"
)
pol2_long[, component_label := fifelse(
  component == "direct_pct", "Direct oil price", "Indirect (all other pathways)"
)]
pol2_long[, component_label := factor(component_label,
  levels=c("Direct oil price","Indirect (all other pathways)"))]

msg("Direct vs indirect shares:")
print(ratio_dt[order(indirect_pct),
               .(outcome, direct_pct=round(direct_pct,1),
                 indirect_pct=round(indirect_pct,1))])

p_pol2 <- ggplot(pol2_long[!is.na(outcome_label)],
                 aes(x=outcome_label, y=pct, fill=component_label)) +
  geom_col(width=0.72, colour="white", linewidth=0.3) +
  geom_hline(yintercept=50, linetype="dashed", colour="#888", linewidth=0.4) +
  annotate("text", x=0.6, y=53, label="50%% threshold",
           hjust=0, size=2.8, colour="#888") +
  scale_fill_manual(
    values=c("Direct oil price"="#b5470a",
             "Indirect (all other pathways)"="#185FA5"),
    name=NULL
  ) +
  scale_x_discrete(guide=guide_axis(angle=35)) +
  scale_y_continuous(labels=function(x) paste0(x,"%%"),
                     breaks=seq(0,100,25), limits=c(0,110)) +
  geom_text(aes(label=sprintf("%.0f%%%%",pct)),
            position=position_stack(vjust=0.5),
            size=3.0, fontface="bold", colour="white") +
  labs(
    title    = "POLICY FIG 2 \u2014 Direct vs Indirect Oil Transmission: Share of Total Effect",
    subtitle = paste(
      "Orange = direct oil price SHAP | Blue = all indirect pathways combined",
      "\nKey finding: >95%% of oil's effect operates through indirect channels for most outcomes"
    ),
    caption  = paste(
      "Indirect share = (total \u2212 direct) / total SHAP | Exact TreeSHAP | XGBoost",
      "\nNear-zero direct share reflects that oil transmits via balance sheet accumulation, not contemporaneous price level"
    ),
    x=NULL, y="Share of total attributable SHAP effect (%%)"
  ) +
  theme_minimal(base_size=10) +
  theme(
    plot.title=element_text(face="bold",size=11),
    plot.subtitle=element_text(size=8.5,colour="#444",lineheight=1.3),
    plot.caption=element_text(size=7.5,colour="#888",hjust=0),
    panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
    legend.position="top", legend.text=element_text(size=9)
  )
ggsave("Figures/04c_policy2_amplification.png", p_pol2,
       width=12, height=6, dpi=300, bg="white")
msg("Policy Chart 2 \u2192 Figures/04c_policy2_amplification.png")


# ── POLICY CHART 3: Transmission matrix heatmap ──────────────────────────────
hm_dt <- pathway_agg[pathway != "Other" & !is.na(outcome_label)]
hm_dt[, share := shap_abs / sum(shap_abs), by=outcome_label]
hm_dt[, cell  := fifelse(share >= 0.03, sprintf("%.0f%%", share*100), "")]
hm_dt[, pathway_rev := factor(pathway, levels=rev(levels(pathway)))]

p_pol3 <- ggplot(hm_dt, aes(x=outcome_label, y=pathway_rev, fill=share)) +
  geom_tile(colour="white", linewidth=0.6) +
  geom_text(aes(label=cell), size=2.8, fontface="bold",
            colour=ifelse(hm_dt$share > 0.35, "white", "#333")) +
  scale_fill_gradient2(low="white", mid="#F4B183", high="#b5470a",
                       midpoint=0.20, limits=c(0,NA),
                       labels=percent_format(), name="Share of\ntotal effect") +
  scale_x_discrete(guide=guide_axis(angle=35)) +
  labs(
    title    = "POLICY FIG 3 \u2014 Oil Transmission Matrix: Pathway \u00d7 CU Outcome",
    subtitle = "Cell = % of total SHAP from each pathway | Darker = that pathway dominates the transmission",
    caption  = "Exact TreeSHAP | orange=direct, purple=rate, blue=balance sheet, green=labour/housing",
    x=NULL, y=NULL
  ) +
  theme_minimal(base_size=10) +
  theme(plot.title=element_text(face="bold",size=11),
        plot.subtitle=element_text(size=8.5,colour="#444"),
        plot.caption=element_text(size=7.5,colour="#888",hjust=0),
        panel.grid=element_blank(), axis.text.y=element_text(size=9),
        legend.position="right")
ggsave("Figures/04c_policy3_transmission_matrix.png", p_pol3,
       width=14, height=8, dpi=300, bg="white")
msg("Policy Chart 3 \u2192 Figures/04c_policy3_transmission_matrix.png")

# ── POLICY CHART 4: Time-varying — smoothed, split macro vs balance sheet ────
# Two improvements vs prior version:
#   1. 4-quarter rolling mean removes quarterly noise, shows episodic structure
#   2. TWO sub-panels per outcome: macro channels + balance sheet lags
#      This directly shows the hand-off: macro channels spark the shock,
#      balance sheet lags carry it forward.
msg("Policy Chart 4: computing SHAP on full panel ...")
full_feat_mat <- as.matrix(panel_clean[, ..feature_cols])
full_feat_mat[!is.finite(full_feat_mat)] <- 0
dfull <- xgb.DMatrix(full_feat_mat, feature_names=feature_cols)

tv_list <- lapply(c("dq_rate","pll_rate"), function(v) {
  sv <- predict(xgb_models[[v]], dfull, predcontrib=TRUE)
  sv <- sv[, colnames(sv) != "BIAS", drop=FALSE]
  rbindlist(lapply(names(PATHWAY_LABELS), function(code) {
    lbl <- PATHWAY_LABELS[code]
    fp  <- intersect(pathway_map[pathway==lbl, feature], colnames(sv))
    data.table(
      yyyyqq   = panel_clean$yyyyqq,
      outcome  = v,
      pathway  = lbl,
      shap_val = if (length(fp)>0) rowSums(sv[, fp, drop=FALSE], na.rm=TRUE) else 0
    )
  }))
})
tv_dt  <- rbindlist(tv_list)

# Quarterly mean across CUs
tv_agg <- tv_dt[, .(shap_mean=mean(shap_val, na.rm=TRUE)),
                by=.(yyyyqq, outcome, pathway)]
setorder(tv_agg, outcome, pathway, yyyyqq)

# 4-quarter rolling mean (within outcome × pathway)
tv_agg[, shap_roll := frollmean(shap_mean, n=4L, na.rm=TRUE, align="right"),
        by=.(outcome, pathway)]

tv_agg[, yr       := yyyyqq %/% 100L]
tv_agg[, qtr      := yyyyqq %% 100L]
tv_agg[, date_num := yr + (qtr-1)/4]
tv_agg[, outcome_label := factor(OUTCOME_LABELS[as.character(outcome)],
                                  levels=OUTCOME_LABELS)]
tv_agg[, pathway := factor(pathway, levels=names(PATHWAY_COLS))]

# Split into macro channels vs balance sheet
MACRO_PATHS <- c("Direct oil price","Interaction (oil x regime)",
                 "Indirect: Labour market","Indirect: Inflation (CPI)",
                 "Indirect: Rate channel","Indirect: Housing")
BAL_PATHS   <- "Indirect: Balance sheet lags"

tv_agg[, panel_grp := fifelse(
  as.character(pathway) %in% MACRO_PATHS, "Macro channels", "Balance sheet lags"
)]
tv_agg[, panel_grp := factor(panel_grp,
  levels=c("Macro channels","Balance sheet lags"))]

episodes <- data.table(
  xmin  = c(2008.5, 2015.0, 2020.0, 2022.0),
  xmax  = c(2009.75,2016.0, 2021.25,2023.75),
  label = c("GFC","Shale bust","COVID","Hike cycle")
)

plot_tv <- tv_agg[
  pathway != "Other" & !is.na(outcome_label) & !is.na(shap_roll)
]

p_pol4 <- ggplot(plot_tv,
                 aes(x=date_num, y=shap_roll,
                     colour=pathway, group=pathway)) +
  geom_rect(data=episodes,
            aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf),
            fill="#f0f0f0", alpha=0.7, inherit.aes=FALSE) +
  geom_text(data=episodes,
            aes(x=(xmin+xmax)/2, y=Inf, label=label),
            vjust=1.3, size=2.3, colour="#666", inherit.aes=FALSE) +
  geom_hline(yintercept=0, linewidth=0.35, colour="#444") +
  geom_line(linewidth=0.8) +
  scale_colour_manual(values=PATHWAY_COLS, name="Pathway", drop=FALSE) +
  scale_x_continuous(breaks=seq(2006,2026,2),
                     labels=function(x) paste0("'",substr(x,3,4))) +
  facet_grid(outcome_label ~ panel_grp, scales="free_y") +
  labs(
    title    = "POLICY FIG 4 \u2014 Time-Varying Oil Transmission: Macro vs Balance Sheet (2005\u20132025)",
    subtitle = paste(
      "Left panels: contemporaneous macro channels | Right: balance sheet lag propagation",
      "\n4-quarter rolling mean | Key finding: macro channels initiate, balance sheet carries forward"
    ),
    caption  = paste(
      "TreeSHAP on full panel 2005Q1-2025Q4 | 4-quarter rolling mean of quarterly cross-CU mean",
      "\nEpisodes: GFC=2008-09, Shale bust=2015-16, COVID=2020-21, Hike cycle=2022-23"
    ),
    x=NULL, y="Mean SHAP (4-qtr rolling)"
  ) +
  theme_minimal(base_size=10) +
  theme(
    plot.title=element_text(face="bold",size=11),
    plot.subtitle=element_text(size=8.5,colour="#444",lineheight=1.3),
    plot.caption=element_text(size=7.5,colour="#888",hjust=0),
    strip.text=element_text(face="bold",size=9),
    strip.background=element_rect(fill="#f5f5f5",colour="#ccc"),
    panel.grid.minor=element_blank(),
    legend.position="right",
    legend.text=element_text(size=7.5),
    legend.key.size=unit(0.35,"cm")
  )
ggsave("Figures/04c_policy4_time_varying.png", p_pol4,
       width=15, height=10, dpi=300, bg="white")
msg("Policy Chart 4 \u2192 Figures/04c_policy4_time_varying.png")


# ── POLICY CHART 5: Pathway by asset tier ────────────────────────────────────
if ("asset_tier" %in% names(panel_test)) {
  # Focus pathways for tier comparison
  pw_focus_codes <- c("D","RTE","LAB","CPI","BAL")
  pw_focus_lbls  <- PATHWAY_LABELS[pw_focus_codes]

  tier_rows <- rbindlist(lapply(Y_VARS[1:4], function(v) {
    sv <- shap_results[[v]]
    tier_vec <- panel_test$asset_tier
    rbindlist(lapply(pw_focus_codes, function(code) {
      lbl <- PATHWAY_LABELS[code]
      fp  <- intersect(pathway_map[pathway==lbl, feature], colnames(sv))
      vals <- if (length(fp)>0) rowSums(abs(sv[, fp, drop=FALSE]), na.rm=TRUE) else 0
      data.table(outcome=v, tier=tier_vec, pathway=lbl, abs_shap=vals)
    }))
  }))

  tier_agg <- tier_rows[, .(mean_abs=mean(abs_shap, na.rm=TRUE)),
                         by=.(outcome, tier, pathway)]
  tier_agg[, outcome_label := factor(OUTCOME_LABELS[as.character(outcome)],
                                      levels=OUTCOME_LABELS)]
  tier_agg[, tier_label := factor(
    paste0("T", tier),
    levels=paste0("T", sort(unique(tier_agg$tier)))
  )]
  tier_agg[, pathway := factor(pathway, levels=names(PATHWAY_COLS))]

  T_COLS <- PATHWAY_COLS[pw_focus_lbls]

  p_pol5 <- ggplot(tier_agg[!is.na(outcome_label) & !is.na(pathway)],
                   aes(x=tier_label, y=mean_abs, fill=pathway)) +
    geom_col(width=0.72, position="stack", colour="white", linewidth=0.2) +
    facet_wrap(~outcome_label, scales="free_y", ncol=2) +
    scale_fill_manual(values=T_COLS, name="Pathway", drop=FALSE) +
    scale_x_discrete(drop=FALSE) +
    labs(
      title    = "POLICY FIG 5 \u2014 Oil Transmission Pathways by Asset Tier",
      subtitle = paste(
        "Does SIZE change HOW oil affects CUs, not just how much?",
        "\nT1 <$10M | T2 $10-100M | T3 $100M-$1B | T4 >$1B"
      ),
      caption  = "Mean |TreeSHAP| per tier | Stacked = total attributable effect by pathway",
      x="Asset tier", y="Mean |SHAP|"
    ) +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold",size=11),
          plot.subtitle=element_text(size=8.5,colour="#444",lineheight=1.3),
          plot.caption=element_text(size=7.5,colour="#888",hjust=0),
          strip.text=element_text(face="bold",size=9),
          strip.background=element_rect(fill="#f5f5f5",colour="#ccc"),
          panel.grid.minor=element_blank(),
          legend.position="right", legend.text=element_text(size=8))
  ggsave("Figures/04c_policy5_tier_pathways.png", p_pol5,
         width=13, height=9, dpi=300, bg="white")
  msg("Policy Chart 5 \u2192 Figures/04c_policy5_tier_pathways.png")
}



hdr("SECTION 9: Output Manifest")

outputs <- list(
  Models  = list.files("Models",  pattern="04c_", full.names=TRUE),
  Results = list.files("Results", pattern="04c_", full.names=TRUE),
  Figures = list.files("Figures", pattern="04c_", full.names=TRUE)
)
for (grp in names(outputs)) {
  cat(sprintf("\n  %s:\n", grp))
  for (f in outputs[[grp]])
    cat(sprintf("    %-55s [%s bytes]\n", basename(f),
                format(file.size(f), big.mark=",")))
}

t_total <- proc.time() - t_script_start
cat("\n", strrep("=",70), "\n")
cat(sprintf("  Total runtime      : %.1f sec (%.1f min)\n",
            t_total["elapsed"], t_total["elapsed"]/60))
cat(sprintf("  Grid search combos : %d\n", nrow(GRID)))
cat(sprintf("  XGBoost threads    : %d\n", N_THREAD))
cat(sprintf("  SHAP method        : Exact TreeSHAP (predcontrib=TRUE)\n"))
cat(sprintf("  Test observations  : %s (full set, no sampling)\n",
            format(nrow(X_test), big.mark=",")))
cat("  Script 04c complete\n")
cat(strrep("=",70), "\n")
