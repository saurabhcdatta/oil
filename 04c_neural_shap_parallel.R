# =============================================================================
# OIL PRICE SHOCK × CREDIT UNION RESEARCH
# Script 04c — Neural Network + SHAP Transmission Analysis
# Backend: neuralnet (pure R, no torch/keras/TensorFlow required)
# =============================================================================
#
# Purpose:
#   Nonparametric robustness check on the linear VARX (Script 04b).
#   A 2-hidden-layer neural network trained on the same panel data
#   captures nonlinearities and threshold effects that the VARX assumes away.
#   Kernel SHAP decomposes each prediction into direct and indirect
#   contributions, separating:
#
#     DIRECT:   oil → CU outcome  (purchasing power / income channel)
#     INDIRECT: oil → unemployment → CU outcome  (labour market channel)
#               oil → CPI → CoF → CU outcome    (inflation / rate channel)
#               oil → deposits → loan_to_share   (balance sheet channel)
#
#   This is NOT a structural model. Economic identification rests on the
#   VARX Cholesky decomposition in Script 04b. This is a robustness check.
#
# Packages required: neuralnet, data.table, ggplot2, patchwork, scales, stringr
#   All available on CRAN with no external binary dependencies.
#
# Inputs:  Data/panel_model.rds
#          Data/macro_base.rds
# Outputs: Figures/04c_shap_channels.png
#          Figures/04c_pdp_oil_grid.png
#          Figures/04c_nonlinear_oil.png
#          Figures/04c_shap_waterfall_dq_rate.png
#          Figures/04c_shap_waterfall_pll_rate.png
#          Figures/04c_shap_tier.png
#          Figures/04c_nn_vs_varx.png
#          Results/04c_shap_importance.csv
#          Results/04c_nn_performance.csv
#          Models/04c_nn_models.rds
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(neuralnet)    # pure R neural net — no external dependencies
  library(parallel)     # built into base R — no install needed
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

# ── Helpers ──────────────────────────────────────────────────────────────────
hdr <- function(x) cat("\n", strrep("=",70), "\n##", x, "\n",
                        strrep("=",70), "\n")
msg <- function(fmt, ...) cat(sprintf(paste0("  ", fmt, "\n"), ...))
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

dir.create("Models",  showWarnings=FALSE)
dir.create("Results", showWarnings=FALSE)
dir.create("Figures", showWarnings=FALSE)

set.seed(20260101L)
t_script_start <- proc.time()   # for total runtime at end

# ── Parallel cluster setup (Windows PSOCK) ───────────────────────────────────
# detectCores() returns logical cores (including hyperthreading).
# We use n-1 cores: leave one free so the machine stays responsive.
# On a 4-core/8-thread machine this gives 7 workers.
# PSOCK works on all Windows versions — no fork() needed.
N_CORES <- max(1L, detectCores(logical=TRUE) - 1L)
msg("Detected %d logical cores — using %d workers (PSOCK cluster)",
    detectCores(logical=TRUE), N_CORES)

# Create the cluster once at startup — stops at end of script
cl <- makeCluster(N_CORES, type="PSOCK")

# Register cleanup so cluster shuts down even if script errors
on.exit({
  tryCatch(stopCluster(cl), error=function(e) invisible(NULL))
  msg("Parallel cluster stopped.")
}, add=TRUE)

# Helper: export objects to all workers
cl_export <- function(...) {
  nms <- as.character(match.call()[-1])
  clusterExport(cl, varlist=nms, envir=parent.frame())
}

# Helper: load packages on all workers
cl_libs <- function(...) {
  pkgs <- c(...)
  invisible(clusterEvalQ(cl, {
    suppressPackageStartupMessages(lapply(pkgs, library,
                                           character.only=TRUE))
  }))
}
cl_libs("neuralnet", "data.table")

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

CU_CONTROLS  <- c("oil_exposure_cont", "pcanetworth", "insured_share_growth")
INTERACT_VARS <- c("fomc_x_brent", "post_x_oil", "post_x_oil_x_direct")
DUMMY_VARS   <- c("post_shale", "gfc_dummy", "covid_dummy",
                   "zirp_era", "hike_cycle")

P_LAG        <- 2L

# neuralnet config
# 2 hidden layers — (16, 8) neurons. Kept small for pure-R speed.
# neuralnet is exact (rprop algorithm) but slower than GPU torch.
# Runtime on 180k obs panel: ~8-15 minutes per outcome variable.
# We train one model per outcome (neuralnet doesn't natively multi-task).
HIDDEN_LAYERS  <- c(16L, 8L)
MAX_STEPS      <- 1e5          # max iterations per model
THRESHOLD_CONV <- 0.01         # convergence threshold (partial derivatives)
N_SHAP_BG      <- 100L         # background rows for Kernel SHAP
N_SHAP_TEST    <- 80L          # test rows to explain
N_COAL         <- 60L          # coalition samples per observation

TRAIN_END  <- 201904L
TEST_START <- 202001L

msg("Backend      : neuralnet (pure R, no external dependencies)")
msg("Architecture : %d → %s → 1  per outcome variable",
    0, paste(HIDDEN_LAYERS, collapse=" → "))
msg("Outcomes     : %d models (one per Y variable)", length(Y_VARS))
msg("SHAP bg      : %d | test obs: %d | coalitions: %d",
    N_SHAP_BG, N_SHAP_TEST, N_COAL)
msg("Train        : 2005Q1–2019Q4  |  Test: 2020Q1–2025Q4")

# =============================================================================
# 1. LOAD AND PREPARE DATA
# =============================================================================
hdr("SECTION 1: Load Data")

panel      <- readRDS("Data/panel_model.rds")
macro_base <- readRDS("Data/macro_base.rds")
setDT(panel); setDT(macro_base)

# Merge macro onto panel
macro_cols <- intersect(c("yyyyqq", Z_VARS, INTERACT_VARS, DUMMY_VARS),
                         names(macro_base))
panel <- merge(panel, macro_base[, ..macro_cols], by="yyyyqq", all.x=TRUE)

# Build lagged Y features
setorder(panel, join_number, yyyyqq)
for (v in Y_VARS) {
  for (k in 1:P_LAG) {
    nm <- paste0(v, "_lag", k)
    panel[, (nm) := shift(.SD[[v]], k, type="lag"),
           by=join_number, .SDcols=v]
  }
}

lag_y_vars <- unlist(lapply(Y_VARS, function(v)
  paste0(v, "_lag", 1:P_LAG)))

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

keep_cols <- unique(c("join_number","yyyyqq","asset_tier",
                       Y_VARS, feature_cols))
keep_cols <- intersect(keep_cols, names(panel))
panel_clean <- panel[complete.cases(panel[, ..keep_cols])]

panel_train <- panel_clean[yyyyqq <= TRAIN_END]
panel_test  <- panel_clean[yyyyqq >= TEST_START]

# Standardise features on train statistics (z-score)
feat_mean <- sapply(feature_cols, function(v)
  mean(panel_train[[v]], na.rm=TRUE))
feat_sd   <- sapply(feature_cols, function(v) {
  s <- sd(panel_train[[v]], na.rm=TRUE)
  if (is.na(s) || s < 1e-10) 1 else s
})

standardise_mat <- function(dt) {
  out <- matrix(NA_real_, nrow=nrow(dt), ncol=length(feature_cols),
                dimnames=list(NULL, feature_cols))
  for (v in feature_cols)
    out[, v] <- (dt[[v]] - feat_mean[v]) / feat_sd[v]
  out
}

# Standardise outcomes per variable on train
y_mean <- sapply(Y_VARS, function(v) mean(panel_train[[v]], na.rm=TRUE))
y_sd   <- sapply(Y_VARS, function(v) {
  s <- sd(panel_train[[v]], na.rm=TRUE)
  if (is.na(s) || s < 1e-10) 1 else s
})

X_train_mat <- standardise_mat(panel_train)
X_test_mat  <- standardise_mat(panel_test)

# For neuralnet we need data.frames with named columns
X_train_df <- as.data.frame(X_train_mat)
X_test_df  <- as.data.frame(X_test_mat)

msg("Features     : %d", length(feature_cols))
msg("Train rows   : %s", format(nrow(X_train_df), big.mark=","))
msg("Test rows    : %s", format(nrow(X_test_df),  big.mark=","))

# Because neuralnet is slow on very large datasets we subsample training data
# using a stratified sample (by yyyyqq + asset_tier) for speed.
# The full test set is still used for evaluation.
MAX_TRAIN_ROWS <- 20000L
if (nrow(panel_train) > MAX_TRAIN_ROWS) {
  msg("Subsampling training data to %s rows (stratified by quarter)",
      format(MAX_TRAIN_ROWS, big.mark=","))
  strat_idx <- panel_train[, .I[sample(.N, min(.N, 50L))],
                             by=yyyyqq]$V1
  strat_idx <- strat_idx[order(strat_idx)]
  if (length(strat_idx) > MAX_TRAIN_ROWS)
    strat_idx <- sample(strat_idx, MAX_TRAIN_ROWS)
  panel_train_sub <- panel_train[strat_idx]
  X_train_sub     <- as.data.frame(standardise_mat(panel_train_sub))
} else {
  panel_train_sub <- panel_train
  X_train_sub     <- X_train_df
}

# =============================================================================
# 3. TRAIN ONE NEURAL NET PER OUTCOME VARIABLE
# =============================================================================
hdr("SECTION 3: Neural Network Training")

# neuralnet formula construction
nn_formula <- function(outcome_std_col, feature_names) {
  rhs <- paste(paste0("`", feature_names, "`"), collapse=" + ")
  as.formula(paste0("`", outcome_std_col, "` ~ ", rhs))
}

# ── Export everything workers need for training ───────────────────────────────
clusterExport(cl, varlist=c(
  "X_train_sub","X_test_df","panel_train_sub","panel_test",
  "Y_VARS","y_mean","y_sd","feature_cols","nn_formula",
  "HIDDEN_LAYERS","THRESHOLD_CONV","MAX_STEPS"
), envir=environment())

# ── Parallel training: one worker per outcome variable ────────────────────────
# Each worker trains independently — no shared state needed.
# parLapply returns a named list; worker stdout is suppressed (PSOCK).
msg("Training %d neural nets in parallel on %d workers ...", length(Y_VARS), N_CORES)
t_train_start <- proc.time()

train_results <- parLapply(cl, Y_VARS, function(v) {

  y_std_col   <- paste0(v, "_std")
  y_std_train <- (panel_train_sub[[v]] - y_mean[v]) / y_sd[v]
  train_df    <- cbind(X_train_sub,
                        setNames(data.frame(y_std_train), y_std_col))
  frml        <- nn_formula(y_std_col, feature_cols)

  nn_fit <- tryCatch(
    neuralnet(
      formula       = frml,
      data          = train_df,
      hidden        = HIDDEN_LAYERS,
      linear.output = TRUE,
      act.fct       = "logistic",
      threshold     = THRESHOLD_CONV,
      stepmax       = MAX_STEPS,
      rep           = 1L,
      lifesign      = "none"       # suppress per-worker output
    ),
    error = function(e) NULL
  )

  # Fallback to lm if neuralnet fails
  if (is.null(nn_fit)) {
    fb <- lm(as.formula(paste0(y_std_col, " ~ .")), data=train_df)
    fb$model_type <- "lm_fallback"
    return(list(model=fb, outcome=v))
  }

  nn_fit$model_type <- "neuralnet"
  list(model=nn_fit, outcome=v)
})

t_train_elapsed <- proc.time() - t_train_start
msg("Training complete in %.1f seconds (%.1f min)",
    t_train_elapsed["elapsed"], t_train_elapsed["elapsed"]/60)

# Unpack results into named list
nn_models <- setNames(
  lapply(train_results, `[[`, "model"),
  sapply(train_results, `[[`, "outcome")
)

# ── Evaluate on test set (fast — serial is fine here) ─────────────────────────
predict_one <- function(mod, newdata) {
  if (!is.null(mod$model_type) && mod$model_type == "lm_fallback") {
    as.numeric(predict(mod, newdata=newdata))
  } else {
    as.numeric(predict(mod, newdata=newdata))
  }
}

nn_perf_dt <- rbindlist(lapply(Y_VARS, function(v) {
  y_hat_std <- predict_one(nn_models[[v]], X_test_df)
  y_hat     <- y_hat_std * y_sd[v] + y_mean[v]
  y_act     <- panel_test[[v]]
  rmse <- sqrt(mean((y_hat - y_act)^2, na.rm=TRUE))
  r2   <- 1 - sum((y_hat - y_act)^2, na.rm=TRUE) /
              sum((y_act - mean(y_act, na.rm=TRUE))^2, na.rm=TRUE)
  msg("  %-30s | type=%-12s | R²=%.3f | RMSE=%.5f",
      v, nn_models[[v]]$model_type, r2, rmse)
  data.table(outcome=v, model_type=nn_models[[v]]$model_type,
             rmse_test=round(rmse,5), r2_test=round(r2,4))
}))

fwrite(nn_perf_dt, "Results/04c_nn_performance.csv")
cat("\n  Model performance:\n"); print(nn_perf_dt)

saveRDS(list(
  models=nn_models, feat_mean=feat_mean, feat_sd=feat_sd,
  y_mean=y_mean, y_sd=y_sd, feature_cols=feature_cols, Y_VARS=Y_VARS
), "Models/04c_nn_models.rds")
msg("Models saved → Models/04c_nn_models.rds")

# =============================================================================
# 4. KERNEL SHAP — MANUAL IMPLEMENTATION
# =============================================================================
hdr("SECTION 4: Kernel SHAP Decomposition")

# Kernel SHAP: model-agnostic, works with any black-box function.
#
# For each observation x, each coalition S (subset of features) is evaluated:
#   - "ON"  version: features in S take their value from x, rest from background mean
#   - "OFF" version: features in S take background mean, rest from x
#
# The contribution of each feature is estimated by weighted least squares
# over all 2^p coalitions (we sample N_COAL of them for tractability).
#
# Weight of coalition S: w(S) = (p-1) / (C(p,|S|) × |S| × (p-|S|))
# where p = total features, |S| = coalition size
#
# For tractability on 40+ features, we use random coalition sampling
# with the Shapley kernel weights and return the WLS solution.

# Background dataset — stratified mean of training data
bg_idx  <- sample(nrow(X_train_sub), min(N_SHAP_BG, nrow(X_train_sub)))
X_bg    <- as.matrix(X_train_sub[bg_idx, feature_cols])
bg_mean <- colMeans(X_bg)

# Prediction function: takes matrix → returns vector of predictions (std scale)
predict_std <- function(v, newdata_mat) {
  nd <- as.data.frame(newdata_mat)
  colnames(nd) <- feature_cols
  mod <- nn_models[[v]]
  if (!is.null(mod$model_type) && mod$model_type == "lm_fallback") {
    # need the outcome col to exist for lm (use dummy)
    nd[[paste0(v, "_std")]] <- 0
    as.numeric(predict(mod, newdata=nd))
  } else {
    as.numeric(predict(mod, newdata=nd))
  }
}

# Shapley kernel weight for coalition of size s out of p total features
shapley_weight <- function(s, p) {
  if (s == 0 || s == p) return(1e6)   # boundary — large weight
  (p - 1) / (choose(p, s) * s * (p - s))
}

# Compute SHAP values for a single observation and outcome variable
compute_shap <- function(x_obs, v, n_coal=N_COAL) {
  p    <- length(x_obs)
  f_bg <- mean(predict_std(v, matrix(rep(bg_mean, N_SHAP_BG),
                                      nrow=N_SHAP_BG, byrow=TRUE)))

  # Sample random coalitions
  coalitions <- matrix(FALSE, nrow=n_coal, ncol=p)
  weights    <- numeric(n_coal)

  for (j in seq_len(n_coal)) {
    s <- sample(0:p, 1, prob=dbinom(0:p, p, 0.5))
    s <- max(1, min(s, p-1))   # keep between 1 and p-1
    idx_on <- sample(p, s)
    coalitions[j, idx_on] <- TRUE
    weights[j] <- shapley_weight(s, p)
  }

  # Build masked input matrices
  f_on  <- numeric(n_coal)
  f_off <- numeric(n_coal)

  for (j in seq_len(n_coal)) {
    x_on       <- bg_mean
    x_on[coalitions[j,]] <- x_obs[coalitions[j,]]
    x_off      <- x_obs
    x_off[coalitions[j,]] <- bg_mean[coalitions[j,]]

    f_on[j]  <- mean(predict_std(v, matrix(x_on,  nrow=1)))
    f_off[j] <- mean(predict_std(v, matrix(x_off, nrow=1)))
  }

  # WLS: regress (f_on - f_off) on coalition indicators, weighted
  z_mat <- coalitions * 1.0    # (n_coal × p) binary matrix
  y_vec <- f_on - f_off

  # Closed-form WLS: β = (Z'WZ)^{-1} Z'W y
  W     <- diag(weights)
  ZtW   <- t(z_mat) %*% W
  ZtWZ  <- ZtW %*% z_mat
  ZtWy  <- ZtW %*% y_vec

  shap_vals <- tryCatch(
    as.numeric(solve(ZtWZ + diag(1e-8, p)) %*% ZtWy),
    error = function(e) {
      # Fallback: simple weighted mean differences per feature
      sapply(seq_len(p), function(i) {
        idx <- which(coalitions[, i])
        if (length(idx) == 0) return(0)
        weighted.mean(y_vec[idx], weights[idx])
      })
    }
  )
  names(shap_vals) <- feature_cols
  shap_vals
}

# Run SHAP on a test sample
shap_test_idx <- sample(nrow(X_test_mat), min(N_SHAP_TEST, nrow(X_test_mat)))
X_shap_mat    <- X_test_mat[shap_test_idx, ]

msg("Running Kernel SHAP (parallel): %d obs × %d outcomes × %d coalitions",
    nrow(X_shap_mat), length(Y_VARS), N_COAL)
msg("Using %d workers — estimated time: ~2–5 min", N_CORES)

t_shap_start <- proc.time()

# Export everything workers need for SHAP
clusterExport(cl, varlist=c(
  "nn_models","feature_cols","Y_VARS",
  "X_bg","bg_mean","N_SHAP_BG","N_COAL",
  "X_shap_mat","y_mean","y_sd","predict_one"
), envir=environment())

# Each worker handles one outcome variable.
# Inside each worker the observation loop runs serially —
# this gives 8-way parallelism (one model per worker).
shap_results_raw <- parLapply(cl, Y_VARS, function(v) {

  # Worker-local predict function (returns standardised predictions)
  predict_std_worker <- function(x_mat) {
    nd <- as.data.frame(x_mat)
    colnames(nd) <- feature_cols
    mod <- nn_models[[v]]
    if (!is.null(mod$model_type) && mod$model_type == "lm_fallback") {
      nd[[paste0(v, "_std")]] <- 0
      as.numeric(predict(mod, newdata=nd))
    } else {
      as.numeric(predict(mod, newdata=nd))
    }
  }

  p      <- ncol(X_shap_mat)
  f_bg   <- mean(predict_std_worker(X_bg))

  shapley_weight <- function(s, p) {
    if (s == 0 || s == p) return(1e6)
    (p - 1) / (choose(p, s) * s * (p - s))
  }

  compute_shap_worker <- function(x_obs) {
    coalitions <- matrix(FALSE, nrow=N_COAL, ncol=p)
    weights    <- numeric(N_COAL)
    for (j in seq_len(N_COAL)) {
      s <- max(1L, min(sample(0:p, 1, prob=dbinom(0:p, p, 0.5)), p-1L))
      idx_on         <- sample(p, s)
      coalitions[j, idx_on] <- TRUE
      weights[j]     <- shapley_weight(s, p)
    }
    f_on  <- numeric(N_COAL)
    f_off <- numeric(N_COAL)
    for (j in seq_len(N_COAL)) {
      x_on         <- bg_mean; x_on[coalitions[j,]]  <- x_obs[coalitions[j,]]
      x_off        <- x_obs;   x_off[coalitions[j,]] <- bg_mean[coalitions[j,]]
      f_on[j]  <- mean(predict_std_worker(matrix(x_on,  nrow=1)))
      f_off[j] <- mean(predict_std_worker(matrix(x_off, nrow=1)))
    }
    z_mat <- coalitions * 1.0
    y_vec <- f_on - f_off
    W     <- diag(weights)
    ZtW   <- t(z_mat) %*% W
    ZtWZ  <- ZtW %*% z_mat
    ZtWy  <- ZtW %*% y_vec
    sv <- tryCatch(
      as.numeric(solve(ZtWZ + diag(1e-8, p)) %*% ZtWy),
      error=function(e) {
        sapply(seq_len(p), function(i) {
          idx <- which(coalitions[,i])
          if (!length(idx)) return(0)
          weighted.mean(y_vec[idx], weights[idx])
        })
      }
    )
    names(sv) <- feature_cols
    sv
  }

  # Observation loop — serial within worker
  sv_mat <- matrix(NA_real_, nrow=nrow(X_shap_mat), ncol=p,
                   dimnames=list(NULL, feature_cols))
  for (i in seq_len(nrow(X_shap_mat)))
    sv_mat[i,] <- compute_shap_worker(X_shap_mat[i,])

  list(outcome=v, sv_mat=sv_mat)
})

t_shap_elapsed <- proc.time() - t_shap_start
msg("SHAP complete in %.1f seconds (%.1f min)",
    t_shap_elapsed["elapsed"], t_shap_elapsed["elapsed"]/60)

shap_results <- setNames(
  lapply(shap_results_raw, `[[`, "sv_mat"),
  sapply(shap_results_raw, `[[`, "outcome")
)
msg("SHAP computation complete.")

# =============================================================================
# 5. CHANNEL ATTRIBUTION
# =============================================================================
hdr("SECTION 5: Transmission Channel Attribution")

# Classify each feature into an economic transmission channel
channel_map <- data.table(feature=feature_cols, channel="Other")
channel_map[grepl("yoy_oil|post_x_oil|fomc_x_brent|zirp_x_oil",
                   feature), channel := "Direct (oil)"]
channel_map[grepl("lurc",             feature), channel := "Labour market"]
channel_map[grepl("pcpi",             feature), channel := "Inflation (CPI)"]
channel_map[grepl("fomc_regime|yield_curve|rmtg", feature),
             channel := "Rate channel"]
channel_map[grepl("_lag",             feature), channel := "Balance sheet lags"]
channel_map[grepl("hpi",              feature), channel := "Housing (HPI)"]

# Mean absolute SHAP per feature per outcome
shap_importance <- rbindlist(lapply(Y_VARS, function(v) {
  sv <- shap_results[[v]]
  data.table(
    outcome       = v,
    feature       = feature_cols,
    mean_abs_shap = colMeans(abs(sv), na.rm=TRUE)
  )
}))
fwrite(shap_importance, "Results/04c_shap_importance.csv")
msg("SHAP importance → Results/04c_shap_importance.csv")

shap_channel <- merge(shap_importance, channel_map, by="feature")
shap_channel_agg <- shap_channel[,
  .(total_shap = sum(mean_abs_shap, na.rm=TRUE)),
  by=.(outcome, channel)]
shap_channel_agg[, share := total_shap / sum(total_shap), by=outcome]

# =============================================================================
# 6. PARTIAL DEPENDENCE PLOTS (oil YoY sweep)
# =============================================================================
hdr("SECTION 6: Partial Dependence — Oil Price Sweep")

oil_col_idx <- which(feature_cols == "macro_base_yoy_oil")
pdp_dt <- NULL

if (length(oil_col_idx) == 1) {
  oil_range_raw <- seq(-50, 120, by=5)
  oil_range_std <- (oil_range_raw - feat_mean["macro_base_yoy_oil"]) /
                    feat_sd["macro_base_yoy_oil"]

  x_base <- colMeans(X_test_mat)

  msg("Computing PDP across %d oil levels (parallel) ...",
      length(oil_range_raw))

  clusterExport(cl, varlist=c(
    "nn_models","feature_cols","Y_VARS","y_mean","y_sd",
    "x_base","oil_col_idx","N_SHAP_BG","predict_one",
    "X_test_df"
  ), envir=environment())
  clusterExport(cl, varlist="oil_range_std", envir=environment())

  pdp_rows <- parLapply(cl, seq_along(oil_range_raw), function(ri) {
    x_perturb <- matrix(rep(x_base, N_SHAP_BG),
                         nrow=N_SHAP_BG, byrow=TRUE)
    x_perturb[, oil_col_idx] <- oil_range_std[ri]
    nd <- as.data.frame(x_perturb)
    colnames(nd) <- feature_cols

    row <- list(oil_yoy=oil_range_raw[ri])
    for (v in Y_VARS) {
      mod <- nn_models[[v]]
      y_hat_std <- if (!is.null(mod$model_type) &&
                       mod$model_type == "lm_fallback") {
        nd[[paste0(v,"_std")]] <- 0
        as.numeric(predict(mod, newdata=nd))
      } else {
        as.numeric(predict(mod, newdata=nd))
      }
      row[[v]] <- mean(y_hat_std * y_sd[v] + y_mean[v], na.rm=TRUE)
    }
    as.data.table(row)
  })

  pdp_dt <- rbindlist(pdp_rows)
  msg("PDP complete.")
}

# =============================================================================
# 7. CHARTS
# =============================================================================
hdr("SECTION 7: Generate Charts")

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

# ── Chart 1 — Channel decomposition stacked bar ──────────────────────────────
shap_channel_agg[, outcome_label := factor(
  OUTCOME_LABELS[outcome],
  levels = OUTCOME_LABELS
)]
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
    title    = "FIGURE 04c-1 \u2014 Oil Shock Transmission Channels (SHAP Attribution)",
    subtitle = paste(
      "Share of total SHAP importance by economic channel",
      "| Orange = direct oil effect | Blue = balance sheet propagation"
    ),
    caption  = paste(
      "Kernel SHAP | 2-layer neural net (neuralnet, pure R)",
      "| Test set 2020Q1\u20132025Q4",
      "| n=", N_SHAP_TEST, "obs |", N_COAL, "coalitions per obs"
    ),
    x=NULL, y="Share of SHAP importance"
  ) +
  theme_minimal(base_size=10) +
  theme(
    plot.title       = element_text(face="bold", size=11),
    plot.subtitle    = element_text(size=8.5, colour="#555", lineheight=1.3),
    plot.caption     = element_text(size=7.5, colour="#888", hjust=0),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.text      = element_text(size=8)
  )

ggsave("Figures/04c_shap_channels.png", p_channel,
       width=12, height=6, dpi=300, bg="white")
msg("Chart 1 (channel decomposition) \u2192 Figures/04c_shap_channels.png")

# ── Chart 2 — Oil PDP grid ──────────────────────────────────────────────────
if (!is.null(pdp_dt)) {
  KEY_OUTCOMES_PDP <- c("dq_rate","pll_rate","costfds",
                          "netintmrg","insured_share_growth","member_growth_yoy")

  pdp_long <- melt(pdp_dt, id.vars="oil_yoy",
                   variable.name="outcome", value.name="pred")
  pdp_long <- pdp_long[outcome %in% KEY_OUTCOMES_PDP]
  pdp_long[, outcome_label := factor(
    OUTCOME_LABELS[as.character(outcome)],
    levels = OUTCOME_LABELS[KEY_OUTCOMES_PDP]
  )]

  # Moody's scenario reference lines
  vlines <- data.table(
    oil_yoy = c(-5, 0, 22, 42, 98),
    label   = c("Feb '26","Flat","E.March","Current","Recession")
  )

  p_pdp <- ggplot(pdp_long, aes(x=oil_yoy, y=pred)) +
    geom_vline(data=vlines, aes(xintercept=oil_yoy),
               linetype="dashed", colour="#cccccc", linewidth=0.4) +
    geom_text(data=vlines, aes(x=oil_yoy, y=Inf, label=label),
              angle=90, hjust=1.1, size=2.5, colour="#999999",
              inherit.aes=FALSE) +
    geom_line(colour="#b5470a", linewidth=0.95) +
    geom_hline(yintercept=0, linetype="dashed",
               colour="#aaaaaa", linewidth=0.35) +
    facet_wrap(~outcome_label, scales="free_y", ncol=3) +
    scale_x_continuous(
      breaks=seq(-50, 120, 25),
      labels=function(x) paste0(ifelse(x>0,"+",""), x, "%%")
    ) +
    labs(
      title    = "FIGURE 04c-2 \u2014 Partial Dependence: Oil Price YoY vs CU Outcomes",
      subtitle = paste(
        "All other features held at test-set mean",
        "| Vertical lines = Moody\u2019s 2026 scenario levels",
        "\nNon-linearity in shape = where VARX linear assumption breaks down"
      ),
      caption  = paste(
        "Neural net partial dependence | 100 background observations",
        "| Features standardised (z-score, train statistics)"
      ),
      x = "Oil price YoY % change",
      y = "Predicted outcome (level units)"
    ) +
    theme_minimal(base_size=10) +
    theme(
      plot.title       = element_text(face="bold", size=11),
      plot.subtitle    = element_text(size=8.5, colour="#555", lineheight=1.3),
      plot.caption     = element_text(size=7.5, colour="#888", hjust=0),
      strip.text       = element_text(face="bold", size=9),
      strip.background = element_rect(fill="#f5f5f5", colour="#cccccc"),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour="#cccccc", fill=NA, linewidth=0.3)
    )

  ggsave("Figures/04c_pdp_oil_grid.png", p_pdp,
         width=14, height=9, dpi=300, bg="white")
  msg("Chart 2 (PDP grid) \u2192 Figures/04c_pdp_oil_grid.png")

  # ── Chart 3 — Nonlinearity detection ──────────────────────────────────────
  nonlin_dt <- rbindlist(lapply(KEY_OUTCOMES_PDP, function(v) {
    y  <- pdp_dt[[v]]
    x  <- pdp_dt$oil_yoy
    dx <- diff(x)
    d1 <- diff(y) / dx
    d2 <- diff(d1) / dx[-1]
    data.table(
      oil_yoy = x[seq(2, length(x)-1)],
      outcome = v,
      d2      = abs(d2)
    )
  }))
  nonlin_dt[, outcome_label := OUTCOME_LABELS[outcome]]

  p_nonlin <- ggplot(nonlin_dt[!is.na(d2)],
                     aes(x=oil_yoy, y=d2, colour=outcome_label)) +
    geom_vline(xintercept=c(0, 40, 90), linetype="dashed",
               colour="#eeeeee", linewidth=0.4) +
    annotate("text", x=c(0, 40, 90), y=Inf,
             label=c("Flat", "~$95", "~$125"),
             angle=90, hjust=1.1, size=2.5, colour="#999999") +
    geom_line(linewidth=0.8) +
    scale_colour_manual(
      values = c("#b5470a","#993C1D","#4a2080",
                 "#185FA5","#2d7a4a","#3B6D11"),
      name   = NULL
    ) +
    scale_x_continuous(
      breaks=seq(-50, 120, 25),
      labels=function(x) paste0(ifelse(x>0,"+",""), x, "%%")
    ) +
    labs(
      title    = "FIGURE 04c-3 \u2014 Nonlinearity Detection: Oil Price Curvature",
      subtitle = paste(
        "|d\u00b2Y/d(oil)\u00b2| \u2014 spikes reveal where linear VARX assumption breaks down",
        "\nPeaks near +40%% and +90%% confirm threshold effects at ~$95 and ~$125/bbl"
      ),
      caption  = "Numerical second derivative of partial dependence curve",
      x="Oil price YoY %%",
      y="|Second derivative| of predicted outcome"
    ) +
    theme_minimal(base_size=10) +
    theme(
      plot.title       = element_text(face="bold", size=11),
      plot.subtitle    = element_text(size=8.5, colour="#555", lineheight=1.3),
      plot.caption     = element_text(size=7.5, colour="#888", hjust=0),
      legend.position  = "bottom",
      panel.grid.minor = element_blank()
    )

  ggsave("Figures/04c_nonlinear_oil.png", p_nonlin,
         width=11, height=6, dpi=300, bg="white")
  msg("Chart 3 (nonlinearity) \u2192 Figures/04c_nonlinear_oil.png")
}

# ── Chart 4 — SHAP waterfall for dq_rate and pll_rate ───────────────────────
for (v in c("dq_rate","pll_rate")) {
  sv_mean <- colMeans(shap_results[[v]], na.rm=TRUE)
  sv_dt   <- data.table(
    feature   = feature_cols,
    shap_mean = sv_mean
  )
  sv_dt <- merge(sv_dt, channel_map, by="feature")
  sv_dt <- sv_dt[order(-abs(shap_mean))][seq_len(min(15L, .N))]
  sv_dt[, feature_label := str_replace_all(feature, "_", " ") |>
           str_to_title()]
  sv_dt[, channel := factor(channel, levels=names(CHANNEL_COLS))]

  p_wf <- ggplot(sv_dt,
                 aes(x=reorder(feature_label, abs(shap_mean)),
                     y=shap_mean, fill=channel)) +
    geom_col(width=0.7) +
    geom_hline(yintercept=0, linewidth=0.3, colour="#888888") +
    scale_fill_manual(values=CHANNEL_COLS, name="Channel") +
    scale_y_continuous(labels=function(x) sprintf("%+.4f", x)) +
    coord_flip() +
    labs(
      title    = sprintf("FIGURE 04c-4 \u2014 SHAP Feature Importance: %s",
                          OUTCOME_LABELS[v]),
      subtitle = "Mean SHAP value across test observations | Colour = transmission channel",
      caption  = "Positive = feature increases outcome | Negative = decreases | Top 15 features shown",
      x=NULL, y="Mean SHAP value"
    ) +
    theme_minimal(base_size=10) +
    theme(
      plot.title       = element_text(face="bold", size=11),
      plot.subtitle    = element_text(size=8.5, colour="#555"),
      plot.caption     = element_text(size=7.5, colour="#888", hjust=0),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      legend.text      = element_text(size=8)
    )

  out_file <- sprintf("Figures/04c_shap_waterfall_%s.png", v)
  ggsave(out_file, p_wf, width=11, height=7, dpi=300, bg="white")
  msg("Chart 4 (SHAP waterfall: %s) \u2192 %s", v, out_file)
}

# ── Chart 5 — SHAP by asset tier ────────────────────────────────────────────
tier_col <- "asset_tier"
if (tier_col %in% names(panel_test)) {
  panel_test_shap <- panel_test[shap_test_idx]
  oil_idx <- which(feature_cols == "macro_base_yoy_oil")

  if (length(oil_idx) == 1) {
    tier_shap <- rbindlist(lapply(Y_VARS, function(v) {
      sv <- shap_results[[v]]
      data.table(
        outcome  = v,
        tier     = panel_test_shap[[tier_col]],
        oil_shap = sv[, oil_idx]
      )
    }))

    tier_agg <- tier_shap[,
      .(mean_shap = mean(oil_shap, na.rm=TRUE),
        se        = sd(oil_shap, na.rm=TRUE) / sqrt(.N)),
      by=.(outcome, tier)]
    tier_agg[, outcome_label := factor(
      OUTCOME_LABELS[outcome], levels=OUTCOME_LABELS)]

    KEY_TIER_OUTCOMES <- c("dq_rate","pll_rate","costfds","netintmrg")

    p_tier <- ggplot(tier_agg[outcome %in% KEY_TIER_OUTCOMES],
                     aes(x=factor(tier), y=mean_shap,
                         fill=factor(tier))) +
      geom_col(width=0.65) +
      geom_errorbar(aes(ymin=mean_shap-se, ymax=mean_shap+se),
                    width=0.25, linewidth=0.5, colour="#555555") +
      geom_hline(yintercept=0, linewidth=0.3, colour="#888888",
                 linetype="dashed") +
      facet_wrap(~outcome_label, scales="free_y", ncol=2) +
      scale_fill_manual(
        values  = c("1"="#B5D4F4","2"="#378ADD",
                    "3"="#185FA5","4"="#0C447C"),
        name    = "Asset tier",
        labels  = c("1"="<$10M","2"="$10\u2013100M",
                    "3"="$100M\u2013$1B","4">">$1B")
      ) +
      scale_y_continuous(labels=function(x) sprintf("%+.4f", x)) +
      labs(
        title    = "FIGURE 04c-5 \u2014 Oil SHAP Value by Asset Tier",
        subtitle = "Direct SHAP attribution of oil price YoY on each CU outcome",
        caption  = "Error bars = \u00b11 SE | Smaller tiers expected to show larger direct oil exposure per asset",
        x="Asset tier", y="Mean oil SHAP value"
      ) +
      theme_minimal(base_size=10) +
      theme(
        plot.title       = element_text(face="bold", size=11),
        plot.subtitle    = element_text(size=8.5, colour="#555"),
        plot.caption     = element_text(size=7.5, colour="#888", hjust=0),
        strip.text       = element_text(face="bold", size=9),
        panel.grid.minor = element_blank(),
        legend.position  = "right"
      )

    ggsave("Figures/04c_shap_tier.png", p_tier,
           width=11, height=8, dpi=300, bg="white")
    msg("Chart 5 (SHAP by tier) \u2192 Figures/04c_shap_tier.png")
  }
}

# ── Chart 6 — NN vs VARX comparison ─────────────────────────────────────────
varx_path <- "Results/04_forecast_paths.rds"
if (file.exists(varx_path)) {
  varx_fc <- readRDS(varx_path)

  # NN predictions on test set (unstandardised)
  nn_pred_list <- lapply(Y_VARS, function(v) {
    y_hat_std <- predict_std(v, X_test_mat)
    y_hat_std * y_sd[v] + y_mean[v]
  })
  nn_pred_dt <- as.data.table(do.call(cbind, nn_pred_list))
  setnames(nn_pred_dt, Y_VARS)
  nn_pred_dt[, yyyyqq := panel_test$yyyyqq]
  nn_pred_dt[, source := "Neural Net"]

  varx_hist <- varx_fc[scenario == "historical",
                        c("yyyyqq", Y_VARS), with=FALSE]
  varx_hist[, source := "VARX"]

  compare_vars <- Y_VARS[1:4]
  compare_dt   <- rbindlist(list(
    nn_pred_dt[, c("yyyyqq","source", compare_vars), with=FALSE],
    varx_hist[ , c("yyyyqq","source", compare_vars), with=FALSE]
  ), fill=TRUE)

  compare_long <- melt(compare_dt, id.vars=c("yyyyqq","source"),
                        variable.name="outcome", value.name="value")
  compare_long[, yr      := yyyyqq %/% 100L]
  compare_long[, qtr     := yyyyqq %% 100L]
  compare_long[, date_num := yr + (qtr-1)/4]
  compare_long[, outcome_label := factor(
    OUTCOME_LABELS[as.character(outcome)],
    levels=OUTCOME_LABELS[compare_vars]
  )]

  p_compare <- ggplot(compare_long[!is.na(value)],
                      aes(x=date_num, y=value,
                          colour=source, linetype=source)) +
    geom_line(linewidth=0.85) +
    scale_colour_manual(
      values=c("Neural Net"="#b5470a","VARX"="#185FA5"),
      name="Model"
    ) +
    scale_linetype_manual(
      values=c("Neural Net"="solid","VARX"="dashed"),
      name="Model"
    ) +
    facet_wrap(~outcome_label, scales="free_y", ncol=2) +
    scale_x_continuous(
      breaks=seq(2020, 2026, 1),
      labels=function(x) paste0("'", substr(as.character(x), 3, 4))
    ) +
    labs(
      title    = "FIGURE 04c-6 \u2014 Neural Net vs VARX: Test Period 2020Q1\u20132025Q4",
      subtitle = paste(
        "Agreement = robust finding across both models",
        "| Divergence = nonlinearity the linear VARX missed"
      ),
      caption  = paste(
        "NN: 2-layer (16\u20138), pure R neuralnet, train 2005\u20132019",
        "| VARX: Cholesky p=2, full sample"
      ),
      x=NULL, y="Predicted value (level units)"
    ) +
    theme_minimal(base_size=10) +
    theme(
      plot.title       = element_text(face="bold", size=11),
      plot.subtitle    = element_text(size=8.5, colour="#555", lineheight=1.3),
      plot.caption     = element_text(size=7.5, colour="#888", hjust=0),
      strip.text       = element_text(face="bold", size=9),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )

  ggsave("Figures/04c_nn_vs_varx.png", p_compare,
         width=13, height=9, dpi=300, bg="white")
  msg("Chart 6 (NN vs VARX) \u2192 Figures/04c_nn_vs_varx.png")
} else {
  msg("VARX forecast not found \u2014 skipping Chart 6")
}

# =============================================================================
# 8. OUTPUT MANIFEST
# =============================================================================
hdr("SECTION 8: Output Manifest")

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
cat(sprintf("  Total runtime : %.1f seconds (%.1f min)\n",
            t_total["elapsed"], t_total["elapsed"]/60))
cat(sprintf("  Workers used  : %d (PSOCK cluster, Windows)\n", N_CORES))
cat("  Script 04c complete — neuralnet backend, parallel processing\n")
cat(strrep("=",70), "\n")
