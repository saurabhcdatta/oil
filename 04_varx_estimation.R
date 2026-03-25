# =============================================================================
# OIL PRICE SHOCK × CREDIT UNION RESEARCH
# Script 04b — VARX Model Estimation
# =============================================================================
#
# Inputs  : Data/panel_model.rds         — from Script 03
#           Data/macro_base.rds          — CCAR 2026 baseline exogenous paths
#           Data/macro_severe.rds        — CCAR 2026 severely adverse paths
#
# Outputs : Models/varx_full.rds         — full-sample VARX object
#           Models/varx_pre2015.rds      — pre-shale subsample
#           Models/varx_post2015.rds     — post-shale subsample
#           Results/04_varx_coef.csv     — coefficient tables (all specs)
#           Results/04_varx_irf.rds      — impulse response functions
#           Results/04_forecast_paths.rds — CCAR scenario projections
#           Figures/04_irf_*.png         — IRF plots
#           Figures/04_forecast_*.png    — scenario fan charts
#
# Model Architecture:
#   Endogenous Y(t): dq_rate, pll_rate, netintmrg, insured_share_growth,
#                    member_growth_yoy, costfds, loan_to_share, pcanetworth
#   Exogenous  Z(t): macro_base_yoy_oil, macro_base_lurc, macro_base_pcpi,
#                    macro_base_yield_curve, macro_base_rmtg, hpi_yoy,
#                    macro_base_fomc_regime
#   Interactions  : fomc_x_brent, post_x_oil,
#                   post_x_oil_x_direct (triple)
#   Dummies       : post_shale, gfc_dummy, covid_dummy, zirp_era, hike_cycle
#   Lag order p   : selected by AIC/BIC in Script 03 (default p=2)
#
# Estimation Strategy:
#   Stage 1  — Aggregate system VARX (vars::VAR with exogenous)
#   Stage 2  — Panel fixed-effects VARX equation-by-equation
#              with CU + quarter FE (lfe::felm)
#   Stage 3  — Asset-tier stratified models (4 tiers)
#   Stage 4  — Pre/post-2015 structural split
#   Stage 5  — IRF + FEVD (Cholesky ordering justified below)
#   Stage 6  — CCAR 2026 scenario fan-chart projections
#
# Cholesky ordering (most → least exogenous within Y block):
#   costfds (rate channel, first to reprice)
#   → pll_rate (forward-looking provisions, leads credit loss)
#   → dq_rate (realised delinquency, lags provisions)
#   → netintmrg (NIM compressed by CoF before loan repricing)
#   → loan_to_share (balance sheet demand)
#   → pcanetworth (capital absorbs losses)
#   → member_growth_yoy (membership response, early warning)
#   → insured_share_growth (deposit flow, last behavioural response)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(vars)          # VAR(), irf(), fevd()
  library(lfe)           # felm() — panel FE with clustering
  library(plm)           # pdim(), pdata.frame()
  library(sandwich)      # vcovHC for robust SE
  library(lmtest)        # coeftest
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

# ── Helpers ──────────────────────────────────────────────────────────────────
hdr <- function(x) cat("\n", strrep("=", 70), "\n##", x, "\n",
                        strrep("=", 70), "\n")
msg <- function(fmt, ...) cat(sprintf(paste0("  ", fmt, "\n"), ...))
ts  <- function() format(Sys.time(), "%H:%M:%S")

dir.create("Models",   showWarnings = FALSE)
dir.create("Results",  showWarnings = FALSE)
dir.create("Figures",  showWarnings = FALSE)

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
hdr("SECTION 0: Configuration")

P_LAG       <- 2L          # Lag order from Script 03 AIC/BIC
N_BOOT      <- 500L        # Bootstrap reps for IRF bands
IRF_HORIZON <- 12L         # Quarters ahead for IRF
FCST_HOR    <- 12L         # Forecast horizon (quarters) — 2026Q1–2028Q4
SEED        <- 20260101L
set.seed(SEED)

# Endogenous block — Y(t)
# Column names exactly as they appear in panel_model.rds (from Script 03)
Y_VARS <- c(
  "dq_rate",             # Delinquency rate (backward-looking credit quality)
  "pll_rate",            # PLL / avg loans — forward-looking provisions (NEW)
  "netintmrg",           # Net interest margin
  "insured_share_growth",# YoY growth in insured shares (deposit channel)
  "member_growth_yoy",   # YoY membership growth (early warning indicator) (NEW)
  "costfds",             # Cost of funds (rate channel)
  "loan_to_share",       # Loan-to-share ratio (balance sheet leverage)
  "pcanetworth"          # Net worth / total assets (capital adequacy)
)

# Exogenous block — Z(t)  [plugged in from CCAR at forecast time]
# Column names exactly as they appear in macro_base.rds / macro_severe.rds
Z_VARS <- c(
  "macro_base_yoy_oil",      # Brent crude oil YoY % change (PBRENT)
  "macro_base_lurc",         # Unemployment rate
  "macro_base_pcpi",         # CPI level
  "macro_base_yield_curve",  # 10Y - 3M spread (derived in Script 01)
  "macro_base_rmtg",         # 30Y mortgage rate
  "hpi_yoy",                 # HPI YoY (derived in Script 03)
  "macro_base_fomc_regime"   # 1=hiking, -1=cutting, 0=hold (derived in Script 01)
)

# Interaction terms (pre-computed in Script 03)
INTERACT_VARS <- c(
  "fomc_x_brent",        # fomc_regime × yoy_oil
  "post_x_oil",          # post_shale × yoy_oil
  "post_x_oil_x_direct"  # post_shale × yoy_oil × oil_exposure_bin (triple)
)

# Structural dummies
DUMMY_VARS <- c(
  "post_shale",          # I(yyyyqq >= 201501)
  "gfc_dummy",           # I(200803 <= yyyyqq <= 200904)
  "covid_dummy",         # I(202001 <= yyyyqq <= 202104)
  "zirp_era",            # I(200901 <= yyyyqq <= 201504)
  "hike_cycle"           # I(202201 <= yyyyqq <= 202304)
)

msg("Lag order p    : %d", P_LAG)
msg("IRF horizon    : %d quarters", IRF_HORIZON)
msg("Forecast hor.  : %d quarters", FCST_HOR)
msg("Bootstrap reps : %d", N_BOOT)
msg("Endogenous vars: %s", paste(Y_VARS, collapse=", "))

# =============================================================================
# 1. LOAD DATA
# =============================================================================
hdr("SECTION 1: Load Data")

panel      <- readRDS("Data/panel_model.rds")
macro_base <- readRDS("Data/macro_base.rds")

# Load severely adverse if available
macro_severe <- tryCatch(readRDS("Data/macro_severe.rds"),
                          error=function(e) { msg("macro_severe.rds not found — severe scenario skipped"); NULL })

setDT(panel)
setDT(macro_base)
if (!is.null(macro_severe)) setDT(macro_severe)

# Add scenario tag and unify for multi-scenario forecast
macro_base[,  scenario := "baseline"]
if (!is.null(macro_severe)) macro_severe[, scenario := "severely_adverse"]
macro <- rbindlist(list(macro_base,
                         if (!is.null(macro_severe)) macro_severe else NULL),
                    fill=TRUE)

# ── macro_severe_aligned: rename macro_severe_* cols to macro_base_* ─────────
# The forecast function expects exogenous paths in macro_base_* naming.
# For the severely adverse scenario we rename so the same function works.
if (!is.null(macro_severe)) {
  macro_severe_aligned <- copy(macro_severe)
  sev_cols  <- grep("^macro_severe_", names(macro_severe_aligned), value=TRUE)
  base_cols <- gsub("^macro_severe_", "macro_base_", sev_cols)
  # Only rename cols that don't already exist as macro_base_*
  rename_idx <- !base_cols %in% names(macro_severe_aligned)
  if (any(rename_idx))
    setnames(macro_severe_aligned,
             sev_cols[rename_idx], base_cols[rename_idx])
  msg("macro_severe_aligned: %d cols renamed to macro_base_* naming",
      sum(rename_idx))
} else {
  macro_severe_aligned <- NULL
  msg("macro_severe_aligned: NULL (severe scenario not available)")
}

# Ensure cal_date present on macro (may arrive without it)
Q_MONTH <- c("1"=1L,"2"=4L,"3"=7L,"4"=10L)
for (.nm in c("macro_base","macro_severe","macro")) {
  .obj <- tryCatch(get(.nm), error=function(e) NULL)
  if (!is.null(.obj) && is.data.table(.obj) && !"cal_date" %in% names(.obj))
    .obj[, cal_date := as.Date(paste(year, Q_MONTH[as.character(quarter)],
                                      "01", sep="-"))]
}
rm(.nm, .obj)

msg("Panel obs      : %s", format(nrow(panel), big.mark=","))
msg("Unique CUs     : %s", format(uniqueN(panel$join_number), big.mark=","))
msg("Date range     : %d to %d", min(panel$yyyyqq), max(panel$yyyyqq))

# Verify all required variables exist — soft stop with informative message
all_required <- c(Y_VARS, Z_VARS, INTERACT_VARS, DUMMY_VARS,
                  "join_number", "yyyyqq", "asset_tier")
missing_vars <- setdiff(all_required, names(panel))
if (length(missing_vars) > 0) {
  msg("WARNING — Missing variables from panel_model.rds:")
  for (mv in missing_vars) msg("  MISSING: %s", mv)
  msg("Re-run Scripts 01-03, then retry. Continuing with available vars...")
  # Trim Y_VARS to only available columns
  Y_VARS <- intersect(Y_VARS, names(panel))
  Z_VARS <- intersect(Z_VARS, names(macro_base))
} else {
  msg("All required variables present \u2713")
}
msg("Y_VARS active  : %s", paste(Y_VARS, collapse=", "))
msg("Z_VARS active  : %s", paste(Z_VARS, collapse=", "))

# =============================================================================
# 2. AGGREGATE TIME-SERIES FOR STAGE 1 VARX
# =============================================================================
hdr("SECTION 2: Aggregate Time-Series Construction")

# Asset-weighted aggregate CU system series
# Use total assets as weight for balance-sheet ratios;
# simple mean for growth rates where weighting would be circular

# Asset-weighted aggregate CU system series
# Balance-sheet ratios: weight by assets; growth rates: simple mean
WEIGHTED_VARS <- intersect(
  c("netintmrg","costfds","loan_to_share","pcanetworth","dq_rate","pll_rate"),
  names(panel))
MEAN_VARS     <- intersect(
  c("insured_share_growth","member_growth_yoy"),
  names(panel))

msg("WEIGHTED_VARS  : %s", paste(WEIGHTED_VARS, collapse=", "))
msg("MEAN_VARS      : %s", paste(MEAN_VARS, collapse=", "))

# Resolve asset weight column — try several common names
asset_col <- intersect(c("assets_tot","total_assets","acct_010","assets"),
                        names(panel))[1]
if (is.na(asset_col)) {
  msg("WARNING: no asset column found — using simple mean for all vars")
  asset_col <- NULL
}

agg <- panel[, {
  if (!is.null(asset_col) && asset_col %in% names(.SD)) {
    w <- get(asset_col) / sum(get(asset_col), na.rm=TRUE)
  } else {
    w <- rep(1 / .N, .N)
  }
  wtd <- lapply(WEIGHTED_VARS, function(v) {
    if (v %in% names(.SD)) weighted.mean(.SD[[v]], w, na.rm=TRUE) else NA_real_
  })
  mn <- lapply(MEAN_VARS, function(v) {
    if (v %in% names(.SD)) mean(.SD[[v]], na.rm=TRUE) else NA_real_
  })
  setNames(c(wtd, mn), c(WEIGHTED_VARS, MEAN_VARS))
}, by=yyyyqq]

setorder(agg, yyyyqq)

# ── WEIGHT_COL: used in tier-stratified aggregation (Section 5) ──────────────
WEIGHT_COL <- intersect(c("assets_tot","total_assets","acct_010","assets"),
                         names(panel))[1]
if (is.na(WEIGHT_COL)) {
  msg("WARNING: no asset weight column found — tier aggregation will use simple mean")
  # Create a unit-weight column so the code doesn't break
  panel[, .wt_unit := 1L]
  WEIGHT_COL <- ".wt_unit"
}
msg("WEIGHT_COL     : %s", WEIGHT_COL)

# ── panel_dum: quarter-level dummies pre-aggregated for merging ───────────────
# These are time-series dummies (same value for all CUs in a quarter)
# Pre-extract to avoid repeated panel aggregation in the tier loop
dum_cols_avail <- intersect(DUMMY_VARS, names(panel))
panel_dum <- unique(panel[, c("yyyyqq", dum_cols_avail), with=FALSE])
setorder(panel_dum, yyyyqq)
msg("panel_dum cols : %s", paste(names(panel_dum), collapse=", "))

# ── avail_* : resolve which exogenous columns actually exist ──────────────────
# Used in both tier loop (Section 5) and subsample loop (Section 6)
avail_z          <- intersect(Z_VARS,        names(agg))
avail_int_panel  <- intersect(INTERACT_VARS, names(agg))
avail_dum        <- intersect(DUMMY_VARS,    names(agg))
msg("avail_z        : %s", paste(avail_z,         collapse=", "))
msg("avail_int_panel: %s", paste(avail_int_panel, collapse=", "))
msg("avail_dum      : %s", paste(avail_dum,       collapse=", "))
macro_ts <- macro_base[, intersect(c("yyyyqq", Z_VARS, INTERACT_VARS,
                                       DUMMY_VARS), names(macro_base)),
                         with=FALSE]
agg <- merge(agg, macro_ts, by="yyyyqq", all.x=TRUE)

msg("Aggregate series: %d quarters (%d to %d)",
    nrow(agg), min(agg$yyyyqq), max(agg$yyyyqq))
msg("NA check (Y vars):")
for (v in Y_VARS)
  msg("  %-30s : %d NAs", v, sum(is.na(agg[[v]])))

# =============================================================================
# 3. STAGE 1 — AGGREGATE VARX (vars package)
# =============================================================================
hdr("SECTION 3: Stage 1 — Aggregate VARX")

# ── 3.1 Build matrices ───────────────────────────────────────────────────────
# Drop rows with any NA in endogenous block
agg_clean <- agg[complete.cases(agg[, ..Y_VARS])]

Y_mat <- as.matrix(agg_clean[, ..Y_VARS])
rownames(Y_mat) <- as.character(agg_clean$yyyyqq)

# Exogenous: Z + interactions + dummies
exog_cols <- intersect(c(Z_VARS, INTERACT_VARS, DUMMY_VARS), names(agg_clean))
X_mat <- as.matrix(agg_clean[, ..exog_cols])

msg("Y matrix : %d × %d", nrow(Y_mat), ncol(Y_mat))
msg("X matrix : %d × %d", nrow(X_mat), ncol(X_mat))

# ── 3.2 Lag order selection (verify / cross-check Script 03) ─────────────────
varselect_full <- VARselect(Y_mat, lag.max=6, type="const",
                             exogen=X_mat[-(1:6), ])  # trim for lag alignment

ic_table <- rbind(
  AIC = varselect_full$criteria["AIC(n)", ],
  BIC = varselect_full$criteria["SC(n)", ],
  HQ  = varselect_full$criteria["HQ(n)", ]
)
cat("\n  Information Criteria by Lag Order:\n")
print(round(ic_table, 4))

p_aic <- varselect_full$selection["AIC(n)"]
p_bic <- varselect_full$selection["SC(n)"]
msg("AIC selects p=%d | BIC selects p=%d | Using configured p=%d",
    p_aic, p_bic, P_LAG)

if (P_LAG != p_aic)
  warning(sprintf("Configured P_LAG=%d differs from AIC selection=%d. Review.",
                  P_LAG, p_aic))

# ── 3.3 Estimate VARX ────────────────────────────────────────────────────────
# vars::VAR with exogen argument implements VARX directly
# type="const" adds intercept; seasonal=0 (quarterly dummies already in X_mat)

varx_full <- VAR(
  y      = Y_mat,
  p      = P_LAG,
  type   = "const",
  exogen = X_mat[-(1:P_LAG), ]   # align rows after lag trimming
)

msg("VARX estimated: %d equations, %d obs, p=%d", length(Y_VARS),
    nrow(summary(varx_full)$varresult[[1]]$model), P_LAG)

# ── 3.4 Model diagnostics ────────────────────────────────────────────────────
cat("\n  --- Equation-level R² ---\n")
for (eq_nm in names(varx_full$varresult)) {
  s   <- summary(varx_full$varresult[[eq_nm]])
  adj <- s$adj.r.squared
  msg("  %-30s : Adj-R² = %.3f", eq_nm, adj)
}

# Serial correlation: Portmanteau test
pt <- serial.test(varx_full, lags.pt=10, type="PT.asymptotic")
msg("\n  Portmanteau (lags=10) p-value: %.4f  [H0: no serial corr; want > 0.05]",
    pt$serial$p.value)

# Stability: check eigenvalues (roots should be < 1)
roots <- roots(varx_full)
msg("  Max modulus root: %.4f  [< 1 = stable]", max(roots))
if (max(roots) >= 1)
  warning("VARX is NOT stable — roots on or outside unit circle. Check specification.")

# Normality (Jarque-Bera)
jb <- normality.test(varx_full, multivariate.only=TRUE)
msg("  Multivariate Jarque-Bera p: %.4f  [note: non-normality common in fin. data]",
    jb$jb.mul$JB$p.value)

saveRDS(varx_full, "Models/varx_full.rds")
msg("\n  Saved → Models/varx_full.rds")

# =============================================================================
# 4. STAGE 2 — PANEL FIXED-EFFECTS VARX (equation-by-equation)
# =============================================================================
hdr("SECTION 4: Stage 2 — Panel FE VARX")

# ── Approach ──────────────────────────────────────────────────────────────────
# For each Y variable, estimate:
#   y_it = α_i + δ_t + Σ_k β_k y_{i,t-k} + γ Z_t + ε_it
#
# α_i = CU fixed effect (controls for time-invariant CU heterogeneity)
# δ_t = quarter FE (controls for common macro shocks NOT in Z)
# SE clustered two-way: CU + quarter
#
# Using lfe::felm() for speed on large panel
# Formula: y ~ [X] | join_number + yyyyqq | 0 | join_number + yyyyqq

# ── 4.1 Build lagged Y matrix in panel ───────────────────────────────────────
setorder(panel, join_number, yyyyqq)

for (v in Y_VARS) {
  for (k in 1:P_LAG) {
    lag_nm <- paste0(v, "_lag", k)
    panel[, (lag_nm) := shift(.SD[[v]], k, type="lag"), by=join_number, .SDcols=v]
  }
}

LAG_Y_VARS <- unlist(lapply(Y_VARS, function(v)
  paste0(v, "_lag", 1:P_LAG)))

# ── 4.2 Estimate each equation ───────────────────────────────────────────────
# RHS: lagged Y + Z + interactions + dummies (no CU/quarter FE in formula;
#      those go in the | FE | cluster block of felm)

panel_exog_cols <- intersect(c(Z_VARS, INTERACT_VARS, DUMMY_VARS), names(panel))

fe_models <- list()
fe_coefs  <- list()

for (dep in Y_VARS) {
  lag_terms  <- paste0(dep, "_lag", 1:P_LAG)
  # Cross-lag terms (all other Y lags) — key for VARX
  cross_lags <- setdiff(LAG_Y_VARS, lag_terms)

  rhs <- paste(c(lag_terms, cross_lags, panel_exog_cols), collapse=" + ")

  frml <- as.formula(paste0(
    dep, " ~ ", rhs,
    " | join_number + yyyyqq | 0 | join_number + yyyyqq"
  ))

  # Drop incomplete cases for this equation
  keep_cols <- c(dep, lag_terms, cross_lags, panel_exog_cols,
                 "join_number", "yyyyqq")
  keep_cols <- intersect(keep_cols, names(panel))
  panel_eq  <- panel[complete.cases(panel[, ..keep_cols])]

  fe_models[[dep]] <- tryCatch(
    felm(frml, data=panel_eq),
    error=function(e) { warning(paste("felm failed for", dep, ":", e$message)); NULL }
  )

  if (!is.null(fe_models[[dep]])) {
    cf <- as.data.table(summary(fe_models[[dep]])$coefficients, keep.rownames=TRUE)
    setnames(cf, c("term","estimate","se","t","p"))
    cf[, dep_var := dep]
    fe_coefs[[dep]] <- cf
    msg("  %-30s : Adj-R² = %.3f | N = %s",
        dep,
        summary(fe_models[[dep]])$adj.r.squared,
        format(nrow(fe_models[[dep]]$response), big.mark=","))
  }
}

# Combine coefficient table
coef_table <- rbindlist(fe_coefs)
fwrite(coef_table, "Results/04_varx_coef_panel_fe.csv")
msg("\n  Coefficient table → Results/04_varx_coef_panel_fe.csv")

# =============================================================================
# 5. STAGE 3 — ASSET-TIER STRATIFIED MODELS
# =============================================================================
hdr("SECTION 5: Stage 3 — Asset Tier Stratified VARX")

# Tiers (from Script 01):
#   T1: < $10M   T2: $10–$100M   T3: $100M–$1B   T4: > $1B

tier_models <- list()
tier_agg    <- list()

for (tier in sort(unique(panel$asset_tier))) {

  panel_t <- panel[asset_tier == tier]

  # Aggregate within tier (asset-weighted for ratios, simple mean for growth)
  agg_t <- panel_t[, {
    wt <- get(WEIGHT_COL)
    w  <- wt / sum(wt, na.rm=TRUE)
    wtd <- lapply(WEIGHTED_VARS, function(v) {
      if (v %in% names(.SD)) weighted.mean(.SD[[v]], w, na.rm=TRUE) else NA_real_
    })
    mn <- lapply(MEAN_VARS, function(v) {
      if (v %in% names(.SD)) mean(.SD[[v]], na.rm=TRUE) else NA_real_
    })
    setNames(c(wtd, mn), c(WEIGHTED_VARS, MEAN_VARS))
  }, by=yyyyqq, .SDcols=c(WEIGHTED_VARS, MEAN_VARS, WEIGHT_COL)]

  setorder(agg_t, yyyyqq)
  agg_t <- merge(agg_t, macro_ts,  by="yyyyqq", all.x=TRUE)
  # panel_dum already contains DUMMY_VARS at quarter level — merge once
  agg_t <- merge(agg_t, panel_dum, by="yyyyqq", all.x=TRUE, allow.cartesian=FALSE)
  tier_agg[[as.character(tier)]] <- agg_t

  agg_t_clean <- agg_t[complete.cases(agg_t[, intersect(Y_VARS, names(agg_t)), with=FALSE])]
  if (nrow(agg_t_clean) < 30) {
    msg("Tier %s: insufficient obs (%d) — skipping", tier, nrow(agg_t_clean))
    next
  }

  y_cols_t  <- intersect(Y_VARS, names(agg_t_clean))
  Y_t       <- as.matrix(agg_t_clean[, ..y_cols_t])
  exog_t_cols <- intersect(c(avail_z, avail_int_panel, avail_dum), names(agg_t_clean))
  X_t <- if (length(exog_t_cols) > 0)
    as.matrix(agg_t_clean[, ..exog_t_cols])
  else
    matrix(0, nrow=nrow(agg_t_clean), ncol=1)

  tier_models[[as.character(tier)]] <- tryCatch(
    VAR(Y_t, p=P_LAG, type="const",
        exogen=X_t[-(1:P_LAG), , drop=FALSE]),
    error=function(e) {
      warning(paste("VAR failed for tier", tier, ":", e$message)); NULL
    })

  if (!is.null(tier_models[[as.character(tier)]])) {
    r2s <- sapply(tier_models[[as.character(tier)]]$varresult,
                  function(m) summary(m)$adj.r.squared)
    msg("Tier %-4s | N_CU=%-6s | Avg Adj-R\u00b2=%.3f",
        tier,
        format(uniqueN(panel_t$join_number), big.mark=","),
        mean(r2s, na.rm=TRUE))
  }
}

saveRDS(tier_models, "Models/varx_tiers.rds")
msg("\nTier models → Models/varx_tiers.rds")

# =============================================================================
# 6. STAGE 4 — PRE/POST-2015 STRUCTURAL SPLIT
# =============================================================================
hdr("SECTION 6: Stage 4 — Pre/Post-2015 Structural Split")

estimate_subsample <- function(data_full, macro_data, label, yyyyqq_filter) {

  sub  <- data_full[yyyyqq_filter]
  mac  <- macro_data[yyyyqq_filter]
  sub  <- merge(sub, mac, by="yyyyqq", all.x=TRUE)
  sub  <- sub[complete.cases(sub[, ..Y_VARS])]

  if (nrow(sub) < P_LAG + ncol(Y_mat) + 5) {
    msg("  %s: too few obs (%d) — skipping", label, nrow(sub)); return(NULL)
  }

  Y_s <- as.matrix(sub[, ..Y_VARS])
  exog_s_cols <- intersect(c(Z_VARS, INTERACT_VARS, DUMMY_VARS), names(sub))
  # Remove post_shale (collinear within subsample)
  exog_s_cols <- setdiff(exog_s_cols, "post_shale")
  X_s <- as.matrix(sub[, ..exog_s_cols])

  m <- tryCatch(
    VAR(Y_s, p=P_LAG, type="const",
        exogen=X_s[-(1:P_LAG), ]),
    error=function(e) { warning(paste(label, "failed:", e$message)); NULL }
  )

  if (!is.null(m)) {
    r2s <- sapply(m$varresult, function(eq) summary(eq)$adj.r.squared)
    msg("  %-15s | Obs=%-5d | Avg Adj-R²=%.3f | Max root=%.4f",
        label, nrow(Y_s), mean(r2s, na.rm=TRUE), max(roots(m)))
  }
  m
}

varx_pre  <- estimate_subsample(agg, macro_ts, "Pre-2015",
                                 agg$yyyyqq < 201501L)
varx_post <- estimate_subsample(agg, macro_ts, "Post-2015",
                                 agg$yyyyqq >= 201501L)

saveRDS(varx_pre,  "Models/varx_pre2015.rds")
saveRDS(varx_post, "Models/varx_post2015.rds")
msg("\nSubsample models saved → Models/varx_pre2015.rds | varx_post2015.rds")

# ── 6.1 Coefficient comparison table: pre vs post ────────────────────────────
extract_coefs <- function(var_model, label, focus_term="yoy_oil") {
  if (is.null(var_model)) return(NULL)
  lapply(names(var_model$varresult), function(eq_nm) {
    cf <- coef(var_model$varresult[[eq_nm]])
    cf_dt <- data.table(
      dep_var = eq_nm,
      term    = names(cf),
      est     = cf,
      period  = label
    )
    cf_dt[str_detect(term, focus_term)]
  }) |> rbindlist()
}

oil_coefs_compare <- rbindlist(list(
  extract_coefs(varx_full, "Full Sample"),
  extract_coefs(varx_pre,  "Pre-2015"),
  extract_coefs(varx_post, "Post-2015")
))

fwrite(oil_coefs_compare, "Results/04_oil_coef_comparison.csv")
msg("Oil coefficient comparison → Results/04_oil_coef_comparison.csv")

# =============================================================================
# 7. STAGE 5 — IMPULSE RESPONSE FUNCTIONS
# =============================================================================
hdr("SECTION 7: Stage 5 — Impulse Response Functions")

# ── 7.1 Cholesky ordering ─────────────────────────────────────────────────────
# 9-variable system ordered most → least exogenous within Y block.
# Theoretical justification:
#   pll_rate      — provisions respond first (management forward-looking signal)
#   dq_rate       — realised delinquency lags pll_rate by 1-2 quarters
#   costfds       — rate channel: reprices quickly with Fed hikes
#   netintmrg     — NIM compressed after CoF rises, before loans reprice
#   loan_to_share — balance sheet leverage: demand-driven, slower
#   pcanetworth   — capital: absorbs accumulated losses, most persistent
#   insured_share_growth — deposit flows: member income/behavioural response
#   member_growth_yoy    — membership: structural, slow-moving
#   cert_share           — deposit mix: certificate migration is last to adjust

chol_order <- intersect(
  c("pll_rate","dq_rate","costfds","netintmrg","loan_to_share",
    "pcanetworth","insured_share_growth","member_growth_yoy","cert_share"),
  colnames(Y_mat))

msg("Cholesky order (%d vars): %s", length(chol_order),
    paste(chol_order, collapse=" → "))

# Reorder Y matrix columns to match Cholesky ordering
chol_idx   <- match(chol_order, colnames(Y_mat))
chol_idx   <- chol_idx[!is.na(chol_idx)]
Y_chol     <- Y_mat[, chol_idx]

varx_chol <- VAR(Y_chol, p=P_LAG, type="const",
                  exogen=X_mat[-(1:P_LAG), , drop=FALSE])

# ── 7.2 Compute IRFs ─────────────────────────────────────────────────────────
# Shock 1: pll_rate — the forward-looking credit quality shock
#   Captures: when management raises provisions, how does the rest of the
#   CU system respond over 12 quarters? Expect dq_rate to follow 1-2Q later.
# Shock 2: costfds — the rate channel shock
#   Captures: when funding costs rise (Fed hike transmission), how does NIM,
#   loan demand, membership, and deposit mix respond?

irf_pll <- irf(
  varx_chol,
  impulse  = if ("pll_rate" %in% chol_order) "pll_rate" else chol_order[1],
  response = chol_order,
  n.ahead  = IRF_HORIZON,
  boot     = TRUE, ci=0.90, runs=N_BOOT, ortho=TRUE
)

irf_cof <- irf(
  varx_chol,
  impulse  = if ("costfds" %in% chol_order) "costfds" else chol_order[1],
  response = chol_order,
  n.ahead  = IRF_HORIZON,
  boot     = TRUE, ci=0.90, runs=N_BOOT, ortho=TRUE
)

# Keep delinquency shock for comparison with earlier version
irf_delq <- irf(
  varx_chol,
  impulse  = if ("dq_rate" %in% chol_order) "dq_rate" else chol_order[2],
  response = chol_order,
  n.ahead  = IRF_HORIZON,
  boot     = TRUE, ci=0.90, runs=N_BOOT, ortho=TRUE
)

saveRDS(list(irf_pll=irf_pll, irf_cof=irf_cof, irf_delq=irf_delq),
        "Results/04_varx_irf.rds")
msg("IRFs saved \u2192 Results/04_varx_irf.rds")

# ── 7.3 Plot IRFs ─────────────────────────────────────────────────────────────
plot_irf <- function(irf_obj, impulse_label, outfile) {

  # Extract IRF data into tidy data.table
  extract_one <- function(series_nm) {
    irf_dt <- data.table(
      horizon = 0:IRF_HORIZON,
      irf     = irf_obj$irf[[series_nm]],
      lower   = irf_obj$Lower[[series_nm]],
      upper   = irf_obj$Upper[[series_nm]],
      response = series_nm
    )
    irf_dt
  }

  all_irf <- rbindlist(lapply(names(irf_obj$irf), extract_one))

  # Clean labels
  all_irf[, response_lbl := str_replace_all(response, "_", " ") |>
                              str_to_title()]

  p <- ggplot(all_irf, aes(x=horizon)) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="#2196F3", alpha=0.15) +
    geom_line(aes(y=irf), color="#1565C0", linewidth=0.8) +
    geom_hline(yintercept=0, linetype="dashed", color="#666", linewidth=0.4) +
    facet_wrap(~response_lbl, scales="free_y", ncol=3) +
    scale_x_continuous(breaks=seq(0, IRF_HORIZON, 4),
                       labels=function(x) paste0("Q", x)) +
    labs(
      title    = sprintf("Impulse Response: %s shock → CU system", impulse_label),
      subtitle = sprintf("Cholesky-identified; 90%% bootstrap CI (%d runs); p=%d",
                         N_BOOT, P_LAG),
      x        = "Quarters after shock",
      y        = "Response (level units of dependent var)"
    ) +
    theme_minimal(base_size=10) +
    theme(
      plot.title    = element_text(face="bold", size=11),
      strip.text    = element_text(face="bold", size=9),
      panel.grid.minor = element_blank()
    )

  ggsave(outfile, p, width=12, height=7, dpi=150)
  msg("  IRF chart → %s", outfile)
  invisible(p)
}

plot_irf(irf_pll,  "PLL Rate (Forward-Looking Provisions)", "Figures/04_irf_pll_rate.png")
plot_irf(irf_cof,  "Cost of Funds (Rate Channel)",          "Figures/04_irf_cost_of_funds.png")
plot_irf(irf_delq, "Delinquency Rate (Credit Channel)",     "Figures/04_irf_delinquency.png")

# ── 7.4 FEVD (Forecast Error Variance Decomposition) ─────────────────────────
fevd_obj  <- fevd(varx_chol, n.ahead=IRF_HORIZON)

# Convert to data.table for export
fevd_list <- lapply(names(fevd_obj), function(resp) {
  dt <- as.data.table(fevd_obj[[resp]])
  dt[, `:=`(response=resp, horizon=1:.N)]
  melt(dt, id.vars=c("response","horizon"),
       variable.name="source", value.name="share")
}) |> rbindlist()

fwrite(fevd_list, "Results/04_fevd.csv")
msg("FEVD table → Results/04_fevd.csv")

# FEVD heatmap at horizon 8 (2 years)
fevd_h8 <- fevd_list[horizon==8]
fevd_h8[, response_lbl := str_replace_all(response, "_", " ") |> str_to_title()]
fevd_h8[, source_lbl   := str_replace_all(as.character(source), "_", " ") |>
                            str_to_title()]

p_fevd <- ggplot(fevd_h8, aes(x=source_lbl, y=response_lbl, fill=share)) +
  geom_tile(color="white", linewidth=0.3) +
  geom_text(aes(label=sprintf("%.0f%%", share*100)),
            size=3, fontface="bold") +
  scale_fill_gradient2(low="white", mid="#BBDEFB", high="#1565C0",
                       midpoint=0.25, limits=c(0,1),
                       labels=percent_format()) +
  labs(
    title    = "Forecast Error Variance Decomposition — Horizon 8Q",
    subtitle = "Share of forecast variance explained by each structural shock",
    x="Source shock", y="Response variable", fill="Share"
  ) +
  theme_minimal(base_size=10) +
  theme(
    axis.text.x  = element_text(angle=35, hjust=1),
    plot.title   = element_text(face="bold"),
    legend.position = "right"
  )

ggsave("Figures/04_fevd_h8.png", p_fevd, width=9, height=6, dpi=150)
msg("FEVD heatmap → Figures/04_fevd_h8.png")

# =============================================================================
# 8. STAGE 6 — CCAR 2026 SCENARIO PROJECTIONS
# =============================================================================
hdr("SECTION 8: Stage 6 — CCAR 2026 Scenario Projections")

# ── 8.1 Scenario paths ───────────────────────────────────────────────────────
# macro_base and macro_severe_aligned were loaded in Section 1.
# macro_severe_aligned has its macro_severe_* columns renamed to macro_base_*
# so the same forecast function works for both scenarios.

last_hist_yyyyqq <- max(agg_clean$yyyyqq)
init_rows        <- tail(agg_clean, P_LAG)

scenario_paths <- list(
  baseline         = macro_base[yyyyqq > last_hist_yyyyqq],
  severely_adverse = if (!is.null(macro_severe_aligned))
                       macro_severe_aligned[yyyyqq > last_hist_yyyyqq]
                     else NULL
)
scenario_paths <- Filter(function(x) !is.null(x) && nrow(x) > 0,
                          scenario_paths)
msg("CCAR scenarios available: %s", paste(names(scenario_paths), collapse=", "))

# ── 8.2 Multi-step forecast function ─────────────────────────────────────────
forecast_varx <- function(var_model, macro_scenario_path, init_Y, p) {
  # var_model         : vars::VAR object
  # macro_scenario_path: data.table with exogenous paths (FCST_HOR rows)
  # init_Y            : matrix (p × k) of initial Y values
  # p                 : lag order

  k     <- ncol(init_Y)
  h     <- nrow(macro_scenario_path)
  fcst  <- matrix(NA_real_, nrow=h, ncol=k,
                  dimnames=list(NULL, colnames(init_Y)))

  # Build companion-form coefficient matrices from var_model
  A_list <- lapply(1:p, function(lag) {
    sapply(var_model$varresult, function(eq) {
      cf    <- coef(eq)
      terms <- paste0(names(var_model$varresult), ".l", lag)
      cf[terms]
    })
  })
  # Intercept vector
  const <- sapply(var_model$varresult, function(eq) coef(eq)["const"])

  # Exogenous coefficient matrix (one column per exog variable, one row per equation)
  exog_nms   <- names(var_model$varresult[[1]]$model)
  exog_nms   <- exog_nms[!str_detect(exog_nms, paste0("^(", paste(colnames(init_Y), collapse="|"), ")")) &
                           exog_nms != "(Intercept)"]
  B_exog <- sapply(var_model$varresult, function(eq) {
    cf <- coef(eq)
    cf[intersect(exog_nms, names(cf))]
  })

  Y_hist <- init_Y  # last p rows for lag construction

  for (t in 1:h) {
    z_t <- unlist(macro_scenario_path[t, intersect(exog_nms, names(macro_scenario_path)), with=FALSE])
    y_hat <- const
    for (lag in 1:p) {
      idx   <- nrow(Y_hist) - lag + 1
      if (idx < 1) break
      y_hat <- y_hat + A_list[[lag]] %*% Y_hist[idx, ]
    }
    if (length(z_t) > 0 && nrow(B_exog) > 0)
      y_hat <- y_hat + B_exog %*% z_t[rownames(B_exog)]

    fcst[t, ]  <- y_hat
    Y_hist     <- rbind(Y_hist, y_hat)
  }
  fcst
}

# ── 8.3 Run projections for each scenario ────────────────────────────────────
y_cols_fcst <- intersect(Y_VARS, colnames(Y_chol))
init_Y      <- as.matrix(tail(agg_clean[, ..y_cols_fcst], P_LAG))

forecast_results <- list()

for (scen in names(scenario_paths)) {
  scen_path <- scenario_paths[[scen]]
  setorder(scen_path, yyyyqq)

  if (nrow(scen_path) == 0) {
    msg("  No forward path for scenario '%s' — skipping", scen); next
  }

  fcst_mat <- tryCatch(
    forecast_varx(varx_full, scen_path, init_Y, P_LAG),
    error=function(e) { warning(paste("Forecast failed for", scen, ":", e$message)); NULL }
  )

  if (!is.null(fcst_mat)) {
    fcst_dt <- as.data.table(fcst_mat)
    fcst_dt[, `:=`(yyyyqq=scen_path$yyyyqq[1:.N], scenario=scen)]
    forecast_results[[scen]] <- fcst_dt
    msg("  Scenario %-20s : %d quarters projected", scen, nrow(fcst_dt))
  }
}

all_forecasts <- rbindlist(forecast_results)

# Attach historical actuals
hist_long <- agg_clean[yyyyqq >= 200501L,
                         c("yyyyqq", Y_VARS), with=FALSE]
hist_long[, scenario := "historical"]

full_series <- rbindlist(list(hist_long, all_forecasts), fill=TRUE)
saveRDS(full_series, "Results/04_forecast_paths.rds")
msg("Forecast paths → Results/04_forecast_paths.rds")

# ── 8.4 Fan chart plots ───────────────────────────────────────────────────────
scenario_colors <- c(
  "historical"      = "#333333",
  "baseline"        = "#1565C0",
  "adverse"         = "#F57C00",
  "severely_adverse"= "#C62828"
)

plot_fan <- function(dep_var, dep_label) {
  plot_dt <- full_series[!is.na(get(dep_var)),
                          .(yyyyqq, scenario, value=get(dep_var))]

  # Convert yyyyqq to approximate date
  plot_dt[, yr  := yyyyqq %/% 100L]
  plot_dt[, qtr := yyyyqq %% 100L]
  plot_dt[, date_num := yr + (qtr - 1) / 4]

  ggplot(plot_dt, aes(x=date_num, y=value,
                       color=scenario, linetype=scenario)) +
    geom_vline(xintercept=last_hist_yyyyqq %/% 100 +
                 (last_hist_yyyyqq %% 100 - 1) / 4,
               linetype="dotted", color="#666", linewidth=0.5) +
    geom_line(linewidth=0.8, na.rm=TRUE) +
    annotate("text", x=2026.25, y=Inf, label="CCAR 2026\nProjection",
             vjust=1.5, hjust=0, size=3, color="#666") +
    scale_color_manual(values=scenario_colors, name="Scenario") +
    scale_linetype_manual(
      values=c(historical="solid", baseline="solid",
               adverse="dashed", severely_adverse="dotdash"),
      name="Scenario") +
    scale_x_continuous(breaks=seq(2006, 2029, 2),
                        labels=function(x) paste0("'", substr(x, 3, 4))) +
    labs(
      title    = dep_label,
      subtitle = "Historical + CCAR 2026 scenario projections",
      x=NULL, y=dep_label
    ) +
    theme_minimal(base_size=9) +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size=8),
      plot.title       = element_text(face="bold", size=10),
      panel.grid.minor = element_blank()
    )
}

fan_plots <- Filter(Negate(is.null), list(
  plot_fan("pll_rate",             "PLL Rate (% of Avg Loans)"),
  plot_fan("dq_rate",              "Delinquency Rate (%)"),
  plot_fan("netintmrg",            "Net Interest Margin (%)"),
  plot_fan("insured_share_growth", "Insured Share Growth (YoY %)"),
  plot_fan("member_growth_yoy",    "Membership Growth (YoY %)"),
  plot_fan("costfds",              "Cost of Funds (%)"),
  plot_fan("loan_to_share",        "Loan-to-Share Ratio"),
  plot_fan("pcanetworth",          "Net Worth Ratio (%)")
))

p_fan_combined <- wrap_plots(fan_plots, ncol=2) +
  plot_annotation(
    title    = "CU System — CCAR 2026 Scenario Fan Charts",
    subtitle = sprintf("VARX(p=%d) | Cholesky-identified | Full sample 2005Q1–2025Q4", P_LAG),
    theme    = theme(plot.title=element_text(face="bold", size=13))
  )

ggsave("Figures/04_fan_charts_all.png", p_fan_combined,
       width=14, height=10, dpi=150)
msg("Fan charts → Figures/04_fan_charts_all.png")

# =============================================================================
# 9. HEADLINE COEFFICIENT EXTRACT — The Three-Way Interaction
# =============================================================================
hdr("SECTION 9: Headline Finding — Triple Interaction Extract")

# Pull post_shale × yoy_oil × direct_exposure from panel FE models
# This is the paper's headline — only meaningful in post-shale + hiking regime

triple_coefs <- coef_table[term %like% "post_x_oil_x_direct",
                             .(dep_var, term, estimate, se, p)]
setorder(triple_coefs, dep_var)

cat("\n  Triple Interaction: post_shale × yoy_oil × direct_exposure\n")
cat("  (Interpretation: oil price rise impact on CUs with direct energy exposure\n")
cat("   in post-shale era vs. pre-shale era)\n\n")

triple_coefs[, sig := fcase(
  p < 0.01, "***",
  p < 0.05, "**",
  p < 0.10, "*",
  default  , ""
)]

print(triple_coefs[, .(dep_var, estimate=round(estimate,5),
                         se=round(se,5), p=round(p,4), sig)])

# Also extract fomc × oil (two-way — rate channel activation)
fomc_oil_coefs <- coef_table[term %like% "fomc_x_brent",
                               .(dep_var, term, estimate, se, p)]
fomc_oil_coefs[, sig := fcase(p<0.01,"***", p<0.05,"**", p<0.10,"*", default,"")]

cat("\n  Two-Way Interaction: fomc_regime × yoy_oil (rate channel)\n\n")
print(fomc_oil_coefs[, .(dep_var, estimate=round(estimate,5),
                           se=round(se,5), p=round(p,4), sig)])

fwrite(rbindlist(list(triple_coefs[, type:="triple"],
                      fomc_oil_coefs[, type:="fomc_oil"]), fill=TRUE),
       "Results/04_headline_coefs.csv")
msg("\nHeadline coefficients → Results/04_headline_coefs.csv")

# =============================================================================
# 10. SESSION & OUTPUT MANIFEST
# =============================================================================
hdr("SECTION 10: Output Manifest")

outputs <- list(
  Models  = list.files("Models",  pattern="\\.rds$", full.names=TRUE),
  Results = list.files("Results", pattern="04_",     full.names=TRUE),
  Figures = list.files("Figures", pattern="04_",     full.names=TRUE)
)

for (grp in names(outputs)) {
  cat(sprintf("\n  %s:\n", grp))
  for (f in outputs[[grp]])
    cat(sprintf("    %-55s  [%s]\n", basename(f),
                format(file.size(f), big.mark=",")))
}

cat("\n", strrep("=", 70), "\n")
cat("  Script 04b complete —", ts(), "\n")
cat("  Next: Script 05 — Results tables & publication charts\n")
cat(strrep("=", 70), "\n")
