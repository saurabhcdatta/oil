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
  library(fixest)        # feols() — panel FE with clustering
  library(zoo)           # na.approx for NA interpolation
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
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[1])) a else b

# ── fill_nas: interpolate NAs in all specified columns of a data.table ────────
# Used wherever a matrix must be NA-free before passing to vars::VAR/VARselect.
# Applies zoo::na.approx (linear) then fills any residual NAs with column mean.
fill_nas <- function(dt, vars) {
  dt <- copy(dt)
  for (v in intersect(vars, names(dt))) {
    n_na <- sum(is.na(dt[[v]]))
    if (n_na > 0) {
      dt[, (v) := zoo::na.approx(get(v), na.rm=FALSE)]
      still_na <- sum(is.na(dt[[v]]))
      if (still_na > 0)
        dt[is.na(get(v)), (v) := mean(dt[[v]], na.rm=TRUE)]
    }
  }
  dt
}

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
# vars::VAR and VARselect require ZERO NAs in Y.
# Strategy:
#   Step 1 — Start with rows where CORE non-sparse vars are all present
#            (excludes pcanetworth which has 51 NAs in recent quarters)
#   Step 2 — Fill remaining NAs in every Y variable via linear interpolation
#            (handles pll_rate 1 NA, insured_share_growth 4, member_growth_yoy 4)
#   Step 3 — Drop any row still NA after interpolation (tail end effect)
#   Step 4 — Verify zero NAs before passing to vars functions

CORE_Y_VARS <- intersect(
  c("dq_rate","netintmrg","costfds","loan_to_share"),  # densest cols only
  names(agg))

agg_clean <- agg[complete.cases(agg[, ..CORE_Y_VARS])]
setorder(agg_clean, yyyyqq)

# Interpolate ALL Y variables — handles any remaining NAs
agg_clean <- fill_nas(agg_clean, Y_VARS)

# Also fill any NAs in X (exog) columns
exog_cols_raw <- intersect(c(Z_VARS, INTERACT_VARS, DUMMY_VARS), names(agg_clean))
agg_clean <- fill_nas(agg_clean, exog_cols_raw)

# Final complete-case drop (should remove nothing if interpolation worked)
y_na_check  <- !complete.cases(agg_clean[, ..Y_VARS])
agg_clean   <- agg_clean[!y_na_check]
msg("  Rows dropped after final NA check: %d", sum(y_na_check))

# Build Y and X matrices — guaranteed zero NA
Y_mat <- as.matrix(agg_clean[, ..Y_VARS])
rownames(Y_mat) <- as.character(agg_clean$yyyyqq)

exog_cols <- intersect(c(Z_VARS, INTERACT_VARS, DUMMY_VARS), names(agg_clean))
X_mat <- as.matrix(agg_clean[, ..exog_cols])

msg("Y matrix : %d \u00d7 %d | NAs: %d", nrow(Y_mat), ncol(Y_mat), sum(is.na(Y_mat)))
msg("X matrix : %d \u00d7 %d | NAs: %d", nrow(X_mat), ncol(X_mat), sum(is.na(X_mat)))

if (sum(is.na(Y_mat)) > 0)
  stop("Y_mat still contains NAs after interpolation — check agg construction")

# ── 3.2 Lag order selection ───────────────────────────────────────────────────
# VARselect requires nrow(Y) == nrow(exogen) — pass full matrices, no pre-trim.
# vars handles alignment internally for each candidate lag order.
LAG_MAX_SELECT <- min(6L, floor(nrow(Y_mat) / (ncol(Y_mat) + ncol(X_mat) + 1L)))
LAG_MAX_SELECT <- max(LAG_MAX_SELECT, 2L)   # always test at least p=1,2
msg("VARselect lag.max = %d (capped by obs/params ratio)", LAG_MAX_SELECT)

varselect_full <- tryCatch(
  VARselect(Y_mat, lag.max=LAG_MAX_SELECT, type="const", exogen=X_mat),
  error=function(e) {
    msg("VARselect with exogen failed: %s", e$message)
    msg("Retrying without exogen (IC only — exog still used in VAR())")
    tryCatch(
      VARselect(Y_mat, lag.max=LAG_MAX_SELECT, type="const"),
      error=function(e2) {
        msg("VARselect failed entirely: %s — using default P_LAG=%d", e2$message, P_LAG)
        NULL
      }
    )
  }
)

# Guard against NULL varselect (both attempts failed) — use configured P_LAG
if (!is.null(varselect_full) && !is.null(varselect_full$criteria)) {
  ic_table <- tryCatch(
    rbind(
      AIC = varselect_full$criteria["AIC(n)", ],
      BIC = varselect_full$criteria["SC(n)",  ],
      HQ  = varselect_full$criteria["HQ(n)",  ]
    ),
    error=function(e) NULL
  )
  if (!is.null(ic_table) && is.numeric(ic_table)) {
    cat("\n  Information Criteria by Lag Order:\n")
    print(round(ic_table, 4))
  }
  p_aic <- tryCatch(as.integer(varselect_full$selection["AIC(n)"]), error=function(e) P_LAG)
  p_bic <- tryCatch(as.integer(varselect_full$selection["SC(n)"]),  error=function(e) P_LAG)
  msg("AIC selects p=%d | BIC selects p=%d | Using configured p=%d", p_aic, p_bic, P_LAG)
  if (P_LAG != p_aic)
    warning(sprintf("Configured P_LAG=%d differs from AIC p=%d. Review.", P_LAG, p_aic))
} else {
  msg("VARselect unavailable — proceeding with configured P_LAG=%d", P_LAG)
}

# ── 3.3 Estimate VARX ────────────────────────────────────────────────────────
# vars::VAR with exogen argument implements VARX directly
# type="const" adds intercept; seasonal=0 (quarterly dummies already in X_mat)

varx_full <- VAR(
  y      = Y_mat,
  p      = P_LAG,
  type   = "const",
  exogen = X_mat   # vars::VAR handles lag alignment internally — do NOT pre-trim
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
# Use fixest::feols (faster, more robust than lfe::felm for large panels)
# RHS: own lags (p) + cross-lags (p) + exog Z + interactions + dummies
# CU FE + Quarter FE absorbed; SE clustered two-way at CU + quarter level
#
# NOTE: macro_base_yoy_oil is COLLINEAR with yyyyqq FE (same value all CUs).
# Interaction terms (oil_x_brent, fomc_x_brent) vary cross-sectionally and ARE
# identified. Pure time-series exog vars will be dropped by feols — expected.

panel_exog_cols <- intersect(
  c(Z_VARS, INTERACT_VARS, DUMMY_VARS), names(panel))

# Remove pure time-series macro vars that will be absorbed by quarter FE
# Keep only cross-sectionally varying interaction terms and dummies
panel_exog_id <- intersect(
  c("oil_x_brent","spillover_x_brent","fomc_x_brent",
    "post_x_oil","post_x_oil_x_direct","zirp_x_oil",
    "post_shale","gfc_dummy","covid_dummy","zirp_era","hike_cycle"),
  names(panel))

fe_models <- list()
fe_coefs  <- list()

for (dep in Y_VARS) {
  own_lags   <- paste0(dep, "_lag", 1:P_LAG)
  cross_lags <- setdiff(LAG_Y_VARS, own_lags)

  rhs_terms  <- c(
    intersect(own_lags,    names(panel)),
    intersect(cross_lags,  names(panel)),
    intersect(panel_exog_id, names(panel))
  )

  if (length(rhs_terms) == 0) next

  rhs    <- paste(rhs_terms, collapse=" + ")
  frml   <- as.formula(paste0(dep, " ~ ", rhs, " | join_number + yyyyqq"))

  keep_cols <- unique(c(dep, rhs_terms, "join_number","yyyyqq"))
  keep_cols <- intersect(keep_cols, names(panel))
  panel_eq  <- panel[complete.cases(panel[, ..keep_cols])]

  if (nrow(panel_eq) < 1000) {
    msg("  SKIP %s: only %d complete obs", dep, nrow(panel_eq)); next
  }

  fe_models[[dep]] <- tryCatch(
    fixest::feols(frml, data=panel_eq,
                  cluster=~join_number, warn=FALSE, notes=FALSE),
    error=function(e) { warning(paste("feols failed for", dep, ":", e$message)); NULL }
  )

  if (!is.null(fe_models[[dep]])) {
    ct <- tryCatch(fixest::coeftable(fe_models[[dep]]), error=function(e) NULL)
    if (!is.null(ct)) {
      cf <- as.data.table(ct, keep.rownames=TRUE)
      setnames(cf, c("term","estimate","se","t","p"))
      cf[, dep_var := dep]
      fe_coefs[[dep]] <- cf
      r2v <- tryCatch(as.numeric(fixest::r2(fe_models[[dep]], "wr2"))[[1L]],
                      error=function(e) NA_real_)
      msg("  %-30s : within-R\u00b2=%.4f | N=%s",
          dep, r2v %||% NA_real_,
          format(nobs(fe_models[[dep]]), big.mark=","))
    }
  }
}

# Combine coefficient table — guard against empty list
fe_coefs_valid <- Filter(function(x) !is.null(x) && nrow(x) > 0, fe_coefs)
if (length(fe_coefs_valid) > 0) {
  coef_table <- rbindlist(fe_coefs_valid, fill=TRUE)
} else {
  msg("  WARNING: no panel FE equations estimated — coef_table empty")
  coef_table <- data.table(dep_var=character(), term=character(),
                            estimate=numeric(), se=numeric(), t=numeric(), p=numeric())
}
fwrite(coef_table, "Results/04_varx_coef_panel_fe.csv")
msg("\n  Coefficient table \u2192 Results/04_varx_coef_panel_fe.csv (%d rows)", nrow(coef_table))

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
  agg_t <- merge(agg_t, panel_dum, by="yyyyqq", all.x=TRUE, allow.cartesian=FALSE)
  tier_agg[[as.character(tier)]] <- agg_t

  # Fill NAs via interpolation before building Y/X matrices
  y_cols_tier  <- intersect(Y_VARS, names(agg_t))
  exog_tier    <- intersect(c(avail_z, avail_int_panel, avail_dum), names(agg_t))
  agg_t_clean  <- fill_nas(agg_t, c(y_cols_tier, exog_tier))
  agg_t_clean  <- agg_t_clean[complete.cases(agg_t_clean[, ..y_cols_tier])]

  if (nrow(agg_t_clean) < 30) {
    msg("Tier %s: insufficient obs (%d) — skipping", tier, nrow(agg_t_clean))
    next
  }

  y_cols_t  <- intersect(Y_VARS, names(agg_t_clean))
  Y_t       <- as.matrix(agg_t_clean[, ..y_cols_t])
  exog_t_cols <- intersect(c(avail_z, avail_int_panel, avail_dum), names(agg_t_clean))
  if (length(exog_t_cols) > 0) {
    X_t <- as.matrix(agg_t_clean[, ..exog_t_cols])
  } else {
    X_t <- matrix(0, nrow=nrow(agg_t_clean), ncol=1)
  }

  tier_models[[as.character(tier)]] <- tryCatch(
    VAR(Y_t, p=P_LAG, type="const",
        exogen=X_t),   # vars::VAR handles lag alignment internally
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

estimate_subsample <- function(data_full, macro_data, label,
                               yyyyqq_from, yyyyqq_to=Inf) {
  # Filter by yyyyqq range — avoids logical vector length mismatch
  sub <- data_full[yyyyqq >= yyyyqq_from & yyyyqq <= yyyyqq_to]
  mac <- macro_data[yyyyqq >= yyyyqq_from & yyyyqq <= yyyyqq_to]

  if (nrow(sub) == 0) {
    msg("  %s: no rows in range — skipping", label); return(NULL)
  }

  sub <- merge(sub, mac, by="yyyyqq", all.x=TRUE)

  # Interpolate NAs in Y and X before building matrices
  all_exog <- intersect(c(Z_VARS, INTERACT_VARS, DUMMY_VARS), names(sub))
  sub <- fill_nas(sub, c(Y_VARS, all_exog))

  sub <- sub[complete.cases(sub[, intersect(Y_VARS, names(sub)), with=FALSE])]

  if (nrow(sub) < P_LAG + ncol(Y_mat) + 5) {
    msg("  %s: too few obs (%d) — skipping", label, nrow(sub)); return(NULL)
  }

  y_cols_s    <- intersect(Y_VARS, names(sub))
  Y_s         <- as.matrix(sub[, ..y_cols_s])
  exog_s_cols <- setdiff(
    intersect(c(Z_VARS, INTERACT_VARS, DUMMY_VARS), names(sub)),
    "post_shale"   # collinear within a single era subsample
  )
  X_s <- if (length(exog_s_cols) > 0) as.matrix(sub[, ..exog_s_cols]) else matrix(rep(1, nrow(sub)), ncol=1)

  m <- tryCatch(
    VAR(Y_s, p=P_LAG, type="const", exogen=X_s),
    error=function(e) { warning(paste(label, "failed:", e$message)); NULL }
  )

  if (!is.null(m)) {
    r2s <- sapply(m$varresult, function(eq) summary(eq)$adj.r.squared)
    msg("  %-15s | Obs=%-5d | Avg Adj-R\u00b2=%.3f | Max root=%.4f",
        label, nrow(Y_s), mean(r2s, na.rm=TRUE), max(roots(m)))
  }
  m
}

varx_pre  <- estimate_subsample(agg, macro_ts, "Pre-2015",
                                 yyyyqq_from=200501L, yyyyqq_to=201404L)
varx_post <- estimate_subsample(agg, macro_ts, "Post-2015",
                                 yyyyqq_from=201501L)

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
                  exogen=X_mat)   # vars::VAR handles lag alignment internally

# ── 7.2 Compute IRFs ─────────────────────────────────────────────────────────
# Shock 1: pll_rate — the forward-looking credit quality shock
#   Captures: when management raises provisions, how does the rest of the
#   CU system respond over 12 quarters? Expect dq_rate to follow 1-2Q later.
# Shock 2: costfds — the rate channel shock
#   Captures: when funding costs rise (Fed hike transmission), how does NIM,
#   loan demand, membership, and deposit mix respond?

# ── Compute ALL 8 IRFs — one per endogenous variable in Cholesky order ───────
# Each shows: "when THIS variable is shocked, how does the WHOLE system respond?"
# This gives the complete picture of dynamic transmission in the CU system.

safe_irf <- function(model, impulse_var, response_vars, horizon, boot, ci, runs) {
  if (!impulse_var %in% colnames(model$y)) {
    msg("  SKIP IRF: impulse '%s' not in model", impulse_var)
    return(NULL)
  }
  tryCatch(
    irf(model, impulse=impulse_var, response=response_vars,
        n.ahead=horizon, boot=boot, ci=ci, runs=runs, ortho=TRUE),
    error=function(e) {
      msg("  IRF failed for impulse '%s': %s", impulse_var, e$message)
      NULL
    }
  )
}

irf_pll    <- safe_irf(varx_chol, "pll_rate",             chol_order, IRF_HORIZON, TRUE, 0.90, N_BOOT)
irf_delq   <- safe_irf(varx_chol, "dq_rate",              chol_order, IRF_HORIZON, TRUE, 0.90, N_BOOT)
irf_cof    <- safe_irf(varx_chol, "costfds",              chol_order, IRF_HORIZON, TRUE, 0.90, N_BOOT)
irf_nim    <- safe_irf(varx_chol, "netintmrg",            chol_order, IRF_HORIZON, TRUE, 0.90, N_BOOT)
irf_dep    <- safe_irf(varx_chol, "insured_share_growth", chol_order, IRF_HORIZON, TRUE, 0.90, N_BOOT)
irf_mem    <- safe_irf(varx_chol, "member_growth_yoy",    chol_order, IRF_HORIZON, TRUE, 0.90, N_BOOT)
irf_lts    <- safe_irf(varx_chol, "loan_to_share",        chol_order, IRF_HORIZON, TRUE, 0.90, N_BOOT)
irf_nw     <- safe_irf(varx_chol, "pcanetworth",          chol_order, IRF_HORIZON, TRUE, 0.90, N_BOOT)

saveRDS(
  list(irf_pll  = irf_pll,  irf_delq = irf_delq, irf_cof = irf_cof,
       irf_nim  = irf_nim,  irf_dep  = irf_dep,  irf_mem = irf_mem,
       irf_lts  = irf_lts,  irf_nw   = irf_nw),
  "Results/04_varx_irf.rds"
)
msg("All 8 IRFs saved \u2192 Results/04_varx_irf.rds")

# ── 7.3 Plot IRFs ─────────────────────────────────────────────────────────────
plot_irf <- function(irf_obj, impulse_label, outfile) {

  # ── Extract IRF matrices into tidy data.table ─────────────────────────────
  safe_vec <- function(x) {
    if (is.matrix(x)) as.numeric(x[, 1L]) else as.numeric(x)
  }

  extract_one <- function(series_nm) {
    irf_vals <- irf_obj$irf[[series_nm]]
    lo_vals  <- irf_obj$Lower[[series_nm]]
    up_vals  <- irf_obj$Upper[[series_nm]]
    if (is.null(irf_vals)) return(NULL)
    n_rows <- if (is.matrix(irf_vals)) nrow(irf_vals) else length(irf_vals)
    data.table(
      horizon  = seq(0L, n_rows - 1L),
      irf      = safe_vec(irf_vals),
      lower    = if (!is.null(lo_vals)) safe_vec(lo_vals) else NA_real_,
      upper    = if (!is.null(up_vals)) safe_vec(up_vals) else NA_real_,
      response = series_nm
    )
  }

  all_irf <- rbindlist(
    Filter(Negate(is.null), lapply(names(irf_obj$irf), extract_one))
  )
  if (nrow(all_irf) == 0) {
    msg("  plot_irf: no data extracted for %s — skipping", impulse_label)
    return(invisible(NULL))
  }

  # ── Statistical significance classification ───────────────────────────────
  # The 90% bootstrap CI is the standard significance test for IRFs.
  # Interpretation:
  #   sig_pos: lower > 0  → CI entirely above zero → significantly POSITIVE
  #   sig_neg: upper < 0  → CI entirely below zero → significantly NEGATIVE
  #   insig  : CI crosses zero → effect not distinguishable from zero
  #
  # This is equivalent to a two-sided test at α = 0.10.
  # We visualise this in three ways:
  #   1. Shaded ribbon: green (pos) / red (neg) / grey (insig)
  #   2. Significance bar: coloured tick strip along x-axis bottom
  #   3. Annotation: first and last significant quarter per panel

  has_ci <- all_irf[, any(!is.na(lower) & !is.na(upper))]

  if (has_ci) {
    all_irf[, sig := fcase(
      !is.na(lower) & !is.na(upper) & lower > 0,  "pos",   # sig positive
      !is.na(lower) & !is.na(upper) & upper < 0,  "neg",   # sig negative
      !is.na(lower) & !is.na(upper),              "insig", # CI crosses 0
      default = "no_ci"
    )]

    # Ribbon fill colour per horizon point — use sig classification
    RIBBON_COLS <- c(
      pos   = "#1B7837",   # green  — significantly positive
      neg   = "#C62828",   # red    — significantly negative
      insig = "#90A4AE",   # grey   — not significant
      no_ci = "#BBDEFB"    # pale blue — no CI available
    )
    RIBBON_ALPHA <- c(pos=0.22, neg=0.22, insig=0.10, no_ci=0.10)

    # Significance bar data: strip along bottom of each panel
    # y position is dynamically placed at panel minimum
    sig_bar <- all_irf[sig != "no_ci", .(
      horizon, sig, response,
      response_lbl = str_replace_all(response, "_", " ") |> str_to_title()
    )]
    sig_bar[, bar_col := fcase(
      sig == "pos",   "#1B7837",
      sig == "neg",   "#C62828",
      sig == "insig", "#CFD8DC",
      default         = "#CFD8DC"
    )]

    # First / last significant quarter per panel (for annotation)
    sig_windows <- all_irf[sig %in% c("pos","neg"),
      .(first_sig = min(horizon), last_sig = max(horizon),
        direction = sig[which.max(abs(irf))]),
      by=response]
  } else {
    all_irf[, sig := "no_ci"]
  }

  # Clean labels
  all_irf[, response_lbl := str_replace_all(response, "_", " ") |>
                              str_to_title()]

  # ── Build plot ────────────────────────────────────────────────────────────
  p <- ggplot(all_irf, aes(x=horizon))

  # Layer 1: CI ribbon — split by significance colour
  # inherit.aes=FALSE means parent x=horizon is dropped → must supply x explicitly
  if (has_ci) {
    for (sig_type in c("insig","pos","neg")) {
      dt_sub <- all_irf[sig == sig_type & !is.na(lower) & !is.na(upper)]
      if (nrow(dt_sub) > 0) {
        p <- p + geom_ribbon(
          data    = dt_sub,
          aes(x=horizon, ymin=lower, ymax=upper),   # x required when inherit.aes=FALSE
          fill    = RIBBON_COLS[sig_type],
          alpha   = RIBBON_ALPHA[sig_type],
          inherit.aes = FALSE
        )
      }
    }
  }

  # Layer 2: Zero reference line
  p <- p + geom_hline(yintercept=0, linetype="dashed",
                       colour="#555555", linewidth=0.4)

  # Layer 3: IRF point estimate — thick+coloured where significant, thin+grey where not
  # Must supply x=horizon explicitly for data= subsets (same reason as ribbon)
  if (has_ci) {
    p <- p +
      geom_line(data=all_irf[sig == "insig"],
                aes(x=horizon, y=irf, group=response),
                colour="#607D8B", linewidth=0.6, inherit.aes=FALSE) +
      geom_line(data=all_irf[sig == "pos"],
                aes(x=horizon, y=irf, group=response),
                colour="#1B7837", linewidth=1.0, inherit.aes=FALSE) +
      geom_line(data=all_irf[sig == "neg"],
                aes(x=horizon, y=irf, group=response),
                colour="#C62828", linewidth=1.0, inherit.aes=FALSE) +
      geom_line(data=all_irf[sig == "no_ci"],
                aes(x=horizon, y=irf, group=response),
                colour="#1565C0", linewidth=0.8, inherit.aes=FALSE)
  } else {
    p <- p + geom_line(aes(y=irf), colour="#1565C0", linewidth=0.8)
  }

  # Layer 4: Significance bar at bottom of each panel via geom_rug
  if (has_ci && nrow(sig_bar) > 0) {
    p <- p + geom_rug(
      data        = sig_bar,
      aes(x=horizon, colour=sig),
      sides       = "b",
      linewidth   = 1.2,
      length      = unit(0.04, "npc"),
      inherit.aes = FALSE
    ) +
    scale_colour_manual(
      values = c(pos="darkgreen", neg="darkred", insig="grey60"),
      labels = c(pos="Sig. positive (90% CI)",
                 neg="Sig. negative (90% CI)",
                 insig="Not significant"),
      name   = NULL,
      guide  = guide_legend(override.aes=list(linewidth=2))
    )
  }

  # Layer 5: Significance window label  (label.size → linewidth from ggplot2 3.5.0)
  if (has_ci && exists("sig_windows") && nrow(sig_windows) > 0) {
    sig_windows[, response_lbl := str_replace_all(response, "_", " ") |>
                                   str_to_title()]
    sig_windows[, label    := sprintf("Q%d\u2013Q%d", first_sig, last_sig)]
    sig_windows[, ann_col  := fifelse(direction=="pos", "#1B7837", "#C62828")]
    p <- p + geom_label(
      data        = sig_windows,
      aes(x=first_sig, y=Inf, label=label),
      colour      = sig_windows$ann_col,
      fill        = "white",
      size        = 2.5,
      fontface    = "bold",
      hjust       = 0,
      vjust       = 1.4,
      linewidth   = 0.2,          # replaces deprecated label.size
      inherit.aes = FALSE
    )
  }

  # Facet + scales + labels
  p <- p +
    facet_wrap(~response_lbl, scales="free_y", ncol=3) +
    scale_x_continuous(
      breaks = seq(0, max(all_irf$horizon), 4),
      labels = function(x) paste0("Q", x)
    ) +
    labs(
      title    = sprintf("Impulse Response: %s shock \u2192 CU system",
                          impulse_label),
      subtitle = sprintf(
        paste("Cholesky-identified | 90%%%% bootstrap CI (%d runs) | p=%d",
              "| Green = sig. positive | Red = sig. negative | Grey = not sig."),
        N_BOOT, P_LAG
      ),
      caption  = paste(
        "Significance: 90%% CI band does not cross zero.",
        "Coloured tick marks at bottom = significant quarters.",
        "QX\u2013QY label = window of statistical significance."
      ),
      x = "Quarters after shock",
      y = "Response (level units of dependent var)"
    ) +
    theme_minimal(base_size=10) +
    theme(
      plot.title       = element_text(face="bold", size=11),
      plot.subtitle    = element_text(size=8, colour="#444444",
                                       lineheight=1.3),
      plot.caption     = element_text(size=7.5, colour="#888888", hjust=0),
      strip.text       = element_text(face="bold", size=9),
      strip.background = element_rect(fill="#f5f5f5", colour="#cccccc"),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour="#dddddd", fill=NA,
                                       linewidth=0.3),
      legend.position  = "bottom",
      legend.text      = element_text(size=8)
    )

  ggsave(outfile, p, width=13, height=8, dpi=200, bg="white")
  msg("  IRF chart \u2192 %s", outfile)
  invisible(p)
}

plot_irf(irf_pll,  "PLL Rate (Forward-Looking Provisions)",    "Figures/04_irf_pll_rate.png")
plot_irf(irf_delq, "Delinquency Rate (Realised Credit Stress)", "Figures/04_irf_delinquency.png")
plot_irf(irf_cof,  "Cost of Funds (Rate Channel)",             "Figures/04_irf_cost_of_funds.png")
plot_irf(irf_nim,  "Net Interest Margin (NIM Compression)",    "Figures/04_irf_net_int_margin.png")
plot_irf(irf_dep,  "Insured Share Growth (Deposit Channel)",   "Figures/04_irf_deposit_growth.png")
plot_irf(irf_mem,  "Membership Growth (Early Warning Signal)", "Figures/04_irf_membership.png")
plot_irf(irf_lts,  "Loan-to-Share Ratio (Balance Sheet)",      "Figures/04_irf_loan_to_share.png")
plot_irf(irf_nw,   "Net Worth Ratio (Capital Adequacy)",       "Figures/04_irf_net_worth.png")

msg("All 8 IRF charts saved to Figures/")

# ── Combined 2×4 IRF summary grid for paper ──────────────────────────────────
# Shows headline response (dq_rate) to each of the 8 structural shocks
# Answers: "Which shock matters most for credit quality?"
hdr("IRF Summary Grid: dq_rate response to all 8 shocks")

irf_list_all <- list(
  list(obj=irf_pll,  label="PLL Rate shock"),
  list(obj=irf_delq, label="Delinquency shock"),
  list(obj=irf_cof,  label="Cost of Funds shock"),
  list(obj=irf_nim,  label="NIM shock"),
  list(obj=irf_dep,  label="Deposit Growth shock"),
  list(obj=irf_mem,  label="Membership shock"),
  list(obj=irf_lts,  label="Loan-to-Share shock"),
  list(obj=irf_nw,   label="Net Worth shock")
)

# For each shock, extract dq_rate response
dq_resp_list <- lapply(irf_list_all, function(x) {
  if (is.null(x$obj) || is.null(x$obj$irf[["dq_rate"]])) return(NULL)
  vals <- x$obj$irf[["dq_rate"]]
  lo   <- x$obj$Lower[["dq_rate"]]
  up   <- x$obj$Upper[["dq_rate"]]
  n    <- if (is.matrix(vals)) nrow(vals) else length(vals)
  data.table(
    horizon  = seq(0L, n-1L),
    irf      = if (is.matrix(vals)) as.numeric(vals[,1]) else as.numeric(vals),
    lower    = if (!is.null(lo) && is.matrix(lo)) as.numeric(lo[,1]) else NA_real_,
    upper    = if (!is.null(up) && is.matrix(up)) as.numeric(up[,1]) else NA_real_,
    impulse  = x$label
  )
})
dq_resp_dt <- rbindlist(Filter(Negate(is.null), dq_resp_list))

if (nrow(dq_resp_dt) > 0) {
  has_ci_grid <- dq_resp_dt[, any(!is.na(lower) & !is.na(upper))]
  p_grid <- ggplot(dq_resp_dt, aes(x=horizon)) +
    { if (has_ci_grid)
        geom_ribbon(aes(ymin=lower, ymax=upper), fill="#d73027", alpha=0.12)
      else list() } +
    geom_line(aes(y=irf), color="#d73027", linewidth=0.8) +
    geom_hline(yintercept=0, linetype="dashed", color="#888", linewidth=0.35) +
    facet_wrap(~impulse, scales="free_y", ncol=4) +
    scale_x_continuous(breaks=seq(0, IRF_HORIZON, 4),
                       labels=function(x) paste0("Q", x)) +
    labs(
      title    = "IRF Summary — Delinquency Rate Response to Each Structural Shock",
      subtitle = paste("Which shock matters most for credit quality?",
                       "| Cholesky-identified | 90% bootstrap CI"),
      caption  = paste("VARX(p=2) | Full sample 2005Q1-2025Q4",
                       "| Red shading = 90% CI band"),
      x="Quarters after shock", y="Response of Delinquency Rate"
    ) +
    theme_minimal(base_size=10) +
    theme(
      plot.title       = element_text(face="bold", size=11),
      strip.text       = element_text(face="bold", size=8.5),
      panel.grid.minor = element_blank()
    )
  ggsave("Figures/04_irf_dq_summary_grid.png", p_grid,
         width=16, height=8, dpi=150, bg="white")
  msg("  Combined IRF grid \u2192 Figures/04_irf_dq_summary_grid.png")
}

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
# 7.5 OIL PRICE SHOCK SIMULATION — The Key Chart for the Paper
# =============================================================================
# oil (macro_base_yoy_oil) is EXOGENOUS — not in the Y endogenous block.
# Standard Cholesky IRF cannot be computed from an exogenous variable.
#
# Method: Conditional forecast simulation
#   Baseline path: oil YoY = 0% every quarter (flat/neutral)
#   Shock path:    oil YoY = +30% at Q1 only (one-time shock), then 0%
#                  oil YoY = +30% sustained for 4 quarters
#                  oil YoY = -30% at Q1 only (oil bust shock)
#
# The DIFFERENCE between the shocked path and the baseline path gives the
# dynamic causal effect of the oil shock on each CU outcome — this is the
# VARX equivalent of an IRF for an exogenous variable.
#
# +30% YoY is approximately 1 standard deviation of historical PBRENT YoY.
# Results interpretable as: "effect of a +1 SD oil price increase"
# =============================================================================
hdr("SECTION 7.5: Oil Shock Simulation — CU Response to Oil Price Change")

# ── Guard: Section 7.5 requires varx_full from Section 3 ─────────────────────
if (!exists("varx_full") || is.null(varx_full)) {
  msg("WARNING: varx_full not found — attempting to load from Models/varx_full.rds")
  varx_full <- tryCatch(readRDS("Models/varx_full.rds"),
                         error=function(e) {
                           msg("ERROR: could not load varx_full: %s", e$message)
                           msg("Section 7.5 skipped — re-run Section 3 first")
                           NULL
                         })
}

if (is.null(varx_full)) {
  msg("Section 7.5 skipped — varx_full unavailable")
} else {

# ── Build z_var_names_sim: exogenous variable names as they appear in varx_full
# These are the column names of the exogenous block in the estimated VAR,
# excluding endogenous lags and the intercept.
# Must be built from varx_full (not Z_VARS directly) because some vars may
# have been dropped during estimation due to collinearity or missing data.
z_var_names_sim <- tryCatch({
  # Get all term names from the first equation of varx_full
  all_terms <- names(varx_full$varresult[[1]]$model)
  # Remove endogenous lag terms (pattern: varname.lN) and intercept
  lag_pattern <- paste0("^(", paste(colnames(Y_mat), collapse="|"), ")\\.l")
  all_terms[!grepl(lag_pattern, all_terms) & all_terms != "(Intercept)"]
}, error = function(e) {
  msg("Could not extract z_var_names_sim from varx_full: %s", e$message)
  msg("Falling back to Z_VARS + INTERACT_VARS + DUMMY_VARS")
  intersect(c(Z_VARS, INTERACT_VARS, DUMMY_VARS), names(agg_clean))
})

msg("VARX exogenous vars (%d): %s", length(z_var_names_sim),
    paste(z_var_names_sim, collapse=", "))
# =============================================================================
# FULL MACRO SIMULATION — All Z variables move together
# =============================================================================
# The previous version held all macro vars at historical means and only moved
# oil. This was wrong — oil price rises cause:
#   unemployment to FALL  (energy sector hiring, multiplier effects)
#   CPI to RISE           (energy costs pass through to all goods)
#   Fed to HIKE           (responds to inflation with rate increases)
#   mortgage rates to RISE (follow Fed funds rate)
#   yield curve to FLATTEN (short rates rise faster than long rates)
#   HPI to soften         (higher rates cool housing)
#
# These co-movements are not optional — they are the transmission mechanism.
# Holding them constant while shocking oil gives the WRONG answer because it
# ignores the Fed's reaction function and inflation channel entirely.
#
# Solution: define a full macro scenario matrix for each Moody's path,
# consistent with how the macro variables historically co-moved with oil.
# We use the CCAR 2026 macro paths (macro_base / macro_severe) as the
# primary source, and construct the Moody's-consistent paths by scaling
# the CCAR severe scenario by the ratio of Moody's oil shock to CCAR oil shock.
# =============================================================================

simulate_full_macro <- function(var_model, init_Y, macro_z_path, p) {
  # macro_z_path: data.table with one row per projection quarter
  #               columns = all Z variable names as they appear in var_model
  tryCatch(
    forecast_varx(var_model, macro_z_path, init_Y, p),
    error=function(e) { msg("  Simulation failed: %s", e$message); NULL }
  )
}

# ── Historical macro means (baseline anchor) ─────────────────────────────────
# Use recent 2-year mean (2024-2025) as the starting point — more relevant than
# 5-year mean which includes ZIRP era and distorts base levels
recent_agg <- agg_clean[yyyyqq >= 202401L]

macro_hist_mean <- sapply(z_var_names_sim, function(v) {
  if (v %in% names(recent_agg)) mean(recent_agg[[v]], na.rm=TRUE) else 0
})

# ── MARKET-CALIBRATED ELASTICITIES (March 25, 2026) ──────────────────────────
# All calibrated against CURRENT market conditions — not historical averages.
# Key sources:
#   Fed Funds: March 2026 FOMC hold; dot plot = 1 cut H2 2026, 1 cut 2027
#              CME FedWatch: ~0% probability of hike in any 2026 scenario
#              Fed explicitly "looking through" energy inflation → FOMC = HOLD
#   CPI: Feb 2026 = 2.4% YoY; tariff inflation already elevated baseline
#        Oil pass-through elasticity still ~0.04 but starts from higher base
#   Unemployment: Feb 2026 = 4.4% and RISING (labor market softening)
#                 Oil employment boost smaller than historically — energy sector
#                 hiring offset by trade/tariff headwinds; reduce elasticity
#   Mortgage rate: 6.22-6.45% currently (Freddie Mac March 19, 2026)
#                  Fed cutting NOT hiking → mortgage rates drift DOWN, not up
#                  Long end rising on geopolitical inflation fears, but slowly
#                  Net elasticity ≈ near zero or mildly negative
#   Yield curve: Currently +0.3pp (re-steepened after Fed cuts)
#                Oil spike → long rates (10Y) rise on inflation premium faster
#                than short rates (3M, anchored by hold/cut Fed) → STEEPENS
#                This is the OPPOSITE of the old assumption (which had it flatten)
#   HPI: Keep positive — oil-region housing demand still valid
#
# Starting (anchor) values as of 2026Q1:
ANCHOR <- list(
  macro_base_lurc         =  4.4,    # BLS Feb 2026 unemployment
  macro_base_pcpi         =  316.0,  # CPI level approx (BLS Feb 2026)
  macro_base_yield_curve  =  0.30,   # 10Y-3M spread, approx March 2026
  macro_base_rmtg         =  6.30,   # Freddie Mac 30Y FRM, March 2026
  hpi_yoy                 =  3.5,    # HPI YoY approx (housing slowing)
  macro_base_fomc_regime  =  0L,     # HOLD — neither hiking nor cutting yet
  macro_base_yoy_oil      =  0.0     # overridden by Moody's path
)

elast <- list(
  # Unemployment: oil rise has SMALLER effect now — labor market already
  # softening from tariff/trade headwinds; energy hiring partially offset
  # Reduced from -0.015 to -0.008 per 1% oil YoY
  macro_base_lurc         = -0.008,

  # CPI: standard oil pass-through; tariff inflation already in baseline
  # so marginal oil effect is the same but starts from elevated level
  macro_base_pcpi         =  0.040,

  # Yield curve: STEEPENS under oil shock — long rates rise faster than
  # short rates which are anchored by Fed hold/cut. POSITIVE elasticity.
  # (Old assumption of flattening was wrong for current regime)
  macro_base_yield_curve  =  0.008,

  # Mortgage rate: Fed cutting → short-rate anchor falling
  # Long end rising on inflation premium but slowly; net effect ≈ flat
  # In Feb baseline (oil fall) rates drift down; in Recession rates tick up
  # Elasticity near-zero: -0.005 (mild inverse — Fed cuts outweigh long-end rise)
  macro_base_rmtg         = -0.005,

  # HPI: oil-region housing demand positive but smaller with higher rates
  hpi_yoy                 =  0.040
)

# ── Build full macro path for each Moody's scenario ──────────────────────────
build_macro_path <- function(oil_yoy_proj, scenario_name) {
  h  <- length(oil_yoy_proj)
  dt <- data.table(matrix(0, nrow=h, ncol=length(z_var_names_sim),
                           dimnames=list(NULL, z_var_names_sim)))

  # 2-quarter smoothed oil (avoids over-reaction to single-quarter spike)
  oil_smooth <- frollmean(oil_yoy_proj, 2, na.rm=TRUE, align="right")
  oil_smooth[is.na(oil_smooth)] <- oil_yoy_proj[is.na(oil_smooth)]

  for (v in z_var_names_sim) {

    base_val <- ANCHOR[[v]] %||% (macro_hist_mean[v] %||% 0)
    if (is.na(base_val)) base_val <- 0

    if (v == "macro_base_yoy_oil") {
      dt[, (v) := oil_yoy_proj]

    } else if (v == "macro_base_fomc_regime") {
      # ── MARKET-CALIBRATED FED REACTION FUNCTION ──────────────────────────
      # Based on March 2026 FOMC statement and dot plot:
      #   - Fed is on HOLD for all of 2026 (0) in baseline/current
      #   - 1 cut expected H2 2026 → by Q3/Q4 shifts to -1 (cutting)
      #   - Even at $125/bbl oil, Fed is "looking through" energy inflation
      #     (Powell explicitly said this; CME shows 0% hike probability)
      #   - ONLY if oil COLLAPSES (<-20% YoY, i.e. bust) does Fed cut sooner
      #   - HIKE (+1) is ruled out entirely in all 2026 scenarios
      fomc_path <- ifelse(
        oil_yoy_proj < -20,         # Oil bust → faster cuts
        -1L,                         # Cutting
        0L                           # Hold (baseline) — NOT hiking
      )
      # After Q4 2026, allow 1 cut regardless (dot plot guidance)
      fomc_path[8:h] <- pmin(fomc_path[8:h], -1L)
      dt[, (v) := as.integer(fomc_path)]

    } else if (v %in% names(elast)) {
      path <- base_val + elast[[v]] * oil_smooth

      # Realistic bounds anchored to current market levels
      if (v == "macro_base_lurc")
        path <- pmax(pmin(path, 5.5), 3.5)   # 3.5-5.5% range
      if (v == "macro_base_rmtg")
        path <- pmax(pmin(path, 7.5), 5.5)   # 5.5-7.5% range
      if (v == "macro_base_yield_curve")
        path <- pmax(pmin(path, 1.5), -0.5)  # -0.5 to +1.5pp range
      if (v == "hpi_yoy")
        path <- pmax(pmin(path, 8.0), -3.0)  # bounded HPI growth

      dt[, (v) := path]

    } else {
      dt[, (v) := base_val]
    }
  }

  # Recompute interaction terms consistently with new paths
  if (all(c("macro_base_fomc_regime","macro_base_yoy_oil","fomc_x_brent") %in% names(dt)))
    dt[, fomc_x_brent := macro_base_fomc_regime * macro_base_yoy_oil]
  if (all(c("post_shale","macro_base_yoy_oil","post_x_oil") %in% names(dt)))
    dt[, post_x_oil := post_shale * macro_base_yoy_oil]

  # Structural dummies for projection period
  if ("post_shale"  %in% names(dt)) dt[, post_shale  := 1L]
  if ("gfc_dummy"   %in% names(dt)) dt[, gfc_dummy   := 0L]
  if ("covid_dummy" %in% names(dt)) dt[, covid_dummy  := 0L]
  if ("zirp_era"    %in% names(dt)) dt[, zirp_era    := 0L]
  # hike_cycle = 0 for all 2026 scenarios (Fed not hiking)
  if ("hike_cycle"  %in% names(dt)) dt[, hike_cycle  := 0L]

  dt[, yyyyqq := moodys_yyyyqq[PROJ_START_IDX:(PROJ_START_IDX + h - 1L)]]
  dt
}

# Print what macro paths look like for Current Baseline and Recession
cat("\n  Macro co-movement under Current Baseline scenario:\n")
tmp <- build_macro_path(oil_proj_current, "Current Baseline")
cat(sprintf("    Oil YoY     : %s\n", paste(sprintf("%+.0f%%", tmp$macro_base_yoy_oil), collapse=", ")))
if ("macro_base_lurc" %in% names(tmp))
  cat(sprintf("    Unemployment: %s\n", paste(sprintf("%.1f", tmp$macro_base_lurc), collapse=", ")))
if ("macro_base_pcpi" %in% names(tmp))
  cat(sprintf("    CPI         : %s\n", paste(sprintf("%.1f", tmp$macro_base_pcpi), collapse=", ")))
if ("macro_base_fomc_regime" %in% names(tmp))
  cat(sprintf("    FOMC regime : %s\n", paste(tmp$macro_base_fomc_regime, collapse=", ")))
if ("macro_base_rmtg" %in% names(tmp))
  cat(sprintf("    Mortgage    : %s\n", paste(sprintf("%.1f", tmp$macro_base_rmtg), collapse=", ")))

# ── Build all 5 full macro paths ─────────────────────────────────────────────
macro_path_flat    <- build_macro_path(oil_proj_flat,    "Flat oil (counterfactual)")
macro_path_feb26   <- build_macro_path(oil_proj_feb26,   "Feb '26 Baseline")
macro_path_march   <- build_macro_path(oil_proj_march,   "Early March Baseline")
macro_path_current <- build_macro_path(oil_proj_current, "Current Baseline")
macro_path_recess  <- build_macro_path(oil_proj_recess,  "Recession Scenario")

# ── Run full-macro simulations ────────────────────────────────────────────────
sim_flat    <- simulate_full_macro(varx_full, init_Y_sim, macro_path_flat,    P_LAG)
sim_feb26   <- simulate_full_macro(varx_full, init_Y_sim, macro_path_feb26,   P_LAG)
sim_march   <- simulate_full_macro(varx_full, init_Y_sim, macro_path_march,   P_LAG)
sim_current <- simulate_full_macro(varx_full, init_Y_sim, macro_path_current, P_LAG)
sim_recess  <- simulate_full_macro(varx_full, init_Y_sim, macro_path_recess,  P_LAG)

# =============================================================================
# MOODY'S ANALYTICS PBRENT LEVEL PATHS
# Source: "Higher For Longer Oil Prices" (Moody's Analytics, March 2026)
# Digitised from chart — quarterly average Brent crude $/bbl
# Historical anchor 2025Q1-Q4 common to all scenarios
# =============================================================================
moodys_yyyyqq <- c(
  202501L,202502L,202503L,202504L,   # 2025 — historical
  202601L,202602L,202603L,202604L,   # 2026 — projection
  202701L,202702L,202703L,202704L,   # 2027 — projection
  202801L,202802L,202803L,202804L    # 2028 — projection
)

# Level paths ($/bbl) — read from chart at quarterly resolution
pbrent_feb26  <- c(75, 67, 65, 63,   64,  64,  64,  64,   65, 66, 67, 68,   69, 69, 70, 70)
pbrent_march  <- c(75, 67, 65, 63,   70,  77,  77,  73,   71, 70, 70, 69,   69, 69, 70, 70)
pbrent_current<- c(75, 67, 65, 63,   80,  95,  93,  83,   76, 73, 72, 71,   70, 70, 70, 70)
pbrent_recess <- c(75, 67, 65, 63,   90, 125, 110,  90,   80, 75, 72, 71,   71, 70, 70, 70)

# Convert level paths → YoY % change
# YoY = (Q_t - Q_{t-4}) / Q_{t-4} * 100
# For 2026 quarters we compare vs 2025 (available); 2027 vs 2026; 2028 vs 2027
level_to_yoy <- function(levels) {
  n  <- length(levels)
  yy <- rep(NA_real_, n)
  for (i in 5:n) yy[i] <- (levels[i] - levels[i-4]) / levels[i-4] * 100
  yy
}

yoy_feb26   <- level_to_yoy(pbrent_feb26)
yoy_march   <- level_to_yoy(pbrent_march)
yoy_current <- level_to_yoy(pbrent_current)
yoy_recess  <- level_to_yoy(pbrent_recess)

# Projection window only (2026Q1 onwards = indices 5:16 of the 16-quarter path)
# The VARX simulation starts at the first projection quarter
PROJ_START_IDX <- 5L   # 2026Q1
n_proj <- SIM_HORIZON - PROJ_START_IDX + 1L   # 12 projection quarters

extract_proj_yoy <- function(yoy_vec) {
  yoy_vec[PROJ_START_IDX:SIM_HORIZON]
}

# Simulation uses 12-quarter projection window (2026Q1-2028Q4)
oil_proj_feb26   <- extract_proj_yoy(yoy_feb26)
oil_proj_march   <- extract_proj_yoy(yoy_march)
oil_proj_current <- extract_proj_yoy(yoy_current)
oil_proj_recess  <- extract_proj_yoy(yoy_recess)

# Print YoY paths for verification
cat("\n  Moody's PBRENT YoY paths (2026Q1-2028Q4):\n")
cat(sprintf("  %-22s : %s\n", "Feb 2026 Baseline",
    paste(sprintf("%+.0f%%", oil_proj_feb26), collapse=", ")))
cat(sprintf("  %-22s : %s\n", "Early March Baseline",
    paste(sprintf("%+.0f%%", oil_proj_march), collapse=", ")))
cat(sprintf("  %-22s : %s\n", "Current Baseline",
    paste(sprintf("%+.0f%%", oil_proj_current), collapse=", ")))
cat(sprintf("  %-22s : %s\n", "Recession Scenario",
    paste(sprintf("%+.0f%%", oil_proj_recess), collapse=", ")))

# ── Run simulations for each Moody's scenario + a flat-oil baseline ──────────
# Flat baseline: oil YoY = 0% (neutral — to compute delta / marginal effect)
oil_proj_flat <- rep(0, n_proj)

SIM_HOR_PROJ <- n_proj   # 12 quarters

sim_flat    <- simulate_oil_path(varx_full, init_Y_sim, oil_proj_flat,    agg_clean, z_var_names_sim, SIM_HOR_PROJ)
sim_feb26   <- simulate_oil_path(varx_full, init_Y_sim, oil_proj_feb26,   agg_clean, z_var_names_sim, SIM_HOR_PROJ)
sim_march   <- simulate_oil_path(varx_full, init_Y_sim, oil_proj_march,   agg_clean, z_var_names_sim, SIM_HOR_PROJ)
sim_current <- simulate_oil_path(varx_full, init_Y_sim, oil_proj_current, agg_clean, z_var_names_sim, SIM_HOR_PROJ)
sim_recess  <- simulate_oil_path(varx_full, init_Y_sim, oil_proj_recess,  agg_clean, z_var_names_sim, SIM_HOR_PROJ)

if (!is.null(sim_flat) && !is.null(sim_current)) {

  # ── Compute DELTA = Moody's scenario minus flat-oil baseline ─────────────
  # Delta = pure marginal effect of oil price path vs neutral flat-oil world
  build_delta <- function(sim_shocked, sim_base, label) {
    if (is.null(sim_shocked)) return(NULL)
    delta <- sim_shocked - sim_base
    dt    <- as.data.table(delta)
    setnames(dt, colnames(delta))
    dt[, horizon  := seq_len(.N) - 1L]
    dt[, scenario := label]
    # Attach calendar quarter labels
    dt[, yyyyqq   := moodys_yyyyqq[PROJ_START_IDX:(PROJ_START_IDX + .N - 1L)]]
    melt(dt, id.vars=c("horizon","scenario","yyyyqq"),
         variable.name="outcome", value.name="delta")
  }

  delta_feb26   <- build_delta(sim_feb26,   sim_flat, "Feb '26 Baseline")
  delta_march   <- build_delta(sim_march,   sim_flat, "Early March Baseline")
  delta_current <- build_delta(sim_current, sim_flat, "Current Baseline")
  delta_recess  <- build_delta(sim_recess,  sim_flat, "Recession Scenario")

  all_deltas <- rbindlist(
    Filter(Negate(is.null),
           list(delta_feb26, delta_march, delta_current, delta_recess))
  )

  saveRDS(all_deltas, "Results/04_oil_shock_simulation.rds")
  msg("  Moody's scenario simulation saved \u2192 Results/04_oil_shock_simulation.rds")

  # ── CHART 7.5 — CU response to Moody's oil paths ─────────────────────────
  key_outcomes <- intersect(
    c("dq_rate","pll_rate","costfds","netintmrg",
      "insured_share_growth","member_growth_yoy"),
    unique(all_deltas$outcome)
  )

  outcome_labels <- c(
    dq_rate              = "Delinquency Rate\n(backward-looking credit)",
    pll_rate             = "PLL Rate\n(forward-looking provisions)",
    costfds              = "Cost of Funds\n(rate channel)",
    netintmrg            = "Net Interest Margin\n(profitability)",
    insured_share_growth = "Deposit Growth\n(buffer hypothesis)",
    member_growth_yoy    = "Membership Growth\n(early warning)"
  )

  # Match Moody's chart colours exactly
  SCEN_COLS <- c(
    "Feb '26 Baseline"     = "#4a90d9",   # blue
    "Early March Baseline" = "#2d7a4a",   # green
    "Current Baseline"     = "#e07b54",   # salmon/pink
    "Recession Scenario"   = "#4a2080"    # purple
  )
  SCEN_LT <- c(
    "Feb '26 Baseline"     = "dotted",
    "Early March Baseline" = "dashed",
    "Current Baseline"     = "solid",
    "Recession Scenario"   = "solid"
  )
  SCEN_SIZE <- c(
    "Feb '26 Baseline"     = 0.7,
    "Early March Baseline" = 0.8,
    "Current Baseline"     = 1.0,
    "Recession Scenario"   = 1.1
  )

  plot_dt_key <- all_deltas[outcome %in% key_outcomes]
  plot_dt_key[, outcome_label := outcome_labels[as.character(outcome)]]
  plot_dt_key[is.na(outcome_label), outcome_label := as.character(outcome)]
  plot_dt_key[, scenario := factor(scenario, levels=names(SCEN_COLS))]

  # Quarter labels for x-axis: Q1 2026 = horizon 0
  q_labels <- c("2026\nQ1","Q2","Q3","Q4",
                 "2027\nQ1","Q2","Q3","Q4",
                 "2028\nQ1","Q2","Q3","Q4")

  p_oil_shock <- ggplot(plot_dt_key[!is.na(delta)],
                        aes(x=horizon, y=delta,
                            colour=scenario, linetype=scenario)) +
    geom_hline(yintercept=0, linewidth=0.5, colour="#aaaaaa",
               linetype="dashed") +
    geom_line(linewidth=0.9) +
    geom_point(data=plot_dt_key[!is.na(delta) & horizon %% 4 == 0],
               size=2, show.legend=FALSE) +
    scale_colour_manual(values=SCEN_COLS, name="Moody's Scenario") +
    scale_linetype_manual(values=SCEN_LT,  name="Moody's Scenario") +
    scale_x_continuous(
      breaks = 0:(n_proj-1),
      labels = q_labels[1:n_proj]
    ) +
    facet_wrap(~outcome_label, scales="free_y", ncol=3) +
    labs(
      title    = "FIGURE 7.5 \u2014 CU System Response to Moody\u2019s Oil Price Scenarios",
      subtitle = paste(
        "Delta vs flat-oil baseline (YoY=0%%) | Colour/style matches Moody\u2019s Analytics chart",
        "\nHigher oil \u2192 higher deposits & CoF, but eventually higher delinquency as inflation bites"
      ),
      caption  = paste(
        "Source: Moody\u2019s Analytics \u2018Higher For Longer Oil Prices\u2019 (March 2026);",
        "VARX conditional forecast simulation",
        "\nDelta = Moody\u2019s scenario path minus flat-oil (YoY=0%%) counterfactual",
        "| All other macro vars at historical mean | VARX(p=2)"
      ),
      x = NULL,
      y = "Change from flat-oil baseline (level units)"
    ) +
    theme_minimal(base_size=10) +
    theme(
      plot.title       = element_text(face="bold", size=12, margin=margin(b=4)),
      plot.subtitle    = element_text(size=9, colour="#444444",
                                       lineheight=1.3, margin=margin(b=8)),
      plot.caption     = element_text(size=7.5, colour="#888888",
                                       hjust=0, lineheight=1.3),
      strip.text       = element_text(face="bold", size=9, lineheight=1.2),
      strip.background = element_rect(fill="#f5f5f5", colour="#cccccc"),
      panel.grid.major = element_line(colour="#eeeeee", linewidth=0.3),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour="#cccccc", fill=NA, linewidth=0.4),
      axis.text.x      = element_text(size=7, lineheight=1.2),
      legend.position  = "bottom",
      legend.title     = element_text(face="bold", size=9),
      legend.text      = element_text(size=8.5),
      legend.key.width = unit(1.5, "cm"),
      plot.margin      = margin(10, 12, 8, 10)
    )

  ggsave("Figures/04_oil_shock_response.png", p_oil_shock,
         width=14, height=9, dpi=300, bg="white")
  msg("  \u2713 Moody\u2019s scenario CU response chart \u2192 Figures/04_oil_shock_response.png")

  # ── Also save the oil price level paths as a reference chart ─────────────
  pbrent_dt <- rbindlist(list(
    data.table(yyyyqq=moodys_yyyyqq, pbrent=pbrent_feb26,   scenario="Feb '26 Baseline"),
    data.table(yyyyqq=moodys_yyyyqq, pbrent=pbrent_march,   scenario="Early March Baseline"),
    data.table(yyyyqq=moodys_yyyyqq, pbrent=pbrent_current, scenario="Current Baseline"),
    data.table(yyyyqq=moodys_yyyyqq, pbrent=pbrent_recess,  scenario="Recession Scenario")
  ))
  pbrent_dt[, yr  := yyyyqq %/% 100L]
  pbrent_dt[, qtr := yyyyqq %% 100L]
  pbrent_dt[, date_num := yr + (qtr - 1) / 4]
  pbrent_dt[, scenario := factor(scenario, levels=names(SCEN_COLS))]

  p_oil_paths <- ggplot(pbrent_dt, aes(x=date_num, y=pbrent,
                                        colour=scenario, linetype=scenario)) +
    geom_vline(xintercept=2026, linetype="dashed",
               colour="#888888", linewidth=0.5) +
    annotate("text", x=2026.05, y=130, label="Projection \u2192",
             hjust=0, size=3, colour="#888888", fontface="italic") +
    geom_line(linewidth=1.0) +
    scale_colour_manual(values=SCEN_COLS, name="Scenario") +
    scale_linetype_manual(values=SCEN_LT,  name="Scenario") +
    scale_x_continuous(breaks=2025:2029,
                       labels=paste0("'", 25:29)) +
    scale_y_continuous(labels=dollar_format(prefix="$", suffix="/bbl"),
                       breaks=seq(40, 140, 20)) +
    labs(
      title    = "Moody\u2019s Analytics PBRENT Scenarios (March 2026)",
      subtitle = "Higher For Longer Oil Prices | Source: Moody\u2019s Analytics",
      caption  = "Digitised from Moody\u2019s Analytics chart. Quarterly average Brent crude $/bbl.",
      x=NULL, y="Brent Crude ($/bbl)"
    ) +
    theme_minimal(base_size=10) +
    theme(
      plot.title       = element_text(face="bold", size=11),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
  ggsave("Figures/04_moodys_oil_paths.png", p_oil_paths,
         width=11, height=5, dpi=300, bg="white")
  msg("  \u2713 Moody\u2019s oil paths chart \u2192 Figures/04_moodys_oil_paths.png")

  # ── Peak response summary table ───────────────────────────────────────────
  peak_resp <- all_deltas[outcome %in% key_outcomes,
    .(peak_delta   = delta[which.max(abs(delta))],
      peak_horizon = horizon[which.max(abs(delta))]),
    by=.(scenario, outcome)]
  peak_resp[, outcome_label := outcome_labels[as.character(outcome)]]
  peak_resp[is.na(outcome_label), outcome_label := as.character(outcome)]
  peak_resp[, quarter_label := q_labels[peak_horizon + 1L]]

  cat("\n  Peak CU Response to Moody's Oil Scenarios (vs flat-oil baseline):\n")
  cat("  ", strrep("-", 90), "\n", sep="")
  cat(sprintf("  %-30s %-25s %12s %10s\n",
              "Outcome","Scenario","Peak delta","At quarter"))
  cat("  ", strrep("-", 90), "\n", sep="")
  for (i in seq_len(nrow(peak_resp))) {
    r <- peak_resp[i]
    cat(sprintf("  %-30s %-25s %+12.4f %10s\n",
                r$outcome_label, r$scenario,
                r$peak_delta, r$quarter_label))
  }
  fwrite(peak_resp, "Results/04_moodys_peak_responses.csv")
  msg("  Peak responses \u2192 Results/04_moodys_peak_responses.csv")

  # ── CHART 7.5b — Macro co-movement paths by scenario ─────────────────────
  # Shows how ALL macro variables differ across the 4 Moody's scenarios
  # This is the key methodological transparency chart

  macro_paths_all <- rbindlist(list(
    cbind(macro_path_feb26,   scenario="Feb '26 Baseline"),
    cbind(macro_path_march,   scenario="Early March Baseline"),
    cbind(macro_path_current, scenario="Current Baseline"),
    cbind(macro_path_recess,  scenario="Recession Scenario")
  ), fill=TRUE)

  macro_paths_all[, yr  := yyyyqq %/% 100L]
  macro_paths_all[, qtr := yyyyqq %% 100L]
  macro_paths_all[, date_num := yr + (qtr - 1) / 4]
  macro_paths_all[, scenario := factor(scenario, levels=names(SCEN_COLS))]

  macro_plot_vars <- list(
    list(col="macro_base_yoy_oil",     lab="Oil Price YoY (%)",        fmt="%+.0f%%"),
    list(col="macro_base_lurc",        lab="Unemployment Rate (%)",     fmt="%.1f%%"),
    list(col="macro_base_pcpi",        lab="CPI Level (index)",         fmt="%.1f"),
    list(col="macro_base_rmtg",        lab="30Y Mortgage Rate (%)",     fmt="%.1f%%"),
    list(col="macro_base_yield_curve", lab="Yield Curve (10Y-3M, pp)",  fmt="%.2f"),
    list(col="macro_base_fomc_regime", lab="FOMC Regime (+1=hike)",     fmt="%.0f")
  )

  macro_panel_list <- lapply(macro_plot_vars, function(v) {
    if (!v$col %in% names(macro_paths_all)) return(NULL)
    ggplot(macro_paths_all[!is.na(get(v$col))],
           aes(x=date_num, y=get(v$col),
               colour=scenario, linetype=scenario)) +
      geom_hline(yintercept=0, linewidth=0.3, colour="#cccccc") +
      geom_line(linewidth=0.85) +
      scale_colour_manual(values=SCEN_COLS, name=NULL) +
      scale_linetype_manual(values=SCEN_LT, name=NULL) +
      scale_x_continuous(breaks=2026:2029,
                         labels=paste0("'", 26:29)) +
      labs(title=v$lab, x=NULL, y=NULL) +
      theme_minimal(base_size=9) +
      theme(
        plot.title       = element_text(face="bold", size=9),
        legend.position  = "none",
        panel.grid.minor = element_blank(),
        panel.border     = element_rect(colour="#dddddd", fill=NA, linewidth=0.3)
      )
  })
  macro_panel_list <- Filter(Negate(is.null), macro_panel_list)

  if (length(macro_panel_list) >= 2) {
    p_macro_paths <- wrap_plots(macro_panel_list, ncol=3) +
      plot_annotation(
        title    = "FIGURE 7.5b \u2014 Macro Co-Movement Paths by Moody\u2019s Scenario",
        subtitle = paste(
          "All 6 VARX exogenous variables shown | Recession scenario: oil spikes \u2192 inflation rises",
          "\n\u2192 Fed hikes \u2192 mortgage rates rise \u2192 yield curve flattens \u2192 unemployment falls then rises"
        ),
        caption  = paste(
          "Macro co-movements derived from empirical oil-macro elasticities",
          "(Hamilton 2009; Kilian 2014)",
          "| Fed reaction function: hike if oil YoY > 20%%, cut if < -20%%"
        ),
        theme    = theme(
          plot.title    = element_text(face="bold", size=11),
          plot.subtitle = element_text(size=8.5, colour="#555", lineheight=1.3),
          plot.caption  = element_text(size=7.5, colour="#888")
        )
      )
    ggsave("Figures/04_macro_scenario_paths.png", p_macro_paths,
           width=13, height=8, dpi=300, bg="white")
    msg("  \u2713 Macro co-movement chart \u2192 Figures/04_macro_scenario_paths.png")
  }

} else {
  msg("  Oil shock simulation skipped — baseline simulation failed")
  msg("  Check that varx_full was estimated successfully in Section 3")
}

}  # end varx_full guard (Section 7.5)

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

  k   <- ncol(init_Y)
  h   <- nrow(macro_scenario_path)
  nms <- colnames(init_Y)
  fcst <- matrix(NA_real_, nrow=h, ncol=k, dimnames=list(NULL, nms))

  # ── Companion-form AR coefficient matrices ───────────────────────────────
  # vars::VAR stores each equation's coefficients separately.
  # For lag l, the (k × k) matrix A_l has:
  #   row i = equation for Y_i
  #   col j = coefficient on Y_j at lag l
  # sapply(...) produces a (k × k) matrix with equations as COLUMNS — transpose needed.
  A_list <- lapply(1:p, function(lag) {
    M <- sapply(var_model$varresult, function(eq) {
      cf    <- coef(eq)
      terms <- paste0(nms, ".l", lag)
      # Return k coefficients in the same order as nms
      vapply(terms, function(t) {
        if (t %in% names(cf)) cf[[t]] else 0
      }, numeric(1))
    })
    # M is k×k: rows=predictors, cols=equations → transpose to equations×predictors
    t(M)
  })

  # Intercept vector (length k)
  const <- vapply(var_model$varresult, function(eq) {
    cf <- coef(eq)
    if ("const" %in% names(cf)) cf[["const"]] else 0
  }, numeric(1))

  # ── Exogenous coefficient matrix B (k × n_exog) ──────────────────────────
  exog_nms <- names(var_model$varresult[[1]]$model)
  # Remove endogenous lags and intercept
  lag_pattern <- paste0("^(", paste(nms, collapse="|"), ")\\.l")
  exog_nms <- exog_nms[!grepl(lag_pattern, exog_nms) & exog_nms != "(Intercept)"]

  if (length(exog_nms) > 0) {
    B_exog <- t(sapply(var_model$varresult, function(eq) {
      cf <- coef(eq)
      vapply(exog_nms, function(nm) {
        if (nm %in% names(cf)) cf[[nm]] else 0
      }, numeric(1))
    }))
    # B_exog: k × n_exog after transpose
  } else {
    B_exog <- matrix(0, nrow=k, ncol=0)
  }

  Y_hist <- init_Y

  for (t in 1:h) {
    y_hat <- matrix(const, ncol=1)   # k × 1

    for (lag in 1:p) {
      idx <- nrow(Y_hist) - lag + 1
      if (idx < 1) break
      y_lag <- matrix(Y_hist[idx, ], ncol=1)   # k × 1 — force column vector
      y_hat <- y_hat + A_list[[lag]] %*% y_lag
    }

    # Add exogenous contribution
    if (length(exog_nms) > 0 && ncol(B_exog) > 0) {
      z_nms <- intersect(exog_nms, names(macro_scenario_path))
      if (length(z_nms) > 0) {
        z_t   <- unlist(macro_scenario_path[t, ..z_nms])
        z_mat <- matrix(z_t, ncol=1)
        B_sub <- B_exog[, match(z_nms, exog_nms), drop=FALSE]
        y_hat <- y_hat + B_sub %*% z_mat
      }
    }

    fcst[t, ] <- as.numeric(y_hat)
    Y_hist    <- rbind(Y_hist, t(y_hat))
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

if (length(forecast_results) > 0) {
  all_forecasts <- rbindlist(forecast_results, fill=TRUE)
} else {
  msg("  WARNING: no scenario forecasts produced")
  all_forecasts <- data.table()
}

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

if (length(fan_plots) >= 2) {
  p_fan_combined <- wrap_plots(fan_plots, ncol=2) +
    plot_annotation(
      title    = "CU System \u2014 CCAR 2026 Scenario Fan Charts",
      subtitle = sprintf("VARX(p=%d) | Cholesky-identified | Full sample 2005Q1\u20132025Q4", P_LAG),
      theme    = theme(plot.title=element_text(face="bold", size=13))
    )
  ggsave("Figures/04_fan_charts_all.png", p_fan_combined,
         width=14, height=10, dpi=150, bg="white")
  msg("Fan charts \u2192 Figures/04_fan_charts_all.png")
} else {
  msg("Fan charts skipped — fewer than 2 dep vars have projection data")
}

# =============================================================================
# 9. HEADLINE COEFFICIENT EXTRACT — The Three-Way Interaction
# =============================================================================
hdr("SECTION 9: Headline Finding — Triple Interaction Extract")

if (!exists("coef_table") || nrow(coef_table) == 0) {
  msg("  coef_table empty — Section 9 skipped (panel FE section must run first)")
} else {

# Pull post_shale × yoy_oil × direct_exposure from panel FE models
triple_coefs <- coef_table[term %like% "post_x_oil_x_direct",
                             .(dep_var, term, estimate, se, p)]
setorder(triple_coefs, dep_var)

cat("\n  Triple Interaction: post_shale \u00d7 yoy_oil \u00d7 direct_exposure\n")
cat("  (Interpretation: oil price rise impact on CUs with direct energy exposure\n")
cat("   in post-shale era vs. pre-shale era)\n\n")

triple_coefs[, sig := fcase(
  p < 0.01, "***",
  p < 0.05, "**",
  p < 0.10, "*",
  default  = ""
)]

if (nrow(triple_coefs) > 0)
  print(triple_coefs[, .(dep_var, estimate=round(estimate,5),
                           se=round(se,5), p=round(p,4), sig)])
else
  msg("  No post_x_oil_x_direct terms found in coef_table")

# Also extract fomc × oil (two-way — rate channel activation)
fomc_oil_coefs <- coef_table[term %like% "fomc_x_brent",
                               .(dep_var, term, estimate, se, p)]
fomc_oil_coefs[, sig := fcase(p<0.01,"***", p<0.05,"**", p<0.10,"*", default="")]

cat("\n  Two-Way Interaction: fomc_regime \u00d7 yoy_oil (rate channel)\n\n")
if (nrow(fomc_oil_coefs) > 0)
  print(fomc_oil_coefs[, .(dep_var, estimate=round(estimate,5),
                             se=round(se,5), p=round(p,4), sig)])
else
  msg("  No fomc_x_brent terms found in coef_table")

out_coefs <- rbindlist(list(
  if (nrow(triple_coefs)   > 0) triple_coefs[,   type := "triple"]   else NULL,
  if (nrow(fomc_oil_coefs) > 0) fomc_oil_coefs[, type := "fomc_oil"] else NULL
), fill=TRUE)

if (nrow(out_coefs) > 0) {
  fwrite(out_coefs, "Results/04_headline_coefs.csv")
  msg("\nHeadline coefficients \u2192 Results/04_headline_coefs.csv")
}

}  # end coef_table guard

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
