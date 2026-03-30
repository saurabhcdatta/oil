# =============================================================================
# SCRIPT 01 — Data Ingestion, Cleaning, and Variable Construction
# Oil Price Shock × Credit Union Financial Performance — v2
# Author: Saurabh C. Datta, Ph.D. | NCUA Office of Chief Economist
# Version: v2.2-DIAGNOSTIC | Built: 2026-03-30
# Working dir: S:/Projects/Oil_Price_Shock_2026/
# =============================================================================
# Inputs
#   Data/call_report.rds                  NCUA call report panel (CU x quarter)
#   Data/FRB_Baseline_2026.xlsx           FRB CCAR 2026 Baseline
#   Data/FRB_Severely_Adverse_2026.xlsx   FRB CCAR 2026 Severely Adverse
#   Data/oil_exposure.rds                 Built by 01b_oil_exposure_v2.R
#
# Outputs
#   Data/call_clean.rds      Cleaned call report with exposure + derived vars
#   Data/macro_base.rds      Baseline macro (wide, macro_base_ prefix)
#   Data/macro_severe.rds    Severely adverse macro (wide, macro_severe_ prefix)
#   Data/panel_base.rds      CU x quarter panel merged with baseline macro
#   Data/panel_severe.rds    CU x quarter panel merged with severely adverse
#   Results/01_outlier_log.csv         Winsorization summary per variable
#   Results/01_variable_summary.csv    Post-cleaning distributional stats
#   Results/01_checkpoint_log.csv      All checkpoint pass/warn/fail records
#
# Identifiers : join_number, year, quarter
# yyyyqq      : integer format year*100+quarter (e.g. 200501, 201504)
# Start date  : 2005Q1
# Run order   : 01b must run first to produce Data/oil_exposure.rds
#
# BACKWARD REASONING ANCHORS (from v1 confirmed findings):
#   - AR(1) approx 0.59 -> lagged outcomes constructed here
#   - Deposit growth = most exposed (3x SHAP) -> insured_share_growth critical
#   - pll_rate = pll / avg(lns_tot) -> confirmed construction from v1 Script 04d
#   - 2015Q1 structural break -> post2015 flag aligned to date here
#   - dq_rate, netintmrg, pcanetworth, costfds may be pre-computed in call_report
# =============================================================================

rm(list = ls())
gc()

# ============================================================================
# SECTION 0: LOGGING INFRASTRUCTURE
# ============================================================================

CHKPT_LOG <- data.frame(
  checkpoint = character(), status = character(), detail = character(),
  stringsAsFactors = FALSE
)

SEP  <- paste(rep("=", 72), collapse = "")
SEP2 <- paste(rep("-", 72), collapse = "")
ts   <- function() format(Sys.time(), "[%H:%M:%S]")

log_checkpoint <- function(label, status, detail = "") {
  icon <- switch(status,
    PASS = "PASS", WARN = "WARN", FAIL = "FAIL", INFO = "INFO"
  )
  cat(sprintf("  %s [%s] %s\n", ts(), icon, label))
  if (nchar(detail) > 0) cat(sprintf("         => %s\n", detail))
  CHKPT_LOG <<- rbind(CHKPT_LOG, data.frame(
    checkpoint = label, status = status, detail = detail,
    stringsAsFactors = FALSE
  ))
  if (status == "FAIL") {
    cat("\n", SEP, "\n  FATAL: Script halted at: ", label, "\n", SEP, "\n", sep = "")
    stop(paste("FATAL checkpoint failure:", label, "|", detail))
  }
}

section <- function(n, title) {
  cat("\n", SEP, "\n", sep = "")
  cat(sprintf("  SECTION %s: %s\n", n, title))
  cat(SEP2, "\n", sep = "")
}

msg <- function(...) cat(sprintf(...), "\n")

describe_var <- function(x, varname, tag = "") {
  n     <- sum(!is.na(x))
  nmiss <- mean(is.na(x))
  if (n == 0) {
    cat(sprintf("  %-28s  ALL MISSING\n", varname))
    return(invisible(NULL))
  }
  qs <- quantile(x, c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99), na.rm = TRUE)
  cat(sprintf("  %-28s  [%s]\n", varname, tag))
  cat(sprintf("    N=%-8s  miss=%5.1f%%  mean=%10.5f  sd=%10.5f\n",
              format(n, big.mark = ","), nmiss * 100,
              mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)))
  cat(sprintf("    p01=%9.5f  p10=%9.5f  p25=%9.5f  p50=%9.5f\n",
              qs[1], qs[2], qs[3], qs[4]))
  cat(sprintf("    p75=%9.5f  p90=%9.5f  p99=%9.5f\n", qs[5], qs[6], qs[7]))
  cat(sprintf("    min=%9.5f  max=%9.5f\n",
              min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
  cat(sprintf("    %s\n", paste(rep("-", 66), collapse = "")))
  invisible(list(n = n, nmiss = nmiss,
                 mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
                 p1 = qs[[1]], p99 = qs[[7]],
                 min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE)))
}

winsorize_audit <- function(x, p_lo = 0.01, p_hi = 0.99) {
  qs   <- quantile(x, c(p_lo, p_hi), na.rm = TRUE)
  lo   <- qs[1]; hi <- qs[2]
  n_lo <- sum(!is.na(x) & x < lo, na.rm = TRUE)
  n_hi <- sum(!is.na(x) & x > hi, na.rm = TRUE)
  list(
    x   = pmax(pmin(x, hi), lo),
    lo  = lo, hi = hi,
    n_lo = n_lo, n_hi = n_hi,
    pct  = round((n_lo + n_hi) / max(sum(!is.na(x)), 1) * 100, 3)
  )
}

# YoY growth by CU (lag-4 within join_number; data must be sorted CU + time)
cu_yoy <- function(dt, col) {
  dt[, .(val = {
    x    <- get(col)
    x_l4 <- shift(x, 4L)
    fifelse(
      !is.na(x) & !is.na(x_l4) & is.finite(x_l4) & x_l4 > 0,
      (x - x_l4) / x_l4 * 100,
      NA_real_
    )
  }), by = join_number]$val
}

# ============================================================================
# SCRIPT HEADER
# ============================================================================
cat(SEP, "\n")
cat("  SCRIPT 01 -- Data Ingestion, Cleaning, and Variable Construction\n")
cat("  Oil Price Shock x Credit Union Financial Performance -- v2\n")
cat("  Author  : Saurabh C. Datta, Ph.D. | NCUA Office of Chief Economist\n")
cat("  Version : v2.2-DIAGNOSTIC\n")
cat(sprintf("  Started : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(SEP, "\n\n")

# ============================================================================
# SECTION 1: LIBRARIES AND CONFIGURATION
# ============================================================================
section("1", "Libraries and Configuration")

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(lubridate)
  library(stringr)
})

for (pkg in c("data.table", "readxl", "lubridate", "stringr")) {
  log_checkpoint(paste("Package:", pkg), "PASS",
                 paste("v", packageVersion(pkg)))
}

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
START_YEAR   <- 2005L
WIN_LOW      <- 0.01
WIN_HIGH     <- 0.99
MIN_CU_COUNT <- 1000L
MIN_QTR      <- 40L

# Asset tier scheme: 8-tier assets_cat2 matching OCE_combined Stata file
# Breaks and labels defined in Section 5 where they are used
# (Removed old 4-tier scheme; replaced with 8-tier assets_cat2)

# Preliminary oil states (superseded by oil_exposure.rds from 01b,
# but used as fallback if 01b has not yet been run)
OIL_STATES_PRELIM <- c("TX", "ND", "LA", "AK", "WY", "OK",
                        "NM", "CO", "WV", "PA", "MT")

# Quarter-to-month mapping
Q_MONTH <- c("1" = 1L, "2" = 4L, "3" = 7L, "4" = 10L)

DATA_DIR    <- "Data/"
RESULTS_DIR <- "Results/"

for (d in c(DATA_DIR, RESULTS_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat(sprintf("\n  START_YEAR  : %d\n", START_YEAR))
cat(sprintf("  Winsorize   : [%.0f%%, %.0f%%]\n", WIN_LOW * 100, WIN_HIGH * 100))
cat("  Asset tiers : 8-tier assets_cat2 (see Section 5)\n")
cat(sprintf("  Min CUs     : %d\n", MIN_CU_COUNT))
cat(sprintf("  Min qtrs    : %d\n\n", MIN_QTR))

log_checkpoint("CK-01: Config loaded", "PASS",
               sprintf("START_YEAR=%d, WIN=[%.0f,%.0f]%%",
                       START_YEAR, WIN_LOW * 100, WIN_HIGH * 100))

# ============================================================================
# SECTION 2: CALL REPORT -- LOAD AND INITIAL CHECKS
# ============================================================================
section("2", "Call Report -- Load and Initial Checks")

RAW_FILE <- "Data/call_report.rds"
cat(sprintf("  File: %s\n\n", RAW_FILE))

if (!file.exists(RAW_FILE)) {
  log_checkpoint("CK-02: call_report.rds exists", "FAIL",
                 paste("Not found:", RAW_FILE))
}
log_checkpoint("CK-02: call_report.rds exists", "PASS", RAW_FILE)

t_load <- system.time({
  cr <- setDT(readRDS(RAW_FILE))
  setnames(cr, tolower(names(cr)))   # harmonise to lowercase
})

cat(sprintf("  Load time   : %.1f sec\n", t_load["elapsed"]))
cat(sprintf("  Raw rows    : %s\n", format(nrow(cr), big.mark = ",")))
cat(sprintf("  Raw columns : %d\n", ncol(cr)))
cat(sprintf("  Object size : %.1f MB\n\n", object.size(cr) / 1e6))

cat("  Column inventory (first 50):\n")
cat(sprintf("    %s\n", paste(head(names(cr), 50), collapse = ", ")))
if (ncol(cr) > 50)
  cat(sprintf("    ... and %d more columns\n", ncol(cr) - 50))

if (nrow(cr) < 10000L) {
  log_checkpoint("CK-03: Raw row count", "WARN",
                 sprintf("Only %s rows -- expected >10,000",
                         format(nrow(cr), big.mark = ",")))
} else {
  log_checkpoint("CK-03: Raw row count", "PASS",
                 sprintf("%s rows", format(nrow(cr), big.mark = ",")))
}

# ============================================================================
# SECTION 3: IDENTIFIER STANDARDIZATION AND YEAR FILTER
# ============================================================================
section("3", "Identifier Standardization and Year Filter")

REQ_IDS      <- c("join_number", "year", "quarter")
missing_ids  <- setdiff(REQ_IDS, names(cr))
if (length(missing_ids) > 0) {
  log_checkpoint("CK-04: Required id columns", "FAIL",
                 paste("Missing:", paste(missing_ids, collapse = ", ")))
}
log_checkpoint("CK-04: Required id columns", "PASS",
               paste(REQ_IDS, collapse = ", "))

cr[, `:=`(year = as.integer(year), quarter = as.integer(quarter))]

cat(sprintf("  year range (raw)    : %d -- %d\n",
            min(cr$year, na.rm = TRUE), max(cr$year, na.rm = TRUE)))
cat(sprintf("  quarter values      : %s\n",
            paste(sort(unique(cr$quarter)), collapse = ", ")))

bad_qtr <- cr[!quarter %in% 1:4, .N]
if (bad_qtr > 0) {
  log_checkpoint("CK-05: Quarter values valid", "WARN",
                 sprintf("%d rows have quarter not in 1:4", bad_qtr))
} else {
  log_checkpoint("CK-05: Quarter values valid", "PASS", "All in {1,2,3,4}")
}

# Start-year filter
n_before  <- nrow(cr)
cr        <- cr[year >= START_YEAR]
n_dropped <- n_before - nrow(cr)

cat(sprintf("\n  Rows before %dQ1 filter : %s\n",
            START_YEAR, format(n_before, big.mark = ",")))
cat(sprintf("  Rows dropped           : %s (%.1f%%)\n",
            format(n_dropped, big.mark = ","), n_dropped / n_before * 100))
cat(sprintf("  Rows retained          : %s\n",
            format(nrow(cr), big.mark = ",")))

log_checkpoint("CK-06: Start-year filter", "PASS",
               sprintf("%s rows retained (>= %dQ1)",
                       format(nrow(cr), big.mark = ","), START_YEAR))

# Construct yyyyqq (integer) and cal_date
cr[, `:=`(
  yyyyqq   = year * 100L + quarter,
  cal_date = as.Date(paste(year, Q_MONTH[as.character(quarter)], "01", sep = "-"))
)]

cat(sprintf("\n  yyyyqq format   : integer (e.g. %d)\n", cr$yyyyqq[1]))
cat(sprintf("  cal_date range  : %s  to  %s\n",
            format(min(cr$cal_date)), format(max(cr$cal_date))))
cat(sprintf("  Unique quarters : %d\n", uniqueN(cr$yyyyqq)))
cat(sprintf("  Unique CUs      : %s\n\n",
            format(uniqueN(cr$join_number), big.mark = ",")))

n_cu  <- uniqueN(cr$join_number)
n_qtr <- uniqueN(cr$yyyyqq)

if (n_cu < MIN_CU_COUNT) {
  log_checkpoint("CK-07: CU count", "WARN",
                 sprintf("Only %d CUs -- expected >= %d", n_cu, MIN_CU_COUNT))
} else {
  log_checkpoint("CK-07: CU count", "PASS",
                 sprintf("%s unique CUs", format(n_cu, big.mark = ",")))
}
if (n_qtr < MIN_QTR) {
  log_checkpoint("CK-08: Quarter count", "WARN",
                 sprintf("Only %d quarters -- expected >= %d", n_qtr, MIN_QTR))
} else {
  log_checkpoint("CK-08: Quarter count", "PASS",
                 sprintf("%d quarters (%s to %s)",
                         n_qtr, format(min(cr$cal_date)), format(max(cr$cal_date))))
}

# Duplicate check
dups <- cr[, .N, by = .(join_number, year, quarter)][N > 1, .N]
if (dups > 0) {
  cat(sprintf("  WARNING: %d duplicate CU-quarter keys -- keeping first\n", dups))
  cr <- unique(cr, by = c("join_number", "year", "quarter"))
  log_checkpoint("CK-09: Duplicates", "WARN",
                 sprintf("%d duplicate keys removed", dups))
} else {
  log_checkpoint("CK-09: Duplicates", "PASS",
                 "No duplicates on join_number x year x quarter")
}

# State column detection
state_col <- intersect(c("reporting_state", "state_code", "state"), names(cr))[1]
if (is.na(state_col)) {
  log_checkpoint("CK-10: State column", "WARN",
                 "No state column found -- oil-state classification will use fallback")
} else {
  n_states <- uniqueN(cr[[state_col]])
  cat(sprintf("\n  State column: %s  |  %d unique values\n",
              state_col, n_states))
  if (n_states < 48L) {
    log_checkpoint("CK-10: State column", "WARN",
                   sprintf("Only %d states -- expected ~50+DC+territories", n_states))
  } else {
    log_checkpoint("CK-10: State column", "PASS",
                   sprintf("%s, %d unique values", state_col, n_states))
  }
}

# Sort before any lag operations
setorderv(cr, c("join_number", "year", "quarter"))

# ============================================================================
# SECTION 4: DIRECT-RATIO VARIABLE AUDIT
# ============================================================================
section("4", "Direct-Ratio Variable Audit")

cat("  Checking for pre-computed ratio columns in call_report.rds...\n\n")

DIRECT_RATIOS <- c("netintmrg", "networth", "networthalt", "pcanetworth",
                   "costfds", "roa", "dq_rate", "chg_tot_lns_ratio")
DIRECT_META   <- list(
  netintmrg         = "Net interest margin",
  networth          = "Net worth ($000s)",
  networthalt       = "Alternative net worth measure",
  pcanetworth       = "Net worth ratio (PCA)",
  costfds           = "Cost of funds",
  roa               = "Return on assets",
  dq_rate           = "Delinquency rate",
  chg_tot_lns_ratio = "Charge-off ratio"
)

cat(sprintf("  %-22s  %-6s  %-8s  %s\n", "Column", "Found", "Non-NA%", "Description"))
cat(sprintf("  %s\n", paste(rep("-", 65), collapse = "")))

direct_present <- character(0)
for (v in DIRECT_RATIOS) {
  found  <- v %in% names(cr)
  pct_ok <- if (found) round(mean(!is.na(cr[[v]])) * 100, 1) else NA
  cat(sprintf("  %-22s  %-6s  %s  %s\n",
              v,
              if (found) "YES" else "NO",
              if (found) sprintf("%5.1f%%", pct_ok) else "  ---",
              DIRECT_META[[v]]))
  if (found) direct_present <- c(direct_present, v)
}
cat(sprintf("\n  Direct ratios confirmed: %s\n",
            if (length(direct_present) > 0)
              paste(direct_present, collapse = ", ")
            else "NONE"))
log_checkpoint("CK-11: Direct ratio audit", "PASS",
               sprintf("%d of %d pre-computed ratios present",
                       length(direct_present), length(DIRECT_RATIOS)))

if (length(direct_present) > 0) {
  cat("\n  Distributions of pre-computed ratio variables:\n\n")
  for (v in direct_present) describe_var(cr[[v]], v, "pre-computed, pre-winsor")
}

# ============================================================================
# SECTION 5: ASSET TIER ASSIGNMENT  (assets_cat2 -- 8 tiers)
# ============================================================================
section("5", "Asset Tier Assignment (assets_cat2 -- 8 tiers)")

# Tier scheme matches assets_cat2 from OCE_combined Stata file exactly.
# Units in call report: $thousands.
# Breaks (in $thousands):  0 | 10,000 | 50,000 | 100,000 | 500,000 |
#                          1,000,000 | 5,000,000 | 10,000,000 | Inf
ASSET_BREAKS2 <- c(0, 10e3, 50e3, 100e3, 500e3, 1e6, 5e6, 10e6, Inf)
ASSET_LABELS2 <- c(
  "T1_under10M",        # < $10M
  "T2_10to50M",         # $10M up to $50M
  "T3_50to100M",        # $50M through $100M
  "T4_100to500M",       # Over $100M through $500M
  "T5_500Mto1B",        # Over $500M through $1B
  "T6_1Bto5B",          # Over $1B through $5B
  "T7_5Bto10B",         # Over $5B through $10B
  "T8_over10B"          # Over $10B
)

cat("  8-tier scheme matching assets_cat2 in OCE_combined Stata file:\n")
cat("    T1_under10M   : assets < $10M              (<      10,000)\n")
cat("    T2_10to50M    : $10M up to $50M            (10,000 --  50,000)\n")
cat("    T3_50to100M   : $50M through $100M         (50,000 -- 100,000)\n")
cat("    T4_100to500M  : Over $100M through $500M   (100,000 -- 500,000)\n")
cat("    T5_500Mto1B   : Over $500M through $1B     (500,000 -- 1,000,000)\n")
cat("    T6_1Bto5B     : Over $1B through $5B       (1,000,000 -- 5,000,000)\n")
cat("    T7_5Bto10B    : Over $5B through $10B      (5,000,000 -- 10,000,000)\n")
cat("    T8_over10B    : Over $10B                  (>= 10,000,000)\n\n")

# Strategy: prefer assets_cat2 if already present in data (pre-computed in Stata)
# Fall back to constructing from assets_tot / acct_010

if ("assets_cat2" %in% names(cr)) {
  # ── Path A: use pre-computed assets_cat2 directly ──────────────────────────
  cat("  assets_cat2 column FOUND in call_report.rds -- using pre-computed tiers\n\n")

  # Audit the pre-computed values
  cat("  Pre-computed assets_cat2 distribution:\n")
  cat(sprintf("  %-24s  %10s  %8s  %8s\n",
              "Category", "N (CU-qtrs)", "Pct", "Cum."))
  cat(sprintf("  %s\n", paste(rep("-", 56), collapse = "")))
  cat_dist <- cr[, .N, by = assets_cat2][order(assets_cat2)]
  cum_pct  <- 0
  for (i in seq_len(nrow(cat_dist))) {
    pct_i <- cat_dist$N[i] / nrow(cr) * 100
    cum_pct <- cum_pct + pct_i
    cat(sprintf("  %-24s  %10s  %7.1f%%  %7.1f%%\n",
                as.character(cat_dist$assets_cat2[i]),
                format(cat_dist$N[i], big.mark = ","),
                pct_i, cum_pct))
  }

  # Map existing labels to canonical T1-T8 codes for regression use
  # assets_cat2 Stata labels -> T-codes
  label_map <- c(
    "Assets < 10 million"        = "T1_under10M",
    "10M up to 50M"              = "T2_10to50M",
    "50M through 100M"           = "T3_50to100M",
    "Over 100M through 500M"     = "T4_100to500M",
    "Over 500M through 1B"       = "T5_500Mto1B",
    "Over 1B through 5B"         = "T6_1Bto5B",
    "Over 5B through 10B"        = "T7_5Bto10B",
    "Over 10B"                   = "T8_over10B"
  )

  cr[, asset_tier := factor(
    label_map[as.character(assets_cat2)],
    levels = ASSET_LABELS2
  )]

  n_unmapped <- cr[is.na(asset_tier) & !is.na(assets_cat2), .N]
  if (n_unmapped > 0) {
    unmapped_vals <- unique(cr[is.na(asset_tier) & !is.na(assets_cat2), assets_cat2])
    cat(sprintf("\n  WARNING: %d rows with unrecognised assets_cat2 labels:\n", n_unmapped))
    cat(sprintf("    %s\n", paste(head(unmapped_vals, 10), collapse = ", ")))
    cat("  These will be NA in asset_tier -- verify label_map above\n")
    log_checkpoint("CK-12a: assets_cat2 label mapping", "WARN",
                   sprintf("%d rows unmapped -- check label_map", n_unmapped))
  } else {
    log_checkpoint("CK-12a: assets_cat2 label mapping", "PASS",
                   "All assets_cat2 labels mapped to T1-T8 codes")
  }

} else {
  # ── Path B: construct from raw asset column ─────────────────────────────────
  cat("  assets_cat2 NOT in call_report.rds -- constructing from raw asset column\n\n")

  asset_col <- intersect(c("assets_tot", "acct_010"), names(cr))[1]

  if (!is.na(asset_col)) {
    cat(sprintf("  Asset column: %s\n\n", asset_col))
    describe_var(cr[[asset_col]], asset_col, "raw")

    # Zero/negative assets are data errors -- set tier to NA
    n_zero_assets <- cr[get(asset_col) <= 0, .N]
    if (n_zero_assets > 0)
      cat(sprintf("  Zero/negative assets: %d rows -> tier set to NA\n", n_zero_assets))

    cr[, asset_tier := factor(
      fifelse(
        get(asset_col) > 0,
        ASSET_LABELS2[findInterval(get(asset_col), ASSET_BREAKS2,
                                   rightmost.closed = FALSE)],
        NA_character_
      ),
      levels = ASSET_LABELS2
    )]

    # Also construct assets_cat2 to match Stata labels (for downstream merges)
    cat2_labels <- c(
      "Assets < 10 million",
      "10M up to 50M",
      "50M through 100M",
      "Over 100M through 500M",
      "Over 500M through 1B",
      "Over 1B through 5B",
      "Over 5B through 10B",
      "Over 10B"
    )
    cr[, assets_cat2 := factor(
      fifelse(
        get(asset_col) > 0,
        cat2_labels[findInterval(get(asset_col), ASSET_BREAKS2,
                                  rightmost.closed = FALSE)],
        NA_character_
      ),
      levels = cat2_labels
    )]
    cat("  Also constructed assets_cat2 (Stata-label format) for compatibility\n")
    log_checkpoint("CK-12b: Asset column found", "INFO",
                   sprintf("Constructed from %s", asset_col))

  } else {
    cat("  WARNING: Neither assets_cat2, assets_tot nor acct_010 found\n")
    cr[, asset_tier  := factor(NA_character_, levels = ASSET_LABELS2)]
    cr[, assets_cat2 := factor(NA_character_)]
    log_checkpoint("CK-12b: Asset column found", "WARN",
                   "No asset column -- tiers set to NA")
  }
}

# ── Final tier distribution (regardless of path taken) ───────────────────────
tier_dist <- cr[, .N, by = asset_tier][order(asset_tier)]
tier_dist[, pct := round(N / nrow(cr) * 100, 1)]
n_na_tier <- cr[is.na(asset_tier), .N]

cat("\n  Final asset_tier distribution (T-codes for regression):\n")
cat(sprintf("  %-22s  %10s  %8s  %8s  %s\n",
            "Tier", "N (CU-qtrs)", "Pct", "Cum.", "Label"))
cat(sprintf("  %s\n", paste(rep("-", 74), collapse = "")))

tier_labels_display <- c(
  T1_under10M  = "< $10M",
  T2_10to50M   = "$10M-$50M",
  T3_50to100M  = "$50M-$100M",
  T4_100to500M = "$100M-$500M",
  T5_500Mto1B  = "$500M-$1B",
  T6_1Bto5B    = "$1B-$5B",
  T7_5Bto10B   = "$5B-$10B",
  T8_over10B   = "> $10B"
)

cum_pct <- 0
for (i in seq_len(nrow(tier_dist))) {
  tier_nm <- as.character(tier_dist$asset_tier[i])
  if (!is.na(tier_nm)) {
    cum_pct <- cum_pct + tier_dist$pct[i]
    cat(sprintf("  %-22s  %10s  %7.1f%%  %7.1f%%  %s\n",
                tier_nm,
                format(tier_dist$N[i], big.mark = ","),
                tier_dist$pct[i],
                cum_pct,
                tier_labels_display[tier_nm]))
  }
}
if (n_na_tier > 0)
  cat(sprintf("  %-22s  %10s  %7.1f%%  (zero/NA assets)\n",
              "NA",
              format(n_na_tier, big.mark = ","),
              n_na_tier / nrow(cr) * 100))

# ── Validation: cross-check against Stata counts where visible ───────────────
# From OCE_combined Stata tab assets_cat2 (616,971 total, 2005Q1+):
#   T1_under10M  : 239,275  (38.78%)
#   T2_10to50M   : 192,470  (31.20%)
#   T3_50to100M  :  65,604  (10.63%)
#   T4_100to500M :  86,153  (13.96%)
#   T5_500Mto1B  :  17,363   (2.81%)
#   T6_1Bto5B    :  14,618   (2.37%)
#   T7_5Bto10B   :   1,075   (0.17%)
#   T8_over10B   :     413   (0.07%)
#   Total        : 616,971
cat("\n  Validation against Stata tab assets_cat2 (from OCE_combined Stata file):\n")
stata_ref <- data.table(
  tier    = ASSET_LABELS2,
  n_stata = c(239275L, 192470L, 65604L, 86153L, 17363L, 14618L, 1075L, 413L),
  pct_stata = c(38.78, 31.20, 10.63, 13.96, 2.81, 2.37, 0.17, 0.07)
)
cat(sprintf("  %-22s  %10s  %8s  %10s  %8s  %s\n",
            "Tier", "N (data)", "Pct", "N (Stata)", "Pct", "Match"))
cat(sprintf("  %s\n", paste(rep("-", 80), collapse = "")))
for (i in seq_len(nrow(stata_ref))) {
  tn    <- stata_ref$tier[i]
  row_d <- tier_dist[asset_tier == tn]
  n_d   <- if (nrow(row_d) > 0) row_d$N[1] else 0L
  pct_d <- if (nrow(row_d) > 0) row_d$pct[1] else 0
  diff  <- abs(pct_d - stata_ref$pct_stata[i])
  flag  <- if (diff > 2.0) "  <- DIFF" else if (diff > 0.5) "  <- NOTE" else "  OK"
  cat(sprintf("  %-22s  %10s  %7.1f%%  %10s  %7.2f%%%s\n",
              tn,
              format(n_d, big.mark = ","),
              pct_d,
              format(stata_ref$n_stata[i], big.mark = ","),
              stata_ref$pct_stata[i],
              flag))
}

# Overall size check
n_within_2pct <- sum(abs(tier_dist$pct -
  stata_ref[match(tier_dist$asset_tier, tier_dist$asset_tier), pct_stata]) < 2,
  na.rm = TRUE)

log_checkpoint("CK-12: Asset tiers (assets_cat2)", "PASS",
               sprintf("8 tiers constructed; %d NA-tier rows; Stata ref validation printed",
                       n_na_tier))

# ============================================================================
# SECTION 6: CONSTRUCTED OUTCOME VARIABLES
# ============================================================================
section("6", "Constructed Outcome Variables")

cat("  Variables computed here: YoY growth rates, pll_rate, level ratios\n")
cat("  All YoY rates use lag-4 within join_number\n")
cat("  Data is pre-sorted by (join_number, year, quarter)\n\n")

OUTLIER_RECORDS <- list()

# Helper: apply winsorize and print the audit line
apply_winsor <- function(dt, col) {
  x   <- dt[[col]]
  idx <- !is.na(x)
  w   <- winsorize_audit(x[idx])
  set(dt, which(idx), col, w$x)
  cat(sprintf("  Winsorize %-26s: lo=%9.5f  hi=%9.5f  clipped lo=%d  hi=%d  (%.3f%% of valid)\n",
              col, w$lo, w$hi, w$n_lo, w$n_hi, w$pct))
  invisible(w)
}

# --------------------------------------------------------------------------
# 6a. insured_share_growth -- YoY change in insured shares
# CRITICAL: #1 outcome by SHAP importance (3x next-highest, v1 finding)
# --------------------------------------------------------------------------
cat("  [6a] insured_share_growth\n")
cat("       Source  : insured_tot  |  YoY % change (lag-4 within CU)\n")
cat("       CRITICAL: #1 SHAP importance outcome in v1 (3x next-highest)\n\n")

if ("insured_tot" %in% names(cr)) {
  describe_var(cr[["insured_tot"]], "insured_tot", "raw level")
  cr[, insured_share_growth := cu_yoy(cr, "insured_tot")]
  n_valid_pre <- sum(!is.na(cr$insured_share_growth))
  cat(sprintf("  Valid obs (non-NA) after YoY construction: %s\n",
              format(n_valid_pre, big.mark = ",")))
  describe_var(cr$insured_share_growth, "insured_share_growth", "pre-winsor")
  w <- apply_winsor(cr, "insured_share_growth")
  OUTLIER_RECORDS[["insured_share_growth"]] <- data.table(
    variable = "insured_share_growth", method = "winsor_1_99",
    lo_boundary = w$lo, hi_boundary = w$hi,
    n_clipped_lo = w$n_lo, n_clipped_hi = w$n_hi, pct_affected = w$pct
  )
  log_checkpoint("CK-13: insured_share_growth", "PASS",
                 sprintf("N=%s valid; clipped=%d",
                         format(n_valid_pre, big.mark = ","), w$n_lo + w$n_hi))
} else {
  cat("  SKIP: 'insured_tot' not found in call_report.rds\n\n")
  cr[, insured_share_growth := NA_real_]
  log_checkpoint("CK-13: insured_share_growth", "WARN",
                 "'insured_tot' not found")
}

# --------------------------------------------------------------------------
# 6b. dep_growth_yoy -- YoY change in total deposits (acct_018)
# --------------------------------------------------------------------------
cat("\n  [6b] dep_growth_yoy\n")
cat("       Source  : acct_018  |  YoY % change\n\n")

if ("acct_018" %in% names(cr)) {
  cr[, dep_growth_yoy := cu_yoy(cr, "acct_018")]
  describe_var(cr$dep_growth_yoy, "dep_growth_yoy", "pre-winsor")
  w <- apply_winsor(cr, "dep_growth_yoy")
  OUTLIER_RECORDS[["dep_growth_yoy"]] <- data.table(
    variable = "dep_growth_yoy", method = "winsor_1_99",
    lo_boundary = w$lo, hi_boundary = w$hi,
    n_clipped_lo = w$n_lo, n_clipped_hi = w$n_hi, pct_affected = w$pct
  )
  log_checkpoint("CK-14: dep_growth_yoy", "PASS",
                 sprintf("clipped=%d", w$n_lo + w$n_hi))
} else {
  cr[, dep_growth_yoy := NA_real_]
  log_checkpoint("CK-14: dep_growth_yoy", "WARN", "acct_018 not found")
}

# --------------------------------------------------------------------------
# 6c. cert_growth_yoy -- YoY change in certificate shares
# --------------------------------------------------------------------------
cat("\n  [6c] cert_growth_yoy\n")
cat("       Source  : dep_shrcert  |  YoY % change\n\n")

if ("dep_shrcert" %in% names(cr)) {
  cr[, cert_growth_yoy := cu_yoy(cr, "dep_shrcert")]
  describe_var(cr$cert_growth_yoy, "cert_growth_yoy", "pre-winsor")
  w <- apply_winsor(cr, "cert_growth_yoy")
  OUTLIER_RECORDS[["cert_growth_yoy"]] <- data.table(
    variable = "cert_growth_yoy", method = "winsor_1_99",
    lo_boundary = w$lo, hi_boundary = w$hi,
    n_clipped_lo = w$n_lo, n_clipped_hi = w$n_hi, pct_affected = w$pct
  )
  log_checkpoint("CK-15: cert_growth_yoy", "PASS",
                 sprintf("clipped=%d", w$n_lo + w$n_hi))
} else {
  cr[, cert_growth_yoy := NA_real_]
  log_checkpoint("CK-15: cert_growth_yoy", "INFO", "dep_shrcert not found -- skipped")
}

# --------------------------------------------------------------------------
# 6d. member_growth_yoy -- YoY change in membership
# Economic note: M&A events create structural breaks (sudden large
# positives/negatives); winsorisation at 1/99 removes M&A noise.
# --------------------------------------------------------------------------
cat("\n  [6d] member_growth_yoy\n")
cat("       Source  : 'members' (fallback: acct_083, number_of_members, tot_members)\n")
cat("       Note    : M&A events create outliers; winsorisation removes M&A noise\n\n")

mem_col <- intersect(c("members", "acct_083", "number_of_members",
                        "tot_members"), names(cr))[1]
if (!is.na(mem_col)) {
  cat(sprintf("  Column used: %s\n", mem_col))
  describe_var(cr[[mem_col]], mem_col, "raw level")
  cr[, member_growth_yoy := cu_yoy(cr, mem_col)]
  describe_var(cr$member_growth_yoy, "member_growth_yoy", "pre-winsor")
  w <- apply_winsor(cr, "member_growth_yoy")
  OUTLIER_RECORDS[["member_growth_yoy"]] <- data.table(
    variable = "member_growth_yoy", method = "winsor_1_99",
    lo_boundary = w$lo, hi_boundary = w$hi,
    n_clipped_lo = w$n_lo, n_clipped_hi = w$n_hi, pct_affected = w$pct
  )
  log_checkpoint("CK-16: member_growth_yoy", "PASS",
                 sprintf("source=%s; clipped=%d", mem_col, w$n_lo + w$n_hi))
} else {
  cat("  SKIP: no membership column found\n")
  cr[, member_growth_yoy := NA_real_]
  log_checkpoint("CK-16: member_growth_yoy", "WARN", "No membership column found")
}

# --------------------------------------------------------------------------
# 6e. cert_share = dep_shrcert / acct_018  (proportion 0-1, bounded)
# --------------------------------------------------------------------------
cat("\n  [6e] cert_share\n")
cat("       Formula : dep_shrcert / acct_018  |  Bounded [0, 1]\n")
cat("       Logic   : Hard bounds, not winsorisation (outside [0,1] = data error)\n\n")

if (all(c("dep_shrcert", "acct_018") %in% names(cr))) {
  n_zero_dep <- cr[acct_018 <= 0, .N]
  cat(sprintf("  acct_018 <= 0 (zero denominator): %s rows -> NA\n",
              format(n_zero_dep, big.mark = ",")))
  cr[, cert_share := fifelse(
    !is.na(acct_018) & is.finite(acct_018) & acct_018 > 0,
    dep_shrcert / acct_018,
    NA_real_
  )]
  n_gt1 <- cr[!is.na(cert_share) & cert_share > 1, .N]
  n_lt0 <- cr[!is.na(cert_share) & cert_share < 0, .N]
  cr[!is.na(cert_share), cert_share := pmin(pmax(cert_share, 0), 1)]
  cat(sprintf("  Hard-bounded at [0,1]: %d above 1, %d below 0\n", n_gt1, n_lt0))
  describe_var(cr$cert_share, "cert_share", "post-bounding")
  OUTLIER_RECORDS[["cert_share"]] <- data.table(
    variable = "cert_share", method = "hard_bound_0_1",
    lo_boundary = 0, hi_boundary = 1,
    n_clipped_lo = n_lt0, n_clipped_hi = n_gt1, pct_affected = NA_real_
  )
  log_checkpoint("CK-17: cert_share", "PASS",
                 sprintf("hard-bounded: %d obs capped", n_gt1 + n_lt0))
} else {
  cr[, cert_share := NA_real_]
  log_checkpoint("CK-17: cert_share", "INFO", "dep_shrcert or acct_018 missing")
}

# --------------------------------------------------------------------------
# 6f. loan_to_share = lns_tot / dep_tot  (bounded 0-2)
# dep_tot preferred; acct_018 fallback
# --------------------------------------------------------------------------
cat("\n  [6f] loan_to_share\n")
cat("       Formula : lns_tot / dep_tot  (fallback denominator: acct_018)\n")
cat("       Bounded : [0, 2]  -- ratios above 2 or below 0 are data errors\n\n")

dep_denom <- intersect(c("dep_tot", "acct_018"), names(cr))[1]
if (!is.na(dep_denom) && "lns_tot" %in% names(cr)) {
  n_zero_dd <- cr[get(dep_denom) <= 0, .N]
  cat(sprintf("  Denominator: %s  |  %d zero/neg rows -> NA\n",
              dep_denom, n_zero_dd))
  cr[, loan_to_share := fifelse(
    !is.na(get(dep_denom)) & is.finite(get(dep_denom)) & get(dep_denom) > 0,
    lns_tot / get(dep_denom),
    NA_real_
  )]
  n_hi_lts <- cr[!is.na(loan_to_share) & loan_to_share > 2, .N]
  n_lo_lts <- cr[!is.na(loan_to_share) & loan_to_share < 0, .N]
  cr[!is.na(loan_to_share), loan_to_share := pmin(pmax(loan_to_share, 0), 2)]
  cat(sprintf("  Hard-bounded at [0,2]: %d above 2, %d below 0\n",
              n_hi_lts, n_lo_lts))
  describe_var(cr$loan_to_share, "loan_to_share", "post-bounding")
  OUTLIER_RECORDS[["loan_to_share"]] <- data.table(
    variable = "loan_to_share", method = "hard_bound_0_2",
    lo_boundary = 0, hi_boundary = 2,
    n_clipped_lo = n_lo_lts, n_clipped_hi = n_hi_lts, pct_affected = NA_real_
  )
  log_checkpoint("CK-18: loan_to_share", "PASS",
                 sprintf("denom=%s; capped=%d", dep_denom, n_hi_lts + n_lo_lts))
} else {
  cr[, loan_to_share := NA_real_]
  log_checkpoint("CK-18: loan_to_share", "WARN",
                 "lns_tot or deposit denominator not found")
}

# --------------------------------------------------------------------------
# 6g. nim_spread = yldavgloans - costfds
# --------------------------------------------------------------------------
cat("\n  [6g] nim_spread\n")
cat("       Formula : yldavgloans - costfds\n\n")

if (all(c("yldavgloans", "costfds") %in% names(cr))) {
  cr[, nim_spread := fifelse(
    !is.na(yldavgloans) & !is.na(costfds) &
      is.finite(yldavgloans) & is.finite(costfds),
    yldavgloans - costfds,
    NA_real_
  )]
  describe_var(cr$nim_spread, "nim_spread", "computed")
  log_checkpoint("CK-19: nim_spread", "PASS", "yldavgloans - costfds")
} else {
  miss_nim <- setdiff(c("yldavgloans", "costfds"), names(cr))
  cr[, nim_spread := NA_real_]
  log_checkpoint("CK-19: nim_spread", "INFO",
                 sprintf("Missing: %s", paste(miss_nim, collapse = ", ")))
}

# --------------------------------------------------------------------------
# 6h. pll_rate = pll / avg(lns_tot) * 100
# CRITICAL construction confirmed in v1 Script 04d
# Bounded [-2%, +5%] THEN winsorized within bounds
# --------------------------------------------------------------------------
cat("\n  [6h] pll_rate -- Provision for Loan Loss Rate\n")
cat("       Formula  : pll / [(lns_tot_t + lns_tot_{t-1}) / 2]  x 100\n")
cat("       Expressed: % of average loans (quarterly, not annualized)\n")
cat("       CRITICAL : avg loan denominator confirmed from v1 Script 04d\n")
cat("       Step 1   : Hard bounds [-2%, +5%]  (outside = data error)\n")
cat("       Step 2   : Winsorize within bounds at 1/99\n\n")
cat("       Economic note:\n")
cat("         Rising oil  => energy-sector stress => management raises provisions proactively\n")
cat("         Falling oil => regional unemployment => PLL spikes in oil-state CUs\n")
cat("         Post-2015   => pll_rate diverges from dq_rate (forward- vs backward-looking)\n\n")

if ("pll" %in% names(cr) && "lns_tot" %in% names(cr)) {
  describe_var(cr[["pll"]],     "pll",     "raw level")
  describe_var(cr[["lns_tot"]], "lns_tot", "raw level")

  cr[, lns_avg := (lns_tot + shift(lns_tot, 1L)) / 2, by = join_number]

  n_avg_na   <- cr[is.na(lns_avg), .N]
  n_avg_zero <- cr[!is.na(lns_avg) & lns_avg <= 0, .N]
  n_neg_pll  <- sum(cr$pll < 0, na.rm = TRUE)

  cat(sprintf("  lns_avg NA (first obs per CU)    : %s rows\n",
              format(n_avg_na, big.mark = ",")))
  cat(sprintf("  lns_avg <= 0 (zero avg loans)    : %s rows -> NA\n",
              format(n_avg_zero, big.mark = ",")))
  cat(sprintf("  Negative pll (recoveries > prov) : %d obs (valid, not an error)\n\n",
              n_neg_pll))

  cr[, pll_rate := fifelse(
    !is.na(pll) & !is.na(lns_avg) & is.finite(lns_avg) & lns_avg > 0,
    pll / lns_avg * 100,
    NA_real_
  )]

  # Step 1: hard bounds
  n_below_minus2 <- cr[!is.na(pll_rate) & pll_rate < -2, .N]
  n_above_5      <- cr[!is.na(pll_rate) & pll_rate >  5, .N]
  cr[!is.na(pll_rate), pll_rate := pmin(pmax(pll_rate, -2), 5)]
  cat(sprintf("  Hard bounds [-2%%, +5%%]:\n"))
  cat(sprintf("    Below -2%%  : %d obs hard-capped\n", n_below_minus2))
  cat(sprintf("    Above +5%%  : %d obs hard-capped\n\n", n_above_5))

  # Step 2: winsorize within bounds
  describe_var(cr$pll_rate, "pll_rate", "pre-winsor (post hard-bound)")
  w <- apply_winsor(cr, "pll_rate")

  cat(sprintf("\n  pll_rate FINAL:\n"))
  cat(sprintf("    N valid : %s\n", format(sum(!is.na(cr$pll_rate)), big.mark = ",")))
  cat(sprintf("    Mean    : %.4f%%\n", mean(cr$pll_rate, na.rm = TRUE)))
  cat(sprintf("    p99     : %.4f%%\n", quantile(cr$pll_rate, 0.99, na.rm = TRUE)))

  OUTLIER_RECORDS[["pll_rate"]] <- data.table(
    variable = "pll_rate", method = "hard_bound_then_winsor_1_99",
    lo_boundary = w$lo, hi_boundary = w$hi,
    n_clipped_lo = w$n_lo, n_clipped_hi = w$n_hi, pct_affected = w$pct,
    n_hardbound_lo = n_below_minus2, n_hardbound_hi = n_above_5
  )

  cr[, lns_avg := NULL]  # drop intermediate column

  log_checkpoint("CK-20: pll_rate", "PASS",
                 sprintf("hard-capped=%d; winsorized=%d",
                         n_below_minus2 + n_above_5, w$n_lo + w$n_hi))
} else {
  miss_pll <- setdiff(c("pll", "lns_tot"), names(cr))
  cat(sprintf("  SKIP: missing columns: %s\n", paste(miss_pll, collapse = ", ")))
  cr[, pll_rate := NA_real_]
  log_checkpoint("CK-20: pll_rate", "WARN",
                 sprintf("Missing: %s", paste(miss_pll, collapse = ", ")))
}

# --------------------------------------------------------------------------
# 6i. pll_per_loan = pll / lns_tot_n  (per-loan count, size-neutral)
# --------------------------------------------------------------------------
cat("\n  [6i] pll_per_loan\n")
cat("       Formula : pll / lns_tot_n  ($ provision per loan count)\n")
cat("       Use     : Size-neutral complement to pll_rate\n\n")

if ("pll" %in% names(cr) && "lns_tot_n" %in% names(cr)) {
  n_zero_n <- cr[lns_tot_n <= 0, .N]
  cat(sprintf("  lns_tot_n <= 0: %d rows -> NA\n", n_zero_n))
  cr[, pll_per_loan := fifelse(
    !is.na(pll) & !is.na(lns_tot_n) & is.finite(lns_tot_n) & lns_tot_n > 0,
    pll / lns_tot_n,
    NA_real_
  )]
  describe_var(cr$pll_per_loan, "pll_per_loan", "pre-winsor")
  w <- apply_winsor(cr, "pll_per_loan")
  OUTLIER_RECORDS[["pll_per_loan"]] <- data.table(
    variable = "pll_per_loan", method = "winsor_1_99",
    lo_boundary = w$lo, hi_boundary = w$hi,
    n_clipped_lo = w$n_lo, n_clipped_hi = w$n_hi, pct_affected = w$pct
  )
  log_checkpoint("CK-21: pll_per_loan", "PASS",
                 sprintf("clipped=%d", w$n_lo + w$n_hi))
} else {
  cr[, pll_per_loan := NA_real_]
  log_checkpoint("CK-21: pll_per_loan", "INFO", "pll or lns_tot_n not found")
}

# --------------------------------------------------------------------------
# 6j. Winsorize pre-computed direct ratio variables
# --------------------------------------------------------------------------
cat("\n  [6j] Winsorizing pre-computed direct ratio variables\n\n")

WINSOR_DIRECT <- intersect(
  c("netintmrg", "pcanetworth", "costfds", "roa",
    "dq_rate", "chg_tot_lns_ratio"),
  names(cr)
)

cat(sprintf("  %-22s  %8s  %8s  %6s  %6s  %8s\n",
            "Variable", "p1 (lo)", "p99 (hi)", "N < p1", "N > p99", "% affect"))
cat(sprintf("  %s\n", paste(rep("-", 68), collapse = "")))

for (v in WINSOR_DIRECT) {
  if (all(is.na(cr[[v]]))) {
    cat(sprintf("  %-22s  ALL MISSING\n", v))
    next
  }
  x_pre <- cr[[v]]
  w     <- winsorize_audit(x_pre)
  cr[, (v) := w$x]
  cat(sprintf("  %-22s  %8.5f  %8.5f  %6d  %6d  %7.3f%%\n",
              v, w$lo, w$hi, w$n_lo, w$n_hi, w$pct))
  OUTLIER_RECORDS[[paste0(v, "_winsor")]] <- data.table(
    variable = v, method = "winsor_1_99",
    lo_boundary = w$lo, hi_boundary = w$hi,
    n_clipped_lo = w$n_lo, n_clipped_hi = w$n_hi, pct_affected = w$pct
  )
}

log_checkpoint("CK-22: Direct ratio winsorization", "PASS",
               sprintf("%d variables processed", length(WINSOR_DIRECT)))

# ============================================================================
# SECTION 7: TIME FLAGS (post2015, first_obs, cecl)
# ============================================================================
section("7", "Time Flags")

cr[, post2015 := as.integer(cal_date >= as.Date("2015-01-01"))]

cat(sprintf("  post2015 == 0 (pre-shale)  : %s CU-quarters\n",
            format(cr[post2015 == 0, .N], big.mark = ",")))
cat(sprintf("  post2015 == 1 (post-shale) : %s CU-quarters\n",
            format(cr[post2015 == 1, .N], big.mark = ",")))
cat("  ANCHOR: 2015Q1 structural break reverses all transmission signs (v1)\n\n")

if (cr[post2015 == 0, .N] == 0 || cr[post2015 == 1, .N] == 0) {
  log_checkpoint("CK-23a: 2015Q1 splits data", "WARN",
                 "One side of the 2015Q1 split is empty -- check date range")
} else {
  log_checkpoint("CK-23a: 2015Q1 splits data", "PASS",
                 sprintf("pre=%s, post=%s",
                         format(cr[post2015==0,.N], big.mark=","),
                         format(cr[post2015==1,.N], big.mark=",")))
}

# First-obs flag (lag unavailable for AR models)
lag_check_vars <- intersect(c("insured_share_growth", "dq_rate"), names(cr))
if (length(lag_check_vars) > 0) {
  for (v in lag_check_vars) {
    lag_nm <- paste0(v, "_lag1")
    cr[, (lag_nm) := shift(get(v), 1L, type = "lag"), by = join_number]
  }
  cr[, first_obs_flag := as.integer(
    Reduce(`|`, lapply(paste0(lag_check_vars, "_lag1"),
                       function(nm) is.na(cr[[nm]])))
  )]
} else {
  cr[, first_obs_flag := 0L]
}

n_first <- cr[first_obs_flag == 1, .N]
cat(sprintf("  first_obs_flag == 1 : %s rows (no lag -- exclude from AR models)\n",
            format(n_first, big.mark = ",")))

# CECL flags
cr[, cecl_flag := as.integer(
  (cal_date >= as.Date("2020-01-01") & cal_date <= as.Date("2020-03-31")) |
    (cal_date >= as.Date("2023-01-01") & cal_date <= as.Date("2023-03-31"))
)]
n_cecl <- cr[cecl_flag == 1, .N]
cat(sprintf("  cecl_flag == 1      : %s rows (2020Q1 large-CU and 2023Q1 FICU CECL)\n",
            format(n_cecl, big.mark = ",")))
cat("  Note: Exclude cecl_flag==1 from PLL/charge-off regressions\n")

log_checkpoint("CK-23: Time flags", "PASS",
               sprintf("post2015=%d/%d; first_obs=%d; cecl=%d",
                       cr[post2015==1,.N], nrow(cr), n_first, n_cecl))

# ============================================================================
# SECTION 8: PRELIMINARY OIL-STATE FLAG
# ============================================================================
section("8", "Preliminary Oil-State Flag")

cat("  This is a BINARY PRELIMINARY flag only.\n")
cat("  It is SUPERSEDED by oil_exposure.rds from 01b_oil_exposure_v2.R\n")
cat("  which provides QCEW-based continuous exposure weights.\n\n")
cat(sprintf("  States: %s\n\n", paste(OIL_STATES_PRELIM, collapse = ", ")))

if (!is.na(state_col)) {
  cr[, oil_state_prelim := toupper(get(state_col)) %in% OIL_STATES_PRELIM]
  n_oil <- cr[oil_state_prelim == TRUE, .N]
  n_non <- cr[oil_state_prelim == FALSE, .N]
  cat(sprintf("  oil_state_prelim == TRUE  : %s CU-quarters\n",
              format(n_oil, big.mark = ",")))
  cat(sprintf("  oil_state_prelim == FALSE : %s CU-quarters\n",
              format(n_non, big.mark = ",")))
  log_checkpoint("CK-24: Preliminary oil flag", "PASS",
                 sprintf("oil=%s, non-oil=%s",
                         format(n_oil, big.mark = ","),
                         format(n_non, big.mark = ",")))
} else {
  cr[, oil_state_prelim := NA]
  log_checkpoint("CK-24: Preliminary oil flag", "WARN",
                 "No state column -- oil_state_prelim set to NA")
}

# ============================================================================
# SECTION 9: SAVE CALL_CLEAN (FIRST SAVE -- BEFORE MACRO MERGE)
# ============================================================================
section("9", "Save call_clean.rds (CU panel, no macro yet)")

setorderv(cr, c("join_number", "year", "quarter"))
saveRDS(cr, "Data/call_clean.rds")
cat(sprintf("  Saved: Data/call_clean.rds  |  %s rows x %d cols  |  %.1f MB\n",
            format(nrow(cr), big.mark = ","), ncol(cr),
            file.size("Data/call_clean.rds") / 1e6))
log_checkpoint("CK-25: call_clean.rds first save", "PASS",
               sprintf("%s rows, %d cols", format(nrow(cr), big.mark = ","), ncol(cr)))

# ============================================================================
# SECTION 10: FRB CCAR SCENARIO LOADER
# ============================================================================
section("10", "FRB CCAR Scenario Loader")

cat("  Files : FRB_Baseline_2026.xlsx  and  FRB_Severely_Adverse_2026.xlsx\n")
cat("  Format: Date column in YYYY.Q format (e.g. 2005.1 = 2005Q1)\n\n")

parse_frb_date <- function(x, path = "") {
  x_num  <- suppressWarnings(as.numeric(as.character(x)))
  yr_p   <- floor(x_num)
  qn_p   <- round((x_num - yr_p) * 10)
  valid  <- !is.na(x_num) & yr_p >= 1900 & yr_p <= 2100 &
    qn_p >= 1 & qn_p <= 4

  if (mean(valid, na.rm = TRUE) > 0.8) {
    mo   <- c("1"="01","2"="04","3"="07","4"="10")[as.character(qn_p)]
    return(list(
      dates = as.Date(paste(yr_p, mo, "01", sep = "-")),
      yr = yr_p, qn = qn_p
    ))
  }
  # Excel serial date fallback
  if (is.numeric(x)) {
    d <- suppressWarnings(as.Date(x, origin = "1899-12-30"))
    if (mean(!is.na(d)) > 0.8)
      return(list(dates = d, yr = year(d), qn = quarter(d)))
  }
  stop(paste("Cannot parse Date column in", path,
             "\n  Sample:", paste(head(as.character(x), 5), collapse = ", ")))
}

load_frb <- function(path, prefix) {
  cat(sprintf("\n  Loading: %s  [prefix: %s]\n", basename(path), prefix))

  if (!file.exists(path)) {
    log_checkpoint(paste("FRB file:", basename(path)), "FAIL",
                   paste("Not found:", path))
  }

  raw <- tryCatch(
    as.data.table(read_excel(path, .name_repair = "unique")),
    error = function(e) stop("Cannot read ", path, ": ", e$message)
  )

  raw <- raw[!apply(raw, 1, function(r)
    all(is.na(r) | !suppressWarnings(!is.na(as.numeric(r[!is.na(r)])))))]

  cat(sprintf("  Dims after cleaning: %d rows x %d cols\n", nrow(raw), ncol(raw)))

  date_idx <- which(tolower(names(raw)) == "date")[1]
  if (is.na(date_idx)) {
    date_idx <- which(vapply(raw, function(x) {
      x2 <- suppressWarnings(as.numeric(as.character(x)))
      mean(!is.na(x2) & x2 > 1990 & x2 < 2100) > 0.4
    }, logical(1)))[1]
  }
  if (is.na(date_idx))
    log_checkpoint(paste("Date col in", basename(path)), "FAIL", "Not found")

  date_col <- names(raw)[date_idx]
  cat(sprintf("  Date column: '%s'\n", date_col))

  parsed  <- parse_frb_date(raw[[date_col]], path)
  dv_num  <- suppressWarnings(as.numeric(as.character(raw[[date_col]])))
  yr_src  <- floor(dv_num)
  qn_src  <- round((dv_num - yr_src) * 10)
  use_src <- !is.na(dv_num) & yr_src >= 1900 & qn_src >= 1 & qn_src <= 4

  raw[, date_parsed := parsed$dates]
  raw[, year    := fifelse(use_src, as.integer(yr_src), year(parsed$dates))]
  raw[, quarter := fifelse(use_src, as.integer(qn_src), quarter(parsed$dates))]
  raw[, yyyyqq  := year * 100L + quarter]
  raw <- raw[!is.na(date_parsed) & year >= 2004]

  id_cols    <- c(date_col, "date_parsed", "year", "quarter", "yyyyqq")
  macro_cols <- setdiff(names(raw), id_cols)
  macro_cols <- macro_cols[vapply(raw[, ..macro_cols], is.numeric, logical(1))]
  macro_cols <- macro_cols[!grepl("^\\.\\.", macro_cols)]
  cat(sprintf("  Macro columns identified: %d\n", length(macro_cols)))

  wide <- raw[, c(.(year = first(year), quarter = first(quarter)),
                  lapply(.SD, function(x) mean(as.numeric(x), na.rm = TRUE))),
              by = yyyyqq, .SDcols = macro_cols]
  wide <- wide[order(yyyyqq)]

  new_nms  <- paste0(prefix, tolower(macro_cols))
  dupe_nm  <- duplicated(new_nms)
  if (any(dupe_nm)) new_nms[dupe_nm] <- paste0(new_nms[dupe_nm], "_v2")
  setnames(wide, macro_cols, new_nms)

  key     <- paste0(prefix, c("pbrent","lurc","pcpi","pcpixfe","rmtg",
                               "phpi","uypsav","ypds","gdps","rff"))
  found   <- intersect(key, names(wide))
  missing <- setdiff(key, names(wide))
  cat(sprintf("  Key vars found   : %s\n", paste(found, collapse = ", ")))
  if (length(missing) > 0)
    cat(sprintf("  Key vars MISSING : %s\n", paste(missing, collapse = ", ")))
  cat(sprintf("  Quarters: %dQ%d -- %dQ%d  (%d rows)\n",
              wide[1, year], wide[1, quarter],
              wide[.N, year], wide[.N, quarter], nrow(wide)))

  wide[]
}

macro_base   <- load_frb("Data/FRB_Baseline_2026.xlsx",         "macro_base_")
macro_severe <- load_frb("Data/FRB_Severely_Adverse_2026.xlsx", "macro_severe_")

saveRDS(macro_base,   "Data/macro_base.rds")
saveRDS(macro_severe, "Data/macro_severe.rds")
cat(sprintf("\n  Saved: macro_base.rds (%d cols)  |  macro_severe.rds (%d cols)\n",
            ncol(macro_base), ncol(macro_severe)))
log_checkpoint("CK-26: FRB scenarios loaded", "PASS",
               sprintf("base=%d cols, severe=%d cols",
                       ncol(macro_base), ncol(macro_severe)))

# ============================================================================
# SECTION 11: DERIVED MACRO VARIABLES
# ============================================================================
section("11", "Derived Macro Variables")

cat("  Computing: oil transforms, yield curve, real rate, FOMC regime, CPI YoY\n")
cat("  Applied to BOTH macro_base_ and macro_severe_ scenarios\n\n")

`%||%` <- function(a, b) if (!is.null(a)) a else b

add_macro_derived <- function(dt, pfx) {
  gc <- function(nm) {
    col <- paste0(pfx, nm)
    if (col %in% names(dt)) dt[[col]] else NULL
  }
  setorderv(dt, "yyyyqq")

  # Oil price transforms
  pb <- gc("pbrent")
  if (!is.null(pb)) {
    dt[, (paste0(pfx,"yoy_oil"))    := (pb - shift(pb,4)) / shift(pb,4) * 100]
    dt[, (paste0(pfx,"qoq_oil"))    := (pb - shift(pb,1)) / shift(pb,1) * 100]
    dt[, (paste0(pfx,"oil_pos"))    := pmax(get(paste0(pfx,"yoy_oil")), 0, na.rm=TRUE)]
    dt[, (paste0(pfx,"oil_neg"))    := pmin(get(paste0(pfx,"yoy_oil")), 0, na.rm=TRUE)]
    for (k in 1:4)
      dt[, (paste0(pfx,"yoy_oil_lag",k)) := shift(get(paste0(pfx,"yoy_oil")), k)]
    dt[, (paste0(pfx,"oil_rsd4"))   := frollapply(get(paste0(pfx,"qoq_oil")),4,sd,na.rm=TRUE,align="right")]
    dt[, (paste0(pfx,"oil_rmean8")) := frollmean(pb, 8, na.rm=TRUE, align="right")]
    dt[, (paste0(pfx,"oil_cyc"))    := pb - get(paste0(pfx,"oil_rmean8"))]
    yoy_col <- paste0(pfx, "yoy_oil")
    cat(sprintf("  [%s] pbrent YoY oil: mean=%+.1f%%  min=%+.1f%%  max=%+.1f%%\n",
                pfx,
                mean(dt[[yoy_col]], na.rm=TRUE),
                min(dt[[yoy_col]], na.rm=TRUE),
                max(dt[[yoy_col]], na.rm=TRUE)))
    cat(sprintf("  [%s] Oil transforms created: YoY, QoQ, pos/neg, lags 1-4, rsd4, cyc\n", pfx))
  } else {
    cat(sprintf("  WARN [%s] pbrent not found -- oil transforms skipped\n", pfx))
  }

  # Yield curve
  gs10 <- gc("rs10y") %||% gc("gs10")
  gs3m <- gc("rs3m")  %||% gc("gs3m")
  if (!is.null(gs10) && !is.null(gs3m)) {
    dt[, (paste0(pfx,"yield_curve"))     := gs10 - gs3m]
    dt[, (paste0(pfx,"yield_curve_inv")) := as.integer(gs10 - gs3m < 0)]
    yc_inv <- sum(dt[[paste0(pfx,"yield_curve_inv")]], na.rm=TRUE)
    cat(sprintf("  [%s] yield_curve (inverted in %d quarters)\n", pfx, yc_inv))
  }

  # Real rate
  rff  <- gc("rff")  %||% gc("fedfunds")
  pcpi <- gc("pcpi") %||% gc("cpi")
  if (!is.null(rff) && !is.null(pcpi)) {
    dt[, (paste0(pfx,"real_rate")) := rff - pcpi]
    cat(sprintf("  [%s] real_rate (rff - pcpi)\n", pfx))
  }

  # FOMC regime + hike run
  if (!is.null(rff)) {
    chg    <- c(NA, diff(rff))
    regime <- fifelse(chg >  0.10,  1L,
               fifelse(chg < -0.10, -1L, 0L))
    run <- integer(length(regime))
    for (i in seq_along(regime)) {
      if (is.na(regime[i]) || regime[i] == 0L) {
        run[i] <- 0L
      } else if (i == 1 || is.na(regime[i-1]) || regime[i] != regime[i-1]) {
        run[i] <- regime[i]
      } else {
        run[i] <- run[i-1] + regime[i]
      }
    }
    dt[, (paste0(pfx,"fomc_regime")) := regime]
    dt[, (paste0(pfx,"hike_run"))    := run]
    n_hike <- sum(regime ==  1L, na.rm=TRUE)
    n_cut  <- sum(regime == -1L, na.rm=TRUE)
    cat(sprintf("  [%s] fomc_regime: %d hike qtrs, %d cut qtrs\n", pfx, n_hike, n_cut))
  }

  # CPI YoY
  if (!is.null(pcpi)) {
    dt[, (paste0(pfx,"cpi_yoy")) := (pcpi - shift(pcpi,4)) / shift(pcpi,4) * 100]
    cat(sprintf("  [%s] cpi_yoy\n", pfx))
  }

  invisible(dt)
}

add_macro_derived(macro_base,   "macro_base_")
cat("\n")
add_macro_derived(macro_severe, "macro_severe_")

saveRDS(macro_base,   "Data/macro_base.rds")
saveRDS(macro_severe, "Data/macro_severe.rds")
cat("\n  Saved: macro_base.rds and macro_severe.rds (with derived vars)\n")
log_checkpoint("CK-27: Derived macro vars", "PASS",
               "Oil transforms, yield curve, FOMC regime for both scenarios")

# ============================================================================
# SECTION 12: PANEL ASSEMBLY (CU x MACRO)
# ============================================================================
section("12", "Panel Assembly -- Call Report x Macro Scenarios")

cat("  Merge key: yyyyqq (integer: year*100+quarter)\n\n")

macro_merge_fn <- function(macro_dt) {
  drop <- intersect(c("year", "quarter"), names(macro_dt))
  macro_dt[, !drop, with = FALSE]
}

panel_base   <- merge(cr, macro_merge_fn(macro_base),   by = "yyyyqq", all.x = TRUE)
panel_severe <- merge(cr, macro_merge_fn(macro_severe),  by = "yyyyqq", all.x = TRUE)

for (nm in c("macro_base_pbrent", "macro_severe_pbrent")) {
  pnl <- if (grepl("base", nm)) panel_base else panel_severe
  if (nm %in% names(pnl)) {
    v    <- pnl[[nm]]
    n_ok <- sum(!is.na(v))
    cat(sprintf("  %-28s  coverage: %s / %s (%.1f%%)  min=%5.1f  max=%5.1f\n",
                nm,
                format(n_ok, big.mark = ","),
                format(nrow(pnl), big.mark = ","),
                n_ok / nrow(pnl) * 100,
                min(v, na.rm=TRUE), max(v, na.rm=TRUE)))
  }
}

cat(sprintf("\n  panel_base   : %s rows x %d cols\n",
            format(nrow(panel_base), big.mark = ","), ncol(panel_base)))
cat(sprintf("  panel_severe : %s rows x %d cols\n",
            format(nrow(panel_severe), big.mark = ","), ncol(panel_severe)))

log_checkpoint("CK-28: Panel assembly", "PASS",
               sprintf("base=%s rows, severe=%s rows",
                       format(nrow(panel_base), big.mark = ","),
                       format(nrow(panel_severe), big.mark = ",")))

# ============================================================================
# SECTION 13: OIL EXPOSURE MERGE (from 01b_oil_exposure_v2.R)
# ============================================================================
section("13", "Oil Exposure Merge")

cat("  Source: Data/oil_exposure.rds  (built by 01b_oil_exposure_v2.R)\n")
cat("  Join  : state_code x yyyyqq\n\n")

merge_exposure <- function(pnl, exp_dt, sc) {
  if (is.na(sc)) {
    cat("  WARNING: no state column -- skipping exposure merge\n")
    return(pnl)
  }
  exp_data_cols <- setdiff(names(exp_dt), c("state_code", "yyyyqq"))
  stale <- intersect(exp_data_cols, names(pnl))
  if (length(stale)) pnl[, (stale) := NULL]

  pnl[,    (sc) := as.character(get(sc))]
  exp_dt[, state_code := as.character(state_code)]
  pnl[,    yyyyqq := as.integer(yyyyqq)]
  exp_dt[, yyyyqq := as.integer(yyyyqq)]

  out <- merge(pnl, exp_dt,
               by.x = c(sc, "yyyyqq"),
               by.y = c("state_code", "yyyyqq"),
               all.x = TRUE)

  fill_zero <- c("oil_exposure_bin", "oil_exposure_cont",
                 "oil_exposure_bin_1pct", "oil_exposure_bin_3pct",
                 "spillover_exposure", "spillover_exposure_wtd")
  for (col in intersect(fill_zero, names(out)))
    out[is.na(get(col)), (col) := 0]

  if ("oil_exposure_bin" %in% names(out))
    out[, oil_exposure_idx := oil_exposure_bin]

  out[]
}

if (file.exists("Data/oil_exposure.rds")) {
  exp_dt   <- setDT(readRDS("Data/oil_exposure.rds"))
  exp_keep <- intersect(
    c("state_code","yyyyqq","mining_emp_share","oil_exposure_cont",
      "oil_exposure_bin","oil_exposure_bin_1pct","oil_exposure_bin_3pct",
      "oil_exposure_tier","oil_exposure_smooth","oil_bartik_iv",
      "spillover_exposure","spillover_exposure_wtd","cu_group"),
    names(exp_dt)
  )
  exp_merge_dt <- exp_dt[, ..exp_keep]

  cat(sprintf("  oil_exposure.rds: %s rows x %d cols\n",
              format(nrow(exp_merge_dt), big.mark = ","), ncol(exp_merge_dt)))
  cat(sprintf("  Exposure cols: %s\n\n",
              paste(setdiff(exp_keep, c("state_code","yyyyqq")), collapse = ", ")))

  cr           <- merge_exposure(cr,           exp_merge_dt, state_col)
  panel_base   <- merge_exposure(panel_base,   exp_merge_dt, state_col)
  panel_severe <- merge_exposure(panel_severe, exp_merge_dt, state_col)

  cat(sprintf("  cr           : %s rows x %d cols\n",
              format(nrow(cr),           big.mark = ","), ncol(cr)))
  cat(sprintf("  panel_base   : %s rows x %d cols\n",
              format(nrow(panel_base),   big.mark = ","), ncol(panel_base)))
  cat(sprintf("  panel_severe : %s rows x %d cols\n",
              format(nrow(panel_severe), big.mark = ","), ncol(panel_severe)))

  if ("cu_group" %in% names(panel_base)) {
    cat("\n  CU group distribution (panel_base):\n")
    grp <- panel_base[, .N, by = cu_group][order(cu_group)]
    for (i in seq_len(nrow(grp))) {
      cat(sprintf("    %-20s : %s (%.1f%%)\n",
                  as.character(grp$cu_group[i]),
                  format(grp$N[i], big.mark = ","),
                  grp$N[i] / nrow(panel_base) * 100))
    }
  }
  log_checkpoint("CK-29: Oil exposure merged", "PASS",
                 sprintf("From oil_exposure.rds -- %d exposure cols",
                         length(exp_keep) - 2))

} else {
  cat("  WARNING: Data/oil_exposure.rds NOT FOUND\n")
  cat("  Run 01b_oil_exposure_v2.R first, then re-run this script\n")
  cat("  Fallback: Provisional binary flag from hard-coded state list\n\n")

  if (!is.na(state_col)) {
    for (pnl_nm in c("cr", "panel_base", "panel_severe")) {
      pnl <- get(pnl_nm)
      pnl[, oil_exposure_idx  := as.integer(toupper(get(state_col)) %in% OIL_STATES_PRELIM)]
      pnl[, oil_exposure_cont := as.numeric(oil_exposure_idx)]
      pnl[, spillover_exposure := 0]
      pnl[, cu_group := fifelse(oil_exposure_idx == 1L, "Direct", "Indirect")]
      assign(pnl_nm, pnl)
    }
    cat("  Provisional oil_exposure_idx set from hard-coded state list\n")
  }
  log_checkpoint("CK-29: Oil exposure merged", "WARN",
                 "oil_exposure.rds not found -- provisional binary flag used")
}

# ============================================================================
# SECTION 14: INTERACTION TERMS
# ============================================================================
section("14", "Interaction Terms")

cat("  Direct channel  : oil_x_brent         = oil_exposure_cont x macro_*_yoy_oil\n")
cat("  Binary direct   : oil_x_brent_bin      = oil_exposure_bin  x macro_*_yoy_oil\n")
cat("  Bartik IV       : bartik_x_brent       = oil_bartik_iv     x macro_*_yoy_oil\n")
cat("  Indirect channel: spillover_x_brent    = spillover_exposure x macro_*_yoy_oil\n")
cat("  FOMC channel    : fomc_x_brent         = fomc_regime        x macro_*_yoy_oil\n\n")

add_interactions <- function(pnl, oil_yoy_col) {
  if (!oil_yoy_col %in% names(pnl)) {
    cat(sprintf("  SKIP interactions: '%s' not in panel\n", oil_yoy_col))
    return(invisible(pnl))
  }
  oy <- pnl[[oil_yoy_col]]
  cat(sprintf("  Oil YoY col: %s  |  mean=%+.2f%%  min=%+.2f%%  max=%+.2f%%\n",
              oil_yoy_col,
              mean(oy, na.rm=TRUE), min(oy, na.rm=TRUE), max(oy, na.rm=TRUE)))

  if ("oil_exposure_cont" %in% names(pnl)) {
    pnl[, oil_x_brent := oil_exposure_cont * oy]
    cat("    oil_x_brent created\n")
  }
  if ("oil_exposure_bin" %in% names(pnl)) {
    pnl[, oil_x_brent_bin := oil_exposure_bin * oy]
    cat("    oil_x_brent_bin created\n")
  }
  if ("oil_bartik_iv" %in% names(pnl)) {
    pnl[, bartik_x_brent := oil_bartik_iv * oy]
    cat("    bartik_x_brent created\n")
  }
  if ("spillover_exposure" %in% names(pnl)) {
    pnl[, spillover_x_brent := spillover_exposure * oy]
    cat("    spillover_x_brent created\n")
  }
  if ("spillover_exposure_wtd" %in% names(pnl)) {
    pnl[, spillover_wtd_x_brent := spillover_exposure_wtd * oy]
    cat("    spillover_wtd_x_brent created\n")
  }
  fomc_col <- sub("yoy_oil", "fomc_regime", oil_yoy_col)
  if (fomc_col %in% names(pnl)) {
    pnl[, fomc_x_brent := get(fomc_col) * oy]
    cat("    fomc_x_brent created\n")
  }
  invisible(pnl)
}

cat("  --- Baseline scenario ---\n")
add_interactions(panel_base,   "macro_base_yoy_oil")
cat("\n  --- Severely adverse scenario ---\n")
add_interactions(panel_severe, "macro_severe_yoy_oil")

log_checkpoint("CK-30: Interaction terms", "PASS",
               "oil_x_brent, spillover_x_brent, fomc_x_brent for both panels")

# ============================================================================
# SECTION 15: SAVE FINAL PANELS
# ============================================================================
section("15", "Save Final Panels")

setorderv(cr,           c("join_number","year","quarter"))
setorderv(panel_base,   c("join_number","year","quarter"))
setorderv(panel_severe, c("join_number","year","quarter"))

saveRDS(cr,           "Data/call_clean.rds")
saveRDS(panel_base,   "Data/panel_base.rds")
saveRDS(panel_severe, "Data/panel_severe.rds")

cat(sprintf("  %-20s  %s rows x %d cols  %.1f MB\n", "call_clean.rds",
            format(nrow(cr),big.mark=","), ncol(cr),
            file.size("Data/call_clean.rds")/1e6))
cat(sprintf("  %-20s  %s rows x %d cols  %.1f MB\n", "panel_base.rds",
            format(nrow(panel_base),big.mark=","), ncol(panel_base),
            file.size("Data/panel_base.rds")/1e6))
cat(sprintf("  %-20s  %s rows x %d cols  %.1f MB\n", "panel_severe.rds",
            format(nrow(panel_severe),big.mark=","), ncol(panel_severe),
            file.size("Data/panel_severe.rds")/1e6))

log_checkpoint("CK-31: Final panels saved", "PASS",
               "call_clean.rds, panel_base.rds, panel_severe.rds")

# ============================================================================
# SECTION 16: DATA QUALITY REPORT
# ============================================================================
section("16", "Data Quality Report")

cat("  CALL REPORT SUMMARY\n")
cat(sprintf("  %-22s : %s\n", "CU-quarter obs",
            format(nrow(cr), big.mark = ",")))
cat(sprintf("  %-22s : %s\n", "Unique CUs",
            format(uniqueN(cr$join_number), big.mark = ",")))
cat(sprintf("  %-22s : %dQ%d -- %dQ%d\n", "Quarter range",
            min(cr$year), cr[which.min(yyyyqq), quarter],
            max(cr$year), cr[which.max(yyyyqq), quarter]))

if ("asset_tier" %in% names(cr)) {
  cat("\n  Asset tiers:\n")
  tier_tbl <- cr[, .N, by = asset_tier][order(asset_tier)]
  for (i in seq_len(nrow(tier_tbl))) {
    cat(sprintf("    %-20s : %s (%.1f%%)\n",
                as.character(tier_tbl$asset_tier[i]),
                format(tier_tbl$N[i], big.mark = ","),
                tier_tbl$N[i] / nrow(cr) * 100))
  }
}

if ("cu_group" %in% names(cr)) {
  cat("\n  CU exposure groups:\n")
  grp_tbl <- cr[, .N, by = cu_group][order(cu_group)]
  for (i in seq_len(nrow(grp_tbl))) {
    cat(sprintf("    %-20s : %s (%.1f%%)\n",
                as.character(grp_tbl$cu_group[i]),
                format(grp_tbl$N[i], big.mark = ","),
                grp_tbl$N[i] / nrow(cr) * 100))
  }
}

# PBRENT range
cat("\n  FRB MACRO PBRENT ($/bbl)\n")
for (nm in c("macro_base_pbrent", "macro_severe_pbrent")) {
  pnl <- if (grepl("base", nm)) panel_base else panel_severe
  if (nm %in% names(pnl)) {
    v <- pnl[[nm]]
    cat(sprintf("  %-30s : min=%5.1f  max=%5.1f  latest=%5.1f\n",
                nm,
                min(v, na.rm = TRUE),
                max(v, na.rm = TRUE),
                v[max(which(!is.na(v)))]))
  }
}

# Missingness table
cat("\n  KEY VARIABLE MISSINGNESS (panel_base)\n")
check_vars <- c(
  "join_number","year","quarter",
  "netintmrg","networth","pcanetworth","costfds","roa",
  "dq_rate","chg_tot_lns_ratio",
  "insured_tot","dep_shrcert","acct_018","members",
  "insured_share_growth","cert_share","loan_to_share",
  "dep_growth_yoy","member_growth_yoy","pll_rate","pll_per_loan","nim_spread",
  "macro_base_pbrent","macro_base_lurc","macro_base_pcpi",
  "macro_base_rmtg","macro_base_phpi","macro_base_yoy_oil",
  "macro_base_yield_curve","macro_base_fomc_regime",
  "oil_exposure_cont","oil_exposure_bin","spillover_exposure",
  "oil_x_brent","spillover_x_brent","fomc_x_brent","oil_bartik_iv"
)
fv  <- intersect(check_vars, names(panel_base))
pct <- sapply(panel_base[, ..fv], function(x) round(mean(is.na(x)) * 100, 1))
miss_tbl <- data.table(variable = names(pct), pct_missing = pct)[order(-pct_missing)]

cat(sprintf("\n  %-34s  %s\n", "Variable", "% Missing"))
cat(sprintf("  %s\n", paste(rep("-", 46), collapse = "")))
for (i in seq_len(nrow(miss_tbl))) {
  flag <- if (miss_tbl$pct_missing[i] > 20) "  <- HIGH"
          else if (miss_tbl$pct_missing[i] > 5) "  <- NOTE"
          else ""
  cat(sprintf("  %-34s  %5.1f%%%s\n",
              miss_tbl$variable[i], miss_tbl$pct_missing[i], flag))
}

log_checkpoint("CK-32: Data quality report", "PASS",
               sprintf("%d variables checked", nrow(miss_tbl)))

# ============================================================================
# SECTION 17: SAVE AUDIT FILES
# ============================================================================
section("17", "Save Audit Files")

outlier_dt <- rbindlist(OUTLIER_RECORDS, fill = TRUE)
fwrite(outlier_dt, file.path(RESULTS_DIR, "01_outlier_log.csv"))
cat(sprintf("  Saved: Results/01_outlier_log.csv  (%d variable records)\n",
            nrow(outlier_dt)))

ALL_OUT_VARS <- intersect(
  c("insured_share_growth","dep_growth_yoy","cert_growth_yoy",
    "member_growth_yoy","pll_rate","pll_per_loan","cert_share",
    "loan_to_share","nim_spread","netintmrg","pcanetworth",
    "costfds","roa","dq_rate","chg_tot_lns_ratio"),
  names(panel_base)
)

var_summ <- rbindlist(lapply(ALL_OUT_VARS, function(v) {
  x  <- panel_base[[v]]
  qs <- tryCatch(
    quantile(x, c(0.01,.10,.25,.50,.75,.90,.99), na.rm = TRUE),
    error = function(e) rep(NA_real_, 7)
  )
  data.table(
    variable = v,
    n        = sum(!is.na(x)),
    pct_miss = round(mean(is.na(x)) * 100, 2),
    mean     = round(mean(x, na.rm=TRUE), 6),
    sd       = round(sd(x, na.rm=TRUE), 6),
    p01 = round(qs[1],6), p10 = round(qs[2],6), p25 = round(qs[3],6),
    p50 = round(qs[4],6), p75 = round(qs[5],6), p90 = round(qs[6],6),
    p99 = round(qs[7],6)
  )
}))

fwrite(var_summ, file.path(RESULTS_DIR, "01_variable_summary.csv"))
cat(sprintf("  Saved: Results/01_variable_summary.csv  (%d variables)\n",
            nrow(var_summ)))

log_checkpoint("CK-33: Audit files saved", "PASS",
               "01_outlier_log.csv, 01_variable_summary.csv")

# ============================================================================
# SECTION 18: CHECKPOINT SUMMARY
# ============================================================================
section("18", "Checkpoint Summary")

n_pass <- sum(CHKPT_LOG$status == "PASS")
n_warn <- sum(CHKPT_LOG$status == "WARN")
n_fail <- sum(CHKPT_LOG$status == "FAIL")
n_info <- sum(CHKPT_LOG$status == "INFO")

cat(sprintf("  Total: %d  |  PASS: %d  |  WARN: %d  |  FAIL: %d  |  INFO: %d\n",
            nrow(CHKPT_LOG), n_pass, n_warn, n_fail, n_info))

if (n_warn > 0) {
  cat("\n  Warnings requiring attention:\n")
  for (i in which(CHKPT_LOG$status == "WARN")) {
    cat(sprintf("    ! %s\n", CHKPT_LOG$checkpoint[i]))
    if (nchar(CHKPT_LOG$detail[i]) > 0)
      cat(sprintf("      => %s\n", CHKPT_LOG$detail[i]))
  }
}

fwrite(as.data.table(CHKPT_LOG),
       file.path(RESULTS_DIR, "01_checkpoint_log.csv"))
cat(sprintf("\n  Checkpoint log saved: Results/01_checkpoint_log.csv\n"))

# ============================================================================
# CLOSING BANNER
# ============================================================================
cat("\n", SEP, "\n", sep = "")
cat("  SCRIPT 01 COMPLETE\n")
cat(sprintf("  Finished : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("\n  Data outputs:\n")
cat("    Data/call_clean.rds      CU panel + derived vars + oil exposure\n")
cat("    Data/macro_base.rds      Baseline macro (macro_base_ prefix)\n")
cat("    Data/macro_severe.rds    Severely adverse (macro_severe_ prefix)\n")
cat("    Data/panel_base.rds      Full panel: CU x macro_base + interactions\n")
cat("    Data/panel_severe.rds    Full panel: CU x macro_severe\n")
cat("\n  Audit outputs:\n")
cat("    Results/01_outlier_log.csv       Winsorization summary per variable\n")
cat("    Results/01_variable_summary.csv  Post-cleaning distributional stats\n")
cat("    Results/01_checkpoint_log.csv    All checkpoint pass/warn/fail records\n")
cat("\n  CU-level derived variables:\n")
cat("    insured_share_growth  YoY change in insured shares  [#1 SHAP outcome, v1]\n")
cat("    dep_growth_yoy        YoY change in total deposits (acct_018)\n")
cat("    cert_growth_yoy       YoY change in certificate shares\n")
cat("    member_growth_yoy     YoY change in membership\n")
cat("    cert_share            Certificate share of deposits (0-1, bounded)\n")
cat("    loan_to_share         Loan-to-share ratio (0-2, bounded)\n")
cat("    nim_spread            yldavgloans - costfds\n")
cat("    pll_rate              pll/avg(lns_tot) x100 (bounded [-2,5], winsorised)\n")
cat("    pll_per_loan          pll/lns_tot_n ($ per loan count)\n")
cat("\n  Direct/indirect effect interactions (both scenarios):\n")
cat("    oil_x_brent           oil_exposure_cont x yoy_oil\n")
cat("    oil_x_brent_bin       oil_exposure_bin  x yoy_oil\n")
cat("    bartik_x_brent        oil_bartik_iv     x yoy_oil\n")
cat("    spillover_x_brent     spillover_exposure x yoy_oil\n")
cat("    fomc_x_brent          fomc_regime        x yoy_oil\n")
cat(sprintf("\n  Checkpoints: %d PASS | %d WARN | %d FAIL\n", n_pass, n_warn, n_fail))
cat(SEP, "\n", sep = "")
