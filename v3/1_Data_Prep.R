# =============================================================================
# SCRIPT 01 — Data Ingestion, Cleaning, and Variable Construction
# Oil Price Shock × Credit Union Financial Performance
# Author: Saurabh C. Datta, Ph.D. | NCUA Office of Chief Economist
# Version: v2.1-DIAGNOSTIC | Built: 2026-03-30
# Working dir: S:/Projects/Oil_Price_Shock_2026/
# =============================================================================
# PURPOSE:
#   Load NCUA Form 5300 call report data, filter to FICUs, construct all
#   8 outcome variables from raw accounts, handle outliers transparently,
#   assign asset tiers, and save clean panel for downstream scripts.
#
# OUTPUT: Data/01_panel_clean.rds
#         Results/01_outlier_log.csv      — every observation winsorized
#         Results/01_variable_summary.csv — pre/post-winsorization stats
#         Results/01_checkpoint_log.csv   — all checkpoint pass/fail records
#
# DIAGNOSTIC PHILOSOPHY:
#   Every transformation is logged to screen. Checkpoints halt on critical
#   failures and warn on soft failures. Outlier treatment is fully auditable.
#   The log is designed to be self-contained — paste it into the project record.
#
# BACKWARD REASONING ANCHORS (from v1 confirmed findings):
#   - AR(1) ≈ 0.59 → lagged outcomes constructed here
#   - Deposit growth = most exposed (3x SHAP) → insured_share_growth critical
#   - pll_rate needs avg loan denominator → lns_tot lagged here
#   - NIM = acct_116 / acct_010 → zero-denominator handled explicitly
#   - 2015Q1 structural break → post2015 flag aligned to date here
# =============================================================================

rm(list = ls())
gc()

# ============================================================================
# SECTION 0: LOGGING INFRASTRUCTURE
# ============================================================================

# Checkpoint log (accumulates across entire script)
CHKPT_LOG <- data.frame(
  checkpoint = character(),
  status     = character(),
  detail     = character(),
  stringsAsFactors = FALSE
)

# Separator lines
SEP  <- paste(rep("=", 72), collapse = "")
SEP2 <- paste(rep("-", 72), collapse = "")

# Timestamp
ts <- function() format(Sys.time(), "[%H:%M:%S]")

# Log a checkpoint — prints to screen AND accumulates in CHKPT_LOG
# status: "PASS", "WARN", "FAIL", "INFO"
log_checkpoint <- function(label, status, detail = "") {
  icon <- switch(status,
    PASS = "✓ PASS",
    WARN = "! WARN",
    FAIL = "✗ FAIL",
    INFO = "  INFO"
  )
  cat(sprintf("  %s %s %s\n", ts(), icon, label))
  if (nchar(detail) > 0) cat(sprintf("         └─ %s\n", detail))
  CHKPT_LOG <<- rbind(CHKPT_LOG, data.frame(
    checkpoint = label,
    status     = status,
    detail     = detail,
    stringsAsFactors = FALSE
  ))
  if (status == "FAIL") {
    cat("\n", SEP, "\n")
    cat("  FATAL: Script halted at checkpoint:", label, "\n")
    cat(SEP, "\n")
    stop(paste("FATAL checkpoint failure:", label, "|", detail))
  }
}

# Section banner
section <- function(n, title) {
  cat("\n", SEP, "\n", sep = "")
  cat(sprintf("  SECTION %s: %s\n", n, title))
  cat(SEP2, "\n", sep = "")
}

# Summary stats row — prints a clean table row
stat_row <- function(varname, n, nmiss, mean, sd, p1, p10, p50, p90, p99,
                     mn, mx, label = "") {
  cat(sprintf(
    "  %-28s  N=%7d  miss=%5.1f%%  mean=%9.4f  sd=%8.4f\n",
    varname, n, nmiss * 100, mean, sd
  ))
  cat(sprintf(
    "  %-28s  p1=%8.4f  p10=%7.4f  p50=%7.4f  p90=%7.4f  p99=%8.4f\n",
    "", p1, p10, p50, p90, p99
  ))
  cat(sprintf(
    "  %-28s  min=%8.4f  max=%7.4f%s\n",
    "", mn, mx, if (nchar(label) > 0) paste0("  [", label, "]") else ""
  ))
  cat(sprintf("  %s\n", paste(rep("-", 70), collapse = "")))
}

# Compute and print summary stats for a vector; return named list
describe_var <- function(x, varname, label = "") {
  n      <- sum(!is.na(x))
  nmiss  <- mean(is.na(x))
  if (n == 0) {
    cat(sprintf("  %-28s  ALL MISSING\n", varname))
    return(invisible(NULL))
  }
  qs <- quantile(x, probs = c(0.01, 0.10, 0.50, 0.90, 0.99), na.rm = TRUE)
  stat_row(
    varname = varname, n = n, nmiss = nmiss,
    mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
    p1 = qs[1], p10 = qs[2], p50 = qs[3], p90 = qs[4], p99 = qs[5],
    mn = min(x, na.rm = TRUE), mx = max(x, na.rm = TRUE),
    label = label
  )
  invisible(list(n = n, nmiss = nmiss, mean = mean(x, na.rm = TRUE),
                 sd = sd(x, na.rm = TRUE), p1 = qs[[1]], p99 = qs[[5]],
                 min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE)))
}

# ============================================================================
# SCRIPT HEADER
# ============================================================================
cat(SEP, "\n")
cat("  SCRIPT 01 — Data Ingestion, Cleaning, and Variable Construction\n")
cat("  Oil Price Shock x Credit Union Financial Performance — v2\n")
cat("  Author  : Saurabh C. Datta, Ph.D. | NCUA Office of Chief Economist\n")
cat("  Version : v2.1-DIAGNOSTIC\n")
cat(sprintf("  Started : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("  Working : S:/Projects/Oil_Price_Shock_2026/\n")
cat(SEP, "\n\n")

# ============================================================================
# SECTION 1: LIBRARIES AND CONFIGURATION
# ============================================================================
section("1", "Libraries and Configuration")

suppressPackageStartupMessages({
  library(data.table)
  library(haven)
  library(stringr)
  library(lubridate)
})

# Confirm packages loaded
for (pkg in c("data.table", "haven", "stringr", "lubridate")) {
  log_checkpoint(paste("Package:", pkg), "PASS",
                 paste("v", packageVersion(pkg)))
}

# ---------------------------------------------------------------------------
# Configuration — edit these paths as needed
# ---------------------------------------------------------------------------
DATA_DIR    <- "Data/"
RESULTS_DIR <- "Results/"

# Primary raw file — set to actual path for production
# Supported: .dta (Stata), .rds, .csv
RAW_FILE    <- "Data/form5300_raw.dta"

# Winsorization thresholds — document in project record if changed
WIN_LOW  <- 0.01   # 1st percentile lower bound
WIN_HIGH <- 0.99   # 99th percentile upper bound

# Minimum cell count thresholds for validity checks
MIN_CU_COUNT   <- 1000   # expect at least 1,000 FICUs
MIN_QTR_COUNT  <- 40     # expect at least 40 quarters (10 years)

cat(sprintf("\n  RAW_FILE    : %s\n", RAW_FILE))
cat(sprintf("  DATA_DIR    : %s\n", DATA_DIR))
cat(sprintf("  RESULTS_DIR : %s\n", RESULTS_DIR))
cat(sprintf("  Winsorize   : [%.0f%%, %.0f%%]\n", WIN_LOW * 100, WIN_HIGH * 100))
cat(sprintf("  Min CU threshold  : %d\n", MIN_CU_COUNT))
cat(sprintf("  Min Qtr threshold : %d\n\n", MIN_QTR_COUNT))

# Create output directories
for (d in c(DATA_DIR, RESULTS_DIR)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    cat(sprintf("  Created directory: %s\n", d))
  }
}

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# yyyyqq parse (BULLETPROOF — validates format before split)
parse_q <- function(x) {
  x      <- as.character(x)
  result <- rep(as.Date(NA), length(x))
  valid  <- grepl("^[0-9]{4}\\.[1-4]$", x)
  if (any(valid)) {
    p             <- str_split_fixed(x[valid], "\\.", 2)
    yr            <- as.integer(p[, 1])
    qtr           <- as.integer(p[, 2])
    mo            <- (qtr - 1L) * 3L + 1L
    result[valid] <- as.Date(sprintf("%04d-%02d-01", yr, mo))
  }
  result
}

# Winsorize and return list: winsorized vector + audit info
winsorize_audit <- function(x, varname,
                            p_low  = WIN_LOW,
                            p_high = WIN_HIGH) {
  qs   <- quantile(x, probs = c(p_low, p_high), na.rm = TRUE)
  lo   <- qs[1]; hi <- qs[2]
  n_lo <- sum(!is.na(x) & x < lo)
  n_hi <- sum(!is.na(x) & x > hi)
  x_w  <- pmax(pmin(x, hi), lo)
  list(
    x       = x_w,
    lo      = lo,
    hi      = hi,
    n_lo    = n_lo,
    n_hi    = n_hi,
    n_total = n_lo + n_hi,
    pct_affected = round((n_lo + n_hi) / sum(!is.na(x)) * 100, 3)
  )
}

log_checkpoint("Helper functions defined", "PASS", "parse_q, winsorize_audit")

# ============================================================================
# SECTION 2: LOAD RAW DATA
# ============================================================================
section("2", "Load Raw Form 5300 Data")

cat(sprintf("  Attempting to load: %s\n\n", RAW_FILE))

if (!file.exists(RAW_FILE)) {
  log_checkpoint("Raw file exists", "FAIL",
                 paste("File not found:", RAW_FILE))
}
log_checkpoint("Raw file exists", "PASS", RAW_FILE)

ext <- tolower(tools::file_ext(RAW_FILE))
cat(sprintf("  Detected file type: .%s\n", ext))

t_load <- system.time({
  if (ext == "dta") {
    raw <- as.data.table(haven::read_dta(RAW_FILE))
    # Strip Stata labels IMMEDIATELY — CRITICAL (pitfall #4)
    # haven_labelled objects break as.double() without this step
    raw <- haven::zap_labels(raw)
    raw <- haven::zap_formats(raw)
    cat("  Applied: haven::zap_labels() + zap_formats() [Stata label strip]\n")
  } else if (ext == "rds") {
    raw <- readRDS(RAW_FILE)
    if (!is.data.table(raw)) setDT(raw)
  } else if (ext == "csv") {
    raw <- fread(RAW_FILE)
  } else {
    log_checkpoint("Supported file type", "FAIL",
                   paste("Unsupported extension:", ext))
  }
})

cat(sprintf("\n  Load time   : %.1f seconds\n", t_load["elapsed"]))
cat(sprintf("  Rows loaded : %s\n", format(nrow(raw), big.mark = ",")))
cat(sprintf("  Columns     : %d\n", ncol(raw)))
cat(sprintf("  Object size : %.1f MB\n", object.size(raw) / 1e6))

# CK-01: Raw file has meaningful data
if (nrow(raw) < 10000) {
  log_checkpoint("CK-01: Raw row count", "WARN",
                 sprintf("Only %d rows loaded — expected >10,000", nrow(raw)))
} else {
  log_checkpoint("CK-01: Raw row count", "PASS",
                 sprintf("%s rows", format(nrow(raw), big.mark = ",")))
}

# ============================================================================
# SECTION 3: FILTER TO FEDERALLY INSURED CREDIT UNIONS (FICUs)
# ============================================================================
section("3", "Filter to FICUs and Type Coercion")

# --- 3a. cu_type filter ---
cat("  Filter: cu_type == 1 (Federally Insured Credit Unions)\n\n")

if ("cu_type" %in% names(raw)) {
  type_dist <- table(raw$cu_type, useNA = "always")
  cat("  cu_type distribution in raw file:\n")
  print(type_dist)
  cat("\n")

  n_before <- nrow(raw)
  panel    <- raw[cu_type == 1]
  n_after  <- nrow(panel)
  n_dropped <- n_before - n_after

  cat(sprintf("  Rows before filter : %s\n", format(n_before, big.mark = ",")))
  cat(sprintf("  Rows after filter  : %s\n", format(n_after,  big.mark = ",")))
  cat(sprintf("  Rows dropped       : %s (%.1f%%)\n",
              format(n_dropped, big.mark = ","),
              n_dropped / n_before * 100))

  if (n_after < MIN_CU_COUNT) {
    log_checkpoint("CK-02: FICU row count after filter", "FAIL",
                   sprintf("Only %d rows after cu_type==1 filter", n_after))
  } else {
    log_checkpoint("CK-02: FICU row count after filter", "PASS",
                   sprintf("%s rows retained", format(n_after, big.mark = ",")))
  }
} else {
  log_checkpoint("CK-02: cu_type column present", "WARN",
                 "cu_type not found — using all rows as FICUs")
  panel <- copy(raw)
}

# --- 3b. Force numeric columns to double ---
cat("\n  Type coercion: forcing numeric/labelled columns to double...\n")
n_coerced <- 0L
n_failed  <- 0L

for (cn in names(panel)) {
  x <- panel[[cn]]
  if (is.numeric(x) || inherits(x, "haven_labelled")) {
    ok <- tryCatch({
      set(panel, j = cn, value = as.double(x))
      TRUE
    }, error = function(e) FALSE)
    if (ok) n_coerced <- n_coerced + 1L else n_failed <- n_failed + 1L
  }
}

cat(sprintf("  Columns coerced to double : %d\n", n_coerced))
cat(sprintf("  Coercion failures         : %d\n", n_failed))
log_checkpoint("CK-03: Type coercion", if (n_failed == 0) "PASS" else "WARN",
               sprintf("%d coerced, %d failures", n_coerced, n_failed))

# --- 3c. Drop non-atomic columns ---
bad_cols <- names(panel)[sapply(names(panel), function(cn) !is.atomic(panel[[cn]]))]
if (length(bad_cols) > 0) {
  cat(sprintf("  Dropping %d non-atomic columns: %s\n",
              length(bad_cols), paste(bad_cols, collapse = ", ")))
  panel[, (bad_cols) := NULL]
  log_checkpoint("CK-04: Non-atomic columns", "WARN",
                 sprintf("Dropped: %s", paste(bad_cols, collapse = ", ")))
} else {
  log_checkpoint("CK-04: Non-atomic columns", "PASS", "None found")
}

# ============================================================================
# SECTION 4: IDENTIFIER STANDARDIZATION
# ============================================================================
section("4", "Identifier Standardization")

# --- 4a. join_number ---
if (!"join_number" %in% names(panel)) {
  log_checkpoint("CK-05: join_number column", "FAIL",
                 "Column not found — cannot identify CUs")
}
n_cu <- uniqueN(panel$join_number)
cat(sprintf("  Unique CU identifiers (join_number): %s\n",
            format(n_cu, big.mark = ",")))
log_checkpoint("CK-05: join_number column", "PASS",
               sprintf("%s unique CUs", format(n_cu, big.mark = ",")))

# --- 4b. yyyyqq date parsing ---
if (!"yyyyqq" %in% names(panel)) {
  log_checkpoint("CK-06: yyyyqq column", "FAIL",
                 "Quarter identifier not found")
}

panel[, yyyyqq := as.character(yyyyqq)]
cat(sprintf("\n  yyyyqq sample values: %s\n",
            paste(head(unique(panel$yyyyqq), 6), collapse = ", ")))

n_before_date <- nrow(panel)
panel[, date := parse_q(yyyyqq)]

n_bad_dates <- sum(is.na(panel$date))
n_valid_fmt <- sum(grepl("^[0-9]{4}\\.[1-4]$", panel$yyyyqq))

cat(sprintf("  Rows with valid yyyyqq format  : %s\n",
            format(n_valid_fmt, big.mark = ",")))
cat(sprintf("  Rows with unparseable yyyyqq   : %s\n",
            format(n_bad_dates, big.mark = ",")))

if (n_bad_dates > 0) {
  bad_samples <- head(unique(panel[is.na(date), yyyyqq]), 5)
  cat(sprintf("  Bad yyyyqq samples             : %s\n",
              paste(bad_samples, collapse = ", ")))
  panel <- panel[!is.na(date)]
  log_checkpoint("CK-06: yyyyqq parsing", "WARN",
                 sprintf("%d rows dropped for bad yyyyqq format", n_bad_dates))
} else {
  log_checkpoint("CK-06: yyyyqq parsing", "PASS", "All rows parse cleanly")
}

# --- 4c. Date range ---
date_min <- min(panel$date)
date_max <- max(panel$date)
n_qtrs   <- uniqueN(panel$yyyyqq)

cat(sprintf("\n  Date range : %s  to  %s\n", format(date_min), format(date_max)))
cat(sprintf("  Quarters   : %d\n", n_qtrs))
cat(sprintf("  Years span : %.1f\n", as.numeric(date_max - date_min) / 365.25))

if (n_qtrs < MIN_QTR_COUNT) {
  log_checkpoint("CK-07: Date range", "WARN",
                 sprintf("Only %d quarters — expected >=%d", n_qtrs, MIN_QTR_COUNT))
} else {
  log_checkpoint("CK-07: Date range", "PASS",
                 sprintf("%d quarters (%s to %s)", n_qtrs,
                         format(date_min), format(date_max)))
}

# --- 4d. reporting_state ---
if ("reporting_state" %in% names(panel)) {
  n_states <- uniqueN(panel$reporting_state)
  cat(sprintf("\n  Unique reporting states: %d\n", n_states))
  cat(sprintf("  State distribution (top 10 by N):\n"))
  state_ct <- panel[, .N, by = reporting_state][order(-N)][1:min(10, .N)]
  for (i in seq_len(nrow(state_ct))) {
    cat(sprintf("    %-4s  %s\n", state_ct$reporting_state[i],
                format(state_ct$N[i], big.mark = ",")))
  }
  log_checkpoint("CK-08: reporting_state", "PASS",
                 sprintf("%d unique states", n_states))
} else {
  log_checkpoint("CK-08: reporting_state", "WARN",
                 "reporting_state column not found — state-level analysis will fail")
}

# Sort panel for all lag operations
setorder(panel, join_number, date)

# ============================================================================
# SECTION 5: RAW ACCOUNT COLUMN AUDIT
# ============================================================================
section("5", "Raw Account Column Audit")

cat("  Checking for required Form 5300 account columns...\n\n")

# Required accounts and what they power
ACCT_MANIFEST <- list(
  acct_041a = "Loans 60+ days delinquent → dq_rate numerator",
  acct_025  = "Total loans outstanding → dq_rate denominator, loan_to_share numerator",
  acct_116  = "Net interest income → netintmrg numerator",
  acct_010  = "Average assets → netintmrg denominator, pcanetworth denominator",
  acct_018  = "Total shares/deposits → insured_share_growth, loan_to_share fallback",
  acct_083  = "Number of members → member_growth_yoy",
  acct_131  = "Interest expense → costfds numerator",
  acct_674  = "Provision for loan losses → pll_rate numerator",
  acct_998  = "Net worth → pcanetworth numerator",
  dep_tot   = "Total deposits (preferred) → loan_to_share denominator"
)

cat(sprintf("  %-12s  %-6s  %-10s  %s\n", "Account", "Found", "Non-NA%", "Powers"))
cat(sprintf("  %s\n", paste(rep("-", 70), collapse = "")))

acct_status <- list()
for (acct in names(ACCT_MANIFEST)) {
  found   <- acct %in% names(panel)
  pct_ok  <- if (found) round(mean(!is.na(panel[[acct]])) * 100, 1) else NA
  status  <- if (!found) "MISSING" else if (pct_ok < 50) "LOW" else "OK"
  cat(sprintf("  %-12s  %-6s  %s  %s\n",
              acct,
              if (found) "YES" else "NO",
              if (found) sprintf("%5.1f%%", pct_ok) else "  ---  ",
              ACCT_MANIFEST[[acct]]))
  acct_status[[acct]] <- list(found = found, pct_ok = pct_ok, status = status)
}

# Critical accounts — fail if missing
CRITICAL_ACCTS <- c("acct_025", "acct_010", "acct_018")
for (acct in CRITICAL_ACCTS) {
  if (!acct_status[[acct]]$found) {
    log_checkpoint(paste("CK-09: Critical account", acct), "FAIL",
                   "Required for multiple outcomes — cannot proceed")
  }
}
log_checkpoint("CK-09: Critical account columns", "PASS",
               paste(CRITICAL_ACCTS, collapse = ", "))

# ============================================================================
# SECTION 6: OUTCOME VARIABLE CONSTRUCTION
# ============================================================================
section("6", "Outcome Variable Construction")

# We track construction status for each variable
OUTCOMES <- c("dq_rate", "pll_rate", "netintmrg", "insured_share_growth",
              "member_growth_yoy", "costfds", "loan_to_share", "pcanetworth")

# Outcome metadata for interpretation
OUTCOME_META <- list(
  dq_rate              = list(label = "Delinquency Rate",           formula = "acct_041a / acct_025",         expected_range = c(0, 0.30)),
  pll_rate             = list(label = "PLL Rate",                    formula = "acct_674 / avg(acct_025)",      expected_range = c(0, 0.05)),
  netintmrg            = list(label = "Net Interest Margin",         formula = "acct_116 / acct_010",          expected_range = c(0, 0.10)),
  insured_share_growth = list(label = "Deposit Growth (QoQ%)",       formula = "Δacct_018 / acct_018_lag",     expected_range = c(-0.50, 0.50)),
  member_growth_yoy    = list(label = "Membership Growth (YoY%)",    formula = "Δacct_083 / acct_083_lag4",    expected_range = c(-0.30, 0.30)),
  costfds              = list(label = "Cost of Funds",               formula = "acct_131 / acct_018 * 4",     expected_range = c(0, 0.15)),
  loan_to_share        = list(label = "Loan-to-Share Ratio",         formula = "acct_025 / dep_tot",           expected_range = c(0, 2.00)),
  pcanetworth          = list(label = "Net Worth Ratio",             formula = "acct_998 / acct_010",          expected_range = c(-0.05, 0.30))
)

# Raw account descriptives before construction
cat("  Raw denominator/driver variables before construction:\n\n")
for (acct in c("acct_025", "acct_010", "acct_018", "acct_083")) {
  if (acct %in% names(panel)) {
    describe_var(panel[[acct]], acct)
  }
}

# --------------------------------------------------------------------------
# 6a. dq_rate = acct_041a / acct_025
# --------------------------------------------------------------------------
cat("\n  [6a] dq_rate — Delinquency Rate\n")
cat("       Formula : acct_041a / acct_025\n")
cat("       Logic   : Loans 60+ days past due / Total loans; zero denominator → NA\n\n")

if (all(c("acct_041a", "acct_025") %in% names(panel))) {
  n_zero_denom <- panel[acct_025 <= 0, .N]
  n_neg_num    <- panel[acct_041a < 0, .N]
  cat(sprintf("  acct_025 <= 0 (zero denominator): %s rows → set to NA\n",
              format(n_zero_denom, big.mark = ",")))
  cat(sprintf("  acct_041a < 0  (negative delinq) : %s rows\n",
              format(n_neg_num, big.mark = ",")))

  panel[, dq_rate := fifelse(acct_025 > 0, acct_041a / acct_025, NA_real_)]

  # Economic sanity: dq_rate > 1 is impossible (more delinquent than loans)
  n_gt1 <- panel[!is.na(dq_rate) & dq_rate > 1, .N]
  if (n_gt1 > 0) {
    cat(sprintf("  WARNING: %d obs have dq_rate > 1.0 (economically impossible)\n",
                n_gt1))
    cat("  Action: Setting dq_rate > 1 to NA before winsorization\n")
    panel[dq_rate > 1, dq_rate := NA_real_]
  }

  cat("\n  Pre-winsorization stats:\n")
  describe_var(panel$dq_rate, "dq_rate", "pre-winsor")
  log_checkpoint("CK-10: dq_rate constructed", "PASS",
                 sprintf("zero-denom: %d, impossible values fixed: %d",
                         n_zero_denom, n_gt1))
} else {
  missing_accts <- setdiff(c("acct_041a", "acct_025"), names(panel))
  panel[, dq_rate := NA_real_]
  log_checkpoint("CK-10: dq_rate constructed", "WARN",
                 paste("Missing accounts:", paste(missing_accts, collapse = ", ")))
}

# --------------------------------------------------------------------------
# 6b. netintmrg = acct_116 / acct_010
# --------------------------------------------------------------------------
cat("\n  [6b] netintmrg — Net Interest Margin\n")
cat("       Formula : acct_116 / acct_010\n")
cat("       Logic   : Net interest income / Average assets; zero denominator → NA\n\n")

if (all(c("acct_116", "acct_010") %in% names(panel))) {
  n_zero_denom <- panel[acct_010 <= 0, .N]
  cat(sprintf("  acct_010 <= 0 (zero/neg assets): %s rows → set to NA\n",
              format(n_zero_denom, big.mark = ",")))

  panel[, netintmrg := fifelse(acct_010 > 0, acct_116 / acct_010, NA_real_)]

  # NIM should be bounded: [-0.05, 0.20] under normal conditions
  n_extreme <- panel[!is.na(netintmrg) & (netintmrg < -0.05 | netintmrg > 0.25), .N]
  if (n_extreme > 0) {
    cat(sprintf("  Extreme NIM values (outside [-5%%, +25%%]): %d obs (flagged)\n",
                n_extreme))
  }

  cat("\n  Pre-winsorization stats:\n")
  describe_var(panel$netintmrg, "netintmrg", "pre-winsor")
  log_checkpoint("CK-11: netintmrg constructed", "PASS",
                 sprintf("zero-denom: %d", n_zero_denom))
} else {
  missing_accts <- setdiff(c("acct_116", "acct_010"), names(panel))
  panel[, netintmrg := NA_real_]
  log_checkpoint("CK-11: netintmrg constructed", "WARN",
                 paste("Missing:", paste(missing_accts, collapse = ", ")))
}

# --------------------------------------------------------------------------
# 6c. insured_share_growth = Δacct_018 / acct_018_lag (QoQ %)
# CRITICAL: Deposit growth is the single most exposed outcome (3x SHAP, v1)
# --------------------------------------------------------------------------
cat("\n  [6c] insured_share_growth — Deposit Growth (QoQ %)\n")
cat("       Formula : (acct_018_t - acct_018_{t-1}) / acct_018_{t-1} * 100\n")
cat("       NOTE    : THIS IS THE #1 OUTCOME BY SHAP IMPORTANCE (v1 finding)\n")
cat("       Logic   : QoQ growth rate; lag = NA for first obs → NA result\n\n")

if ("acct_018" %in% names(panel)) {
  panel[, acct_018_lag := shift(acct_018, 1L, type = "lag"), by = join_number]

  n_first_obs    <- panel[is.na(acct_018_lag), .N]
  n_zero_lag     <- panel[!is.na(acct_018_lag) & acct_018_lag <= 0, .N]
  n_neg_deposits <- panel[acct_018 < 0, .N]

  cat(sprintf("  First-obs NA (no lag)           : %s rows\n",
              format(n_first_obs, big.mark = ",")))
  cat(sprintf("  Zero/negative lagged deposits   : %s rows → NA\n",
              format(n_zero_lag, big.mark = ",")))
  cat(sprintf("  Negative current deposits       : %s rows\n",
              format(n_neg_deposits, big.mark = ",")))

  panel[, insured_share_growth := fifelse(
    !is.na(acct_018_lag) & acct_018_lag > 0,
    (acct_018 - acct_018_lag) / acct_018_lag * 100,
    NA_real_
  )]

  # Flag implausible growth rates (>200% or < -90% QoQ are data issues)
  n_implaus_hi <- panel[!is.na(insured_share_growth) & insured_share_growth > 200, .N]
  n_implaus_lo <- panel[!is.na(insured_share_growth) & insured_share_growth < -90, .N]
  if (n_implaus_hi + n_implaus_lo > 0) {
    cat(sprintf("  Implausible growth >200%%        : %d obs (will be winsorized)\n",
                n_implaus_hi))
    cat(sprintf("  Implausible growth < -90%%       : %d obs (will be winsorized)\n",
                n_implaus_lo))
  }

  cat("\n  Pre-winsorization stats:\n")
  describe_var(panel$insured_share_growth, "insured_share_growth", "pre-winsor")
  log_checkpoint("CK-12: insured_share_growth constructed", "PASS",
                 sprintf("first-obs NAs: %d, zero-lag dropped: %d",
                         n_first_obs, n_zero_lag))
} else {
  panel[, insured_share_growth := NA_real_]
  log_checkpoint("CK-12: insured_share_growth constructed", "WARN",
                 "acct_018 missing — deposit growth will be all NA")
}

# --------------------------------------------------------------------------
# 6d. member_growth_yoy = Δacct_083 / acct_083_lag4 (YoY %)
# --------------------------------------------------------------------------
cat("\n  [6d] member_growth_yoy — Membership Growth (YoY %)\n")
cat("       Formula : (acct_083_t - acct_083_{t-4}) / acct_083_{t-4} * 100\n")
cat("       Logic   : YoY rate; lag-4 = NA for first 4 obs per CU\n\n")

if ("acct_083" %in% names(panel)) {
  panel[, acct_083_lag4 := shift(acct_083, 4L, type = "lag"), by = join_number]

  n_first_obs4 <- panel[is.na(acct_083_lag4), .N]
  n_zero_lag4  <- panel[!is.na(acct_083_lag4) & acct_083_lag4 <= 0, .N]

  cat(sprintf("  First-4-obs NA (no 4-lag)       : %s rows\n",
              format(n_first_obs4, big.mark = ",")))
  cat(sprintf("  Zero lagged members             : %s rows → NA\n",
              format(n_zero_lag4, big.mark = ",")))

  panel[, member_growth_yoy := fifelse(
    !is.na(acct_083_lag4) & acct_083_lag4 > 0,
    (acct_083 - acct_083_lag4) / acct_083_lag4 * 100,
    NA_real_
  )]

  cat("\n  Pre-winsorization stats:\n")
  describe_var(panel$member_growth_yoy, "member_growth_yoy", "pre-winsor")
  log_checkpoint("CK-13: member_growth_yoy constructed", "PASS",
                 sprintf("first-4-obs NAs: %d", n_first_obs4))
} else {
  panel[, member_growth_yoy := NA_real_]
  log_checkpoint("CK-13: member_growth_yoy constructed", "WARN",
                 "acct_083 missing")
}

# --------------------------------------------------------------------------
# 6e. costfds = acct_131 / acct_018 * 4 (annualized rate)
# --------------------------------------------------------------------------
cat("\n  [6e] costfds — Cost of Funds (annualized)\n")
cat("       Formula : acct_131 / acct_018 * 4\n")
cat("       Logic   : Interest expense / Deposits × 4 = annualized rate\n\n")

if (all(c("acct_131", "acct_018") %in% names(panel))) {
  n_zero_dep <- panel[acct_018 <= 0, .N]
  n_neg_exp  <- panel[acct_131 < 0, .N]
  cat(sprintf("  Zero/neg deposits (denominator) : %s rows → NA\n",
              format(n_zero_dep, big.mark = ",")))
  cat(sprintf("  Negative interest expense       : %s rows\n",
              format(n_neg_exp, big.mark = ",")))

  panel[, costfds := fifelse(acct_018 > 0, acct_131 / acct_018 * 4, NA_real_)]
  cat("\n  Pre-winsorization stats:\n")
  describe_var(panel$costfds, "costfds", "pre-winsor")
  log_checkpoint("CK-14: costfds constructed", "PASS",
                 sprintf("zero-denom: %d", n_zero_dep))
} else if ("costfds" %in% names(panel)) {
  panel[, costfds := as.double(costfds)]
  cat("  Using pre-computed costfds column from raw data\n")
  describe_var(panel$costfds, "costfds", "pre-computed")
  log_checkpoint("CK-14: costfds constructed", "INFO",
                 "Used pre-computed column (acct_131/018 not available)")
} else {
  panel[, costfds := NA_real_]
  log_checkpoint("CK-14: costfds constructed", "WARN", "All accounts missing")
}

# --------------------------------------------------------------------------
# 6f. loan_to_share = acct_025 / dep_tot  (fallback: acct_025 / acct_018)
# --------------------------------------------------------------------------
cat("\n  [6f] loan_to_share — Loan-to-Share Ratio\n")
cat("       Formula : acct_025 / dep_tot  [fallback: acct_025 / acct_018]\n")
cat("       Logic   : dep_tot preferred; acct_018 fallback if dep_tot missing\n\n")

if (all(c("acct_025", "dep_tot") %in% names(panel))) {
  n_zero_dt <- panel[dep_tot <= 0, .N]
  cat(sprintf("  Using dep_tot (preferred path)\n"))
  cat(sprintf("  dep_tot <= 0                    : %s rows → NA\n",
              format(n_zero_dt, big.mark = ",")))
  panel[, loan_to_share := fifelse(dep_tot > 0, acct_025 / dep_tot, NA_real_)]
  denom_used <- "dep_tot"
} else if (all(c("acct_025", "acct_018") %in% names(panel))) {
  n_zero_18 <- panel[acct_018 <= 0, .N]
  cat(sprintf("  dep_tot not found — using acct_018 FALLBACK\n"))
  cat(sprintf("  acct_018 <= 0                   : %s rows → NA\n",
              format(n_zero_18, big.mark = ",")))
  panel[, loan_to_share := fifelse(acct_018 > 0, acct_025 / acct_018, NA_real_)]
  denom_used <- "acct_018 (FALLBACK)"
  log_checkpoint("CK-15a: loan_to_share denominator", "WARN",
                 "dep_tot not found — using acct_018 fallback")
} else {
  panel[, loan_to_share := NA_real_]
  denom_used <- "NONE"
}

# Sanity: loan-to-share > 3 is economically unusual
n_extreme_lts <- panel[!is.na(loan_to_share) & loan_to_share > 3, .N]
cat(sprintf("  Denominator used                : %s\n", denom_used))
cat(sprintf("  loan_to_share > 3.0 (unusual)   : %d obs\n", n_extreme_lts))
cat("\n  Pre-winsorization stats:\n")
describe_var(panel$loan_to_share, "loan_to_share", "pre-winsor")
log_checkpoint("CK-15: loan_to_share constructed", "PASS",
               sprintf("denom: %s", denom_used))

# --------------------------------------------------------------------------
# 6g. pcanetworth = acct_998 / acct_010
# --------------------------------------------------------------------------
cat("\n  [6g] pcanetworth — Net Worth Ratio\n")
cat("       Formula : acct_998 / acct_010\n")
cat("       Logic   : Net worth / Average assets; NCUA PCA thresholds at 6%/7%/10%\n\n")

if (all(c("acct_998", "acct_010") %in% names(panel))) {
  n_zero_assets <- panel[acct_010 <= 0, .N]
  n_neg_nw      <- panel[acct_998 < 0, .N]
  cat(sprintf("  Zero/neg assets (denominator)   : %s rows → NA\n",
              format(n_zero_assets, big.mark = ",")))
  cat(sprintf("  Negative net worth (undercap'd) : %s rows\n",
              format(n_neg_nw, big.mark = ",")))

  panel[, pcanetworth := fifelse(acct_010 > 0, acct_998 / acct_010, NA_real_)]

  # NCUA PCA context: flag critically undercapitalized (<2%)
  n_critical <- panel[!is.na(pcanetworth) & pcanetworth < 0.02, .N]
  if (n_critical > 0) {
    cat(sprintf("  Critically undercapitalized (<2%%): %d CU-quarters (informational)\n",
                n_critical))
  }

  cat("\n  Pre-winsorization stats:\n")
  describe_var(panel$pcanetworth, "pcanetworth", "pre-winsor")
  log_checkpoint("CK-16: pcanetworth constructed", "PASS",
                 sprintf("zero-denom: %d, neg NW: %d", n_zero_assets, n_neg_nw))
} else if ("pcanetworth" %in% names(panel)) {
  panel[, pcanetworth := as.double(pcanetworth)]
  cat("  Using pre-computed pcanetworth column\n")
  describe_var(panel$pcanetworth, "pcanetworth", "pre-computed")
  log_checkpoint("CK-16: pcanetworth constructed", "INFO",
                 "Used pre-computed column")
} else {
  panel[, pcanetworth := NA_real_]
  log_checkpoint("CK-16: pcanetworth constructed", "WARN", "All accounts missing")
}

# --------------------------------------------------------------------------
# 6h. pll_rate = acct_674 / avg(acct_025)
# CRITICAL: avg loan denominator confirmed in v1 Script 04d
# --------------------------------------------------------------------------
cat("\n  [6h] pll_rate — PLL Rate\n")
cat("       Formula : acct_674 / [(acct_025_t + acct_025_{t-1}) / 2]\n")
cat("       Logic   : Provision / Average loans; confirmed construction from v1 04d\n\n")

if ("acct_025" %in% names(panel)) {
  panel[, lns_tot_lag := shift(acct_025, 1L, type = "lag"), by = join_number]
  panel[, lns_tot_avg := (acct_025 + lns_tot_lag) / 2]

  n_avg_na <- panel[is.na(lns_tot_avg), .N]
  n_avg_zero <- panel[!is.na(lns_tot_avg) & lns_tot_avg <= 0, .N]
  cat(sprintf("  lns_tot_avg NA (first obs)      : %s rows\n",
              format(n_avg_na, big.mark = ",")))
  cat(sprintf("  lns_tot_avg <= 0                : %s rows → NA\n",
              format(n_avg_zero, big.mark = ",")))
} else {
  panel[, lns_tot_avg := NA_real_]
  n_avg_na <- nrow(panel); n_avg_zero <- 0
}

if ("acct_674" %in% names(panel)) {
  n_neg_pll <- panel[acct_674 < 0, .N]
  cat(sprintf("  Negative PLL (recoveries > prov): %d obs\n", n_neg_pll))
  panel[, pll_rate := fifelse(!is.na(lns_tot_avg) & lns_tot_avg > 0,
                              acct_674 / lns_tot_avg,
                              NA_real_)]
  cat("\n  Pre-winsorization stats:\n")
  describe_var(panel$pll_rate, "pll_rate", "pre-winsor")
  log_checkpoint("CK-17: pll_rate constructed", "PASS",
                 sprintf("avg-loan NAs: %d, zero-denom: %d", n_avg_na, n_avg_zero))
} else if ("pll_rate" %in% names(panel)) {
  panel[, pll_rate := as.double(pll_rate)]
  cat("  Using pre-computed pll_rate column\n")
  describe_var(panel$pll_rate, "pll_rate", "pre-computed")
  log_checkpoint("CK-17: pll_rate constructed", "INFO", "Used pre-computed column")
} else {
  panel[, pll_rate := NA_real_]
  log_checkpoint("CK-17: pll_rate constructed", "WARN",
                 "acct_674 missing and no pre-computed column")
}

# --------------------------------------------------------------------------
# 6i. chg_totlns_ratio — charge-off column name audit
# --------------------------------------------------------------------------
cat("\n  [6i] chg_totlns_ratio — Charge-off Ratio (column audit only)\n")
if ("chg_totlns_ratio" %in% names(panel)) {
  cat("  Found: chg_totlns_ratio (exact v1 name — CORRECT)\n")
  log_checkpoint("CK-18: chg_totlns_ratio name", "PASS", "Exact v1 column name found")
} else if ("chg_tot_lns_ratio" %in% names(panel)) {
  cat("  Found variant: chg_tot_lns_ratio → renaming to chg_totlns_ratio\n")
  setnames(panel, "chg_tot_lns_ratio", "chg_totlns_ratio")
  log_checkpoint("CK-18: chg_totlns_ratio name", "WARN",
                 "Renamed from chg_tot_lns_ratio (common mis-naming from v1)")
} else {
  panel[, chg_totlns_ratio := NA_real_]
  log_checkpoint("CK-18: chg_totlns_ratio name", "WARN",
                 "Column not found — set to NA (not a primary outcome variable)")
}

# ============================================================================
# SECTION 7: PRE-WINSORIZATION MISSINGNESS REPORT
# ============================================================================
section("7", "Pre-Winsorization Missingness Report")

cat("  Outcome missingness BEFORE winsorization:\n\n")
cat(sprintf("  %-28s  %8s  %8s  %8s\n", "Variable", "N (valid)", "N (miss)", "% Miss"))
cat(sprintf("  %s\n", paste(rep("-", 62), collapse = "")))

present_outcomes <- intersect(OUTCOMES, names(panel))
for (v in present_outcomes) {
  n_valid <- sum(!is.na(panel[[v]]))
  n_miss  <- sum(is.na(panel[[v]]))
  pct_m   <- n_miss / nrow(panel) * 100
  flag    <- if (pct_m > 20) "  ← HIGH" else if (pct_m > 5) "  ← NOTE" else ""
  cat(sprintf("  %-28s  %8s  %8s  %7.2f%%%s\n",
              v,
              format(n_valid, big.mark = ","),
              format(n_miss,  big.mark = ","),
              pct_m, flag))
}

# CK: every outcome must be constructible
n_all_na <- sum(sapply(present_outcomes, function(v) all(is.na(panel[[v]]))))
if (n_all_na > 0) {
  log_checkpoint("CK-19: Outcome construction", "WARN",
                 sprintf("%d outcomes are entirely NA", n_all_na))
} else {
  log_checkpoint("CK-19: Outcome construction", "PASS",
                 "All outcomes have at least some valid values")
}

# ============================================================================
# SECTION 8: OUTLIER TREATMENT (WINSORIZATION)
# ============================================================================
section("8", "Outlier Treatment — Winsorization")

cat(sprintf("  Method      : Winsorization at [%.0f%%, %.0f%%] percentiles\n",
            WIN_LOW * 100, WIN_HIGH * 100))
cat("  Scope       : Applied to each outcome variable across full panel\n")
cat("  Rationale   : Preserves distributional shape; avoids deletion of valid extremes\n")
cat("                Outlier observations are clipped to boundary, not dropped\n")
cat("  Auditability: Every observation affected is logged to Results/01_outlier_log.csv\n\n")

cat(sprintf("  %-28s  %8s  %8s  %8s  %8s  %8s  %8s\n",
            "Variable", "p1 (lo)", "p99 (hi)", "N clipped", "N clipped",
            "% aff.", "Change"))
cat(sprintf("  %-28s  %8s  %8s  %8s  %8s  %8s  %8s\n",
            "", "boundary", "boundary", "below p1", "above p99", "", "in mean"))
cat(sprintf("  %s\n", paste(rep("-", 86), collapse = "")))

OUTLIER_RECORDS <- list()
PRE_MEANS  <- list()
POST_MEANS <- list()

for (v in present_outcomes) {
  x_pre  <- panel[[v]]
  if (all(is.na(x_pre))) {
    cat(sprintf("  %-28s  ALL MISSING — skipped\n", v))
    next
  }

  PRE_MEANS[[v]]  <- mean(x_pre, na.rm = TRUE)
  w               <- winsorize_audit(x_pre, v, WIN_LOW, WIN_HIGH)
  panel[, (v) := w$x]
  POST_MEANS[[v]] <- mean(w$x, na.rm = TRUE)
  delta_mean      <- POST_MEANS[[v]] - PRE_MEANS[[v]]

  cat(sprintf("  %-28s  %8.4f  %8.4f  %8d  %8d  %7.3f%%  %+8.5f\n",
              v,
              w$lo, w$hi,
              w$n_lo, w$n_hi,
              w$pct_affected,
              delta_mean))

  # Build audit records for observations clipped BELOW p1
  if (w$n_lo > 0) {
    idx_lo <- which(!is.na(x_pre) & x_pre < w$lo)
    OUTLIER_RECORDS[[paste0(v, "_lo")]] <- data.table(
      variable     = v,
      direction    = "below_p1",
      join_number  = panel$join_number[idx_lo],
      yyyyqq       = panel$yyyyqq[idx_lo],
      raw_value    = x_pre[idx_lo],
      clipped_to   = w$lo,
      clip_amount  = x_pre[idx_lo] - w$lo
    )
  }

  # Build audit records for observations clipped ABOVE p99
  if (w$n_hi > 0) {
    idx_hi <- which(!is.na(x_pre) & x_pre > w$hi)
    OUTLIER_RECORDS[[paste0(v, "_hi")]] <- data.table(
      variable     = v,
      direction    = "above_p99",
      join_number  = panel$join_number[idx_hi],
      yyyyqq       = panel$yyyyqq[idx_hi],
      raw_value    = x_pre[idx_hi],
      clipped_to   = w$hi,
      clip_amount  = x_pre[idx_hi] - w$hi
    )
  }
}

# Save outlier audit log
outlier_log <- rbindlist(OUTLIER_RECORDS, fill = TRUE)
fwrite(outlier_log, file.path(RESULTS_DIR, "01_outlier_log.csv"))
cat(sprintf("\n  Total obs winsorized        : %s\n",
            format(nrow(outlier_log), big.mark = ",")))
cat(sprintf("  Outlier log saved to        : Results/01_outlier_log.csv\n"))
cat(sprintf("  Columns in log              : %s\n",
            paste(names(outlier_log), collapse = ", ")))

log_checkpoint("CK-20: Winsorization complete", "PASS",
               sprintf("%s obs winsorized across %d outcomes",
                       format(nrow(outlier_log), big.mark = ","),
                       length(present_outcomes)))

# ============================================================================
# SECTION 9: POST-WINSORIZATION SUMMARY STATISTICS
# ============================================================================
section("9", "Post-Winsorization Summary Statistics")

cat("  Full distribution for each outcome AFTER winsorization:\n\n")

VAR_SUMMARY <- list()

for (v in present_outcomes) {
  x <- panel[[v]]
  if (all(is.na(x))) next

  n      <- sum(!is.na(x))
  nmiss  <- mean(is.na(x))
  qs     <- quantile(x, probs = c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99),
                     na.rm = TRUE)
  mn_val <- mean(x, na.rm = TRUE)
  sd_val <- sd(x, na.rm = TRUE)
  sk_val <- (mn_val - qs["50%"]) / sd_val   # simple skewness proxy

  # Check against economic expected range
  meta       <- OUTCOME_META[[v]]
  exp_lo     <- meta$expected_range[1]
  exp_hi     <- meta$expected_range[2]
  pct_in_range <- mean(x >= exp_lo & x <= exp_hi, na.rm = TRUE) * 100

  cat(sprintf("  %s — %s\n", v, meta$label))
  cat(sprintf("    N=%-8s  miss=%.1f%%  mean=%8.5f  sd=%8.5f  skew=%+.3f\n",
              format(n, big.mark = ","), nmiss * 100, mn_val, sd_val, sk_val))
  cat(sprintf("    p01=%8.5f  p10=%8.5f  p25=%8.5f  p50=%8.5f\n",
              qs[1], qs[2], qs[3], qs[4]))
  cat(sprintf("    p75=%8.5f  p90=%8.5f  p99=%8.5f\n",
              qs[5], qs[6], qs[7]))
  cat(sprintf("    Expected range: [%.4f, %.4f]  |  %.1f%% of obs within range\n",
              exp_lo, exp_hi, pct_in_range))
  if (pct_in_range < 80) {
    cat(sprintf("    *** WARNING: Only %.1f%% within expected range — check construction\n",
                pct_in_range))
  }
  cat("\n")

  VAR_SUMMARY[[v]] <- data.table(
    variable = v, label = meta$label,
    n = n, pct_miss = round(nmiss * 100, 2),
    mean = round(mn_val, 6), sd = round(sd_val, 6),
    p01 = round(qs[1], 6), p10 = round(qs[2], 6),
    p25 = round(qs[3], 6), p50 = round(qs[4], 6),
    p75 = round(qs[5], 6), p90 = round(qs[6], 6),
    p99 = round(qs[7], 6),
    expected_lo = exp_lo, expected_hi = exp_hi,
    pct_in_expected_range = round(pct_in_range, 2)
  )
}

var_summary_dt <- rbindlist(VAR_SUMMARY, fill = TRUE)
fwrite(var_summary_dt, file.path(RESULTS_DIR, "01_variable_summary.csv"))
cat(sprintf("  Saved: Results/01_variable_summary.csv\n"))

log_checkpoint("CK-21: Post-winsorization stats", "PASS",
               sprintf("%d outcome variables summarized", nrow(var_summary_dt)))

# ============================================================================
# SECTION 10: LAGGED OUTCOMES FOR AR(1) STRUCTURE
# ============================================================================
section("10", "Lagged Outcomes (AR(1) Structure)")

cat("  Constructing one-period lags for each outcome variable.\n")
cat("  ANCHOR: v1 AR(1) coefficient ≈ 0.59 — lags required for Script 04a.\n\n")

cat(sprintf("  %-28s  %8s  %8s  %8s\n",
            "Variable", "Lag name", "N (valid)", "% of parent"))
cat(sprintf("  %s\n", paste(rep("-", 62), collapse = "")))

for (v in present_outcomes) {
  lag_nm <- paste0(v, "_lag1")
  panel[, (lag_nm) := shift(get(v), 1L, type = "lag"), by = join_number]

  n_parent <- sum(!is.na(panel[[v]]))
  n_lag    <- sum(!is.na(panel[[lag_nm]]))
  pct_lag  <- if (n_parent > 0) n_lag / n_parent * 100 else 0

  cat(sprintf("  %-28s  %-28s  %8s  %7.1f%%\n",
              v, lag_nm,
              format(n_lag, big.mark = ","),
              pct_lag))
}

log_checkpoint("CK-22: Lag construction", "PASS",
               sprintf("%d outcome lags constructed", length(present_outcomes)))

# ============================================================================
# SECTION 11: FLAGS — FIRST OBS AND CECL
# ============================================================================
section("11", "Observation Flags (first_obs_flag, cecl_flag)")

# first_obs_flag: any lag unavailable → row should not enter AR models
panel[, first_obs_flag := as.integer(
  is.na(insured_share_growth_lag1) | is.na(dq_rate_lag1)
)]
n_first <- panel[first_obs_flag == 1, .N]
n_usable <- panel[first_obs_flag == 0, .N]

cat(sprintf("  first_obs_flag == 1 (lag missing, exclude from AR models): %s rows (%.1f%%)\n",
            format(n_first, big.mark = ","),
            n_first / nrow(panel) * 100))
cat(sprintf("  first_obs_flag == 0 (usable for AR models)               : %s rows (%.1f%%)\n",
            format(n_usable, big.mark = ","),
            n_usable / nrow(panel) * 100))

# cecl_flag: CECL transition creates accounting break in PLL/charge-offs
# Large CUs: 2020Q1; Small CUs (FICUs): 2023Q1
panel[, cecl_flag := as.integer(
  (date >= as.Date("2020-01-01") & date <= as.Date("2020-03-31")) |
  (date >= as.Date("2023-01-01") & date <= as.Date("2023-03-31"))
)]
n_cecl <- panel[cecl_flag == 1, .N]
cat(sprintf("\n  cecl_flag == 1 (CECL transition quarters):                 %s rows (%.2f%%)\n",
            format(n_cecl, big.mark = ","),
            n_cecl / nrow(panel) * 100))
cat("  CECL quarters flagged: 2020Q1 (large CU adoption), 2023Q1 (FICU adoption)\n")
cat("  Recommendation: Exclude cecl_flag==1 from PLL/charge-off regressions\n")

log_checkpoint("CK-23: Observation flags", "PASS",
               sprintf("first_obs=%d, cecl=%d", n_first, n_cecl))

# ============================================================================
# SECTION 12: ASSET TIERS
# ============================================================================
section("12", "Asset Tier Assignment")

cat("  Tier thresholds (call report units = $thousands):\n")
cat("    T1 : < $50M           (<     50,000)\n")
cat("    T2 : $50M – $100M     (50,000 – 100,000)\n")
cat("    T3 : $100M – $500M    (100,000 – 500,000)\n")
cat("    T4 : $500M – $1B      (500,000 – 1,000,000)\n")
cat("    T5 : > $1B            (>= 1,000,000)\n\n")

ASSET_COL <- if ("acct_010" %in% names(panel)) "acct_010" else
             if ("total_assets" %in% names(panel)) "total_assets" else NA_character_

if (!is.na(ASSET_COL)) {
  cat(sprintf("  Asset column used: %s\n\n", ASSET_COL))

  describe_var(panel[[ASSET_COL]], ASSET_COL, "raw assets")

  panel[, asset_tier := fcase(
    get(ASSET_COL) <  50000,                                 "T1",
    get(ASSET_COL) >= 50000  & get(ASSET_COL) < 100000,     "T2",
    get(ASSET_COL) >= 100000 & get(ASSET_COL) < 500000,     "T3",
    get(ASSET_COL) >= 500000 & get(ASSET_COL) < 1000000,    "T4",
    get(ASSET_COL) >= 1000000,                               "T5",
    default = NA_character_
  )]

  tier_dist <- panel[, .N, by = asset_tier][order(asset_tier)]
  tier_dist[, pct := round(N / nrow(panel) * 100, 1)]

  cat("\n  Asset tier distribution:\n")
  cat(sprintf("  %-6s  %10s  %8s  %s\n", "Tier", "N (CU-qtrs)", "Pct", "Range"))
  cat(sprintf("  %s\n", paste(rep("-", 50), collapse = "")))
  tier_labels <- c(T1 = "<$50M", T2 = "$50-100M", T3 = "$100-500M",
                   T4 = "$500M-$1B", T5 = ">$1B")
  for (i in seq_len(nrow(tier_dist))) {
    tier <- tier_dist$asset_tier[i]
    if (!is.na(tier)) {
      cat(sprintf("  %-6s  %10s  %7.1f%%  %s\n",
                  tier, format(tier_dist$N[i], big.mark = ","),
                  tier_dist$pct[i],
                  tier_labels[tier]))
    }
  }
  n_na_tier <- panel[is.na(asset_tier), .N]
  if (n_na_tier > 0) cat(sprintf("  %-6s  %10s  %7.1f%%  (zero/NA assets)\n",
                                  "NA", format(n_na_tier, big.mark = ","),
                                  n_na_tier / nrow(panel) * 100))

  log_checkpoint("CK-24: Asset tier assignment", "PASS",
                 sprintf("5 tiers assigned, %d NA", n_na_tier))
} else {
  panel[, asset_tier := NA_character_]
  log_checkpoint("CK-24: Asset tier assignment", "WARN",
                 "No asset column found — all tiers set to NA")
}

# ============================================================================
# SECTION 13: TIME VARIABLES AND STRUCTURAL BREAK FLAG
# ============================================================================
section("13", "Time Variables and Structural Break Flag")

panel[, year     := year(date)]
panel[, qtr      := quarter(date)]
panel[, post2015 := as.integer(date >= as.Date("2015-01-01"))]
panel[, t        := as.integer(factor(date))]   # integer time index (1, 2, 3...)

cat("  Variables created:\n")
cat("    year     : calendar year from date\n")
cat("    qtr      : quarter (1–4)\n")
cat("    post2015 : 1 if date >= 2015-01-01  [STRUCTURAL BREAK — v1 finding]\n")
cat("    t        : integer time index (1=earliest quarter)\n\n")

cat(sprintf("  post2015 == 0 (pre-shale revolution) : %s CU-quarters\n",
            format(panel[post2015 == 0, .N], big.mark = ",")))
cat(sprintf("  post2015 == 1 (post-shale revolution): %s CU-quarters\n",
            format(panel[post2015 == 1, .N], big.mark = ",")))
cat(sprintf("  First t value : %d (= %s)\n", min(panel$t), format(min(panel$date))))
cat(sprintf("  Last  t value : %d (= %s)\n", max(panel$t), format(max(panel$date))))

# Confirm 2015Q1 cutoff falls in the data
if (as.Date("2015-01-01") < date_min || as.Date("2015-01-01") > date_max) {
  log_checkpoint("CK-25: 2015Q1 in date range", "WARN",
                 "2015Q1 outside data range — post2015 flag may be all 0 or all 1")
} else {
  log_checkpoint("CK-25: 2015Q1 in date range", "PASS",
                 sprintf("post2015: %d pre, %d post",
                         panel[post2015 == 0, .N],
                         panel[post2015 == 1, .N]))
}

# ============================================================================
# SECTION 14: PANEL BALANCE DIAGNOSTICS
# ============================================================================
section("14", "Panel Balance Diagnostics")

cat("  Checking for unbalanced panel issues...\n\n")

# Observations per CU
obs_per_cu <- panel[, .N, by = join_number]
cat(sprintf("  Obs per CU — min: %d  median: %.0f  mean: %.1f  max: %d\n",
            min(obs_per_cu$N), median(obs_per_cu$N),
            mean(obs_per_cu$N), max(obs_per_cu$N)))

# Short-panel CUs (fewer than 4 quarters — insufficient for lags)
n_short <- obs_per_cu[N < 4, .N]
cat(sprintf("  CUs with < 4 quarters (lag-limited): %d (%.1f%% of CUs)\n",
            n_short, n_short / nrow(obs_per_cu) * 100))

# CUs per quarter
cu_per_qtr <- panel[, .(n_cu = uniqueN(join_number)), by = yyyyqq]
setorder(cu_per_qtr, yyyyqq)
cat(sprintf("\n  CUs per quarter — min: %d  median: %.0f  mean: %.1f  max: %d\n",
            min(cu_per_qtr$n_cu), median(cu_per_qtr$n_cu),
            mean(cu_per_qtr$n_cu), max(cu_per_qtr$n_cu)))

# Quarters with big drops in CU count (potential data issues)
cu_per_qtr[, n_cu_lag := shift(n_cu, 1L)]
cu_per_qtr[, drop_pct := (n_cu - n_cu_lag) / n_cu_lag * 100]
big_drops <- cu_per_qtr[!is.na(drop_pct) & drop_pct < -5]
if (nrow(big_drops) > 0) {
  cat("\n  Quarters with >5% drop in CU count (investigate):\n")
  print(big_drops[, .(yyyyqq, n_cu, n_cu_lag, drop_pct)])
  log_checkpoint("CK-26: Panel CU count stability", "WARN",
                 sprintf("%d quarters with >5%% CU count drop", nrow(big_drops)))
} else {
  log_checkpoint("CK-26: Panel CU count stability", "PASS",
                 "No quarters with >5% drop in CU count")
}

# State coverage
if ("reporting_state" %in% names(panel)) {
  n_states_panel <- uniqueN(panel$reporting_state)
  cat(sprintf("\n  States represented: %d\n", n_states_panel))
  if (n_states_panel < 48) {
    log_checkpoint("CK-27: State coverage", "WARN",
                   sprintf("Only %d states — may affect oil-state classification", n_states_panel))
  } else {
    log_checkpoint("CK-27: State coverage", "PASS",
                   sprintf("%d states", n_states_panel))
  }
}

# ============================================================================
# SECTION 15: FINAL PANEL SUMMARY
# ============================================================================
section("15", "Final Panel Summary")

cat(sprintf("  Total rows             : %s\n", format(nrow(panel), big.mark = ",")))
cat(sprintf("  Total columns          : %d\n", ncol(panel)))
cat(sprintf("  Unique CUs             : %s\n", format(uniqueN(panel$join_number), big.mark = ",")))
cat(sprintf("  Unique quarters        : %d\n", uniqueN(panel$yyyyqq)))
cat(sprintf("  Date range             : %s  to  %s\n", format(date_min), format(date_max)))
cat(sprintf("  Object size            : %.1f MB\n", object.size(panel) / 1e6))
cat(sprintf("  post2015 == 0 rows     : %s\n", format(panel[post2015==0, .N], big.mark=",")))
cat(sprintf("  post2015 == 1 rows     : %s\n", format(panel[post2015==1, .N], big.mark=",")))
cat(sprintf("  first_obs_flag == 1    : %s\n", format(panel[first_obs_flag==1, .N], big.mark=",")))
cat(sprintf("  cecl_flag == 1         : %s\n", format(panel[cecl_flag==1, .N], big.mark=",")))

cat("\n  Outcome variables in final panel:\n")
cat(sprintf("  %-28s  %8s  %8s\n", "Variable", "N (valid)", "% Valid"))
cat(sprintf("  %s\n", paste(rep("-", 52), collapse = "")))
for (v in present_outcomes) {
  n_v   <- sum(!is.na(panel[[v]]))
  pct_v <- n_v / nrow(panel) * 100
  cat(sprintf("  %-28s  %8s  %7.1f%%\n",
              v, format(n_v, big.mark = ","), pct_v))
}

cat("\n  Lag variables in final panel:\n")
cat(sprintf("  %-28s  %8s  %8s\n", "Variable", "N (valid)", "% Valid"))
cat(sprintf("  %s\n", paste(rep("-", 52), collapse = "")))
for (v in present_outcomes) {
  lag_nm <- paste0(v, "_lag1")
  if (lag_nm %in% names(panel)) {
    n_v   <- sum(!is.na(panel[[lag_nm]]))
    pct_v <- n_v / nrow(panel) * 100
    cat(sprintf("  %-28s  %8s  %7.1f%%\n",
                lag_nm, format(n_v, big.mark = ","), pct_v))
  }
}

# ============================================================================
# SECTION 16: CHECKPOINT SUMMARY
# ============================================================================
section("16", "Checkpoint Summary")

n_pass <- sum(CHKPT_LOG$status == "PASS")
n_warn <- sum(CHKPT_LOG$status == "WARN")
n_fail <- sum(CHKPT_LOG$status == "FAIL")
n_info <- sum(CHKPT_LOG$status == "INFO")

cat(sprintf("  Total checkpoints  : %d\n", nrow(CHKPT_LOG)))
cat(sprintf("  PASS               : %d\n", n_pass))
cat(sprintf("  WARN               : %d\n", n_warn))
cat(sprintf("  FAIL               : %d  (script would have halted)\n", n_fail))
cat(sprintf("  INFO               : %d\n", n_info))

if (n_warn > 0) {
  cat("\n  Warnings requiring attention:\n")
  for (i in which(CHKPT_LOG$status == "WARN")) {
    cat(sprintf("    ! %s\n", CHKPT_LOG$checkpoint[i]))
    if (nchar(CHKPT_LOG$detail[i]) > 0)
      cat(sprintf("      └─ %s\n", CHKPT_LOG$detail[i]))
  }
}

fwrite(as.data.table(CHKPT_LOG),
       file.path(RESULTS_DIR, "01_checkpoint_log.csv"))
cat(sprintf("\n  Checkpoint log saved: Results/01_checkpoint_log.csv\n"))

# ============================================================================
# SECTION 17: SAVE
# ============================================================================
section("17", "Save Output")

if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)

t_save <- system.time({
  saveRDS(panel, file.path(DATA_DIR, "01_panel_clean.rds"))
})

fsize <- round(file.size(file.path(DATA_DIR, "01_panel_clean.rds")) / 1e6, 1)
cat(sprintf("  Saved : Data/01_panel_clean.rds\n"))
cat(sprintf("  Size  : %.1f MB\n", fsize))
cat(sprintf("  Time  : %.1f seconds\n", t_save["elapsed"]))

log_checkpoint("CK-28: File saved", "PASS",
               sprintf("Data/01_panel_clean.rds — %.1f MB", fsize))

# ============================================================================
# CLOSING BANNER
# ============================================================================
cat("\n", SEP, "\n", sep = "")
cat("  SCRIPT 01 COMPLETE\n")
cat(sprintf("  Finished  : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("\n  Output files:\n")
cat("    Data/01_panel_clean.rds          — Main panel for Script 02\n")
cat("    Results/01_outlier_log.csv       — Every obs winsorized (CU + quarter + value)\n")
cat("    Results/01_variable_summary.csv  — Pre/post-winsorization stats\n")
cat("    Results/01_checkpoint_log.csv    — All checkpoint pass/warn/fail records\n")
cat(sprintf("\n  Checkpoints: %d PASS | %d WARN | %d FAIL\n", n_pass, n_warn, n_fail))
cat(SEP, "\n", sep = "")
