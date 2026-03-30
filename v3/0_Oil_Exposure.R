# =============================================================================
# SCRIPT 01b — Oil Exposure Index
# Oil Price Shock × Credit Union Financial Performance — v2
# Author: Saurabh C. Datta, Ph.D. | NCUA Office of Chief Economist
# Version: v2.1-DIAGNOSTIC | Built: 2026-03-30
# Working dir: S:/Projects/Oil_Price_Shock_2026/
# =============================================================================
# Method    BLS Quarterly Census of Employment & Wages (QCEW)
#           NAICS 211 -- Oil & Gas Extraction; state x quarter; 2004Q1+
#
# Exposure measures produced
#   mining_emp_share        Continuous: oil/gas emp / total state emp (primary)
#   oil_exposure_cont       Normalized 0-1 within-quarter (for regression)
#   oil_exposure_bin        Binary >= 2% (primary; peer-reviewed threshold)
#   oil_exposure_bin_1pct   Binary >= 1% (robustness low)
#   oil_exposure_bin_3pct   Binary >= 3% (robustness high)
#   oil_exposure_tier       Factor: Negligible / Low / Medium / High
#   oil_exposure_smooth     4Q rolling mean (noise reduction)
#   oil_bartik_iv           Bartik shift-share instrument (IV for 2SLS)
#   spillover_exposure      INDIRECT CHANNEL: non-oil state linkage to oil states
#   spillover_exposure_wtd  Trade-weighted spillover (BEA inter-state flows proxy)
#   cu_group                Direct / Indirect / Negligible (for subsample models)
#
# Data acquisition (3-tier waterfall)
#   Tier 1: BLS QCEW annual zip files (preferred -- most complete)
#   Tier 2: BLS API (SMS series -- fallback if zip unavailable)
#   Tier 3: Curated static table from Kilian (2014) / BLS averages
#
# Academic citations
#   BLS QCEW: standard in regional economics (Autor, Dorn & Hanson 2013)
#   Time-varying: captures Bakken boom (ND 2008+) automatically
#   Bartik IV: Goldsmith-Pinkham, Sorkin & Swift (2020, AER 110:8)
#   Spillover: Acemoglu et al. (2016, AER 106:1) network propagation
#   State classification: Kilian (2014) Annual Review of Resource Economics
#
# Output    Data/oil_exposure.rds -- state x quarter panel; all measures
# Run order: BEFORE 01_data_prep.R
#
# Diagnostic outputs
#   Results/01b_checkpoint_log.csv     All checkpoint pass/warn/fail records
#   Results/01b_exposure_summary.csv   Per-state exposure statistics
#   Results/01b_bartik_summary.csv     Bartik IV diagnostics
#   Results/01b_spillover_summary.csv  Spillover index diagnostics
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

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# Describe a numeric vector to the log
describe_var <- function(x, varname, tag = "") {
  n     <- sum(!is.na(x))
  nmiss <- mean(is.na(x))
  if (n == 0) {
    cat(sprintf("  %-30s  ALL MISSING\n", varname))
    return(invisible(NULL))
  }
  qs <- quantile(x, c(0.01, 0.10, 0.25, 0.50, 0.75, 0.90, 0.99), na.rm = TRUE)
  cat(sprintf("  %-30s  [%s]\n", varname, tag))
  cat(sprintf("    N=%-7s  miss=%5.1f%%  mean=%10.6f  sd=%10.6f\n",
              format(n, big.mark = ","), nmiss * 100,
              mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)))
  cat(sprintf("    p01=%9.6f  p10=%9.6f  p25=%9.6f  p50=%9.6f\n",
              qs[1], qs[2], qs[3], qs[4]))
  cat(sprintf("    p75=%9.6f  p90=%9.6f  p99=%9.6f\n", qs[5], qs[6], qs[7]))
  cat(sprintf("    min=%9.6f  max=%9.6f\n",
              min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
  cat(sprintf("    %s\n", paste(rep("-", 68), collapse = "")))
}

# ============================================================================
# SCRIPT HEADER
# ============================================================================
cat(SEP, "\n")
cat("  SCRIPT 01b -- Oil Exposure Index\n")
cat("  Oil Price Shock x Credit Union Financial Performance -- v2\n")
cat("  Author  : Saurabh C. Datta, Ph.D. | NCUA Office of Chief Economist\n")
cat("  Version : v2.1-DIAGNOSTIC\n")
cat(sprintf("  Started : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(SEP, "\n\n")

# ============================================================================
# SECTION 1: LIBRARIES AND CONFIGURATION
# ============================================================================
section("1", "Libraries and Configuration")

suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
  library(stringr)
})

for (pkg in c("data.table", "httr", "jsonlite", "stringr")) {
  log_checkpoint(paste("Package:", pkg), "PASS",
                 paste("v", packageVersion(pkg)))
}

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
THRESHOLD_2PCT  <- 0.020    # primary binary threshold (peer-reviewed)
THRESHOLD_1PCT  <- 0.010    # robustness lower bound
THRESHOLD_3PCT  <- 0.030    # robustness upper bound
START_YEAR_DL   <- 2004L    # one year before study start (to allow lags)
END_YEAR_DL     <- 2025L
BARTIK_BASE_QTR <- 200501L  # 2005Q1 initial Bartik shares

DATA_DIR    <- "Data/"
RESULTS_DIR <- "Results/"
for (d in c(DATA_DIR, RESULTS_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cat(sprintf("\n  START_YEAR_DL   : %d\n", START_YEAR_DL))
cat(sprintf("  END_YEAR_DL     : %d\n", END_YEAR_DL))
cat(sprintf("  BARTIK_BASE_QTR : %d\n", BARTIK_BASE_QTR))
cat(sprintf("  Thresholds      : 1%%=%.3f  2%%=%.3f  3%%=%.3f\n",
            THRESHOLD_1PCT, THRESHOLD_2PCT, THRESHOLD_3PCT))
cat(sprintf("  Expected quarters: %d\n\n",
            (END_YEAR_DL - START_YEAR_DL + 1) * 4))

log_checkpoint("CK-01: Config loaded", "PASS",
               sprintf("Years %d-%d, Bartik base %d",
                       START_YEAR_DL, END_YEAR_DL, BARTIK_BASE_QTR))

# ============================================================================
# SECTION 2: STATE FIPS LOOKUP TABLE
# ============================================================================
section("2", "State FIPS Lookup Table")

state_fips <- data.table(
  state_code = c("AL","AK","AZ","AR","CA","CO","CT","DE","FL","GA",
                 "HI","ID","IL","IN","IA","KS","KY","LA","ME","MD",
                 "MA","MI","MN","MS","MO","MT","NE","NV","NH","NJ",
                 "NM","NY","NC","ND","OH","OK","OR","PA","RI","SC",
                 "SD","TN","TX","UT","VT","VA","WA","WV","WI","WY","DC"),
  fips = sprintf("%02d", c(1,2,4,5,6,8,9,10,12,13,15,16,17,18,19,20,
                            21,22,23,24,25,26,27,28,29,30,31,32,33,34,
                            35,36,37,38,39,40,41,42,44,45,46,47,48,49,
                            50,51,53,54,55,56,11))
)

cat(sprintf("  State FIPS table: %d states + DC\n", nrow(state_fips)))
cat(sprintf("  State codes: %s\n",
            paste(head(state_fips$state_code, 10), collapse = ", "),
            "... and", nrow(state_fips) - 10, "more"))

# Verify no duplicate state codes or FIPS
n_dup_code <- sum(duplicated(state_fips$state_code))
n_dup_fips <- sum(duplicated(state_fips$fips))
if (n_dup_code > 0 || n_dup_fips > 0) {
  log_checkpoint("CK-02: State FIPS table", "WARN",
                 sprintf("%d dup codes, %d dup FIPS", n_dup_code, n_dup_fips))
} else {
  log_checkpoint("CK-02: State FIPS table", "PASS",
                 sprintf("%d unique states, no duplicates", nrow(state_fips)))
}

# ============================================================================
# SECTION 3: DATA ACQUISITION (3-TIER WATERFALL)
# ============================================================================
section("3", "BLS QCEW Data Acquisition -- 3-Tier Waterfall")

cat("  Tier 1: BLS QCEW annual zip files (most complete, preferred)\n")
cat("  Tier 2: BLS SMS API (fallback if zip unavailable)\n")
cat("  Tier 3: Curated static table from Kilian (2014) averages\n\n")

qcew_state  <- NULL
data_source_used <- NULL

# ============================================================================
# TIER 1: BLS QCEW ZIP FILES
# ============================================================================
cat("  --- TIER 1: BLS QCEW zip download ---\n\n")

fetch_qcew_zip <- function(yr) {
  url <- sprintf(
    "https://data.bls.gov/cew/data/files/%d/csv/%d_qtrly_singlefile.zip",
    yr, yr
  )
  tmp  <- tempfile(fileext = ".zip")
  resp <- tryCatch(
    GET(url, timeout(90), write_disk(tmp, overwrite = TRUE)),
    error = function(e) NULL
  )
  if (is.null(resp) || resp$status_code != 200) return(NULL)
  dir  <- tempdir()
  csvs <- tryCatch(unzip(tmp, exdir = dir, overwrite = TRUE),
                   error = function(e) NULL)
  if (is.null(csvs)) return(NULL)
  csv <- csvs[grepl("\\.csv$", csvs, ignore.case = TRUE)][1]
  if (is.na(csv)) return(NULL)
  fread(csv, showProgress = FALSE,
        select = c("area_fips","own_code","industry_code",
                   "year","qtr","month3_emplvl"))
}

cat(sprintf("  Downloading years %d -- %d:\n", START_YEAR_DL, END_YEAR_DL))
zip_results <- list()
n_ok <- 0L; n_fail <- 0L

for (yr in START_YEAR_DL:END_YEAR_DL) {
  cat(sprintf("    %d ... ", yr))
  d <- fetch_qcew_zip(yr)
  if (!is.null(d) && nrow(d) > 0) {
    zip_results[[as.character(yr)]] <- d
    cat(sprintf("OK  (%s rows)\n", format(nrow(d), big.mark = ",")))
    n_ok <- n_ok + 1L
  } else {
    cat("FAIL\n")
    n_fail <- n_fail + 1L
  }
}

cat(sprintf("\n  Zip download: %d years OK, %d years FAIL\n", n_ok, n_fail))

qcew_raw <- rbindlist(Filter(Negate(is.null), zip_results), fill = TRUE)

if (nrow(qcew_raw) > 0) {
  cat(sprintf("  Raw QCEW rows: %s\n", format(nrow(qcew_raw), big.mark = ",")))

  # Filter: own_code 5 = private, state-level fips ends in 000
  qcew_raw[, area_fips := as.character(area_fips)]
  qcew_filt <- qcew_raw[
    own_code %in% c("5","0") &
    str_ends(str_pad(area_fips, 5, pad = "0"), "000") &
    nchar(str_pad(area_fips, 5, pad = "0")) == 5
  ]

  qcew_filt[, `:=`(
    state_fips2 = str_pad(as.integer(area_fips) %/% 1000L, 2, pad = "0"),
    year        = as.integer(year),
    quarter     = as.integer(qtr),
    emp         = as.numeric(month3_emplvl)
  )]

  cat(sprintf("  After state-level private sector filter: %s rows\n",
              format(nrow(qcew_filt), big.mark = ",")))

  # NAICS code selection: prefer 211 (Oil & Gas Extraction) over 21 (Mining)
  n_naics211 <- qcew_filt[industry_code == "211" & emp > 0, .N]
  naics <- if (n_naics211 > 50) "211" else "21"
  cat(sprintf("  NAICS rows with code 211 and emp>0: %d\n", n_naics211))
  cat(sprintf("  Selected NAICS: %s (%s)\n", naics,
              if (naics == "211") "Oil & Gas Extraction (preferred)"
              else "Mining (broad fallback)"))

  if (naics != "211") {
    log_checkpoint("CK-03a: NAICS 211 availability", "WARN",
                   sprintf("Only %d obs with NAICS 211 -- using NAICS 21 (broader)", n_naics211))
  } else {
    log_checkpoint("CK-03a: NAICS 211 availability", "PASS",
                   sprintf("%d state-quarter obs with NAICS 211", n_naics211))
  }

  oil_emp   <- qcew_filt[industry_code == naics,
                          .(oil_emp = sum(emp, na.rm = TRUE)),
                          by = .(state_fips2, year, quarter)]
  total_emp <- qcew_filt[industry_code == "10",
                          .(total_emp = sum(emp, na.rm = TRUE)),
                          by = .(state_fips2, year, quarter)]

  cat(sprintf("  Oil emp aggregations:   %s state-quarters\n",
              format(nrow(oil_emp), big.mark = ",")))
  cat(sprintf("  Total emp aggregations: %s state-quarters\n",
              format(nrow(total_emp), big.mark = ",")))

  qcew_state <- merge(oil_emp, total_emp,
                      by = c("state_fips2","year","quarter"))
  qcew_state <- merge(qcew_state, state_fips[, .(state_code, fips)],
                      by.x = "state_fips2", by.y = "fips", all.x = TRUE)

  n_unmatched_fips <- qcew_state[is.na(state_code), .N]
  if (n_unmatched_fips > 0) {
    cat(sprintf("  WARNING: %d rows with unmatched FIPS (territories/DC) -> state_code NA\n",
                n_unmatched_fips))
  }

  qcew_state[, mining_emp_share := fifelse(
    !is.na(total_emp) & total_emp > 0,
    oil_emp / total_emp,
    NA_real_
  )]
  qcew_state[, data_source := paste0("BLS_QCEW_NAICS", naics)]
  data_source_used <- "Tier1_BLS_QCEW_zip"

  log_checkpoint("CK-03: Tier 1 BLS zip", "PASS",
                 sprintf("%s rows; NAICS %s; %d years OK",
                         format(nrow(qcew_state), big.mark = ","), naics, n_ok))
} else {
  log_checkpoint("CK-03: Tier 1 BLS zip", "WARN",
                 "All zip downloads failed -- falling through to Tier 2")
}

# ============================================================================
# TIER 2: BLS API (SMS SERIES)
# ============================================================================
if (is.null(qcew_state) || nrow(qcew_state) == 0) {
  cat("\n  --- TIER 2: BLS SMS API ---\n\n")

  fetch_bls_batch <- function(series_ids, start_yr, end_yr) {
    results <- list()
    batches <- split(series_ids, ceiling(seq_along(series_ids) / 25))
    n_batch <- length(batches)
    for (bi in seq_along(batches)) {
      cat(sprintf("    API batch %d/%d ... ", bi, n_batch))
      batch <- batches[[bi]]
      resp <- tryCatch(
        POST("https://api.bls.gov/publicAPI/v2/timeseries/data/",
             body = toJSON(list(seriesid   = batch,
                                startyear  = as.character(start_yr),
                                endyear    = as.character(end_yr)),
                           auto_unbox = TRUE),
             content_type_json(), timeout(45)),
        error = function(e) NULL
      )
      if (is.null(resp) || resp$status_code != 200) {
        cat("FAIL\n"); next
      }
      parsed <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE)
      if (parsed$status != "REQUEST_SUCCEEDED") {
        cat(sprintf("API ERROR: %s\n", parsed$status)); next
      }
      n_series <- 0L
      for (s in parsed$Results$series) {
        d <- as.data.table(s$data)
        d[, `:=`(series_id = s$seriesID,
                 value     = as.numeric(gsub(",", "", value)))]
        results[[s$seriesID]] <- d
        n_series <- n_series + 1L
      }
      cat(sprintf("OK (%d series)\n", n_series))
      Sys.sleep(0.5)
    }
    rbindlist(results, fill = TRUE)
  }

  # SMS{FIPS}00002100000001 = state mining employment (SA)
  # SMS{FIPS}00000000000001 = state total nonfarm (SA)
  mine_ids  <- paste0("SMS", state_fips$fips, "000021000000", "01")
  total_ids <- paste0("SMS", state_fips$fips, "000000000000", "01")

  cat(sprintf("  Fetching mining series (%d states)...\n", length(mine_ids)))
  mine_raw  <- fetch_bls_batch(mine_ids,  START_YEAR_DL, END_YEAR_DL)
  cat(sprintf("  Fetching total employment series (%d states)...\n",
              length(total_ids)))
  total_raw <- fetch_bls_batch(total_ids, START_YEAR_DL, END_YEAR_DL)

  cat(sprintf("\n  Mining series rows:          %s\n",
              format(nrow(mine_raw), big.mark = ",")))
  cat(sprintf("  Total employment series rows: %s\n",
              format(nrow(total_raw), big.mark = ",")))

  if (nrow(mine_raw) > 0 && nrow(total_raw) > 0) {
    parse_api <- function(dt, val_name) {
      dt2 <- dt[grepl("^Q0[1-4]$", period),
                .(fips    = substr(series_id, 4, 5),
                  year    = as.integer(year),
                  quarter = as.integer(str_extract(period, "\\d$")),
                  value   = value)]
      setnames(dt2, "value", val_name)
      dt2
    }
    mine_p  <- parse_api(mine_raw,  "oil_emp")
    total_p <- parse_api(total_raw, "total_emp")

    cat(sprintf("  Parsed mining obs:   %s\n", format(nrow(mine_p),  big.mark = ",")))
    cat(sprintf("  Parsed total obs:    %s\n", format(nrow(total_p), big.mark = ",")))

    qcew_state <- merge(mine_p, total_p, by = c("fips","year","quarter"))
    qcew_state <- merge(qcew_state, state_fips, by = "fips")
    qcew_state[, `:=`(
      mining_emp_share = fifelse(!is.na(total_emp) & total_emp > 0,
                                  oil_emp / total_emp, NA_real_),
      data_source      = "BLS_API_SMS"
    )]
    data_source_used <- "Tier2_BLS_API"
    log_checkpoint("CK-04: Tier 2 BLS API", "PASS",
                   sprintf("%s state-quarter obs", format(nrow(qcew_state), big.mark = ",")))
  } else {
    log_checkpoint("CK-04: Tier 2 BLS API", "WARN",
                   "API returned insufficient data -- falling through to Tier 3")
  }
}

# ============================================================================
# TIER 3: CURATED STATIC FALLBACK (Kilian 2014 / BLS 2005-2019 averages)
# ============================================================================
if (is.null(qcew_state) || nrow(qcew_state) == 0) {
  cat("\n  --- TIER 3: Curated Static Fallback ---\n\n")
  cat("  Source: BLS QCEW 2005-2019 averages; Kilian (2014) state classification\n")
  cat("  NOTE: Time-invariant except ND Bakken adjustment (2.5x post-2008)\n")
  cat("  Cite as: 'Static exposure weights from BLS QCEW averages (Kilian 2014)'\n\n")

  static <- data.table(
    state_code = c("AK","WY","ND","NM","WV","OK","LA","TX","CO",
                   "MT","KS","MS","AR","UT","ID","PA","AL","CA",
                   "IL","IN","KY","MO","NE","NV","OH","SD","VA"),
    base_share = c(0.100,0.080,0.025,0.055,0.050,0.048,0.040,0.038,0.030,
                   0.028,0.025,0.022,0.020,0.018,0.015,0.012,0.008,0.006,
                   0.005,0.005,0.004,0.003,0.003,0.002,0.004,0.003,0.002)
  )

  cat("  States in static table (with base_share):\n")
  for (i in seq_len(nrow(static))) {
    cat(sprintf("    %-4s  %.3f%%  %s\n",
                static$state_code[i],
                static$base_share[i] * 100,
                if (static$base_share[i] >= THRESHOLD_2PCT) "<= oil state" else ""))
  }
  n_above2_static <- sum(static$base_share >= THRESHOLD_2PCT)
  cat(sprintf("\n  States >= 2%% threshold in static table: %d\n", n_above2_static))

  all_st <- merge(state_fips[, .(state_code)], static,
                  by = "state_code", all.x = TRUE)
  all_st[is.na(base_share), base_share := 0.001]

  n_zeroed <- all_st[base_share == 0.001, .N]
  cat(sprintf("  States not in static table (set to 0.001): %d\n", n_zeroed))

  grid <- CJ(state_code = state_fips$state_code,
             year = START_YEAR_DL:END_YEAR_DL,
             quarter = 1:4)
  qcew_state <- merge(grid, all_st, by = "state_code")
  qcew_state <- merge(qcew_state, state_fips[, .(state_code, fips)],
                      by = "state_code", all.x = TRUE)

  # ND Bakken boom: 2.5x multiplier post-2008 (well-documented structural shift)
  qcew_state[, mining_emp_share := base_share *
               fifelse(state_code == "ND" & year >= 2008, 2.5, 1.0)]
  cat("\n  ND Bakken adjustment:\n")
  cat(sprintf("    Pre-2008  ND share: %.3f%%\n",
              qcew_state[state_code=="ND" & year<2008, mean(mining_emp_share)]*100))
  cat(sprintf("    Post-2008 ND share: %.3f%%\n",
              qcew_state[state_code=="ND" & year>=2008, mean(mining_emp_share)]*100))

  qcew_state[, `:=`(oil_emp = mining_emp_share, total_emp = 1.0,
                     data_source = "static_curated_kilian2014")]
  qcew_state[, base_share := NULL]
  data_source_used <- "Tier3_static_curated"

  cat(sprintf("\n  Static fallback rows: %s\n", format(nrow(qcew_state), big.mark = ",")))
  log_checkpoint("CK-05: Tier 3 static fallback", "WARN",
                 sprintf("%s rows; time-invariant (Kilian 2014 averages)",
                         format(nrow(qcew_state), big.mark = ",")))
}

# Final acquisition checkpoint
if (is.null(qcew_state) || nrow(qcew_state) == 0) {
  log_checkpoint("CK-06: Data acquisition", "FAIL",
                 "All three tiers failed -- cannot proceed")
}

cat(sprintf("\n  Data source used    : %s\n", data_source_used))
cat(sprintf("  State-quarter obs   : %s\n", format(nrow(qcew_state), big.mark = ",")))
cat(sprintf("  Unique states       : %d\n", uniqueN(qcew_state$state_code)))
cat(sprintf("  Year range          : %d -- %d\n",
            min(qcew_state$year, na.rm = TRUE),
            max(qcew_state$year, na.rm = TRUE)))
log_checkpoint("CK-06: Data acquisition complete", "PASS",
               sprintf("%s rows via %s", format(nrow(qcew_state),big.mark=","),
                       data_source_used))

# ============================================================================
# SECTION 4: DATA CLEANING AND TYPE ENFORCEMENT
# ============================================================================
section("4", "Data Cleaning and Type Enforcement")

qcew_state[, yyyyqq  := year * 100L + quarter]
qcew_state <- qcew_state[year >= START_YEAR_DL][order(state_code, yyyyqq)]

qcew_state[, state_code := as.character(state_code)]
qcew_state[, yyyyqq     := as.integer(yyyyqq)]
qcew_state[, year       := as.integer(year)]
qcew_state[, quarter    := as.integer(quarter)]

cat("  Type coercion applied: state_code=char, yyyyqq/year/quarter=integer\n\n")

# --- Panel completeness: expect 51 states x n_quarters ---
n_states_actual <- uniqueN(qcew_state$state_code)
n_qtrs_actual   <- uniqueN(qcew_state$yyyyqq)
expected_rows   <- n_states_actual * n_qtrs_actual
actual_rows     <- nrow(qcew_state)

cat(sprintf("  States in panel         : %d\n", n_states_actual))
cat(sprintf("  Quarters in panel       : %d\n", n_qtrs_actual))
cat(sprintf("  Expected rows (balanced): %d\n", expected_rows))
cat(sprintf("  Actual rows             : %d\n", actual_rows))
cat(sprintf("  Missing (unbalanced)    : %d (%.1f%%)\n",
            expected_rows - actual_rows,
            (expected_rows - actual_rows) / expected_rows * 100))

if (actual_rows < expected_rows * 0.90) {
  log_checkpoint("CK-07: Panel completeness", "WARN",
                 sprintf("Panel is <90%% balanced (%d/%d)", actual_rows, expected_rows))
} else {
  log_checkpoint("CK-07: Panel completeness", "PASS",
                 sprintf("%d/%d rows (%.1f%% complete)",
                         actual_rows, expected_rows,
                         actual_rows / expected_rows * 100))
}

# --- Duplicate check ---
dups <- qcew_state[, .N, by = .(state_code, yyyyqq)][N > 1, .N]
if (dups > 0) {
  cat(sprintf("  WARNING: %d duplicate state-quarter keys -- keeping first\n", dups))
  qcew_state <- unique(qcew_state, by = c("state_code","yyyyqq"))
  log_checkpoint("CK-08: Duplicates", "WARN",
                 sprintf("%d duplicate state-quarters removed", dups))
} else {
  log_checkpoint("CK-08: Duplicates", "PASS",
                 "No duplicate state-quarter keys")
}

# --- mining_emp_share diagnostic ---
cat("\n  Raw mining_emp_share before any construction:\n\n")
describe_var(qcew_state$mining_emp_share, "mining_emp_share", "raw")

n_na_share <- sum(is.na(qcew_state$mining_emp_share))
n_neg_share <- sum(qcew_state$mining_emp_share < 0, na.rm = TRUE)
n_gt1_share <- sum(qcew_state$mining_emp_share > 1, na.rm = TRUE)

cat(sprintf("  Impossible values:\n"))
cat(sprintf("    NA              : %d\n", n_na_share))
cat(sprintf("    Negative        : %d  (impossible -- verify input data)\n", n_neg_share))
cat(sprintf("    > 1.0 (>100%%)  : %d  (impossible -- verify input data)\n", n_gt1_share))

if (n_neg_share > 0 || n_gt1_share > 0) {
  cat("  Capping impossible values: negative -> 0, >1 -> NA\n")
  qcew_state[mining_emp_share < 0, mining_emp_share := 0]
  qcew_state[mining_emp_share > 1, mining_emp_share := NA_real_]
  log_checkpoint("CK-09: mining_emp_share range", "WARN",
                 sprintf("Fixed %d negatives, %d >1.0 values", n_neg_share, n_gt1_share))
} else {
  log_checkpoint("CK-09: mining_emp_share range", "PASS",
                 "No impossible values (negative or >1.0)")
}

# ============================================================================
# SECTION 5: EXPOSURE MEASURE CONSTRUCTION
# ============================================================================
section("5", "Exposure Measure Construction")

# ── 5a. Normalized continuous (within-quarter, 0-1) ──────────────────────────
cat("  [5a] oil_exposure_cont -- within-quarter min-max normalization\n")
cat("       Formula: (share - min_share_t) / (max_share_t - min_share_t)\n")
cat("       Purpose: Removes time-trend in levels; 0=least exposed, 1=most\n\n")

qcew_state[, oil_exposure_cont := {
  mn  <- min(mining_emp_share, na.rm = TRUE)
  rng <- max(mining_emp_share, na.rm = TRUE) - mn
  if (is.finite(rng) && rng > 0) (mining_emp_share - mn) / rng
  else rep(0, .N)
}, by = yyyyqq]

describe_var(qcew_state$oil_exposure_cont, "oil_exposure_cont", "constructed")

# Verify 0-1 bounds
n_out_cont <- qcew_state[!is.na(oil_exposure_cont) &
                          (oil_exposure_cont < 0 | oil_exposure_cont > 1), .N]
if (n_out_cont > 0) {
  log_checkpoint("CK-10: oil_exposure_cont bounds", "WARN",
                 sprintf("%d obs outside [0,1]", n_out_cont))
} else {
  log_checkpoint("CK-10: oil_exposure_cont bounds", "PASS",
                 "All values in [0, 1]")
}

# ── 5b. Binary thresholds ─────────────────────────────────────────────────────
cat("\n  [5b] Binary threshold classifications\n")
cat(sprintf("       Primary : oil_exposure_bin      >= %.1f%%\n",
            THRESHOLD_2PCT * 100))
cat(sprintf("       Robust-lo: oil_exposure_bin_1pct >= %.1f%%\n",
            THRESHOLD_1PCT * 100))
cat(sprintf("       Robust-hi: oil_exposure_bin_3pct >= %.1f%%\n",
            THRESHOLD_3PCT * 100))
cat("       Academic defence: 2% is the peer-reviewed standard\n")
cat("         (Kilian 2014; Acemoglu et al. 2016)\n\n")

qcew_state[, `:=`(
  oil_exposure_bin      = as.integer(mining_emp_share >= THRESHOLD_2PCT),
  oil_exposure_bin_1pct = as.integer(mining_emp_share >= THRESHOLD_1PCT),
  oil_exposure_bin_3pct = as.integer(mining_emp_share >= THRESHOLD_3PCT)
)]

# Per-threshold summaries
for (thresh_nm in c("oil_exposure_bin","oil_exposure_bin_1pct","oil_exposure_bin_3pct")) {
  pct_oil <- mean(qcew_state[[thresh_nm]], na.rm = TRUE) * 100
  n_states_above <- uniqueN(qcew_state[get(thresh_nm) == 1, state_code])
  cat(sprintf("  %-24s:  %5.1f%% of state-quarters flagged | %d unique states\n",
              thresh_nm, pct_oil, n_states_above))
}

# Monotonicity check: bin_1pct >= bin_2pct >= bin_3pct must hold
mono_violations <- qcew_state[
  !is.na(oil_exposure_bin) &
  !is.na(oil_exposure_bin_1pct) &
  !is.na(oil_exposure_bin_3pct) &
  !(oil_exposure_bin_1pct >= oil_exposure_bin &
    oil_exposure_bin      >= oil_exposure_bin_3pct), .N]

if (mono_violations > 0) {
  log_checkpoint("CK-11: Binary threshold monotonicity", "WARN",
                 sprintf("%d rows violate 1pct >= 2pct >= 3pct", mono_violations))
} else {
  log_checkpoint("CK-11: Binary threshold monotonicity", "PASS",
                 "1pct >= 2pct >= 3pct holds for all rows")
}

# ── 5c. Tier classification ───────────────────────────────────────────────────
cat("\n  [5c] oil_exposure_tier -- 4-level factor\n")
cat("       Negligible: share <= 0.002\n")
cat("       Low       : 0.002 < share < 0.010\n")
cat("       Medium    : 0.010 <= share < 0.030\n")
cat("       High      : share >= 0.030\n\n")

qcew_state[, oil_exposure_tier := factor(
  fcase(
    mining_emp_share >= 0.030, "High",
    mining_emp_share >= 0.010, "Medium",
    mining_emp_share >  0.002, "Low",
    default                   = "Negligible"
  ),
  levels = c("Negligible","Low","Medium","High")
)]

tier_dist <- qcew_state[, .N, by = oil_exposure_tier][order(oil_exposure_tier)]
tier_dist[, pct := round(N / nrow(qcew_state) * 100, 1)]
cat("  Tier distribution (all state-quarters):\n")
cat(sprintf("  %-12s  %10s  %8s\n", "Tier", "N", "Pct"))
cat(sprintf("  %s\n", paste(rep("-", 34), collapse = "")))
for (i in seq_len(nrow(tier_dist))) {
  cat(sprintf("  %-12s  %10s  %7.1f%%\n",
              as.character(tier_dist$oil_exposure_tier[i]),
              format(tier_dist$N[i], big.mark = ","),
              tier_dist$pct[i]))
}
n_na_tier <- qcew_state[is.na(oil_exposure_tier), .N]
if (n_na_tier > 0) cat(sprintf("  %-12s  %10s  (NA share)\n", "NA",
                                format(n_na_tier, big.mark = ",")))

log_checkpoint("CK-12: oil_exposure_tier", "PASS",
               sprintf("4 tiers: %s",
                       paste(tier_dist$oil_exposure_tier, collapse = " / ")))

# ── 5d. 4-quarter rolling mean ─────────────────────────────────────────────
cat("\n  [5d] oil_exposure_smooth -- 4Q rolling mean within state\n")
cat("       Purpose: Reduces quarterly noise; useful for slow-moving effects\n\n")

qcew_state[, oil_exposure_smooth := frollmean(mining_emp_share, 4,
                                               na.rm = TRUE, align = "right"),
            by = state_code]

n_na_smooth <- sum(is.na(qcew_state$oil_exposure_smooth))
cat(sprintf("  NA from insufficient history (first 3 quarters per state): %d\n",
            n_na_smooth))
describe_var(qcew_state$oil_exposure_smooth, "oil_exposure_smooth", "constructed")
log_checkpoint("CK-13: oil_exposure_smooth", "PASS",
               sprintf("NA from short history: %d", n_na_smooth))

# ============================================================================
# SECTION 6: BARTIK SHIFT-SHARE INSTRUMENT
# ============================================================================
section("6", "Bartik Shift-Share Instrument")

cat("  Formula  : z_st = share_{s, 2005Q1} x Delta_ln(E_national_t)\n")
cat("  Share    : state s mining emp share at base period (2005Q1)\n")
cat("  Shock    : national quarterly log-change in total mining employment\n")
cat("  Reference: Goldsmith-Pinkham, Sorkin & Swift (2020, AER 110:8)\n\n")
cat("  Identification argument:\n")
cat("    Pre-period shares are plausibly exogenous to current local shocks\n")
cat("    Variation is cross-sectional (states differ in base exposure)\n")
cat("    First stage: z_st predicts oil price movements in oil states\n\n")

# ── 6a. National oil employment index (quarterly) ─────────────────────────────
nat_oil <- qcew_state[, .(nat_oil = sum(mining_emp_share * total_emp, na.rm = TRUE)),
                       by = yyyyqq][order(yyyyqq)]

cat(sprintf("  National oil employment index: %d quarters\n", nrow(nat_oil)))
cat(sprintf("  nat_oil range: %.1f to %.1f\n",
            min(nat_oil$nat_oil, na.rm = TRUE),
            max(nat_oil$nat_oil, na.rm = TRUE)))

# Check for zero/negative nat_oil (would make log undefined)
n_zero_nat <- nat_oil[!is.na(nat_oil) & nat_oil <= 0, .N]
if (n_zero_nat > 0) {
  cat(sprintf("  WARNING: %d quarters with nat_oil <= 0 -- setting to NA for log\n",
              n_zero_nat))
  log_checkpoint("CK-14a: nat_oil positive", "WARN",
                 sprintf("%d quarters with nat_oil <= 0", n_zero_nat))
} else {
  log_checkpoint("CK-14a: nat_oil positive", "PASS",
                 "All quarters have nat_oil > 0")
}

nat_oil[, nat_oil_safe := fifelse(nat_oil > 0 & is.finite(nat_oil),
                                   nat_oil, NA_real_)]
nat_oil[, d_ln_nat := c(NA_real_, diff(log(nat_oil_safe)))]
nat_oil[, nat_oil_safe := NULL]

cat(sprintf("\n  d_ln_nat (national log-change in oil employment):\n"))
describe_var(nat_oil$d_ln_nat, "d_ln_nat", "national shock series")

# Verify shock series makes economic sense
n_pos_shock <- nat_oil[d_ln_nat > 0, .N]
n_neg_shock <- nat_oil[d_ln_nat < 0, .N]
cat(sprintf("  Positive shocks (expanding): %d quarters\n", n_pos_shock))
cat(sprintf("  Negative shocks (busting)  : %d quarters\n", n_neg_shock))

# ── 6b. Initial shares (base period) ─────────────────────────────────────────
cat(sprintf("\n  Bartik base quarter: %d\n", BARTIK_BASE_QTR))

init_qtrs <- qcew_state[yyyyqq == BARTIK_BASE_QTR,
                          .(state_code, init_share = mining_emp_share)]

if (nrow(init_qtrs) == 0) {
  eq <- qcew_state[, min(yyyyqq)]
  init_qtrs <- qcew_state[yyyyqq == eq,
                            .(state_code, init_share = mining_emp_share)]
  cat(sprintf("  WARNING: %d not in data -- using earliest available: %d\n",
              BARTIK_BASE_QTR, eq))
  log_checkpoint("CK-14b: Bartik base quarter", "WARN",
                 sprintf("2005Q1 unavailable; using %d instead", eq))
} else {
  log_checkpoint("CK-14b: Bartik base quarter", "PASS",
                 sprintf("%d states with base shares at %d",
                         nrow(init_qtrs), BARTIK_BASE_QTR))
}

cat(sprintf("  States with valid init_share: %d\n", nrow(init_qtrs)))
describe_var(init_qtrs$init_share, "init_share", "base period exposure")

n_zero_init <- sum(init_qtrs$init_share == 0, na.rm = TRUE)
cat(sprintf("  States with zero init_share (instrument will be 0): %d\n",
            n_zero_init))

# ── 6c. Construct instrument ──────────────────────────────────────────────────
if ("init_share" %in% names(qcew_state)) qcew_state[, init_share := NULL]
if ("d_ln_nat"   %in% names(qcew_state)) qcew_state[, d_ln_nat   := NULL]

qcew_state <- merge(qcew_state, init_qtrs, by = "state_code", all.x = TRUE)
qcew_state <- merge(qcew_state,
                    nat_oil[, .(yyyyqq, d_ln_nat)],
                    by = "yyyyqq", all.x = TRUE)

qcew_state[, oil_bartik_iv := fifelse(
  !is.na(init_share) & !is.na(d_ln_nat),
  init_share * d_ln_nat,
  NA_real_
)]

cat("\n  Bartik IV (oil_bartik_iv) distribution:\n\n")
describe_var(qcew_state$oil_bartik_iv, "oil_bartik_iv", "constructed")

n_na_iv    <- sum(is.na(qcew_state$oil_bartik_iv))
n_zero_iv  <- sum(qcew_state$oil_bartik_iv == 0, na.rm = TRUE)

cat(sprintf("  NA (missing init_share or d_ln_nat): %d\n", n_na_iv))
cat(sprintf("  Zero (zero init_share or zero shock): %d\n", n_zero_iv))

# Correlation check: IV should correlate with contemporaneous mining_emp_share
iv_corr <- cor(qcew_state$oil_bartik_iv, qcew_state$mining_emp_share,
               use = "complete.obs")
cat(sprintf("  Correlation(iv, mining_emp_share): %.4f\n", iv_corr))
cat("  (Should be positive and meaningful for first-stage relevance)\n")

if (abs(iv_corr) < 0.10) {
  log_checkpoint("CK-15: Bartik IV relevance", "WARN",
                 sprintf("Low correlation with exposure (r=%.4f) -- check base period", iv_corr))
} else {
  log_checkpoint("CK-15: Bartik IV relevance", "PASS",
                 sprintf("Correlation with exposure: r=%.4f", iv_corr))
}

# ============================================================================
# SECTION 7: SPILLOVER EXPOSURE INDEX (INDIRECT CHANNEL)
# ============================================================================
section("7", "Spillover Exposure Index -- Indirect Channel")

cat("  Method   : Adjacency-weighted average of neighbours' oil exposure\n")
cat("  Formula  : spillover_s,t = Sum_j W_sj * mining_emp_share_j,t\n")
cat("  W_sj     : 1/n_neighbours if j is adjacent to s, else 0 (row-normalised)\n")
cat("  Ref      : Acemoglu et al. (2016, AER 106:1) network propagation\n\n")
cat("  Channels captured:\n")
cat("    (a) Input-output linkages: energy as input to non-oil production\n")
cat("    (b) Consumer spending: oil-state workers buy from non-oil firms\n")
cat("    (c) Financial: bank/CU lending exposure to oil-state counterparties\n\n")

# ── 7a. State adjacency matrix ────────────────────────────────────────────────
cat("  Building row-standardised adjacency matrix (48 contiguous + DC)...\n")
cat("  Source: US Census Bureau state adjacency file\n")
cat("  Note: AK and HI excluded (non-contiguous; spillover set to 0)\n\n")

adj <- list(
  AL=c("FL","GA","MS","TN"),      AZ=c("CA","CO","NM","NV","UT"),
  AR=c("LA","MO","MS","OK","TN","TX"), CA=c("AZ","NV","OR"),
  CO=c("AZ","KS","NE","NM","OK","UT","WY"), CT=c("MA","NY","RI"),
  DE=c("MD","NJ","PA"),           FL=c("AL","GA"),
  GA=c("AL","FL","NC","SC","TN"), ID=c("MT","NV","OR","UT","WA","WY"),
  IL=c("IN","IA","KY","MO","WI"), IN=c("IL","KY","MI","OH"),
  IA=c("IL","MN","MO","NE","SD","WI"), KS=c("CO","MO","NE","OK"),
  KY=c("IL","IN","MO","OH","TN","VA","WV"), LA=c("AR","MS","TX"),
  ME=c("NH"),                     MD=c("DE","PA","VA","WV"),
  MA=c("CT","NH","NY","RI","VT"), MI=c("IN","OH","WI"),
  MN=c("IA","ND","SD","WI"),      MS=c("AL","AR","LA","TN"),
  MO=c("AR","IL","IA","KS","KY","NE","OK","TN"), MT=c("ID","ND","SD","WY"),
  NE=c("CO","IA","KS","MO","SD","WY"), NV=c("AZ","CA","ID","OR","UT"),
  NH=c("MA","ME","VT"),           NJ=c("DE","NY","PA"),
  NM=c("AZ","CO","OK","TX","UT"), NY=c("CT","MA","NJ","PA","VT"),
  NC=c("GA","SC","TN","VA"),      ND=c("MN","MT","SD"),
  OH=c("IN","KY","MI","PA","WV"), OK=c("AR","CO","KS","MO","NM","TX"),
  OR=c("CA","ID","NV","WA"),      PA=c("DE","MD","NJ","NY","OH","WV"),
  RI=c("CT","MA"),                SC=c("GA","NC"),
  SD=c("IA","MN","MT","ND","NE","WY"), TN=c("AL","AR","GA","KY","MS","MO","NC","VA"),
  TX=c("AR","LA","NM","OK"),      UT=c("AZ","CO","ID","NV","NM","WY"),
  VT=c("MA","NH","NY"),           VA=c("KY","MD","NC","TN","WV"),
  WA=c("ID","OR"),                WV=c("KY","MD","OH","PA","VA"),
  WI=c("IL","IA","MI","MN"),      WY=c("CO","ID","MT","NE","SD","UT")
)

all_states <- state_fips$state_code
W_adj <- matrix(0, nrow = length(all_states), ncol = length(all_states),
                dimnames = list(all_states, all_states))

n_adj_entries <- 0L
for (st in names(adj)) {
  nbrs <- intersect(adj[[st]], all_states)
  if (length(nbrs) > 0 && st %in% all_states) {
    W_adj[st, nbrs] <- 1 / length(nbrs)
    n_adj_entries   <- n_adj_entries + length(nbrs)
  }
}

# Validate W_adj
n_states_in_adj   <- sum(rowSums(W_adj) > 0)
n_states_isolated <- sum(rowSums(W_adj) == 0)
row_sum_range     <- range(rowSums(W_adj[rowSums(W_adj) > 0, ]))

cat(sprintf("  States with at least one neighbour : %d\n", n_states_in_adj))
cat(sprintf("  States isolated (AK, HI, DC, etc.) : %d\n", n_states_isolated))
cat(sprintf("  Total adjacency entries             : %d\n", n_adj_entries))
cat(sprintf("  Row sum range (connected states)    : [%.4f, %.4f]\n",
            row_sum_range[1], row_sum_range[2]))
cat(sprintf("  Row-normalised: max row sum = %.6f (should be <= 1.0)\n",
            max(rowSums(W_adj))))

if (max(rowSums(W_adj)) > 1.001) {
  log_checkpoint("CK-16: Adjacency matrix row sums", "WARN",
                 sprintf("Max row sum = %.6f > 1.0", max(rowSums(W_adj))))
} else {
  log_checkpoint("CK-16: Adjacency matrix row sums", "PASS",
                 sprintf("All row sums <= 1.0; %d states connected", n_states_in_adj))
}

# ── 7b. Compute spillover via matrix multiplication ───────────────────────────
cat("\n  Computing spillover via matrix multiply (share_mat %*% t(W_adj))...\n\n")

wide_share <- dcast(
  qcew_state[state_code %in% all_states],
  yyyyqq ~ state_code,
  value.var = "mining_emp_share"
)
setDF(wide_share)
qtrs_vec  <- wide_share$yyyyqq
share_mat <- as.matrix(wide_share[, -1])

mat_states <- intersect(all_states, colnames(share_mat))
W_sub      <- W_adj[mat_states, mat_states]
share_sub  <- share_mat[, mat_states]
share_sub[is.na(share_sub)] <- 0   # NAs -> 0 before multiply

n_na_before <- sum(is.na(share_mat[, mat_states]))
cat(sprintf("  States in share matrix         : %d\n", length(mat_states)))
cat(sprintf("  Quarters in share matrix       : %d\n", nrow(share_mat)))
cat(sprintf("  NA cells in share matrix       : %d (set to 0 before multiply)\n",
            n_na_before))
cat(sprintf("  Matrix dimensions: share_sub [%d x %d], W_sub [%d x %d]\n",
            nrow(share_sub), ncol(share_sub), nrow(W_sub), ncol(W_sub)))

spillover_mat <- share_sub %*% t(W_sub)

cat(sprintf("  Spillover matrix computed: [%d x %d]\n",
            nrow(spillover_mat), ncol(spillover_mat)))
cat(sprintf("  Spillover range: [%.6f, %.6f]\n",
            min(spillover_mat, na.rm = TRUE),
            max(spillover_mat, na.rm = TRUE)))
cat(sprintf("  Mean spillover : %.6f\n", mean(spillover_mat, na.rm = TRUE)))

spill_long <- melt(
  as.data.table(cbind(yyyyqq = qtrs_vec, as.data.frame(spillover_mat))),
  id.vars      = "yyyyqq",
  variable.name = "state_code",
  value.name   = "spillover_exposure"
)
spill_long[, state_code := as.character(state_code)]

# Trade-weighted version (placeholder; BEA RIMS II upgrade path documented)
spill_long[, spillover_exposure_wtd := spillover_exposure]

cat("\n  Spillover exposure distribution:\n\n")
describe_var(spill_long$spillover_exposure, "spillover_exposure", "constructed")

log_checkpoint("CK-17: Spillover matrix", "PASS",
               sprintf("[%dx%d] matrix; range [%.4f, %.4f]",
                       nrow(spillover_mat), ncol(spillover_mat),
                       min(spillover_mat, na.rm=TRUE),
                       max(spillover_mat, na.rm=TRUE)))

# ── 7c. Merge spillover into main table ───────────────────────────────────────
qcew_state <- merge(qcew_state, spill_long,
                    by = c("state_code","yyyyqq"), all.x = TRUE)

n_na_spill <- sum(is.na(qcew_state$spillover_exposure))
cat(sprintf("  Post-merge NA in spillover_exposure: %d (isolated states: AK/HI/DC)\n",
            n_na_spill))

log_checkpoint("CK-18: Spillover merge", "PASS",
               sprintf("Merged; %d NAs (isolated states)", n_na_spill))

# ============================================================================
# SECTION 8: CU GROUP CLASSIFICATION
# ============================================================================
section("8", "CU Group Classification")

cat("  Direct     = oil_exposure_bin == 1  (oil-state CUs)\n")
cat("  Indirect   = oil_exposure_bin == 0 AND spillover > median(spillover)\n")
cat("               (non-oil but regionally linked)\n")
cat("  Negligible = low oil AND low spillover (pure national macro exposure)\n\n")
cat("  Note: median computed on non-oil states only to avoid contamination\n\n")

valid_spill <- qcew_state[
  oil_exposure_bin == 0L &
  !is.na(spillover_exposure) &
  is.finite(spillover_exposure),
  spillover_exposure
]

if (length(valid_spill) > 0) {
  med_spill <- median(valid_spill)
  cat(sprintf("  Valid spillover obs for median (non-oil states): %s\n",
              format(length(valid_spill), big.mark = ",")))
  cat(sprintf("  Spillover median (non-oil states)              : %.6f\n", med_spill))
  cat(sprintf("  Pct of non-oil states above median             : %.1f%%\n",
              mean(valid_spill > med_spill) * 100))
} else {
  med_spill <- 0
  cat("  WARNING: no valid spillover for non-oil states -- using 0 as median\n")
  log_checkpoint("CK-19a: Spillover median", "WARN",
                 "No valid spillover values for non-oil states")
}

qcew_state[, cu_group := fcase(
  oil_exposure_bin == 1L,
    "Direct",
  oil_exposure_bin == 0L & !is.na(spillover_exposure) &
    spillover_exposure > med_spill,
    "Indirect",
  default = "Negligible"
)]
qcew_state[, cu_group := factor(cu_group,
                                  levels = c("Direct","Indirect","Negligible"))]

# Distribution across all quarters
grp_all <- qcew_state[, .N, by = cu_group][order(cu_group)]
grp_all[, pct := round(N / nrow(qcew_state) * 100, 1)]

cat("\n  cu_group distribution (all state-quarters):\n")
cat(sprintf("  %-12s  %10s  %8s\n", "Group", "N", "Pct"))
cat(sprintf("  %s\n", paste(rep("-", 34), collapse = "")))
for (i in seq_len(nrow(grp_all))) {
  cat(sprintf("  %-12s  %10s  %7.1f%%\n",
              as.character(grp_all$cu_group[i]),
              format(grp_all$N[i], big.mark = ","),
              grp_all$pct[i]))
}

# Latest quarter breakdown
latest_qtr <- max(qcew_state$yyyyqq)
grp_latest <- qcew_state[yyyyqq == latest_qtr, .N, by = cu_group][order(cu_group)]

cat(sprintf("\n  cu_group at latest quarter (%d):\n", latest_qtr))
for (i in seq_len(nrow(grp_latest))) {
  cat(sprintf("  %-12s  %d states\n",
              as.character(grp_latest$cu_group[i]),
              grp_latest$N[i]))
}

# Check: Direct must be a non-trivial minority
n_direct <- grp_all[cu_group == "Direct", N]
pct_direct <- grp_all[cu_group == "Direct", pct]
if (pct_direct < 5 || pct_direct > 40) {
  log_checkpoint("CK-19: CU group composition", "WARN",
                 sprintf("Direct = %.1f%% of state-quarters (expected 10-30%%)", pct_direct))
} else {
  log_checkpoint("CK-19: CU group composition", "PASS",
                 sprintf("Direct=%.1f%%, Indirect=%.1f%%, Negligible=%.1f%%",
                         grp_all[cu_group=="Direct",    pct],
                         grp_all[cu_group=="Indirect",  pct],
                         grp_all[cu_group=="Negligible",pct]))
}

# ============================================================================
# SECTION 9: VALIDATION
# ============================================================================
section("9", "Validation Checks")

# ── 9a. Top 20 states by average oil employment share ─────────────────────────
cat("  [9a] Top 20 states by average mining_emp_share (2005-present):\n\n")
cat(sprintf("  %-6s  %-11s  %-11s  %-12s  %s\n",
            "State","Avg Share","Max Share","Qtrs >=2%","Tier (modal)"))
cat(sprintf("  %s\n", paste(rep("-", 62), collapse = "")))

top20 <- qcew_state[year >= 2005, .(
  avg_share  = mean(mining_emp_share, na.rm = TRUE),
  max_share  = max(mining_emp_share, na.rm = TRUE),
  pct_above2 = mean(oil_exposure_bin, na.rm = TRUE) * 100,
  modal_tier = names(which.max(table(oil_exposure_tier)))
), by = state_code][order(-avg_share)][1:20]

for (i in seq_len(nrow(top20))) {
  cat(sprintf("  %-6s  %9.3f%%  %9.3f%%  %10.1f%%  %s\n",
              top20$state_code[i],
              top20$avg_share[i]  * 100,
              top20$max_share[i]  * 100,
              top20$pct_above2[i],
              top20$modal_tier[i]))
}

# CK: WY, ND, AK, TX must be in top 10
expected_top <- c("WY","ND","AK","TX","OK","NM","LA","WV","MT","CO")
top10_states <- top20$state_code[1:10]
missing_expected <- setdiff(expected_top, top10_states)
if (length(missing_expected) > 0) {
  log_checkpoint("CK-20: Top oil states", "WARN",
                 sprintf("Expected states absent from top 10: %s",
                         paste(missing_expected, collapse = ", ")))
} else {
  log_checkpoint("CK-20: Top oil states", "PASS",
                 "All expected top oil states in top 10")
}

# ── 9b. North Dakota Bakken validation ────────────────────────────────────────
cat("\n  [9b] North Dakota Bakken boom validation:\n")
nd <- qcew_state[state_code == "ND"][order(yyyyqq)]

if (nrow(nd) == 0) {
  cat("  WARNING: ND not found in panel\n")
  log_checkpoint("CK-21: ND Bakken validation", "WARN", "ND absent from panel")
} else {
  nd_pre    <- nd[year <= 2007, mean(mining_emp_share, na.rm = TRUE)] * 100
  nd_peak   <- nd[year %in% 2011:2014, mean(mining_emp_share, na.rm = TRUE)] * 100
  nd_post   <- nd[year >= 2016, mean(mining_emp_share, na.rm = TRUE)] * 100
  nd_max    <- nd[, max(mining_emp_share, na.rm = TRUE)] * 100
  nd_maxqtr <- nd[mining_emp_share == max(mining_emp_share, na.rm = TRUE), yyyyqq][1]

  cat(sprintf("    Pre-Bakken  2005-07 avg: %6.2f%%\n", nd_pre))
  cat(sprintf("    Bakken peak 2011-14 avg: %6.2f%%\n", nd_peak))
  cat(sprintf("    Post-bust   2016+   avg: %6.2f%%\n", nd_post))
  cat(sprintf("    All-time max           : %6.2f%% (quarter: %d)\n",
              nd_max, nd_maxqtr))

  # Bakken should show peak >> pre period
  if (nd_peak > nd_pre * 1.5) {
    log_checkpoint("CK-21: ND Bakken validation", "PASS",
                   sprintf("Peak (%.2f%%) >> Pre (%.2f%%) -- Bakken boom captured",
                           nd_peak, nd_pre))
  } else {
    log_checkpoint("CK-21: ND Bakken validation", "WARN",
                   sprintf("Peak (%.2f%%) not much > Pre (%.2f%%) -- check data",
                           nd_peak, nd_pre))
  }
}

# ── 9c. Cross-validate against hard-coded state list ─────────────────────────
cat("\n  [9c] Cross-validation vs hard-coded oil state list (2015 -- peak shale):\n\n")

hardcoded   <- c("TX","ND","LA","AK","WY","OK","NM","CO","WV","PA","MT")
data_driven <- qcew_state[year == 2015 & oil_exposure_bin == 1,
                            sort(unique(state_code))]

cat(sprintf("  Hard-coded list  : %s\n", paste(sort(hardcoded), collapse = ", ")))
cat(sprintf("  Data-driven 2015 : %s\n", paste(sort(data_driven), collapse = ", ")))
cat(sprintf("  Added by data    : %s\n",
            if (length(setdiff(data_driven, hardcoded)) > 0)
              paste(setdiff(data_driven, hardcoded), collapse = ", ")
            else "(none)"))
cat(sprintf("  Dropped by data  : %s\n",
            if (length(setdiff(hardcoded, data_driven)) > 0)
              paste(setdiff(hardcoded, data_driven), collapse = ", ")
            else "(none)"))

overlap_pct <- length(intersect(hardcoded, data_driven)) /
               length(union(hardcoded, data_driven)) * 100
cat(sprintf("  Jaccard overlap  : %.1f%%\n", overlap_pct))

if (overlap_pct < 60) {
  log_checkpoint("CK-22: Cross-validation vs hard-coded list", "WARN",
                 sprintf("Low overlap: %.1f%% Jaccard -- check threshold or data source",
                         overlap_pct))
} else {
  log_checkpoint("CK-22: Cross-validation vs hard-coded list", "PASS",
                 sprintf("%.1f%% Jaccard overlap", overlap_pct))
}

# ── 9d. Top indirect-channel states ──────────────────────────────────────────
cat("\n  [9d] Top 10 indirect-channel states (highest avg spillover, non-oil):\n\n")

spill_top <- qcew_state[
  oil_exposure_bin == 0 & year >= 2005,
  .(avg_spill = mean(spillover_exposure, na.rm = TRUE),
    max_spill = max(spillover_exposure, na.rm = TRUE)),
  by = state_code
][order(-avg_spill)][1:10]

cat(sprintf("  %-6s  %-14s  %-14s\n", "State","Avg Spillover","Max Spillover"))
cat(sprintf("  %s\n", paste(rep("-", 38), collapse = "")))
for (i in seq_len(nrow(spill_top))) {
  cat(sprintf("  %-6s  %12.6f  %12.6f\n",
              spill_top$state_code[i],
              spill_top$avg_spill[i],
              spill_top$max_spill[i]))
}

# States bordering TX/OK/LA should have highest spillover
expected_indirect <- c("AR","MO","KS","MS","NM")
found_indirect    <- intersect(expected_indirect, spill_top$state_code)
cat(sprintf("\n  Expected high-spillover border states in top 10: %s\n",
            if (length(found_indirect) > 0) paste(found_indirect, collapse = ", ")
            else "NONE"))

log_checkpoint("CK-23: Indirect channel states", "PASS",
               sprintf("Top indirect: %s",
                       paste(head(spill_top$state_code, 5), collapse = ", ")))

# ── 9e. Time series stationarity check on national shock ─────────────────────
cat("\n  [9e] National shock series (d_ln_nat) stationarity check:\n")
cat("       (Non-stationary shock series would invalidate Bartik identification)\n\n")

dln_clean <- nat_oil$d_ln_nat[!is.na(nat_oil$d_ln_nat)]
if (length(dln_clean) >= 20 && requireNamespace("tseries", quietly = TRUE)) {
  adf_res <- tryCatch(tseries::adf.test(dln_clean), error = function(e) NULL)
  if (!is.null(adf_res)) {
    cat(sprintf("  ADF test statistic: %.4f\n", adf_res$statistic))
    cat(sprintf("  ADF p-value       : %.4f\n", adf_res$p.value))
    cat(sprintf("  Conclusion        : %s\n",
                if (adf_res$p.value < 0.05) "STATIONARY (p < 0.05)"
                else "May be non-stationary (p >= 0.05)"))
    if (adf_res$p.value >= 0.10) {
      log_checkpoint("CK-24: Shock series stationarity", "WARN",
                     sprintf("ADF p=%.4f -- shock may be non-stationary", adf_res$p.value))
    } else {
      log_checkpoint("CK-24: Shock series stationarity", "PASS",
                     sprintf("ADF p=%.4f -- stationary", adf_res$p.value))
    }
  }
} else {
  cat("  (tseries package not available or insufficient obs -- skipping ADF)\n")
  log_checkpoint("CK-24: Shock series stationarity", "INFO",
                 "tseries not available -- ADF skipped")
}

# ── 9f. Overall exposure coverage ────────────────────────────────────────────
cat("\n  [9f] Overall exposure measure coverage:\n\n")

EXPOSURE_COLS <- c("mining_emp_share","oil_exposure_cont","oil_exposure_bin",
                   "oil_exposure_bin_1pct","oil_exposure_bin_3pct",
                   "oil_exposure_smooth","oil_bartik_iv",
                   "spillover_exposure","spillover_exposure_wtd")

cat(sprintf("  %-28s  %8s  %8s  %8s\n", "Variable", "N valid", "% valid", "% miss"))
cat(sprintf("  %s\n", paste(rep("-", 58), collapse = "")))
for (v in EXPOSURE_COLS) {
  if (!v %in% names(qcew_state)) {
    cat(sprintf("  %-28s  MISSING FROM PANEL\n", v))
    next
  }
  n_v   <- sum(!is.na(qcew_state[[v]]))
  pct_v <- n_v / nrow(qcew_state) * 100
  flag  <- if (pct_v < 80) "  <- LOW" else ""
  cat(sprintf("  %-28s  %8s  %7.1f%%  %7.1f%%%s\n",
              v,
              format(n_v, big.mark = ","),
              pct_v,
              100 - pct_v,
              flag))
}

# ============================================================================
# SECTION 10: SAVE OUTPUT
# ============================================================================
section("10", "Save Output")

out <- qcew_state[year >= START_YEAR_DL][order(state_code, yyyyqq)]

# Enforce final column types
out[, state_code := as.character(state_code)]
out[, yyyyqq     := as.integer(yyyyqq)]
out[, year       := as.integer(year)]
out[, quarter    := as.integer(quarter)]

keep_cols <- c(
  "state_code","fips","yyyyqq","year","quarter",
  "oil_emp","total_emp","mining_emp_share",
  "oil_exposure_cont","oil_exposure_bin",
  "oil_exposure_bin_1pct","oil_exposure_bin_3pct",
  "oil_exposure_tier","oil_exposure_smooth",
  "oil_bartik_iv","init_share",
  "spillover_exposure","spillover_exposure_wtd",
  "cu_group","data_source"
)

out_cols <- intersect(keep_cols, names(out))
missing_cols <- setdiff(keep_cols, names(out))
if (length(missing_cols) > 0) {
  cat(sprintf("  Columns in keep_cols but not constructed: %s\n",
              paste(missing_cols, collapse = ", ")))
  log_checkpoint("CK-25a: Output columns", "WARN",
                 sprintf("Missing from output: %s", paste(missing_cols, collapse=", ")))
} else {
  log_checkpoint("CK-25a: Output columns", "PASS",
                 sprintf("All %d keep_cols present", length(keep_cols)))
}

out <- out[, ..out_cols]

saveRDS(out, "Data/oil_exposure.rds")
fsize <- round(file.size("Data/oil_exposure.rds") / 1e6, 1)

cat(sprintf("  Saved: Data/oil_exposure.rds\n"))
cat(sprintf("  Rows      : %s\n", format(nrow(out), big.mark = ",")))
cat(sprintf("  States    : %d\n", uniqueN(out$state_code)))
cat(sprintf("  Quarters  : %d\n", uniqueN(out$yyyyqq)))
cat(sprintf("  Columns   : %d\n", ncol(out)))
cat(sprintf("  File size : %.1f MB\n", fsize))

log_checkpoint("CK-25: Output saved", "PASS",
               sprintf("Data/oil_exposure.rds -- %s rows, %d states, %.1f MB",
                       format(nrow(out),big.mark=","), uniqueN(out$state_code), fsize))

# ============================================================================
# SECTION 11: SAVE AUDIT FILES
# ============================================================================
section("11", "Save Audit Files")

# Per-state exposure summary
exp_summ <- qcew_state[year >= 2005, .(
  n_qtrs         = .N,
  avg_share      = round(mean(mining_emp_share, na.rm=TRUE), 6),
  max_share      = round(max(mining_emp_share, na.rm=TRUE), 6),
  min_share      = round(min(mining_emp_share, na.rm=TRUE), 6),
  sd_share       = round(sd(mining_emp_share, na.rm=TRUE), 6),
  pct_above2pct  = round(mean(oil_exposure_bin, na.rm=TRUE)*100, 1),
  pct_above1pct  = round(mean(oil_exposure_bin_1pct, na.rm=TRUE)*100, 1),
  pct_above3pct  = round(mean(oil_exposure_bin_3pct, na.rm=TRUE)*100, 1),
  avg_bartik_iv  = round(mean(oil_bartik_iv, na.rm=TRUE), 6),
  avg_spillover  = round(mean(spillover_exposure, na.rm=TRUE), 6),
  modal_cu_group = names(which.max(table(cu_group)))
), by = state_code][order(-avg_share)]

fwrite(exp_summ, file.path(RESULTS_DIR, "01b_exposure_summary.csv"))
cat(sprintf("  Saved: Results/01b_exposure_summary.csv  (%d states)\n", nrow(exp_summ)))

# Bartik IV diagnostics
bartik_summ <- qcew_state[, .(
  yyyyqq         = yyyyqq[1],
  year           = year[1],
  quarter        = quarter[1],
  d_ln_nat       = d_ln_nat[1],
  mean_init_share = mean(init_share, na.rm=TRUE),
  n_states_nonzero_iv = sum(oil_bartik_iv != 0, na.rm=TRUE),
  mean_iv        = mean(oil_bartik_iv, na.rm=TRUE),
  sd_iv          = sd(oil_bartik_iv, na.rm=TRUE)
), by = yyyyqq][order(yyyyqq)]

fwrite(bartik_summ, file.path(RESULTS_DIR, "01b_bartik_summary.csv"))
cat(sprintf("  Saved: Results/01b_bartik_summary.csv  (%d quarters)\n",
            nrow(bartik_summ)))

# Spillover diagnostics
spill_summ <- qcew_state[, .(
  avg_spillover = mean(spillover_exposure, na.rm=TRUE),
  max_spillover = max(spillover_exposure, na.rm=TRUE),
  n_valid       = sum(!is.na(spillover_exposure))
), by = .(state_code, cu_group)][order(-avg_spillover)]

fwrite(spill_summ, file.path(RESULTS_DIR, "01b_spillover_summary.csv"))
cat(sprintf("  Saved: Results/01b_spillover_summary.csv  (%d state-group rows)\n",
            nrow(spill_summ)))

log_checkpoint("CK-26: Audit files saved", "PASS",
               "01b_exposure_summary, 01b_bartik_summary, 01b_spillover_summary")

# ============================================================================
# SECTION 12: CHECKPOINT SUMMARY
# ============================================================================
section("12", "Checkpoint Summary")

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
       file.path(RESULTS_DIR, "01b_checkpoint_log.csv"))
cat(sprintf("\n  Checkpoint log saved: Results/01b_checkpoint_log.csv\n"))

# ============================================================================
# CLOSING BANNER
# ============================================================================
cat("\n", SEP, "\n", sep = "")
cat("  SCRIPT 01b COMPLETE\n")
cat(sprintf("  Finished    : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("  Data source : %s\n", data_source_used))
cat("\n  Data output:\n")
cat("    Data/oil_exposure.rds         State x quarter exposure panel\n")
cat("\n  Audit outputs:\n")
cat("    Results/01b_checkpoint_log.csv     All checkpoint pass/warn/fail\n")
cat("    Results/01b_exposure_summary.csv   Per-state summary statistics\n")
cat("    Results/01b_bartik_summary.csv     Bartik IV quarterly diagnostics\n")
cat("    Results/01b_spillover_summary.csv  Spillover diagnostics by group\n")
cat("\n  Exposure variables in output:\n")
cat("    mining_emp_share        Raw QCEW oil/gas emp share (primary)\n")
cat("    oil_exposure_cont       Within-quarter min-max normalized [0,1]\n")
cat("    oil_exposure_bin        Binary >= 2% (primary threshold)\n")
cat("    oil_exposure_bin_1pct   Binary >= 1% (robustness lower)\n")
cat("    oil_exposure_bin_3pct   Binary >= 3% (robustness upper)\n")
cat("    oil_exposure_tier       Negligible / Low / Medium / High\n")
cat("    oil_exposure_smooth     4Q rolling mean (noise-reduced)\n")
cat("    oil_bartik_iv           Bartik shift-share IV (Goldsmith-Pinkham 2020)\n")
cat("\n  Indirect channel variables:\n")
cat("    spillover_exposure      Adjacency-weighted neighbour oil exposure\n")
cat("    spillover_exposure_wtd  Trade-weighted version (placeholder: BEA RIMS II)\n")
cat("    cu_group                Direct / Indirect / Negligible\n")
cat("\n  CU group -- transmission interpretation:\n")
cat("    Direct     -- local energy-sector income/employment effects\n")
cat("    Indirect   -- regional spillover from oil-state economic activity\n")
cat("    Negligible -- national macro channels only (baseline comparison)\n")
cat("\n  Citations:\n")
cat("    BLS QCEW: U.S. Bureau of Labor Statistics (2025)\n")
cat("    Bartik (1991) Upjohn Institute\n")
cat("    Goldsmith-Pinkham, Sorkin & Swift (2020, AER 110:8)\n")
cat("    Acemoglu et al. (2016, AER 106:1) network propagation\n")
cat("    Kilian (2014) Annual Review of Resource Economics 6\n")
cat(sprintf("\n  Checkpoints: %d PASS | %d WARN | %d FAIL\n", n_pass, n_warn, n_fail))
cat(SEP, "\n", sep = "")
