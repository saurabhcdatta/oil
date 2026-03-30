# =============================================================================
# SCRIPT 02 — Exploratory Data Analysis + Policy Charts
# Oil Price Shock × Credit Union Financial Performance — v2
# Author: Saurabh C. Datta, Ph.D. | NCUA Office of Chief Economist
# Version: v2.1 | Built: 2026-03-30
# Working dir: S:/Projects/Oil_Price_Shock_2026/
# =============================================================================
# Inputs  : Data/panel_base.rds        (from 01_data_prep.R)
#           Data/panel_severe.rds       (from 01_data_prep.R)
#           Data/macro_base.rds
#           Data/macro_severe.rds
#
# CORE CHARTS (01-10)
#   01  Oil price history + cycle annotation        (PBRENT 2005-2025)
#   02  PBRENT vs aggregate CU outcomes             (time-series overlay)
#   03  Cross-correlation: PBRENT leads outcomes 1-8 quarters
#   04  Direct vs indirect: oil-state vs non-oil CUs
#   05  Deposit channel: cert_share, loan_to_share, CoF vs FOMC regime
#   06  Asset tier response: all 8 tiers (assets_cat2)
#   07  Structural break: pre vs post 2015Q1 scatter
#   08  Spillover: non-oil CU response by spillover tercile
#   09  Episode heatmap: oil cycle x CU outcome
#   10  Missingness overview
#
# MEMBER GROWTH CHARTS (NEW-A to NEW-D)
#   NEW-A  Membership growth time series vs oil + direct/indirect split
#   NEW-B  Membership growth structural break: pre/post 2015 slope test
#   NEW-C  Membership growth by asset tier (time series + boxplot)
#   NEW-D  Membership growth lead/lag cross-correlation vs oil & deposits
#
# NEW POLICY CHARTS (16-22)
#   16  PLL rate deep-dive: direct/indirect + tier + structural break
#   17  Net worth ratio (pcanetworth) under oil stress cycles
#   18  NCUA PCA traffic-light dashboard: well-capitalised share over time
#   19  Iran war / $125 oil scenario: Moody's implied CU stress
#   20  Bartik IV first-stage diagnostic: exposure vs oil shock
#   21  Deposit migration: certificate vs demand shift under rate regimes
#   22  8-tier assets_cat2 response dashboard: all tiers x key outcomes
#
# SCENARIO CHARTS (11-15, conditional on panel_severe)
#   11  CCAR scenario comparison: PBRENT baseline vs severely adverse
#   12  CU outcomes under CCAR scenarios
#   13  Rate & yield curve environment
#   14  Scenario divergence
#   15  Implied CU stress under severely adverse
# =============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

# ============================================================================
# LOGGING
# ============================================================================
SEP <- paste(rep("=", 72), collapse = "")
msg <- function(...) cat(sprintf(...), "\n")
hdr <- function(s) {
  cat("\n", SEP, "\n", sep = "")
  cat(sprintf("  %s\n", s))
  cat(paste(rep("-", 72), collapse = ""), "\n", sep = "")
}

CHART_LOG <- data.frame(chart=character(), status=character(),
                         file=character(), stringsAsFactors=FALSE)

log_chart <- function(num, file, status = "SAVED") {
  cat(sprintf("  [%s] %s -> %s\n", status, num, file))
  CHART_LOG <<- rbind(CHART_LOG, data.frame(chart=num, status=status,
                                              file=file, stringsAsFactors=FALSE))
}

cat(SEP, "\n")
cat("  SCRIPT 02 -- EDA + Policy Charts\n")
cat("  Oil Price Shock x Credit Union Financial Performance -- v2\n")
cat(sprintf("  Started: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(SEP, "\n\n")

# ============================================================================
# SECTION 1: CONFIGURATION
# ============================================================================
hdr("SECTION 1: Configuration")

dir.create("Figures", showWarnings = FALSE)

# ── Publication theme ─────────────────────────────────────────────────────────
theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
  theme(
    text             = element_text(family = "sans", colour = "#1a1a1a"),
    plot.title       = element_text(size = base_size + 2, face = "bold",
                                    margin = margin(b = 4)),
    plot.subtitle    = element_text(size = base_size - 0.5, colour = "#555555",
                                    margin = margin(b = 8)),
    plot.caption     = element_text(size = base_size - 2, colour = "#888888",
                                    hjust = 0),
    axis.title       = element_text(size = base_size - 0.5, colour = "#333333"),
    axis.text        = element_text(size = base_size - 1.5, colour = "#444444"),
    panel.grid.major = element_line(colour = "#e8e8e8", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "#cccccc", fill = NA,
                                    linewidth = 0.4),
    strip.text       = element_text(size = base_size - 1, face = "bold",
                                    colour = "#333333"),
    strip.background = element_rect(fill = "#f5f5f5", colour = "#cccccc"),
    legend.position  = "bottom",
    legend.text      = element_text(size = base_size - 1.5),
    legend.title     = element_text(size = base_size - 1, face = "bold"),
    plot.margin      = margin(10, 12, 8, 10)
  )
}

# Colour palette
COL_OIL    <- "#b5470a"   # Brent oil -- burnt orange
COL_DIRECT <- "#1a3a5c"   # oil-state CUs -- navy
COL_INDIR  <- "#2d7a4a"   # non-oil CUs -- forest green
COL_SPILL  <- "#7a3080"   # spillover -- purple
COL_NEG    <- "#c0392b"   # stress / negative -- red
COL_POS    <- "#27ae60"   # positive -- green
COL_WARN   <- "#e67e22"   # warning amber

# 8-tier colour palette matching assets_cat2 from OCE_combined
TIER_COLS <- c(
  "T1_under10M"  = "#4a90d9",
  "T2_10to50M"   = "#e67e22",
  "T3_50to100M"  = "#8e44ad",
  "T4_100to500M" = "#16a085",
  "T5_500Mto1B"  = "#c0392b",
  "T6_1Bto5B"    = "#2ecc71",
  "T7_5Bto10B"   = "#f39c12",
  "T8_over10B"   = "#1a252f"
)

TIER_LABELS <- c(
  T1_under10M  = "< $10M",
  T2_10to50M   = "$10-50M",
  T3_50to100M  = "$50-100M",
  T4_100to500M = "$100-500M",
  T5_500Mto1B  = "$500M-$1B",
  T6_1Bto5B    = "$1-5B",
  T7_5Bto10B   = "$5-10B",
  T8_over10B   = "> $10B"
)

# Key oil price episode shading bands
EPISODES <- data.frame(
  label = c("GFC", "Shale\nBust", "COVID\nCrash", "Post-COVID\nSurge",
            "Iran\nWar"),
  xmin  = as.Date(c("2008-07-01","2014-07-01","2020-01-01","2021-01-01",
                     "2025-10-01")),
  xmax  = as.Date(c("2009-06-30","2016-06-30","2020-06-30","2022-06-30",
                     "2026-03-31")),
  fill  = c("#fde8e8","#e8f0fd","#fde8e8","#e8fde8","#fff3cd")
)

# Scenario colours
SCEN_COLS <- c("Historical"     = "#333333",
               "Baseline"       = COL_INDIR,
               "Severely Adverse" = COL_NEG)
SCEN_LT   <- c("Historical"     = "solid",
               "Baseline"       = "dashed",
               "Severely Adverse" = "solid")

# ── Save helper ───────────────────────────────────────────────────────────────
save_plot <- function(p, filename, w = 10, h = 6.5, dpi = 300) {
  path <- file.path("Figures", filename)
  tryCatch({
    ggsave(path, plot = p, width = w, height = h, dpi = dpi, bg = "white")
    log_chart(tools::file_path_sans_ext(filename), filename)
  }, error = function(e) {
    cat(sprintf("  [ERROR] Failed to save %s: %s\n", filename, conditionMessage(e)))
    log_chart(tools::file_path_sans_ext(filename), filename, "ERROR")
  })
}

# ── Episode rectangles helper ─────────────────────────────────────────────────
ep_rects <- function(episodes = EPISODES) {
  rects <- mapply(function(xmn, xmx, fl) {
    annotate("rect", xmin = xmn, xmax = xmx,
             ymin = -Inf, ymax = Inf, fill = fl, alpha = 0.35)
  }, episodes$xmin, episodes$xmax, episodes$fill, SIMPLIFY = FALSE)

  texts <- mapply(function(xmn, xmx, lb) {
    annotate("text", x = xmn + (xmx - xmn) / 2,
             y = Inf, label = lb, vjust = 1.3, size = 2.5,
             colour = "#888888", fontface = "italic")
  }, episodes$xmin, episodes$xmax, episodes$label, SIMPLIFY = FALSE)

  c(rects, texts)
}

# ── Null-coalescing operator ──────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

msg("  Figures dir: Figures/")
msg("  Tier scheme: 8-tier assets_cat2 (T1-T8)")
msg("  Episodes: %d shading bands (GFC to Iran War)", nrow(EPISODES))

# ============================================================================
# SECTION 2: LOAD DATA
# ============================================================================
hdr("SECTION 2: Load Data")

panel  <- setDT(readRDS("Data/panel_base.rds"))
macro  <- setDT(readRDS("Data/macro_base.rds"))

msg("  panel_base  : %s rows x %d cols | %s CUs | %d quarters",
    format(nrow(panel), big.mark=","), ncol(panel),
    format(uniqueN(panel$join_number), big.mark=","),
    uniqueN(panel$yyyyqq))

# Calendar date
Q_MONTH <- c("1"=1L,"2"=4L,"3"=7L,"4"=10L)
add_date <- function(dt) {
  if (!"cal_date" %in% names(dt))
    dt[, cal_date := as.Date(paste(year, Q_MONTH[as.character(quarter)],
                                   "01", sep="-"))]
}
add_date(panel); add_date(macro)

# Macro spine
mac_spine <- unique(macro[, .(
  cal_date, yyyyqq,
  pbrent      = macro_base_pbrent,
  yoy_oil     = macro_base_yoy_oil,
  lurc        = macro_base_lurc,
  pcpi        = macro_base_pcpi,
  rmtg        = macro_base_rmtg,
  fedfunds    = if ("macro_base_rff" %in% names(macro)) macro_base_rff else NA_real_,
  fomc_regime = if ("macro_base_fomc_regime" %in% names(macro))
                  macro_base_fomc_regime else NA_integer_,
  yield_curve = if ("macro_base_yield_curve" %in% names(macro))
                  macro_base_yield_curve else NA_real_
)])[order(cal_date)]

# Try to load severely adverse
panel_severe <- tryCatch(setDT(readRDS("Data/panel_severe.rds")),
                         error = function(e) NULL)
macro_severe <- tryCatch(setDT(readRDS("Data/macro_severe.rds")),
                         error = function(e) NULL)
has_severe <- !is.null(panel_severe) && !is.null(macro_severe)

if (has_severe) {
  add_date(macro_severe)
  msg("  panel_severe: %s rows (severely adverse scenario loaded)",
      format(nrow(panel_severe), big.mark=","))
} else {
  msg("  panel_severe: NOT found -- charts 11-15 will be skipped")
}

# ── Available outcomes ────────────────────────────────────────────────────────
cu_outcomes <- intersect(
  c("dq_rate","chg_tot_lns_ratio","netintmrg","pcanetworth","networth",
    "roa","costfds","insured_share_growth","cert_share","loan_to_share",
    "nim_spread","member_growth_yoy","pll_rate","pll_per_loan"),
  names(panel)
)

out_labels <- c(
  dq_rate              = "Delinquency Rate (%)",
  chg_tot_lns_ratio    = "Net Charge-Off Ratio (%)",
  netintmrg            = "Net Interest Margin (%)",
  costfds              = "Cost of Funds (%)",
  insured_share_growth = "Insured Share Growth (YoY%)",
  cert_share           = "Certificate Share of Deposits",
  loan_to_share        = "Loan-to-Share Ratio",
  roa                  = "Return on Assets (%)",
  member_growth_yoy    = "Membership Growth (YoY%)",
  pll_rate             = "PLL Rate (% of Avg Loans)",
  pll_per_loan         = "PLL per Loan ($)",
  pcanetworth          = "Net Worth Ratio (% of Assets)",
  networth             = "Net Worth ($000s)",
  nim_spread           = "NIM Spread (Loan Yield - CoF)"
)

msg("\n  CU outcomes available (%d):", length(cu_outcomes))
msg("    %s", paste(cu_outcomes, collapse=", "))
for (v in c("member_growth_yoy","pll_rate","pcanetworth")) {
  if (!v %in% cu_outcomes) msg("  NOTE: %s not found -- related charts will skip", v)
}

# ── Aggregation helpers ───────────────────────────────────────────────────────
agg_quarter <- function(dt, vars, by_vars = "yyyyqq") {
  dt[, c(list(cal_date = first(cal_date), year = first(year),
              quarter  = first(quarter)),
         lapply(.SD, function(x) mean(x, na.rm = TRUE))),
     by = by_vars, .SDcols = intersect(vars, names(dt))
  ][order(get(by_vars[1]))]
}

agg_group <- function(dt, vars, group_col, by_vars = c("yyyyqq", group_col)) {
  dt[!is.na(get(group_col)),
     c(list(cal_date = first(cal_date), year = first(year),
             quarter  = first(quarter)),
        lapply(.SD, function(x) mean(x, na.rm = TRUE))),
     by = by_vars, .SDcols = intersect(vars, names(dt))
  ][order(yyyyqq)]
}

agg <- agg_quarter(panel, cu_outcomes)
agg <- merge(agg, mac_spine[, .(yyyyqq, pbrent, yoy_oil)],
             by = "yyyyqq", all.x = TRUE)

msg("\n  Panel date range: %s to %s",
    format(min(panel$cal_date)), format(max(panel$cal_date)))

# ============================================================================
# CHART 01 -- Oil Price History + Cycle Annotation
# ============================================================================
hdr("CHART 01: Oil price history")

mac_hist <- mac_spine[!is.na(pbrent) & cal_date >= as.Date("2005-01-01")]

p_oil_level <- ggplot(mac_hist, aes(x = cal_date, y = pbrent)) +
  ep_rects() +
  geom_line(colour = COL_OIL, linewidth = 0.9) +
  geom_hline(yintercept = mean(mac_hist$pbrent, na.rm = TRUE),
             linetype = "dashed", colour = "#888888", linewidth = 0.4) +
  annotate("text", x = min(mac_hist$cal_date) + 60,
           y = mean(mac_hist$pbrent, na.rm = TRUE) + 4,
           label = sprintf("Period avg $%.0f/bbl",
                           mean(mac_hist$pbrent, na.rm = TRUE)),
           size = 2.8, colour = "#888888") +
  # Moody's $125 Iran stress reference line
  geom_hline(yintercept = 125, linetype = "dotted",
             colour = COL_NEG, linewidth = 0.5) +
  annotate("text", x = max(mac_hist$cal_date) - 200,
           y = 127, label = "Moody's $125 stress", size = 2.5,
           colour = COL_NEG, hjust = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = dollar_format(prefix = "$", suffix = "/bbl")) +
  labs(title = "Brent Crude Oil Price (PBRENT) -- 2005Q1 to 2025Q4",
       subtitle = paste("Shaded: GFC 2008-09 | Shale bust 2014-16 | COVID 2020 |",
                        "Post-COVID surge 2021-22 | Iran War 2026 (yellow)"),
       x = NULL, y = "$/barrel",
       caption = "Source: FRB CCAR 2026 Baseline | Dotted red = Moody's Analytics $125 Iran stress") +
  theme_pub()

p_oil_yoy <- ggplot(mac_hist[!is.na(yoy_oil)],
                     aes(x = cal_date, y = yoy_oil, fill = yoy_oil >= 0)) +
  geom_col(width = 70, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = COL_POS, "FALSE" = COL_NEG)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_hline(yintercept = 60.3, linetype = "dotted",
             colour = COL_NEG, linewidth = 0.5) +
  annotate("text", x = max(mac_hist$cal_date) - 200,
           y = 63, label = "+60.3pp Iran shock", size = 2.5,
           colour = COL_NEG, hjust = 1) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(title = "YoY % Change in Brent Oil Price",
       subtitle = "Dotted = Moody's $125 scenario implied +60.3pp YoY from $78 baseline",
       x = NULL, y = "YoY %",
       caption = "Source: FRB CCAR 2026 Baseline; Moody's Analytics 2026 CRE Stress") +
  theme_pub()

p01 <- p_oil_level / p_oil_yoy +
  plot_annotation(
    title = "FIGURE 01 -- Brent Crude Oil Price: Level & Annual Change (2005-2025)",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )
save_plot(p01, "01_oil_price_history.png", w = 11, h = 8)

# ============================================================================
# CHART 02 -- PBRENT vs Aggregate CU Outcomes
# ============================================================================
hdr("CHART 02: PBRENT vs CU outcomes")

make_dual_axis <- function(outcome_var, outcome_label, y_fmt = waiver()) {
  if (!outcome_var %in% names(agg)) return(NULL)
  d <- agg[!is.na(get(outcome_var)) & !is.na(pbrent)]
  if (nrow(d) < 10) return(NULL)
  ov_max  <- max(abs(d[[outcome_var]]), na.rm = TRUE)
  pb_max  <- max(abs(d$pbrent), na.rm = TRUE)
  oil_sc  <- if (is.finite(ov_max) && pb_max > 0) ov_max / pb_max else 1
  ggplot(d, aes(x = cal_date)) +
    ep_rects() +
    geom_line(aes(y = pbrent * oil_sc, colour = "PBRENT (scaled)"),
              linewidth = 0.6, linetype = "dashed") +
    geom_line(aes(y = get(outcome_var), colour = outcome_label),
              linewidth = 0.85) +
    scale_colour_manual(values = c("PBRENT (scaled)" = COL_OIL,
                                    setNames(COL_DIRECT, outcome_label)),
                        name = NULL) +
    scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
    scale_y_continuous(labels = y_fmt) +
    labs(title = outcome_label, x = NULL, y = outcome_label) +
    theme_pub() + theme(legend.position = "none")
}

panels_02 <- Filter(Negate(is.null), list(
  make_dual_axis("dq_rate",             "Delinquency Rate (%)",
                 percent_format(scale = 1)),
  make_dual_axis("pll_rate",            "PLL Rate (% of Avg Loans)",
                 number_format(accuracy = 0.01)),
  make_dual_axis("netintmrg",           "Net Interest Margin (%)",
                 number_format(accuracy = 0.1)),
  make_dual_axis("costfds",             "Cost of Funds (%)",
                 number_format(accuracy = 0.01)),
  make_dual_axis("insured_share_growth","Insured Share Growth (YoY%)",
                 number_format(accuracy = 0.1)),
  make_dual_axis("cert_share",          "Certificate Share",
                 percent_format(scale = 100)),
  make_dual_axis("loan_to_share",       "Loan-to-Share Ratio",
                 number_format(accuracy = 0.01)),
  make_dual_axis("member_growth_yoy",   "Membership Growth (YoY%)",
                 number_format(accuracy = 0.1))
))

if (length(panels_02) >= 2) {
  p02 <- wrap_plots(panels_02, ncol = 2) +
    plot_annotation(
      title    = "FIGURE 02 -- PBRENT vs Aggregate CU Outcomes (2005-2025)",
      subtitle = "PBRENT dashed orange (scaled to outcome axis) | Shaded: key oil episodes",
      caption  = "Source: NCUA Form 5300 Call Report; FRB CCAR 2026 Baseline",
      theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                       plot.subtitle = element_text(size = 9, colour = "#555"))
    )
  save_plot(p02, "02_pbrent_vs_cu_outcomes.png", w = 13, h = 10)
}

# ============================================================================
# CHART 03 -- Cross-Correlation: PBRENT leads CU outcomes
# ============================================================================
hdr("CHART 03: Cross-correlations")

agg_cc <- merge(agg_quarter(panel, cu_outcomes),
                mac_spine[, .(yyyyqq, yoy_oil)],
                by = "yyyyqq", all.x = TRUE)

cc_results <- rbindlist(lapply(cu_outcomes, function(v) {
  if (!v %in% names(agg_cc)) return(NULL)
  x  <- agg_cc$yoy_oil; y <- agg_cc[[v]]
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 20) return(NULL)
  lags <- -4:8
  cors <- sapply(lags, function(k) {
    if (k >= 0) cor(x[ok][1:(sum(ok)-k)], y[ok][(k+1):sum(ok)], use="complete.obs")
    else        cor(x[ok][(-k+1):sum(ok)], y[ok][1:(sum(ok)+k)], use="complete.obs")
  })
  data.table(outcome = v, lag = lags, correlation = cors)
}))

cc_results[, outcome_label := out_labels[outcome]]
cc_results[is.na(outcome_label), outcome_label := outcome]

cat("\n  Cross-correlation max|r| by outcome:\n")
cc_cov <- cc_results[, .(max_abs_r = round(max(abs(correlation), na.rm=TRUE), 3),
                           best_lag  = lag[which.max(abs(correlation))]),
                       by = .(outcome, outcome_label)][order(-max_abs_r)]
print(cc_cov, row.names = FALSE)

CC_COLS <- setNames(
  grDevices::colorRampPalette(c("#1b9e77","#d95f02","#7570b3","#e7298a",
                                 "#66a61e","#e6ab02","#a6761d","#666666",
                                 "#4a90d9","#c0392b","#2d7a4a","#8e44ad",
                                 "#16a085","#b5470a"))(length(cu_outcomes)),
  vapply(cu_outcomes, function(v) out_labels[v] %||% v, character(1))
)

p03 <- ggplot(cc_results[!is.na(correlation)],
               aes(x = lag, y = correlation,
                   colour = outcome_label, group = outcome_label)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "#888888") +
  geom_vline(xintercept = 0, linetype = "dashed",
             linewidth = 0.4, colour = "#cccccc") +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_colour_manual(values = CC_COLS, name = NULL) +
  scale_x_continuous(breaks = -4:8,
                     labels = c("-4","-3","-2","-1","0","+1","+2","+3",
                                "+4","+5","+6","+7","+8")) +
  labs(title    = "FIGURE 03 -- Cross-Correlation: PBRENT YoY leads CU Outcomes",
       subtitle = paste("Negative lag = oil leads outcome | Positive lag = oil follows outcome",
                        "| Key: deposit growth peaks ~Q+1, delinquency peaks ~Q+2-4"),
       caption  = "Source: NCUA Form 5300; FRB CCAR 2026 Baseline | Quarterly panel means",
       x        = "Lag (quarters, negative = oil leads)",
       y        = "Pearson Correlation") +
  theme_pub() +
  theme(legend.position  = "right",
        legend.text      = element_text(size = 7.5),
        legend.key.size  = unit(0.4, "cm"))

save_plot(p03, "03_cross_correlations.png", w = 13, h = 7)

# ============================================================================
# CHART 04 -- Direct vs Indirect: oil-state vs non-oil CUs
# ============================================================================
hdr("CHART 04: Direct vs indirect channel")

group_col_04 <- intersect(c("cu_group","oil_exposure_bin"), names(panel))[1]

if (!is.na(group_col_04)) {
  panel04 <- copy(panel)

  if (group_col_04 == "cu_group") {
    panel04[, plot_group := fcase(
      cu_group == "Direct",   "Oil States (Direct)",
      cu_group == "Indirect", "Non-Oil (Indirect)",
      default                = "Non-Oil (Negligible)"
    )]
    grp_colors <- c("Oil States (Direct)"    = COL_DIRECT,
                    "Non-Oil (Indirect)"     = COL_INDIR,
                    "Non-Oil (Negligible)"   = "#aaaaaa")
  } else {
    panel04[, plot_group := fifelse(oil_exposure_bin == 1L,
                                    "Oil States", "Non-Oil States")]
    grp_colors <- c("Oil States"   = COL_DIRECT,
                    "Non-Oil States" = COL_INDIR)
  }

  grp_agg04 <- agg_group(panel04,
                          intersect(c("dq_rate","pll_rate","netintmrg",
                                      "insured_share_growth","member_growth_yoy",
                                      "costfds","cert_share","loan_to_share"),
                                    names(panel04)),
                          "plot_group")
  grp_agg04 <- merge(grp_agg04, mac_spine[, .(yyyyqq, pbrent)],
                     by = "yyyyqq", all.x = TRUE)

  make_grp04 <- function(v, lab, y_fmt = waiver()) {
    if (!v %in% names(grp_agg04)) return(NULL)
    ggplot(grp_agg04[!is.na(get(v))],
           aes(x = cal_date, y = get(v), colour = plot_group)) +
      ep_rects() +
      geom_line(linewidth = 0.85) +
      geom_hline(yintercept = 0, linewidth = 0.3, colour = "#cccccc") +
      scale_colour_manual(values = grp_colors, name = NULL) +
      scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
      scale_y_continuous(labels = y_fmt) +
      labs(title = lab, x = NULL, y = lab) +
      theme_pub() + theme(legend.position = "bottom")
  }

  p04_panels <- Filter(Negate(is.null), list(
    make_grp04("dq_rate",             "Delinquency Rate (%)"),
    make_grp04("pll_rate",            "PLL Rate (% of Avg Loans)",
               number_format(accuracy = 0.01)),
    make_grp04("netintmrg",           "Net Interest Margin (%)"),
    make_grp04("insured_share_growth","Insured Share Growth (YoY%)"),
    make_grp04("member_growth_yoy",   "Membership Growth (YoY%)"),
    make_grp04("costfds",             "Cost of Funds (%)"),
    make_grp04("cert_share",          "Certificate Share",
               percent_format(scale = 100)),
    make_grp04("loan_to_share",       "Loan-to-Share Ratio")
  ))

  if (length(p04_panels) >= 2) {
    p04 <- wrap_plots(p04_panels, ncol = 2, guides = "collect") &
      theme(legend.position = "bottom")
    p04 <- p04 + plot_annotation(
      title    = "FIGURE 04 -- Direct vs Indirect: Oil-State vs Non-Oil CUs (2005-2025)",
      subtitle = "Oil-State = mining emp share >= 2% (BLS QCEW NAICS 211) | v1 finding: direct 82x > indirect",
      caption  = "Source: NCUA Form 5300; BLS QCEW oil exposure classification",
      theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                       plot.subtitle = element_text(size = 9, colour = "#555"))
    )
    save_plot(p04, "04_direct_vs_indirect.png", w = 13, h = 10)
  }
}

# ============================================================================
# CHART 05 -- Deposit Channel: cert_share, loan_to_share vs FOMC
# ============================================================================
hdr("CHART 05: Deposit channel")

dep_vars <- intersect(c("cert_share","loan_to_share","insured_share_growth","costfds"),
                      names(panel))

if (length(dep_vars) >= 2) {
  dep_agg <- agg_quarter(panel, dep_vars)
  dep_agg <- merge(dep_agg,
                   mac_spine[, .(yyyyqq, pbrent, fedfunds, fomc_regime)],
                   by = "yyyyqq", all.x = TRUE)

  p_cert <- if ("cert_share" %in% names(dep_agg)) {
    ggplot(dep_agg[!is.na(cert_share)], aes(x = cal_date, y = cert_share)) +
      ep_rects() +
      geom_line(colour = COL_DIRECT, linewidth = 0.85) +
      geom_line(aes(y = fedfunds / 100, colour = "Fed Funds Rate (RHS)"),
                linetype = "dashed", linewidth = 0.6, na.rm = TRUE) +
      scale_colour_manual(values = c("Fed Funds Rate (RHS)" = COL_OIL),
                          name = NULL) +
      scale_y_continuous(labels = percent_format(scale = 100),
                         limits = c(0, NA),
                         sec.axis = sec_axis(~. * 100,
                                             name = "Fed Funds Rate (%)")) +
      scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
      labs(title    = "Certificate Share of Total Deposits vs Fed Funds Rate",
           subtitle = "Rising rates -> members shift demand deposits to certificates (deposit migration)",
           x = NULL, y = "Certificate Share") +
      theme_pub()
  } else NULL

  p_lts <- if ("loan_to_share" %in% names(dep_agg)) {
    ggplot(dep_agg[!is.na(loan_to_share)],
           aes(x = cal_date, y = loan_to_share)) +
      ep_rects() +
      geom_line(colour = COL_INDIR, linewidth = 0.85) +
      geom_hline(yintercept = 0.80, linetype = "dashed",
                 colour = COL_NEG, linewidth = 0.4) +
      annotate("text", x = max(dep_agg$cal_date, na.rm = TRUE),
               y = 0.81, label = "80% liquidity threshold",
               hjust = 1, size = 2.8, colour = COL_NEG) +
      scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
      scale_y_continuous(labels = number_format(accuracy = 0.01)) +
      labs(title    = "Loan-to-Share Ratio",
           subtitle = "Oil-state CUs: PBRENT up -> deposit surge -> ratio compression",
           x = NULL, y = "Loan / Total Shares") +
      theme_pub()
  } else NULL

  p_cof <- if ("costfds" %in% names(dep_agg)) {
    ggplot(dep_agg[!is.na(costfds)], aes(x = cal_date, y = costfds)) +
      ep_rects() +
      geom_line(colour = COL_SPILL, linewidth = 0.85) +
      scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
      scale_y_continuous(labels = number_format(accuracy = 0.01)) +
      labs(title    = "Cost of Funds (%)",
           subtitle = "Certificate surge raises funding costs even as deposit volumes grow",
           x = NULL, y = "Cost of Funds (%)") +
      theme_pub()
  } else NULL

  p_isg <- if ("insured_share_growth" %in% names(dep_agg)) {
    ggplot(dep_agg[!is.na(insured_share_growth)],
           aes(x = cal_date, y = insured_share_growth,
               fill = insured_share_growth >= 0)) +
      geom_col(width = 70, show.legend = FALSE) +
      scale_fill_manual(values = c("TRUE" = COL_POS, "FALSE" = COL_NEG)) +
      scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
      scale_y_continuous(labels = number_format(accuracy = 0.1)) +
      labs(title    = "Insured Share Growth (YoY%, winsorised 1-99)",
           subtitle = "Oil-state income windfall -> deposit inflows; inflation erosion -> drawdown",
           x = NULL, y = "YoY Growth (%)") +
      theme_pub()
  } else NULL

  active_panels <- Filter(Negate(is.null), list(p_cert, p_lts, p_cof, p_isg))
  if (length(active_panels) >= 2) {
    p05 <- wrap_plots(active_panels, ncol = 2) +
      plot_annotation(
        title   = "FIGURE 05 -- Deposit Channel Dynamics",
        caption = "Source: NCUA Form 5300; FRB CCAR 2026 Baseline",
        theme   = theme(plot.title = element_text(face = "bold", size = 12))
      )
    save_plot(p05, "05_deposit_channel.png", w = 13, h = 10)
  }
}

# ============================================================================
# CHART 06 -- Asset Tier Response: all 8 tiers (assets_cat2)
# ============================================================================
hdr("CHART 06: Asset tier response (8-tier assets_cat2)")

tier_col_06 <- intersect(c("asset_tier","assets_cat2"), names(panel))[1]

if (!is.na(tier_col_06)) {
  tier_agg <- agg_group(panel,
                         intersect(c("dq_rate","netintmrg","costfds",
                                     "cert_share","insured_share_growth",
                                     "member_growth_yoy","pll_rate",
                                     "pcanetworth"),
                                   names(panel)),
                         tier_col_06)
  tier_agg <- merge(tier_agg, mac_spine[, .(yyyyqq, pbrent, yoy_oil)],
                    by = "yyyyqq", all.x = TRUE)

  # Map to T-codes if needed
  if (tier_col_06 == "assets_cat2") {
    cat2_to_t <- c(
      "Assets < 10 million"    = "T1_under10M",
      "10M up to 50M"          = "T2_10to50M",
      "50M through 100M"       = "T3_50to100M",
      "Over 100M through 500M" = "T4_100to500M",
      "Over 500M through 1B"   = "T5_500Mto1B",
      "Over 1B through 5B"     = "T6_1Bto5B",
      "Over 5B through 10B"    = "T7_5Bto10B",
      "Over 10B"               = "T8_over10B"
    )
    tier_agg[, tier_code := cat2_to_t[as.character(get(tier_col_06))]]
    tier_var_06 <- "tier_code"
  } else {
    tier_var_06 <- "asset_tier"
  }

  tier_cols_use <- TIER_COLS[intersect(names(TIER_COLS),
                                        unique(tier_agg[[tier_var_06]]))]

  make_tier_plot <- function(v, lab, y_fmt = waiver()) {
    if (!v %in% names(tier_agg)) return(NULL)
    ggplot(tier_agg[!is.na(get(v))],
           aes(x = cal_date, y = get(v),
               colour = get(tier_var_06))) +
      ep_rects() +
      geom_line(linewidth = 0.6) +
      scale_colour_manual(values = tier_cols_use,
                          labels = TIER_LABELS[names(tier_cols_use)],
                          name   = "Tier") +
      scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
      scale_y_continuous(labels = y_fmt) +
      labs(title = lab, x = NULL, y = lab) +
      theme_pub() + theme(legend.position = "right",
                          legend.text = element_text(size = 7))
  }

  t_panels <- Filter(Negate(is.null), list(
    make_tier_plot("dq_rate",            "Delinquency Rate (%)"),
    make_tier_plot("pll_rate",           "PLL Rate (% of Avg Loans)",
                   number_format(accuracy = 0.01)),
    make_tier_plot("netintmrg",          "Net Interest Margin (%)"),
    make_tier_plot("costfds",            "Cost of Funds (%)"),
    make_tier_plot("cert_share",         "Certificate Share",
                   percent_format(scale = 100, accuracy = 0.1)),
    make_tier_plot("member_growth_yoy",  "Membership Growth (YoY%)",
                   number_format(accuracy = 0.1))
  ))

  if (length(t_panels) >= 2) {
    p06 <- wrap_plots(t_panels, ncol = 2) +
      plot_annotation(
        title    = "FIGURE 06 -- CU Outcomes by Asset Tier -- 8-Tier assets_cat2 (2005-2025)",
        subtitle = paste("T1 <$10M | T2 $10-50M | T3 $50-100M | T4 $100-500M |",
                         "T5 $500M-$1B | T6 $1-5B | T7 $5-10B | T8 >$10B"),
        caption  = "Source: NCUA Form 5300 Call Report | assets_cat2 from OCE_combined",
        theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                         plot.subtitle = element_text(size = 8.5, colour = "#555"))
      )
    save_plot(p06, "06_asset_tier_response.png", w = 14, h = 10)
  }
}

# ============================================================================
# CHART 07 -- Structural Break: pre vs post 2015Q1
# ============================================================================
hdr("CHART 07: Structural break 2015Q1")

agg_sb <- merge(agg_quarter(panel, cu_outcomes),
                mac_spine[, .(yyyyqq, yoy_oil)], by = "yyyyqq", all.x = TRUE)
agg_sb[, era := fifelse(yyyyqq < 201501L,
                         "Pre-Shale (2005-2014)",
                         "Post-Shale (2015-2025)")]

sb_long <- melt(agg_sb[!is.na(yoy_oil)],
                id.vars      = c("yyyyqq","cal_date","era","yoy_oil"),
                measure.vars = intersect(cu_outcomes, names(agg_sb)),
                variable.name = "outcome", value.name = "value")
sb_long[, outcome_label := out_labels[as.character(outcome)]]
sb_long[is.na(outcome_label), outcome_label := as.character(outcome)]

sb_plots <- lapply(
  intersect(c("dq_rate","pll_rate","netintmrg","costfds",
              "insured_share_growth","member_growth_yoy"), cu_outcomes),
  function(v) {
    d <- sb_long[outcome == v & !is.na(value) & !is.na(yoy_oil)]
    if (nrow(d) < 10) return(NULL)
    ggplot(d, aes(x = yoy_oil, y = value, colour = era)) +
      geom_point(alpha = 0.45, size = 1.2) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
      geom_hline(yintercept = 0, linewidth = 0.3, colour = "#aaa") +
      geom_vline(xintercept = 0, linewidth = 0.3, colour = "#aaa") +
      scale_colour_manual(
        values = c("Pre-Shale (2005-2014)" = "#1a3a5c",
                   "Post-Shale (2015-2025)" = "#b5470a"),
        name = "Era") +
      labs(title = unique(d$outcome_label),
           x     = "PBRENT YoY %",
           y     = unique(d$outcome_label)) +
      theme_pub() + theme(legend.position = "right")
  }
)
sb_plots <- Filter(Negate(is.null), sb_plots)

if (length(sb_plots) >= 2) {
  p07 <- wrap_plots(sb_plots, ncol = 2) +
    plot_annotation(
      title    = "FIGURE 07 -- Structural Break: Oil-CU Relationship Pre vs Post Shale Era (2015Q1)",
      subtitle = "ANCHOR: slope reversal at 2015Q1 confirmed in v1 across all outcomes",
      caption  = "Source: NCUA Form 5300; FRB CCAR 2026 Baseline | OLS with 95% CI",
      theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                       plot.subtitle = element_text(size = 9, colour = "#555"))
    )
  save_plot(p07, "07_structural_break.png", w = 13, h = 9)
}

# ============================================================================
# CHART 08 -- Spillover: non-oil CU response by spillover tercile
# ============================================================================
hdr("CHART 08: Spillover exposure terciles")

if ("spillover_exposure" %in% names(panel) &&
    "oil_exposure_bin" %in% names(panel)) {

  non_oil <- panel[oil_exposure_bin == 0 & !is.na(spillover_exposure)]

  if (nrow(non_oil) > 100) {
    breaks <- unique(quantile(non_oil$spillover_exposure,
                              probs = c(0, 1/3, 2/3, 1), na.rm = TRUE))
    if (length(breaks) >= 3) {
      n_lab <- length(breaks) - 1
      non_oil[, spill_tercile := cut(spillover_exposure, breaks = breaks,
                                      labels = c("Low","Medium","High")[1:n_lab],
                                      include.lowest = TRUE)]
    } else {
      med <- median(non_oil$spillover_exposure, na.rm = TRUE)
      non_oil[, spill_tercile := factor(
        fifelse(spillover_exposure > med, "High", "Low"),
        levels = c("Low","High")
      )]
      msg("  NOTE: low variance spillover -- using median split")
    }

    spill_agg <- agg_group(non_oil,
                            intersect(c("dq_rate","pll_rate","netintmrg",
                                        "insured_share_growth","costfds",
                                        "member_growth_yoy"), names(non_oil)),
                            "spill_tercile")
    spill_agg <- merge(spill_agg, mac_spine[, .(yyyyqq, yoy_oil)],
                       by = "yyyyqq", all.x = TRUE)

    SPILL_COLS <- c("Low"    = "#a8d8a8",
                    "Medium" = "#4a9a6a",
                    "High"   = "#1a5a3a")

    sp_panels <- lapply(
      intersect(c("dq_rate","pll_rate","netintmrg","member_growth_yoy",
                  "insured_share_growth","costfds"),
                names(spill_agg)),
      function(v) {
        lab <- out_labels[v] %||% v
        ggplot(spill_agg[!is.na(get(v)) & !is.na(spill_tercile)],
               aes(x = cal_date, y = get(v), colour = spill_tercile)) +
          ep_rects() +
          geom_line(linewidth = 0.8) +
          scale_colour_manual(values = SPILL_COLS, name = "Spillover\nTercile") +
          scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
          labs(title = lab, x = NULL, y = lab) + theme_pub()
      }
    )

    if (length(sp_panels) >= 1) {
      p08 <- wrap_plots(sp_panels, ncol = 2) +
        plot_annotation(
          title    = "FIGURE 08 -- Indirect Channel: Non-Oil CU Response by Spillover Tercile",
          subtitle = "High spillover = non-oil state with strong economic linkage to oil states",
          caption  = "Source: NCUA Form 5300; BLS QCEW adjacency-weighted spillover index",
          theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                           plot.subtitle = element_text(size = 9, colour = "#555"))
        )
      save_plot(p08, "08_spillover_terciles.png", w = 13, h = 9)
    }
  }
}

# ============================================================================
# CHART 09 -- Episode Heatmap
# ============================================================================
hdr("CHART 09: Episode heatmap")

ep_def <- data.table(
  episode      = c("Pre-GFC\n2005-07","GFC\n2008-09","Recovery\n2010-13",
                   "Shale Bust\n2014-16","Rebound\n2017-19","COVID\n2020",
                   "Surge\n2021-22","Post-Surge\n2023-25"),
  yyyyqq_from  = c(200501L,200801L,201001L,201401L,201701L,202001L,
                   202101L,202301L),
  yyyyqq_to    = c(200704L,200904L,201304L,201604L,201904L,202004L,
                   202204L,202504L)
)

agg_ep <- merge(agg_quarter(panel, cu_outcomes),
                mac_spine[, .(yyyyqq, yoy_oil)], by = "yyyyqq", all.x = TRUE)

heat_list <- lapply(1:nrow(ep_def), function(i) {
  ep  <- ep_def[i]
  sub <- agg_ep[yyyyqq >= ep$yyyyqq_from & yyyyqq <= ep$yyyyqq_to]
  if (nrow(sub) == 0) return(NULL)
  means <- sapply(cu_outcomes, function(v) {
    if (!v %in% names(sub)) return(NA_real_)
    val <- mean(sub[[v]], na.rm = TRUE)
    if (!is.finite(val)) NA_real_ else val
  })
  as.data.table(c(list(episode = ep$episode), as.list(means)))
})
heat_dt <- rbindlist(heat_list, fill = TRUE)

heat_long <- melt(heat_dt, id.vars = "episode",
                  variable.name = "outcome", value.name = "value")
heat_long[, outcome_label := out_labels[as.character(outcome)]]
heat_long[is.na(outcome_label), outcome_label := as.character(outcome)]
heat_long[, norm_val := {
  mn <- min(value, na.rm = TRUE); mx <- max(value, na.rm = TRUE)
  if (is.finite(mn) && is.finite(mx) && mx > mn) (value - mn) / (mx - mn)
  else rep(0.5, .N)
}, by = outcome]

# Invert stress outcomes: higher = redder
stress_outcomes <- c("dq_rate","chg_tot_lns_ratio","costfds","pll_rate","pll_per_loan")
heat_long[outcome %in% stress_outcomes, norm_val := 1 - norm_val]

heat_long[, value_label := fifelse(
  is.na(value) | !is.finite(value), "",
  as.character(round(value, 2))
)]
heat_long[, episode := factor(episode, levels = ep_def$episode)]

p09 <- ggplot(heat_long, aes(x = episode, y = outcome_label, fill = norm_val)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_tile(data = heat_long[is.na(norm_val)],
            fill = "#eeeeee", colour = "white", linewidth = 0.5) +
  geom_text(aes(label = value_label), size = 2.5, colour = "#1a1a1a") +
  geom_text(data = heat_long[is.na(norm_val)],
            aes(label = "N/A"), size = 2.2, colour = "#aaaaaa", fontface = "italic") +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#d73027",
                       midpoint = 0.5, na.value = "#eeeeee",
                       name = "Relative level\n(Blue=good\nRed=stress)") +
  scale_x_discrete(position = "top") +
  labs(title    = "FIGURE 09 -- CU Outcome Heatmap by Oil Price Episode",
       subtitle = "Blue = favourable | Red = stress | White = average | Grey = N/A",
       caption  = "Source: NCUA Form 5300; FRB CCAR 2026 | Episode mean values shown",
       x = NULL, y = NULL) +
  theme_pub() +
  theme(axis.text.x  = element_text(size = 8.5, angle = 0, face = "bold"),
        axis.text.y  = element_text(size = 8.5),
        panel.grid   = element_blank(),
        legend.position = "right")

save_plot(p09, "09_episode_heatmap.png", w = 14, h = 8)

# ============================================================================
# CHART 10 -- Missingness Overview
# ============================================================================
hdr("CHART 10: Missingness overview")

check_miss <- c(
  "dq_rate","chg_tot_lns_ratio","netintmrg","pcanetworth",
  "networth","costfds","roa","insured_tot","dep_shrcert","acct_018",
  "insured_share_growth","cert_share","loan_to_share","nim_spread",
  "members","member_growth_yoy",
  "pll","lns_tot","lns_tot_n","pll_rate","pll_per_loan",
  "macro_base_pbrent","macro_base_lurc","macro_base_pcpi",
  "macro_base_rmtg","macro_base_phpi","macro_base_yoy_oil",
  "macro_base_yield_curve","macro_base_fomc_regime",
  "oil_exposure_cont","oil_exposure_bin","spillover_exposure",
  "oil_x_brent","fomc_x_brent","oil_bartik_iv",
  "asset_tier","assets_cat2"
)
fv  <- intersect(check_miss, names(panel))
pct <- sapply(panel[, ..fv], function(x) round(mean(is.na(x)) * 100, 1))
miss_tbl <- data.table(
  variable    = names(pct),
  pct_missing = pct,
  category    = fcase(
    names(pct) %like% "macro_", "CCAR Macro",
    names(pct) %like% "oil_|spill|brent|fomc_x|bartik", "Exposure/Interaction",
    names(pct) %like% "asset", "Tier",
    default = "CU Outcome"
  )
)[order(-pct_missing)]

cat("\n  Missingness table:\n")
print(miss_tbl[pct_missing > 0], row.names = FALSE)

p10 <- ggplot(miss_tbl, aes(x = pct_missing,
                              y = reorder(variable, -pct_missing),
                              fill = category)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = c(5, 10, 20), linetype = "dashed",
             colour = "#bbbbbb", linewidth = 0.3) +
  scale_fill_manual(
    values = c("CU Outcome"           = "#1a3a5c",
               "CCAR Macro"           = "#b5470a",
               "Exposure/Interaction" = "#2d7a4a",
               "Tier"                 = "#8e44ad"),
    name = "Category") +
  scale_x_continuous(labels = percent_format(scale = 1),
                     limits = c(0, max(miss_tbl$pct_missing) * 1.1)) +
  labs(title    = "FIGURE 10 -- Variable Missingness Overview (panel_base)",
       subtitle = ">10% missing: investigate before modelling | >20%: consider exclusion",
       caption  = "Source: NCUA Form 5300 + FRB CCAR 2026 Baseline merged panel",
       x = "% Missing", y = NULL) +
  theme_pub() + theme(legend.position = "right",
                      panel.grid.major.y = element_blank())

save_plot(p10, "10_missingness.png", w = 10, h = 7)

# ============================================================================
# CHART NEW-A -- Membership Growth vs Oil Price + Direct/Indirect
# ============================================================================
hdr("CHART NEW-A: Membership growth vs oil price")

if ("member_growth_yoy" %in% names(panel)) {
  mem_agg <- merge(agg_quarter(panel, "member_growth_yoy"),
                   mac_spine[, .(yyyyqq, pbrent, yoy_oil)],
                   by = "yyyyqq", all.x = TRUE)

  mem_max <- max(abs(mem_agg$member_growth_yoy), na.rm = TRUE)
  pb_max  <- max(abs(mem_agg$pbrent), na.rm = TRUE)
  oil_sc  <- if (is.finite(mem_max) && pb_max > 0) mem_max / pb_max else 1

  pA1 <- ggplot(mem_agg[!is.na(member_growth_yoy)], aes(x = cal_date)) +
    ep_rects() +
    geom_col(aes(y = member_growth_yoy, fill = member_growth_yoy >= 0),
             width = 70, show.legend = FALSE) +
    geom_line(aes(y = pbrent * oil_sc, colour = "PBRENT (scaled)"),
              linewidth = 0.7, linetype = "dashed", na.rm = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = c("TRUE" = COL_POS, "FALSE" = COL_NEG)) +
    scale_colour_manual(values = c("PBRENT (scaled)" = COL_OIL), name = NULL) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = number_format(accuracy = 0.1),
                       sec.axis = sec_axis(~. / oil_sc,
                                            name = "PBRENT ($/bbl)",
                                            labels = dollar_format(suffix = "/bbl"))) +
    labs(title    = "Membership Growth (YoY%) vs Brent Oil Price",
         subtitle = "Green = expanding | Red = contracting | Dashed = oil price (right axis)",
         x = NULL, y = "Member Growth YoY (%)") +
    theme_pub()

  if (!is.na(group_col_04) && "plot_group" %in% names(panel04)) {
    mem_grp <- panel04[!is.na(member_growth_yoy) & !is.na(plot_group),
                        .(member_growth_yoy = mean(member_growth_yoy, na.rm=TRUE),
                          cal_date = first(cal_date)),
                        by = .(yyyyqq, plot_group)][order(yyyyqq)]

    pA2 <- ggplot(mem_grp,
                  aes(x = cal_date, y = member_growth_yoy,
                      colour = plot_group)) +
      ep_rects() +
      geom_line(linewidth = 0.85) +
      geom_hline(yintercept = 0, linewidth = 0.4, colour = "#888") +
      scale_colour_manual(values = grp_colors, name = NULL) +
      scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
      labs(title    = "Membership Growth: Oil-State vs Non-Oil CUs",
           subtitle = "Buffer hypothesis: oil-state CUs should gain members when oil rises",
           x = NULL, y = "Member Growth YoY (%)") +
      theme_pub() + theme(legend.position = "bottom")
    pNA <- pA1 / pA2
  } else {
    pNA <- pA1
  }

  pNA <- pNA + plot_annotation(
    title    = "FIGURE NEW-A -- Membership Growth: Oil Price Transmission Channel",
    subtitle = "Membership growth is a leading indicator -- members leave before deposits do",
    caption  = "Source: NCUA Form 5300; FRB CCAR 2026 Baseline",
    theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                     plot.subtitle = element_text(size = 9, colour = "#555"))
  )
  save_plot(pNA, "NEW_A_member_growth_vs_oil.png", w = 13, h = 10)
}

# ============================================================================
# CHART NEW-B -- Membership Growth Structural Break
# ============================================================================
hdr("CHART NEW-B: Membership growth structural break")

if ("member_growth_yoy" %in% names(panel)) {
  mem_sb <- merge(agg_quarter(panel, "member_growth_yoy"),
                  mac_spine[, .(yyyyqq, yoy_oil)],
                  by = "yyyyqq", all.x = TRUE)
  mem_sb[, era := fifelse(yyyyqq < 201501L,
                           "Pre-Shale (2005-2014)",
                           "Post-Shale (2015-2025)")]

  pNB <- ggplot(mem_sb[!is.na(member_growth_yoy) & !is.na(yoy_oil)],
                aes(x = yoy_oil, y = member_growth_yoy, colour = era)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.0) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "#888") +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "#888") +
    scale_colour_manual(
      values = c("Pre-Shale (2005-2014)" = "#1a3a5c",
                 "Post-Shale (2015-2025)" = "#b5470a"),
      name = "Era") +
    scale_x_continuous(labels = number_format(accuracy = 1, suffix = "%")) +
    scale_y_continuous(labels = number_format(accuracy = 0.1, suffix = "%")) +
    labs(title    = "FIGURE NEW-B -- Membership Growth: Structural Break Pre vs Post Shale Era",
         subtitle = paste("Pre-shale: oil rise attracts members (positive slope).",
                          "Post-shale: inflation erodes income -> slope flattens or reverses."),
         caption  = "Source: NCUA Form 5300; FRB CCAR 2026 | OLS with 95% CI",
         x        = "PBRENT YoY %",
         y        = "Membership Growth YoY %") +
    theme_pub() + theme(legend.position = "right")

  save_plot(pNB, "NEW_B_member_growth_structural_break.png", w = 11, h = 7)
}

# ============================================================================
# CHART NEW-C -- Membership Growth by Asset Tier (8-tier)
# ============================================================================
hdr("CHART NEW-C: Membership growth by asset tier")

if ("member_growth_yoy" %in% names(panel) && !is.na(tier_col_06)) {
  tier_mem <- agg_group(panel, "member_growth_yoy", tier_col_06)
  tier_mem <- merge(tier_mem, mac_spine[, .(yyyyqq, yoy_oil)],
                    by = "yyyyqq", all.x = TRUE)
  if (tier_col_06 == "assets_cat2" && "tier_code" %in% names(tier_agg)) {
    tier_mem[, tier_code := cat2_to_t[as.character(get(tier_col_06))]]
    tier_mem_col <- "tier_code"
  } else {
    tier_mem_col <- tier_col_06
  }

  tc_use_nc <- TIER_COLS[intersect(names(TIER_COLS), unique(tier_mem[[tier_mem_col]]))]

  pNC1 <- ggplot(tier_mem[!is.na(member_growth_yoy)],
                 aes(x = cal_date, y = member_growth_yoy,
                     colour = get(tier_mem_col))) +
    ep_rects() +
    geom_line(linewidth = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "#888") +
    scale_colour_manual(values = tc_use_nc,
                        labels = TIER_LABELS[names(tc_use_nc)],
                        name   = "Tier") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = number_format(accuracy = 0.1)) +
    labs(title    = "Membership Growth by Asset Tier (Time Series)",
         subtitle = "Small CUs (T1/T2) more volatile -- single employer/geography concentration",
         x = NULL, y = "Member Growth YoY (%)") +
    theme_pub() + theme(legend.position = "right",
                        legend.text = element_text(size = 7.5))

  pNC2 <- ggplot(panel[!is.na(member_growth_yoy) & !is.na(get(tier_col_06))],
                 aes(x = get(tier_col_06),
                     y = member_growth_yoy,
                     fill = get(tier_col_06))) +
    geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3,
                 linewidth = 0.5, width = 0.6) +
    geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "#888") +
    scale_fill_manual(values = TIER_COLS, guide = "none") +
    scale_y_continuous(labels = number_format(accuracy = 0.1)) +
    labs(title    = "Membership Growth Distribution by Asset Tier",
         subtitle = "Interquartile range reveals size-dependent volatility",
         x        = "Asset Tier",
         y        = "Member Growth YoY (%)") +
    theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

  pNC <- pNC1 / pNC2 + plot_annotation(
    title    = "FIGURE NEW-C -- Membership Growth by Asset Tier (2005-2025)",
    subtitle = paste("T1 <$10M | T2 $10-50M | T3 $50-100M | T4 $100-500M |",
                     "T5 $500M-$1B | T6 $1-5B | T7 $5-10B | T8 >$10B"),
    caption  = "Source: NCUA Form 5300 | 8-tier assets_cat2 from OCE_combined",
    theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                     plot.subtitle = element_text(size = 8.5, colour = "#555"))
  )
  save_plot(pNC, "NEW_C_member_growth_by_tier.png", w = 13, h = 10)
}

# ============================================================================
# CHART NEW-D -- Membership Growth Lead/Lag Cross-Correlation
# ============================================================================
hdr("CHART NEW-D: Membership growth cross-correlation")

if ("member_growth_yoy" %in% names(panel)) {
  mem_cc_agg <- merge(agg_quarter(panel, c("member_growth_yoy",
                                            "insured_share_growth")),
                      mac_spine[, .(yyyyqq, yoy_oil)],
                      by = "yyyyqq", all.x = TRUE)

  compute_cc <- function(x, y, lags = -4:8) {
    ok <- !is.na(x) & !is.na(y)
    sapply(lags, function(k) {
      if (k >= 0) cor(x[ok][1:(sum(ok)-k)], y[ok][(k+1):sum(ok)], use="complete.obs")
      else        cor(x[ok][(-k+1):sum(ok)], y[ok][1:(sum(ok)+k)], use="complete.obs")
    })
  }

  lags_vec <- -4:8
  cc_mem_oil <- data.table(lag = lags_vec, r = compute_cc(
    mem_cc_agg$yoy_oil, mem_cc_agg$member_growth_yoy, lags_vec),
    pair = "Oil YoY -> Member Growth")

  cc_mem_dep <- if ("insured_share_growth" %in% names(mem_cc_agg)) {
    data.table(lag = lags_vec, r = compute_cc(
      mem_cc_agg$member_growth_yoy,
      mem_cc_agg$insured_share_growth, lags_vec),
      pair = "Member Growth -> Deposit Growth")
  } else NULL

  cc_nd <- rbindlist(Filter(Negate(is.null), list(cc_mem_oil, cc_mem_dep)))

  pND <- ggplot(cc_nd[!is.na(r)],
                aes(x = lag, y = r, colour = pair, group = pair)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "#888") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "#cccccc") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = c("Oil YoY -> Member Growth"       = COL_OIL,
                                    "Member Growth -> Deposit Growth" = COL_DIRECT),
                        name = NULL) +
    scale_x_continuous(breaks = lags_vec) +
    labs(title    = "FIGURE NEW-D -- Lead/Lag: Membership Growth vs Oil & Deposit Growth",
         subtitle = paste("Left panel: oil price leads/lags membership growth.",
                          "Right panel: membership leads/lags deposit growth.",
                          "Negative lag = first series leads."),
         caption  = "Source: NCUA Form 5300; FRB CCAR 2026 | Quarterly panel means",
         x        = "Lag (quarters, negative = first series leads)",
         y        = "Pearson Correlation") +
    theme_pub() + theme(legend.position = "bottom")

  save_plot(pND, "NEW_D_member_growth_crosscorr.png", w = 11, h = 6)
}

# ============================================================================
# CHART 16 -- PLL Rate Deep-Dive
# (forward-looking credit quality: direct/indirect, tier, structural break)
# ============================================================================
hdr("CHART 16: PLL rate deep-dive")

if ("pll_rate" %in% names(panel)) {
  pll_agg <- merge(agg_quarter(panel, "pll_rate"),
                   mac_spine[, .(yyyyqq, yoy_oil, pbrent)],
                   by = "yyyyqq", all.x = TRUE)
  pll_agg[, era := fifelse(yyyyqq < 201501L, "Pre-2015", "Post-2015")]

  p16a <- ggplot(pll_agg[!is.na(pll_rate)],
                 aes(x = cal_date, y = pll_rate)) +
    ep_rects() +
    geom_line(colour = COL_NEG, linewidth = 0.9) +
    geom_hline(yintercept = mean(pll_agg$pll_rate, na.rm = TRUE),
               linetype = "dashed", colour = "#888", linewidth = 0.4) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = number_format(accuracy = 0.01)) +
    labs(title    = "PLL Rate (% of Avg Loans) -- All CUs",
         subtitle = "Provisions are forward-looking: management signals expected credit losses",
         x = NULL, y = "PLL Rate (%)") +
    theme_pub()

  # Scatter: pre vs post 2015 slope
  p16b <- ggplot(pll_agg[!is.na(pll_rate) & !is.na(yoy_oil)],
                 aes(x = yoy_oil, y = pll_rate, colour = era)) +
    geom_point(alpha = 0.55, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
    scale_colour_manual(values = c("Pre-2015" = "#1a3a5c",
                                    "Post-2015" = "#b5470a"),
                        name = "Era") +
    labs(title    = "PLL Rate vs PBRENT YoY (Structural Break)",
         subtitle = "Post-2015: oil rise -> less provisioning (economic resilience)",
         x        = "PBRENT YoY %", y = "PLL Rate (%)") +
    theme_pub() + theme(legend.position = "right")

  # Direct vs indirect PLL rate
  p16c <- NULL
  if (!is.na(group_col_04) && "plot_group" %in% names(panel04)) {
    pll_grp <- panel04[!is.na(pll_rate) & !is.na(plot_group),
                        .(pll_rate = mean(pll_rate, na.rm = TRUE),
                          cal_date = first(cal_date)),
                        by = .(yyyyqq, plot_group)][order(yyyyqq)]
    p16c <- ggplot(pll_grp, aes(x = cal_date, y = pll_rate, colour = plot_group)) +
      ep_rects() +
      geom_line(linewidth = 0.85) +
      scale_colour_manual(values = grp_colors, name = NULL) +
      scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
      labs(title    = "PLL Rate: Direct vs Indirect CUs",
           subtitle = "GFC and COVID spikes diverge: oil-state CUs absorb energy-sector losses",
           x = NULL, y = "PLL Rate (%)") +
      theme_pub() + theme(legend.position = "bottom")
  }

  # By tier
  pll_tier <- if (!is.na(tier_col_06)) {
    tier_pll <- agg_group(panel, "pll_rate", tier_col_06)
    if (tier_col_06 == "assets_cat2") {
      tier_pll[, tier_code := cat2_to_t[as.character(get(tier_col_06))]]
      tc_pll <- "tier_code"
    } else { tc_pll <- tier_col_06 }
    tc_use_pll <- TIER_COLS[intersect(names(TIER_COLS), unique(tier_pll[[tc_pll]]))]
    ggplot(tier_pll[!is.na(pll_rate)],
           aes(x = cal_date, y = pll_rate, colour = get(tc_pll))) +
      ep_rects() +
      geom_line(linewidth = 0.65) +
      scale_colour_manual(values = tc_use_pll, labels = TIER_LABELS[names(tc_use_pll)],
                          name = "Tier") +
      scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
      labs(title    = "PLL Rate by Asset Tier",
           subtitle = "Large CUs (T6-T8) provision earlier due to CECL adoption; small CUs lag",
           x = NULL, y = "PLL Rate (%)") +
      theme_pub() + theme(legend.position = "right", legend.text = element_text(size = 7))
  } else NULL

  panels_16 <- Filter(Negate(is.null), list(p16a, p16b, p16c, pll_tier))
  if (length(panels_16) >= 2) {
    p16 <- wrap_plots(panels_16, ncol = 2) +
      plot_annotation(
        title    = "FIGURE 16 -- PLL Rate Deep-Dive: Forward-Looking Credit Quality",
        subtitle = paste("PLL rate = provision / avg loans | Forward-looking vs dq_rate (backward)",
                         "| CECL adoption 2020 (large CUs) / 2023 (FICUs)"),
        caption  = "Source: NCUA Form 5300 | pll_rate = pll / avg(lns_tot) | Bounded [-2%, +5%]",
        theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                         plot.subtitle = element_text(size = 9, colour = "#555"))
      )
    save_plot(p16, "16_pll_rate_deep_dive.png", w = 13, h = 10)
  }
}

# ============================================================================
# CHART 17 -- Net Worth Ratio (pcanetworth) Under Oil Stress
# ============================================================================
hdr("CHART 17: Net worth ratio under oil stress")

if ("pcanetworth" %in% names(panel)) {
  nw_agg <- merge(agg_quarter(panel, "pcanetworth"),
                  mac_spine[, .(yyyyqq, yoy_oil, pbrent)],
                  by = "yyyyqq", all.x = TRUE)

  p17a <- ggplot(nw_agg[!is.na(pcanetworth)], aes(x = cal_date, y = pcanetworth)) +
    ep_rects() +
    geom_line(colour = COL_DIRECT, linewidth = 0.9) +
    # NCUA PCA thresholds
    geom_hline(yintercept = 0.10, linetype = "dashed",
               colour = COL_POS, linewidth = 0.5) +
    geom_hline(yintercept = 0.07, linetype = "dashed",
               colour = COL_WARN, linewidth = 0.5) +
    geom_hline(yintercept = 0.06, linetype = "dashed",
               colour = COL_NEG, linewidth = 0.5) +
    annotate("text", x = max(nw_agg$cal_date, na.rm = TRUE),
             y = c(0.101, 0.071, 0.061),
             label = c("Well-Capitalised (10%)", "Adequately (7%)", "Undercap'd (6%)"),
             hjust = 1, size = 2.5,
             colour = c(COL_POS, COL_WARN, COL_NEG)) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = percent_format(scale = 100)) +
    labs(title    = "Net Worth Ratio (pcanetworth) vs NCUA PCA Thresholds",
         subtitle = "v1 finding: pcanetworth stable under oil shocks (AR1=0.98, highest persistence)",
         x = NULL, y = "Net Worth Ratio") +
    theme_pub()

  # Direct vs indirect NW ratio
  p17b <- NULL
  if (!is.na(group_col_04) && "plot_group" %in% names(panel04)) {
    nw_grp <- panel04[!is.na(pcanetworth) & !is.na(plot_group),
                       .(pcanetworth = mean(pcanetworth, na.rm = TRUE),
                         cal_date    = first(cal_date)),
                       by = .(yyyyqq, plot_group)][order(yyyyqq)]
    p17b <- ggplot(nw_grp, aes(x = cal_date, y = pcanetworth,
                                colour = plot_group)) +
      ep_rects() +
      geom_line(linewidth = 0.85) +
      geom_hline(yintercept = 0.07, linetype = "dashed",
                 colour = COL_WARN, linewidth = 0.4) +
      scale_colour_manual(values = grp_colors, name = NULL) +
      scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
      scale_y_continuous(labels = percent_format(scale = 100)) +
      labs(title    = "Net Worth Ratio: Direct vs Indirect CUs",
           subtitle = "Post-2015 oil-state CUs are net beneficiaries -- NW ratio should be higher",
           x = NULL, y = "Net Worth Ratio") +
      theme_pub() + theme(legend.position = "bottom")
  }

  p17 <- if (!is.null(p17b)) p17a / p17b else p17a
  p17 <- p17 + plot_annotation(
    title    = "FIGURE 17 -- Net Worth Ratio Under Oil Stress Cycles",
    subtitle = paste("NCUA PCA thresholds: well-capitalised >=10% | adequate >=7% | undercap'd <6%",
                     "| v1: NW ratio most persistent outcome (AR1=0.98)"),
    caption  = "Source: NCUA Form 5300 | pcanetworth = acct_998 / acct_010",
    theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                     plot.subtitle = element_text(size = 9, colour = "#555"))
  )
  save_plot(p17, "17_net_worth_ratio.png", w = 13, h = 9)
}

# ============================================================================
# CHART 18 -- NCUA PCA Traffic-Light Dashboard
# ============================================================================
hdr("CHART 18: NCUA PCA traffic-light dashboard")

if ("pcanetworth" %in% names(panel)) {
  # Classify each CU-quarter by PCA category
  panel_pca <- copy(panel[!is.na(pcanetworth)])
  panel_pca[, pca_class := fcase(
    pcanetworth >= 0.10, "Well-Capitalised (>=10%)",
    pcanetworth >= 0.07, "Adequately Capitalised (7-10%)",
    pcanetworth >= 0.06, "Under-Capitalised (6-7%)",
    pcanetworth >= 0.04, "Significantly Under-Cap'd (4-6%)",
    default             = "Critically Under-Cap'd (<4%)"
  )]
  panel_pca[, pca_class := factor(pca_class, levels = c(
    "Critically Under-Cap'd (<4%)",
    "Significantly Under-Cap'd (4-6%)",
    "Under-Capitalised (6-7%)",
    "Adequately Capitalised (7-10%)",
    "Well-Capitalised (>=10%)"
  ))]

  pca_qtr <- panel_pca[, .(pct = .N), by = .(yyyyqq, cal_date, pca_class)]
  pca_qtr[, total := sum(pct), by = yyyyqq]
  pca_qtr[, pct := pct / total * 100]

  PCA_COLS <- c(
    "Critically Under-Cap'd (<4%)"        = "#c0392b",
    "Significantly Under-Cap'd (4-6%)"    = "#e67e22",
    "Under-Capitalised (6-7%)"            = "#f1c40f",
    "Adequately Capitalised (7-10%)"      = "#2ecc71",
    "Well-Capitalised (>=10%)"            = "#1a3a5c"
  )

  p18 <- ggplot(pca_qtr[!is.na(cal_date)],
                aes(x = cal_date, y = pct, fill = pca_class)) +
    ep_rects() +
    geom_area(position = "stack", alpha = 0.85) +
    scale_fill_manual(values = PCA_COLS, name = "PCA Category") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = percent_format(scale = 1),
                       limits = c(0, 101)) +
    labs(title    = "FIGURE 18 -- NCUA PCA Traffic-Light Dashboard (2005-2025)",
         subtitle = paste("Share of CU-quarters by PCA capitalisation category",
                          "| Red/Orange = supervisory concern"),
         caption  = "Source: NCUA Form 5300 | pcanetworth = net worth / avg assets",
         x = NULL, y = "Share of CU-Quarters (%)") +
    theme_pub() +
    theme(legend.position = "right",
          legend.text     = element_text(size = 8))

  save_plot(p18, "18_pca_traffic_light.png", w = 12, h = 7)
}

# ============================================================================
# CHART 19 -- Iran War / $125 Oil Scenario (Moody's Analytics)
# ============================================================================
hdr("CHART 19: Iran war / $125 oil scenario (Moody's)")

# v1 confirmed findings -- ground projections in these anchors:
# - Direct effect: 0.14287 per 1pp oil YoY
# - Indirect effect: 0.00175 per 1pp oil YoY
# - AR(1) = 0.59 -> half-life 1.3Q -> LR multiplier 2.44x
# - +60.3pp YoY shock (from $78 baseline to $125 Moody's stressed)
IRAN_SHOCK_PP  <- 60.3
BASELINE_BRENT <- 78
STRESS_BRENT   <- 125
AR1            <- 0.59
LR_MULT        <- 1 / (1 - AR1)   # 2.44x

# Confirmed v1 impacts at +60.3pp
iran_impacts <- data.table(
  outcome  = c("insured_share_growth","netintmrg","dq_rate",
                "pll_rate","pcanetworth"),
  label    = c("Deposit Growth (YoY%)","Net Interest Margin",
                "Delinquency Rate","PLL Rate","Net Worth Ratio"),
  Q1_impact = c(-2.01,  +0.11, +0.11, +0.01, +0.13),
  LR_impact = c(-4.91,  +0.27, +0.27, +0.02, +0.32),
  direction = c("neg",  "pos", "neg", "neg", "pos")
)
cat("\n  Iran $125 scenario: confirmed v1 impacts at +60.3pp YoY\n")
print(iran_impacts, row.names = FALSE)

# Impulse-response fan chart (analytical, not stochastic)
iran_quarters <- 0:8
iran_list <- lapply(1:nrow(iran_impacts), function(i) {
  row <- iran_impacts[i]
  # Impulse at Q=1, then AR(1) decay
  impact_path <- row$Q1_impact * AR1^(iran_quarters - 1)
  impact_path[1] <- 0   # pre-shock = 0 deviation
  impact_path[2] <- row$Q1_impact

  # Historical baseline from panel
  hist_val <- if (row$outcome %in% names(agg))
    mean(agg[[row$outcome]], na.rm = TRUE) else 0

  data.table(
    qtr          = iran_quarters,
    deviation    = impact_path,
    cum_deviation = cumsum(impact_path),
    outcome      = row$outcome,
    label        = row$label,
    direction    = row$direction,
    hist_mean    = hist_val
  )
})
iran_dt <- rbindlist(iran_list)

p19 <- ggplot(iran_dt, aes(x = qtr, y = deviation,
                             fill = direction, colour = direction)) +
  geom_col(width = 0.7, alpha = 0.75, show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(aes(y = cum_deviation), linewidth = 0.8,
            linetype = "dashed", show.legend = FALSE) +
  geom_point(aes(y = cum_deviation), size = 1.8, show.legend = FALSE) +
  scale_fill_manual(values   = c("neg" = COL_NEG,  "pos" = COL_POS)) +
  scale_colour_manual(values = c("neg" = COL_NEG,  "pos" = COL_POS)) +
  scale_x_continuous(breaks = 0:8,
                     labels = c("Pre", paste0("Q+", 1:8))) +
  facet_wrap(~label, scales = "free_y", ncol = 2) +
  labs(title    = "FIGURE 19 -- Iran War / $125 Oil Scenario: Projected CU Impact",
       subtitle = sprintf(paste("Moody's Analytics $125/bbl stress = +%.1fpp YoY from $%.0f baseline",
                                "| v1 anchors: direct=0.14287 | AR(1)=0.59 | LR-mult=%.2fx",
                                "| Bars = quarterly impulse | Dashed = cumulative"),
                          IRAN_SHOCK_PP, BASELINE_BRENT, LR_MULT),
       caption  = paste("v1 confirmed findings (3-method validation: VARX, Baron-Kenny, XGBoost SHAP).",
                        "Geopolitical shocks = funding stability risk, NOT credit quality risk."),
       x = "Quarter Relative to Shock", y = "Deviation from Pre-Shock Level") +
  theme_pub() +
  theme(strip.text = element_text(size = 8.5, face = "bold"))

save_plot(p19, "19_iran_war_scenario.png", w = 13, h = 10)

# ============================================================================
# CHART 20 -- Bartik IV First-Stage Diagnostic
# ============================================================================
hdr("CHART 20: Bartik IV first-stage diagnostic")

if (all(c("oil_bartik_iv","mining_emp_share") %in% names(panel)) ||
    "oil_bartik_iv" %in% names(panel)) {

  iv_cols <- intersect(c("oil_bartik_iv","oil_exposure_cont","mining_emp_share"),
                       names(panel))
  if (length(iv_cols) >= 1 && "macro_base_yoy_oil" %in% names(panel)) {
    iv_agg <- merge(
      panel[, c("yyyyqq","cal_date", iv_cols), with = FALSE][
        , lapply(.SD, mean, na.rm = TRUE),
          by = .(yyyyqq, cal_date), .SDcols = iv_cols],
      mac_spine[, .(yyyyqq, yoy_oil)],
      by = "yyyyqq", all.x = TRUE
    )

    iv_plots <- lapply(iv_cols, function(v) {
      d <- iv_agg[!is.na(get(v)) & !is.na(yoy_oil)]
      if (nrow(d) < 10) return(NULL)
      r2 <- round(cor(d[[v]], d$yoy_oil, use = "complete.obs")^2, 3)
      ggplot(d, aes(x = get(v), y = yoy_oil)) +
        geom_point(colour = COL_DIRECT, alpha = 0.6, size = 1.5) +
        geom_smooth(method = "lm", se = TRUE, colour = COL_OIL,
                    linewidth = 0.9) +
        labs(title    = sprintf("%s vs PBRENT YoY (R^2 = %.3f)", v, r2),
             subtitle = sprintf("First-stage relationship | Partial R^2 = %.3f", r2),
             x        = v, y = "PBRENT YoY %") +
        theme_pub()
    })
    iv_plots <- Filter(Negate(is.null), iv_plots)

    if (length(iv_plots) >= 1) {
      p20 <- wrap_plots(iv_plots, ncol = min(2, length(iv_plots))) +
        plot_annotation(
          title    = "FIGURE 20 -- Bartik IV First-Stage Diagnostic",
          subtitle = paste("IV relevance: higher R^2 = stronger first stage.",
                           "Goldsmith-Pinkham et al. (2020) require F-stat > 10."),
          caption  = "Source: NCUA Form 5300; BLS QCEW NAICS 211 | Full F-stat in Script 04a",
          theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                           plot.subtitle = element_text(size = 9, colour = "#555"))
        )
      save_plot(p20, "20_bartik_iv_first_stage.png", w = 12, h = 7)
    }
  }
}

# ============================================================================
# CHART 21 -- Deposit Migration: Certificate vs Demand Shift Under Rate Regimes
# ============================================================================
hdr("CHART 21: Deposit migration under rate regimes")

if (all(c("cert_share","insured_share_growth") %in% names(panel)) &&
    "fomc_regime" %in% names(mac_spine)) {

  dep_regime <- merge(
    agg_quarter(panel, c("cert_share","insured_share_growth",
                         "dep_growth_yoy","costfds")),
    mac_spine[, .(yyyyqq, yoy_oil, fomc_regime)],
    by = "yyyyqq", all.x = TRUE
  )
  dep_regime[, regime_label := factor(
    fcase(fomc_regime ==  1L, "Hiking",
          fomc_regime == -1L, "Cutting",
          default            = "Hold"),
    levels = c("Cutting","Hold","Hiking")
  )]
  dep_regime[, era := fifelse(yyyyqq < 201501L, "Pre-2015", "Post-2015")]

  REGIME_COLS <- c("Hiking"  = COL_NEG,
                   "Hold"    = COL_WARN,
                   "Cutting" = COL_POS)

  # Scatter: cert_share vs yoy_oil by regime
  p21a <- ggplot(dep_regime[!is.na(cert_share) & !is.na(yoy_oil)],
                 aes(x = yoy_oil, y = cert_share,
                     colour = regime_label)) +
    geom_point(alpha = 0.55, size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    scale_colour_manual(values = REGIME_COLS, name = "FOMC Regime") +
    scale_y_continuous(labels = percent_format(scale = 100)) +
    labs(title    = "Certificate Share vs Oil YoY (by FOMC Regime)",
         subtitle = "Hiking + oil rise: dual deposit migration pressure on CUs",
         x        = "PBRENT YoY %", y = "Certificate Share") +
    theme_pub() + theme(legend.position = "right")

  # Boxplot: deposit growth by regime
  p21b <- if ("dep_growth_yoy" %in% names(dep_regime)) {
    ggplot(dep_regime[!is.na(dep_growth_yoy) & !is.na(regime_label)],
           aes(x = regime_label, y = dep_growth_yoy, fill = regime_label)) +
      geom_boxplot(width = 0.6, outlier.size = 0.5, alpha = 0.75) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "#888") +
      scale_fill_manual(values = REGIME_COLS, guide = "none") +
      scale_y_continuous(labels = number_format(accuracy = 0.1)) +
      labs(title    = "Deposit Growth Distribution by FOMC Regime",
           subtitle = "Cutting regime -> deposit inflows | Hiking -> outflows/migration",
           x        = "FOMC Regime", y = "Deposit Growth YoY (%)") +
      theme_pub()
  } else NULL

  # Time series: cert share with regime background
  p21c <- ggplot(dep_regime[!is.na(cert_share)],
                 aes(x = cal_date, y = cert_share, colour = regime_label)) +
    geom_line(linewidth = 0.85) +
    scale_colour_manual(values = REGIME_COLS, name = "FOMC Regime") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = percent_format(scale = 100)) +
    labs(title    = "Certificate Share Over Time by FOMC Regime Colouring",
         subtitle = "Rate hikes (red) coincide with certificate surges",
         x = NULL, y = "Certificate Share") +
    theme_pub() + theme(legend.position = "bottom")

  panels_21 <- Filter(Negate(is.null), list(p21a, p21b, p21c))
  if (length(panels_21) >= 2) {
    p21 <- wrap_plots(panels_21, ncol = 2) +
      plot_annotation(
        title    = "FIGURE 21 -- Deposit Migration: Certificate vs Demand under Rate Regimes",
        subtitle = paste("FOMC hiking + oil rise = dual pressure: volumes shift to certificates",
                         "AND grow (deposit migration + oil-state income inflows)"),
        caption  = "Source: NCUA Form 5300; FRB CCAR 2026 Baseline | fomc_regime from derived macro vars",
        theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                         plot.subtitle = element_text(size = 9, colour = "#555"))
      )
    save_plot(p21, "21_deposit_migration.png", w = 13, h = 10)
  }
}

# ============================================================================
# CHART 22 -- 8-Tier assets_cat2 Response Dashboard
# ============================================================================
hdr("CHART 22: 8-tier assets_cat2 response dashboard")

if (!is.na(tier_col_06) && nrow(tier_agg) > 0) {
  # Long format: all key outcomes by tier
  dash_outcomes <- intersect(c("dq_rate","pll_rate","netintmrg","costfds",
                                "insured_share_growth","pcanetworth"),
                              names(tier_agg))

  if (tier_col_06 == "assets_cat2" && "tier_code" %in% names(tier_agg)) {
    tier_agg_dash <- tier_agg
  } else {
    tier_agg_dash <- tier_agg
    tier_agg_dash[, tier_code := get(tier_col_06)]
  }

  # Compute episode means by tier for heatmap
  tier_ep_list <- lapply(1:nrow(ep_def), function(i) {
    ep  <- ep_def[i]
    sub <- tier_agg_dash[yyyyqq >= ep$yyyyqq_from & yyyyqq <= ep$yyyyqq_to]
    if (nrow(sub) == 0) return(NULL)
    means <- sub[, lapply(.SD, function(x) mean(x, na.rm = TRUE)),
                  by = tier_code, .SDcols = intersect(dash_outcomes, names(sub))]
    means[, episode := ep$episode]
    means
  })
  tier_ep_dt <- rbindlist(Filter(Negate(is.null), tier_ep_list), fill = TRUE)

  if (nrow(tier_ep_dt) > 0 && length(dash_outcomes) >= 1) {
    tier_ep_long <- melt(tier_ep_dt,
                          id.vars = c("tier_code","episode"),
                          variable.name = "outcome", value.name = "val")
    tier_ep_long[, out_label := out_labels[as.character(outcome)]]
    tier_ep_long[is.na(out_label), out_label := as.character(outcome)]
    tier_ep_long[, tier_label := TIER_LABELS[as.character(tier_code)]]
    tier_ep_long[is.na(tier_label), tier_label := as.character(tier_code)]

    # Normalise within outcome
    tier_ep_long[, norm_val := {
      mn <- min(val, na.rm = TRUE); mx <- max(val, na.rm = TRUE)
      if (is.finite(mn) && is.finite(mx) && mx > mn)
        (val - mn) / (mx - mn) else rep(0.5, .N)
    }, by = outcome]
    tier_ep_long[outcome %in% stress_outcomes, norm_val := 1 - norm_val]

    p22 <- ggplot(tier_ep_long[!is.na(tier_label) & !is.na(out_label)],
                  aes(x = episode, y = tier_label, fill = norm_val)) +
      geom_tile(colour = "white", linewidth = 0.4) +
      geom_text(aes(label = round(val, 2)), size = 2.2, colour = "#1a1a1a",
                na.rm = TRUE) +
      scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#d73027",
                           midpoint = 0.5, na.value = "#eeeeee",
                           name = "Relative\n(Blue=good\nRed=stress)") +
      scale_x_discrete(position = "top") +
      facet_wrap(~out_label, ncol = 2) +
      labs(title    = "FIGURE 22 -- 8-Tier Asset Response Dashboard: Episode x Tier Heatmap",
           subtitle = paste("Rows = 8 asset tiers (T1 <$10M to T8 >$10B) |",
                            "Columns = oil price episodes | Values = episode mean"),
           caption  = "Source: NCUA Form 5300 | assets_cat2 from OCE_combined Stata file",
           x = NULL, y = NULL) +
      theme_pub() +
      theme(axis.text.x  = element_text(size = 7, angle = 0, face = "bold"),
            axis.text.y  = element_text(size = 7.5),
            strip.text   = element_text(size = 8, face = "bold"),
            panel.grid   = element_blank(),
            legend.position = "right")

    save_plot(p22, "22_asset_tier_dashboard.png", w = 16, h = 12)
  }
}

# ============================================================================
# CHARTS 11-15 -- CCAR Severely Adverse Scenario (conditional)
# ============================================================================
if (has_severe) {

  hdr("CHARTS 11-15: CCAR Severely Adverse Scenario")

  macro_base_dt   <- macro
  proj_start      <- mac_spine[, max(cal_date, na.rm = TRUE)] -
                     365 * 3   # last 3 years treated as projection window

  out_labels_full <- out_labels   # alias

  # CHART 11 -- PBRENT paths
  hdr("CHART 11: PBRENT scenario paths")
  mac_sev_spine <- unique(macro_severe[, .(
    cal_date, yyyyqq,
    pbrent_sev = if ("macro_severe_pbrent" %in% names(macro_severe))
                   macro_severe_pbrent else NA_real_
  )])

  sev_compare <- merge(
    mac_spine[, .(cal_date, yyyyqq, pbrent_base = pbrent)],
    mac_sev_spine, by = c("cal_date","yyyyqq"), all = TRUE
  )
  sev_long <- melt(sev_compare, id.vars = c("cal_date","yyyyqq"),
                   measure.vars = c("pbrent_base","pbrent_sev"),
                   variable.name = "scenario", value.name = "pbrent")
  sev_long[, scenario := fifelse(scenario == "pbrent_base",
                                  "Baseline", "Severely Adverse")]
  sev_long <- sev_long[!is.na(pbrent)]

  p11 <- ggplot(sev_long, aes(x = cal_date, y = pbrent,
                               colour = scenario, linetype = scenario)) +
    geom_rect(xmin = proj_start, xmax = as.Date("2030-01-01"),
              ymin = -Inf, ymax = Inf,
              fill = "#fafafa", alpha = 0.05, inherit.aes = FALSE) +
    annotate("text", x = proj_start + 60, y = Inf,
             label = "Projection period", vjust = 1.3, size = 2.8,
             colour = "#888") +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = SCEN_COLS, name = "Scenario") +
    scale_linetype_manual(values = SCEN_LT, name = "Scenario") +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
    scale_y_continuous(labels = dollar_format(prefix = "$", suffix = "/bbl")) +
    labs(title    = "FIGURE 11 -- CCAR 2026: PBRENT Baseline vs Severely Adverse",
         subtitle = "Severely adverse tests CU resilience under deep oil price stress",
         caption  = "Source: FRB CCAR 2026 Baseline & Severely Adverse",
         x = NULL, y = "PBRENT ($/bbl)") +
    theme_pub()
  save_plot(p11, "11_pbrent_scenarios.png", w = 11, h = 6)

  # CHART 12 -- CU outcomes under both scenarios
  hdr("CHART 12: CU outcomes under CCAR scenarios")
  if (!is.null(panel_severe)) {
    add_date(panel_severe)
    agg_sev <- merge(
      agg_quarter(panel_severe,
                  intersect(cu_outcomes, names(panel_severe))),
      unique(macro_severe[, .(yyyyqq,
               pbrent_sev = if ("macro_severe_pbrent" %in% names(macro_severe))
                              macro_severe_pbrent else NA_real_)]),
      by = "yyyyqq", all.x = TRUE
    )

    make_scen_panel <- function(v, lab) {
      if (!v %in% names(agg) || !v %in% names(agg_sev)) return(NULL)
      base_d <- agg[!is.na(get(v)), .(cal_date, val = get(v),
                                        scenario = "Baseline")]
      sev_d  <- agg_sev[!is.na(get(v)), .(cal_date, val = get(v),
                                             scenario = "Severely Adverse")]
      both   <- rbindlist(list(base_d, sev_d))
      ggplot(both, aes(x = cal_date, y = val, colour = scenario,
                       linetype = scenario)) +
        geom_line(linewidth = 0.8) +
        scale_colour_manual(values = c("Baseline"         = COL_INDIR,
                                        "Severely Adverse" = COL_NEG),
                            name = NULL) +
        scale_linetype_manual(values = c("Baseline"         = "dashed",
                                          "Severely Adverse" = "solid"),
                              name = NULL) +
        scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
        labs(title = lab, x = NULL, y = lab) +
        theme_pub() + theme(legend.position = "bottom")
    }

    scen_panels <- Filter(Negate(is.null), list(
      make_scen_panel("dq_rate",             "Delinquency Rate (%)"),
      make_scen_panel("netintmrg",           "Net Interest Margin (%)"),
      make_scen_panel("insured_share_growth","Insured Share Growth (YoY%)"),
      make_scen_panel("costfds",             "Cost of Funds (%)")
    ))

    if (length(scen_panels) >= 2) {
      p12 <- wrap_plots(scen_panels, ncol = 2, guides = "collect") &
        theme(legend.position = "bottom")
      p12 <- p12 + plot_annotation(
        title    = "FIGURE 12 -- CU Outcomes: Baseline vs CCAR Severely Adverse",
        caption  = "Source: NCUA Form 5300; FRB CCAR 2026",
        theme    = theme(plot.title = element_text(face = "bold", size = 12))
      )
      save_plot(p12, "12_cu_outcomes_scenarios.png", w = 13, h = 10)
    }
  }

  msg("  Scenario charts 11-12 complete. Charts 13-15 from previous version retained.")

} else {
  msg("  Severely adverse data not available -- charts 11-15 skipped")
}

# ============================================================================
# CLOSING SUMMARY
# ============================================================================
cat("\n", SEP, "\n", sep = "")
cat("  SCRIPT 02 COMPLETE -- FIGURES SAVED TO Figures/\n")
cat(sprintf("  Finished: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

saved_charts <- CHART_LOG[CHART_LOG$status == "SAVED", ]
error_charts <- CHART_LOG[CHART_LOG$status == "ERROR", ]

cat(sprintf("  Charts saved  : %d\n", nrow(saved_charts)))
cat(sprintf("  Charts errored: %d\n\n", nrow(error_charts)))

cat("  FIGURE INDEX:\n")
cat(paste(rep("-", 68), collapse = ""), "\n")
cat(sprintf("  %-10s  %-52s  %s\n", "Chart", "Title", "File"))
cat(paste(rep("-", 68), collapse = ""), "\n")

chart_index <- data.frame(
  num   = c("01","02","03","04","05","06","07","08","09","10",
             "NEW-A","NEW-B","NEW-C","NEW-D",
             "16","17","18","19","20","21","22"),
  title = c("Oil price history + Iran annotation",
             "PBRENT vs aggregate CU outcomes",
             "Cross-correlation: PBRENT leads outcomes 1-8Q",
             "Direct vs indirect: oil-state vs non-oil CUs",
             "Deposit channel: cert_share, LTS, CoF",
             "Asset tier response: 8-tier assets_cat2",
             "Structural break: pre vs post 2015Q1",
             "Spillover: non-oil CU by tercile",
             "Episode heatmap: oil cycle x outcome",
             "Missingness overview",
             "Membership growth vs oil + direct/indirect",
             "Membership growth structural break",
             "Membership growth by 8-tier",
             "Membership growth lead/lag cross-corr",
             "PLL rate deep-dive: direct/tier/break",
             "Net worth ratio under oil stress cycles",
             "NCUA PCA traffic-light dashboard",
             "Iran war $125 Moody's scenario",
             "Bartik IV first-stage diagnostic",
             "Deposit migration under rate regimes",
             "8-tier assets_cat2 response dashboard"),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(chart_index))) {
  cat(sprintf("  %-10s  %-52s\n",
              chart_index$num[i], chart_index$title[i]))
}

if (nrow(error_charts) > 0) {
  cat("\n  ERROR CHARTS (investigate):\n")
  for (i in seq_len(nrow(error_charts))) {
    cat(sprintf("    %-10s  %s\n", error_charts$chart[i], error_charts$file[i]))
  }
}

cat("\n")
figs <- sort(list.files("Figures", pattern = "\\.png$", full.names = FALSE))
cat(sprintf("  PNG files in Figures/ (%d total):\n", length(figs)))
for (f in figs) cat(sprintf("    %s\n", f))
cat(SEP, "\n", sep = "")
