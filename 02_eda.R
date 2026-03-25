# =============================================================================
# OIL PRICE SHOCK × CREDIT UNION RESEARCH
# Script 02 — Exploratory Data Analysis
# =============================================================================
# Input  : Data/panel_base.rds   (from 01_data_prep_v2.R)
#          Data/macro_base.rds
# Output : Figures/                directory of PNG charts
#
# Charts produced:
#  01  Oil price history + cycle annotation          (PBRENT 2005-2025)
#  02  PBRENT vs aggregate CU outcomes — time series overlay
#  03  Cross-correlation: PBRENT leads CU outcomes by 1-8 quarters
#  04  Direct vs indirect: oil-state vs non-oil CU outcomes over time
#  05  Deposit channel: cert_share & loan_to_share vs fomc_regime
#  06  Asset tier response: dq_rate and netintmrg by tier across cycles
#  07  Structural break: pre/post 2015Q1 (shale era)
#  08  Spillover: non-oil CU response by spillover exposure tercile
#  09  Heatmap: oil price cycle episodes × CU outcome variables
#  10  Missingness overview
#  NEW-A  Member growth: YoY membership change vs oil price & by group
#  NEW-B  Member growth structural break: pre/post 2015 slope comparison
#  NEW-C  Member growth by asset tier: small vs large CU membership dynamics
#  NEW-D  Member growth lead/lag: cross-correlation with oil price
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

msg <- function(...) cat(sprintf(...), "\n")
hdr <- function(s)   cat("\n---", s, "---\n")

cat("=================================================================\n")
cat(" OIL SHOCK × CU  |  SCRIPT 02: EDA\n")
cat("=================================================================\n")

# =============================================================================
# CONFIG
# =============================================================================
dir.create("Figures", showWarnings = FALSE)

# ── Publication theme ─────────────────────────────────────────────────────────
theme_pub <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
  theme(
    text              = element_text(family = "sans", colour = "#1a1a1a"),
    plot.title        = element_text(size = base_size + 2, face = "bold",
                                     margin = margin(b = 4)),
    plot.subtitle     = element_text(size = base_size - 0.5, colour = "#555555",
                                     margin = margin(b = 8)),
    plot.caption      = element_text(size = base_size - 2, colour = "#888888",
                                     hjust = 0),
    axis.title        = element_text(size = base_size - 0.5, colour = "#333333"),
    axis.text         = element_text(size = base_size - 1.5, colour = "#444444"),
    panel.grid.major  = element_line(colour = "#e8e8e8", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(colour = "#cccccc", fill = NA,
                                     linewidth = 0.4),
    strip.text        = element_text(size = base_size - 1, face = "bold",
                                     colour = "#333333"),
    strip.background  = element_rect(fill = "#f5f5f5", colour = "#cccccc"),
    legend.position   = "bottom",
    legend.text       = element_text(size = base_size - 1.5),
    legend.title      = element_text(size = base_size - 1, face = "bold"),
    plot.margin       = margin(10, 12, 8, 10)
  )
}

# Colour palettes
COL_OIL    <- "#b5470a"          # Brent oil — burnt orange
COL_DIRECT <- "#1a3a5c"          # oil-state CUs — navy
COL_INDIR  <- "#2d7a4a"          # non-oil CUs — forest green
COL_SPILL  <- "#7a3080"          # spillover — purple
COL_NEG    <- "#c0392b"          # stress / negative — red
COL_POS    <- "#27ae60"          # positive — green
TIER_COLS  <- c("T1_under10M"  = "#4a90d9",
                "T2_10to100M"  = "#e67e22",
                "T3_100Mto1B"  = "#8e44ad",
                "T4_over1B"    = "#16a085")

# Key oil price episode shading bands
EPISODES <- data.frame(
  label = c("GFC","Shale\nBust","COVID\nCrash","Post-COVID\nSurge"),
  xmin  = as.Date(c("2008-07-01","2014-07-01","2020-01-01","2021-01-01")),
  xmax  = as.Date(c("2009-06-30","2016-06-30","2020-06-30","2022-06-30")),
  fill  = c("#fde8e8","#e8f0fd","#fde8e8","#e8fde8")
)

save_plot <- function(p, filename, w = 10, h = 6.5, dpi = 300) {
  path <- file.path("Figures", filename)
  ggsave(path, plot = p, width = w, height = h, dpi = dpi,
         bg = "white")
  msg("  Saved: %s", path)
}

# =============================================================================
# LOAD DATA
# =============================================================================
hdr("Loading data")

panel <- readRDS("Data/panel_base.rds")
macro <- readRDS("Data/macro_base.rds")
setDT(panel); setDT(macro)

# Calendar date for plotting
Q_MONTH <- c("1"=1L,"2"=4L,"3"=7L,"4"=10L)
panel[, cal_date := as.Date(paste(year,
                                   Q_MONTH[as.character(quarter)],
                                   "01", sep="-"))]
macro[, cal_date := as.Date(paste(year,
                                   Q_MONTH[as.character(quarter)],
                                   "01", sep="-"))]

# Macro spine — unique quarters
mac_spine <- unique(macro[, .(cal_date, yyyyqq,
                               pbrent   = macro_base_pbrent,
                               yoy_oil  = macro_base_yoy_oil,
                               lurc     = macro_base_lurc,
                               pcpi     = macro_base_pcpi,
                               rmtg     = macro_base_rmtg,
                               fedfunds = macro_base_rff,
                               fomc_regime = macro_base_fomc_regime)])[order(cal_date)]

msg("  panel: %s rows × %s cols | %s CUs | quarters %sQ%s–%sQ%s",
    format(nrow(panel),big.mark=","), ncol(panel),
    format(uniqueN(panel$join_number),big.mark=","),
    min(panel$year), panel[which.min(yyyyqq), quarter],
    max(panel$year), panel[which.max(yyyyqq), quarter])

# =============================================================================
# HELPER: aggregate CU outcomes to quarter level
# =============================================================================
agg_quarter <- function(dt, vars, by_vars = "yyyyqq") {
  dt[, c(list(cal_date = first(cal_date), year = first(year),
              quarter  = first(quarter)),
         lapply(.SD, function(x) mean(x, na.rm=TRUE))),
     by = by_vars, .SDcols = intersect(vars, names(dt))
  ][order(get(by_vars[1]))]
}

agg_group <- function(dt, vars, group_col, by_vars = c("yyyyqq", group_col)) {
  dt[!is.na(get(group_col)),
     c(list(cal_date = first(cal_date), year = first(year),
             quarter  = first(quarter)),
        lapply(.SD, function(x) mean(x, na.rm=TRUE))),
     by = by_vars, .SDcols = intersect(vars, names(dt))
  ][order(yyyyqq)]
}

# Available outcome vars — includes member_growth_yoy and pll_rate from 01_data_prep.R
cu_outcomes <- intersect(c("dq_rate","chg_tot_lns_ratio","netintmrg",
                             "pcanetworth","networth","roa","costfds",
                             "insured_share_growth","cert_share",
                             "loan_to_share","nim_spread",
                             "member_growth_yoy",
                             "pll_rate","pll_per_loan"),
                          names(panel))

msg("  CU outcomes available: %s", paste(cu_outcomes, collapse=", "))
if (!"member_growth_yoy" %in% cu_outcomes)
  msg("  NOTE: member_growth_yoy not found — Charts NEW-A/B/C/D will be skipped")
if (!"pll_rate" %in% cu_outcomes)
  msg("  NOTE: pll_rate not found — check pll and lns_tot are in call report")

# Episode rectangles helper — returns flat list for ggplot layer addition
ep_rects <- function(episodes = EPISODES) {
  rects <- mapply(function(xmn, xmx, fl) {
    annotate("rect", xmin=xmn, xmax=xmx,
             ymin=-Inf, ymax=Inf, fill=fl, alpha=0.35)
  }, episodes$xmin, episodes$xmax, episodes$fill, SIMPLIFY=FALSE)

  texts <- mapply(function(xmn, xmx, lb) {
    annotate("text", x=xmn + (xmx-xmn)/2,
             y=Inf, label=lb, vjust=1.3, size=2.5,
             colour="#888888", fontface="italic")
  }, episodes$xmin, episodes$xmax, episodes$label, SIMPLIFY=FALSE)

  # Flatten to single list so ggplot + list() works correctly
  c(rects, texts)
}

# =============================================================================
# CHART 01 — PBRENT Oil Price History + Cycle Annotation
# =============================================================================
hdr("Chart 01: Oil price history")

mac_hist <- mac_spine[!is.na(pbrent) & cal_date >= as.Date("2005-01-01")]

p_oil_level <- ggplot(mac_hist, aes(x = cal_date, y = pbrent)) +
  ep_rects() +
  geom_line(colour = COL_OIL, linewidth = 0.9) +
  geom_hline(yintercept = mean(mac_hist$pbrent, na.rm=TRUE),
             linetype = "dashed", colour = "#888888", linewidth = 0.4) +
  annotate("text", x = min(mac_hist$cal_date) + 60,
           y = mean(mac_hist$pbrent, na.rm=TRUE) + 4,
           label = "Period average", size = 2.8, colour = "#888888") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = dollar_format(prefix = "$", suffix = "/bbl")) +
  labs(title    = "Brent Crude Oil Price (PBRENT) — 2005Q1 to 2025Q4",
       subtitle = "Shaded: GFC 2008–09 | Shale bust 2014–16 | COVID 2020 | Post-COVID surge 2021–22",
       x = NULL, y = "$/barrel",
       caption  = "Source: FRB CCAR 2026 Baseline (macro_base_pbrent)") +
  theme_pub()

p_oil_yoy <- ggplot(mac_hist[!is.na(yoy_oil)], aes(x = cal_date, y = yoy_oil,
                                                      fill = yoy_oil >= 0)) +
  geom_col(width = 70, show.legend = FALSE) +
  scale_fill_manual(values = c("TRUE" = COL_POS, "FALSE" = COL_NEG)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  labs(title = "YoY % Change in Brent Oil Price",
       x = NULL, y = "YoY %",
       caption = "Source: FRB CCAR 2026 Baseline") +
  theme_pub()

p01 <- p_oil_level / p_oil_yoy +
  plot_annotation(
    title   = "FIGURE 01 — Brent Crude Oil Price: Level & Annual Change",
    theme   = theme(plot.title = element_text(face="bold", size=12))
  )

save_plot(p01, "01_oil_price_history.png", w=11, h=8)

# =============================================================================
# CHART 02 — PBRENT vs Aggregate CU Outcomes (time-series overlay)
# =============================================================================
hdr("Chart 02: PBRENT vs CU outcomes")

agg <- agg_quarter(panel, cu_outcomes)
agg <- merge(agg, mac_spine[, .(yyyyqq, pbrent, yoy_oil)],
             by="yyyyqq", all.x=TRUE)

make_dual_axis <- function(outcome_var, outcome_label, y_fmt = waiver()) {
  if (!outcome_var %in% names(agg)) return(NULL)

  d <- agg[!is.na(get(outcome_var)) & !is.na(pbrent)]
  ov_max <- max(abs(d[[outcome_var]]), na.rm=TRUE)
  pb_max <- max(abs(d$pbrent),         na.rm=TRUE)
  oil_scale <- if (is.finite(ov_max) && is.finite(pb_max) && pb_max > 0)
                 ov_max / pb_max else 1

  ggplot(d, aes(x = cal_date)) +
    ep_rects() +
    geom_line(aes(y = pbrent * oil_scale, colour = "PBRENT (scaled)"),
              linewidth = 0.6, linetype = "dashed") +
    geom_line(aes(y = get(outcome_var), colour = outcome_label),
              linewidth = 0.85) +
    scale_colour_manual(values = c("PBRENT (scaled)" = COL_OIL,
                                    setNames(COL_DIRECT, outcome_label)),
                        name = NULL) +
    scale_x_date(date_breaks = "3 years", date_labels = "%Y") +
    scale_y_continuous(labels = y_fmt) +
    labs(title   = outcome_label,
         x = NULL, y = outcome_label) +
    theme_pub() +
    theme(legend.position = "none")
}

panels_02 <- list(
  make_dual_axis("dq_rate",            "Delinquency Rate (%)",    percent_format(scale=1)),
  make_dual_axis("pll_rate",           "PLL Rate (% of Avg Loans)", number_format(accuracy=0.01)),
  make_dual_axis("netintmrg",          "Net Interest Margin (%)", number_format(accuracy=0.1)),
  make_dual_axis("costfds",            "Cost of Funds (%)",       number_format(accuracy=0.01)),
  make_dual_axis("insured_share_growth","Insured Share Growth (YoY%)", number_format(accuracy=0.1)),
  make_dual_axis("cert_share",         "Certificate Share of Deposits", percent_format(scale=1)),
  make_dual_axis("loan_to_share",      "Loan-to-Share Ratio",    number_format(accuracy=0.01)),
  make_dual_axis("member_growth_yoy",  "Membership Growth (YoY%)",number_format(accuracy=0.1))
)
panels_02 <- Filter(Negate(is.null), panels_02)

if (length(panels_02) >= 2) {
  p02 <- wrap_plots(panels_02, ncol = 2) +
    plot_annotation(
      title    = "FIGURE 02 — PBRENT vs Aggregate CU Outcomes (2005–2025)",
      subtitle = "PBRENT dashed orange line (scaled to outcome axis) | Shaded: key oil episodes",
      caption  = "Source: NCUA Form 5300 Call Report; FRB CCAR 2026 Baseline",
      theme    = theme(plot.title    = element_text(face="bold", size=12),
                       plot.subtitle = element_text(size=9, colour="#555"))
    )
  save_plot(p02, "02_pbrent_vs_cu_outcomes.png", w=13, h=10)
}

# =============================================================================
# CHART 03 — Cross-Correlation: PBRENT leads CU outcomes
# =============================================================================
hdr("Chart 03: Cross-correlations")

agg_cc <- merge(agg_quarter(panel, cu_outcomes),
                mac_spine[, .(yyyyqq, yoy_oil)],
                by = "yyyyqq", all.x = TRUE)

cc_results <- rbindlist(lapply(cu_outcomes, function(v) {
  if (!v %in% names(agg_cc)) return(NULL)
  x  <- agg_cc$yoy_oil
  y  <- agg_cc[[v]]
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < 20) return(NULL)
  lags <- -4:8  # negative = oil leads by that many quarters
  cors <- sapply(lags, function(k) {
    if (k >= 0) cor(x[ok][1:(sum(ok)-k)], y[ok][(k+1):sum(ok)],
                    use="complete.obs")
    else        cor(x[ok][(-k+1):sum(ok)], y[ok][1:(sum(ok)+k)],
                    use="complete.obs")
  })
  data.table(outcome=v, lag=lags, correlation=cors)
}))

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
cc_results[, outcome_label := out_labels[outcome]]
cc_results[is.na(outcome_label), outcome_label := outcome]

# ── Diagnostic: report which outcomes have valid correlations ─────────────────
cc_coverage <- cc_results[, .(
  n_valid   = sum(!is.na(correlation)),
  max_abs   = max(abs(correlation), na.rm=TRUE)
), by=.(outcome, outcome_label)][order(-max_abs)]

cat("\n  Cross-correlation coverage:\n")
print(cc_coverage, row.names=FALSE)

# ── Colour palette — must support all outcomes (Dark2 only has 8) ─────────────
CC_COLS <- c(
  "Cost of Funds (%)"               = "#1b9e77",
  "Delinquency Rate (%)"            = "#d95f02",
  "Insured Share Growth (YoY%)"     = "#7570b3",
  "Loan-to-Share Ratio"             = "#e7298a",
  "Net Interest Margin (%)"         = "#66a61e",
  "Membership Growth (YoY%)"        = "#e6ab02",
  "PLL Rate (% of Avg Loans)"       = "#a6761d",
  "Return on Assets (%)"            = "#666666",
  "Certificate Share of Deposits"   = "#1f78b4",
  "PLL per Loan ($)"                = "#b2df8a",
  "Net Charge-Off Ratio (%)"        = "#fb9a99",
  "Net Worth Ratio (% of Assets)"   = "#cab2d6",
  "Net Worth ($000s)"               = "#fdbf6f",
  "NIM Spread (Loan Yield - CoF)"   = "#b15928"
)

# Only plot outcomes that have at least some valid correlations
plot_outcomes_03 <- cc_coverage[n_valid > 0 & is.finite(max_abs), outcome_label]

p03 <- ggplot(cc_results[!is.na(correlation) &
                           outcome_label %in% plot_outcomes_03],
              aes(x = lag, y = correlation, colour = outcome_label)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#888") +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype="dashed",
             colour = COL_OIL) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.8) +
  annotate("text", x=0.1, y=Inf, label="\u2190 Oil leads  |  Oil lags \u2192",
           vjust=1.5, hjust=0, size=2.8, colour="#888") +
  scale_colour_manual(
    values = CC_COLS,
    name   = "CU Outcome",
    breaks = plot_outcomes_03   # legend order matches plot
  ) +
  scale_x_continuous(breaks=-4:8,
                     labels=c(paste0("Lag\n",4:1),"0",
                               paste0("Lead\n",1:8))) +
  labs(title    = "FIGURE 03 — Cross-Correlation: PBRENT YoY vs CU Outcomes",
       subtitle = "Positive lag = oil change precedes CU outcome by that many quarters",
       x = "Quarter lag (positive = PBRENT leads)",
       y = "Pearson Correlation",
       caption = "Source: NCUA Form 5300; FRB CCAR 2026 Baseline") +
  theme_pub() +
  theme(legend.position = "right",
        legend.key.height = unit(0.5,"cm"))

save_plot(p03, "03_cross_correlation.png", w=11, h=6.5)

# =============================================================================
# CHART 04 — Direct vs Indirect: Oil-State vs Non-Oil CU outcomes
# =============================================================================
hdr("Chart 04: Direct vs indirect effect")

# ── Step 1: Build a reliable 2-level grouping ─────────────────────────────────
# Primary: oil_exposure_bin (1 = oil state, 0 = non-oil) — always available
# This guarantees BOTH groups appear regardless of cu_group classification
panel04 <- copy(panel)

if ("oil_exposure_bin" %in% names(panel04)) {
  panel04[, plot_group := fifelse(
    !is.na(oil_exposure_bin) & oil_exposure_bin == 1L,
    "Oil-State CUs (Direct)",
    "Non-Oil CUs (Indirect/Spillover)"
  )]
  msg("  Using oil_exposure_bin for Chart 04 grouping")
} else if ("cu_group" %in% names(panel04)) {
  panel04[, plot_group := fcase(
    cu_group == "Direct",    "Oil-State CUs (Direct)",
    cu_group == "Indirect",  "Non-Oil CUs (Indirect/Spillover)",
    default                = "Non-Oil CUs (Indirect/Spillover)"
  )]
  msg("  Using cu_group for Chart 04 grouping")
} else {
  msg("  WARNING: no grouping variable found for Chart 04")
  panel04 <- NULL
}

if (!is.null(panel04)) {
  # Diagnostic
  cat("  Chart 04 group distribution:
")
  print(panel04[, .N, by=plot_group][order(plot_group)])

  # Emergency fallback: diagnose state column then apply correct matching
  n_direct <- panel04[plot_group == "Oil-State CUs (Direct)", .N]
  if (n_direct == 0) {
    msg("  WARNING: 0 Direct CU-qtrs — diagnosing state column...")

    # Find state column and inspect values
    sc_col <- intersect(c("reporting_state","state_code","state",
                           "reporting_state_full","state_full"), names(panel04))[1]

    if (!is.na(sc_col)) {
      sample_vals <- unique(panel04[[sc_col]])[1:10]
      msg("  State column: %s | Sample values: %s", sc_col,
          paste(head(sample_vals,10), collapse=", "))

      # Detect format and build oil-state indicator accordingly
      sv <- as.character(na.omit(unique(panel04[[sc_col]])))

      # FIPS state codes: numeric 1-56
      if (all(suppressWarnings(!is.na(as.numeric(sv[sv!=""]))) )) {
        msg("  Detected: numeric FIPS codes")
        # Oil-state FIPS: TX=48,ND=38,LA=22,AK=2,WY=56,OK=40,NM=35,CO=8,WV=54,PA=42
        OIL_FIPS <- as.character(c(48,38,22,2,56,40,35,8,54,42,30,20))
        panel04[, plot_group := fifelse(
          as.character(get(sc_col)) %in% OIL_FIPS,
          "Oil-State CUs (Direct)",
          "Non-Oil CUs (Indirect/Spillover)")]

      # 2-letter abbreviations
      } else if (all(nchar(sv[sv!=""]) <= 2)) {
        msg("  Detected: 2-letter state abbreviations")
        OIL_ABB <- c("TX","ND","LA","AK","WY","OK","NM","CO","WV","PA","MT","KS")
        panel04[, plot_group := fifelse(
          toupper(get(sc_col)) %in% OIL_ABB,
          "Oil-State CUs (Direct)",
          "Non-Oil CUs (Indirect/Spillover)")]

      # Full state names
      } else {
        msg("  Detected: full state names")
        OIL_NAMES <- c("Texas","North Dakota","Louisiana","Alaska","Wyoming",
                        "Oklahoma","New Mexico","Colorado","West Virginia",
                        "Pennsylvania","Montana","Kansas")
        panel04[, plot_group := fifelse(
          str_to_title(get(sc_col)) %in% OIL_NAMES,
          "Oil-State CUs (Direct)",
          "Non-Oil CUs (Indirect/Spillover)")]
      }

      cat("  Fallback group distribution:\n")
      print(panel04[, .N, by=plot_group])
    } else {
      msg("  ERROR: no state column found — cannot create Direct/Indirect split")
    }
  }


  GRP_COLS04 <- c("Oil-State CUs (Direct)"            = COL_DIRECT,
                   "Non-Oil CUs (Indirect/Spillover)"  = COL_INDIR)
  GRP_LT04   <- c("Oil-State CUs (Direct)"            = "solid",
                   "Non-Oil CUs (Indirect/Spillover)"  = "solid")

  # ── Step 2: Aggregate by group × quarter ───────────────────────────────────
  agg_grp04 <- panel04[!is.na(plot_group),
    c(list(cal_date = first(cal_date), year=first(year), quarter=first(quarter)),
      lapply(.SD, function(x) mean(x, na.rm=TRUE))),
    by = .(yyyyqq, plot_group),
    .SDcols = intersect(cu_outcomes, names(panel04))
  ][order(yyyyqq)]

  # Add PBRENT for reference
  agg_grp04 <- merge(agg_grp04,
                      mac_spine[, .(yyyyqq, pbrent, yoy_oil)],
                      by="yyyyqq", all.x=TRUE)

  # ── Step 3: Plot function ───────────────────────────────────────────────────
  make_grp04 <- function(v, lab, y_fmt=waiver()) {
    if (!v %in% names(agg_grp04)) return(NULL)
    d <- agg_grp04[!is.na(get(v)) & !is.na(plot_group)]
    if (uniqueN(d$plot_group) < 2) {
      msg("  NOTE: only 1 group for %s — both groups may overlap", v)
    }
    if (nrow(d) == 0) return(NULL)

    ggplot(d, aes(x=cal_date, y=get(v),
                   colour=plot_group, linetype=plot_group)) +
      ep_rects() +
      geom_line(linewidth=0.85) +
      scale_colour_manual(values=GRP_COLS04, name=NULL,
                           drop=FALSE) +
      scale_linetype_manual(values=GRP_LT04, name=NULL,
                             drop=FALSE) +
      scale_x_date(date_breaks="3 years", date_labels="%Y") +
      scale_y_continuous(labels=y_fmt) +
      labs(title=lab, x=NULL, y=lab) +
      theme_pub() +
      theme(legend.position="bottom",
            legend.text=element_text(size=7.5))
  }

  p04_panels <- list(
    make_grp04("dq_rate",             "Delinquency Rate (%)"),
    make_grp04("pll_rate",            "PLL Rate (% of Avg Loans)"),
    make_grp04("netintmrg",           "Net Interest Margin (%)"),
    make_grp04("insured_share_growth","Insured Share Growth (YoY%)"),
    make_grp04("member_growth_yoy",   "Membership Growth (YoY%)"),
    make_grp04("costfds",             "Cost of Funds (%)"),
    make_grp04("cert_share",          "Certificate Share",
                percent_format(scale=100)),
    make_grp04("loan_to_share",       "Loan-to-Share Ratio")
  )
  p04_panels <- Filter(Negate(is.null), p04_panels)

  if (length(p04_panels) >= 2) {
    p04 <- wrap_plots(p04_panels, ncol=2,
                       guides="collect") &
      theme(legend.position="bottom")

    p04 <- p04 +
      plot_annotation(
        title    = "FIGURE 04 — Direct vs Indirect: Oil-State vs Non-Oil CUs (2005–2025)",
        subtitle = "Oil-State = mining emp share ≥ 2% (BLS QCEW)  |  Non-Oil = remaining CUs",
        caption  = "Source: NCUA Form 5300 Call Report; BLS QCEW oil exposure classification",
        theme    = theme(plot.title    = element_text(face="bold", size=12),
                         plot.subtitle = element_text(size=9, colour="#555"))
      )
    save_plot(p04, "04_direct_vs_indirect.png", w=13, h=10)
  }
}

# =============================================================================
# CHART 05 — Deposit Channel: cert_share & loan_to_share vs FOMC regime
# =============================================================================
hdr("Chart 05: Deposit channel")

if (all(c("cert_share","loan_to_share","fomc_x_brent") %in% names(panel))) {

  dep_agg <- agg_quarter(panel, c("cert_share","loan_to_share",
                                   "insured_share_growth","costfds"))
  dep_agg <- merge(dep_agg,
                   mac_spine[, .(yyyyqq, pbrent, fedfunds, fomc_regime)],
                   by="yyyyqq", all.x=TRUE)

  # Regime shading
  dep_agg[, regime_label := fcase(
    fomc_regime ==  1L, "Hiking",
    fomc_regime == -1L, "Cutting",
    default           = "Hold"
  )]

  p_cert <- ggplot(dep_agg[!is.na(cert_share)],
                   aes(x=cal_date, y=cert_share)) +
    ep_rects() +
    geom_line(colour=COL_DIRECT, linewidth=0.85) +
    geom_line(aes(y=fedfunds/100, colour="Fed Funds Rate (RHS)"),
              linetype="dashed", linewidth=0.6) +
    scale_colour_manual(values=c("Fed Funds Rate (RHS)"=COL_OIL), name=NULL) +
    scale_y_continuous(labels=percent_format(scale=100),
                       limits=c(0, NA),
                       sec.axis=sec_axis(~.*100, name="Fed Funds Rate (%)")) +
    scale_x_date(date_breaks="2 years", date_labels="%Y") +
    labs(title="Certificate Share of Total Deposits vs Fed Funds Rate",
         subtitle="Rising rates → members shift from demand deposits to certificates (◆ deposit migration)",
         x=NULL, y="Certificate Share") +
    theme_pub()

  p_lts <- ggplot(dep_agg[!is.na(loan_to_share)],
                  aes(x=cal_date, y=loan_to_share)) +
    ep_rects() +
    geom_line(colour=COL_INDIR, linewidth=0.85) +
    geom_hline(yintercept=0.8, linetype="dashed",
               colour=COL_NEG, linewidth=0.4) +
    annotate("text", x=max(dep_agg$cal_date, na.rm=TRUE),
             y=0.81, label="80% liquidity threshold",
             hjust=1, size=2.8, colour=COL_NEG) +
    scale_x_date(date_breaks="2 years", date_labels="%Y") +
    scale_y_continuous(labels=number_format(accuracy=0.01)) +
    labs(title="Loan-to-Share Ratio",
         subtitle="Oil-state CUs: PBRENT↑ → deposit surge → ratio compression",
         x=NULL, y="Loan / Total Shares") +
    theme_pub()

  p_cof <- ggplot(dep_agg[!is.na(costfds)],
                  aes(x=cal_date, y=costfds)) +
    ep_rects() +
    geom_line(colour=COL_SPILL, linewidth=0.85) +
    scale_x_date(date_breaks="2 years", date_labels="%Y") +
    scale_y_continuous(labels=number_format(accuracy=0.01)) +
    labs(title="Cost of Funds (costfds, %)",
         subtitle="CoF squeeze: certificate surge raises funding costs even as deposit volumes grow",
         x=NULL, y="Cost of Funds (%)") +
    theme_pub()

  p_isg <- ggplot(dep_agg[!is.na(insured_share_growth)],
                  aes(x=cal_date, y=insured_share_growth,
                      fill=insured_share_growth>=0)) +
    geom_col(width=70, show.legend=FALSE) +
    scale_fill_manual(values=c("TRUE"=COL_POS,"FALSE"=COL_NEG)) +
    scale_x_date(date_breaks="2 years", date_labels="%Y") +
    scale_y_continuous(labels=number_format(accuracy=0.1),
                       limits=c(
                         min(dep_agg$insured_share_growth, na.rm=TRUE) * 1.05,
                         max(dep_agg$insured_share_growth, na.rm=TRUE) * 1.05
                       )) +
    labs(title="Insured Share Growth (YoY %, winsorised 1-99th pctile)",
         subtitle="Oil-state income windfall → deposit inflows; inflation erosion → drawdown",
         x=NULL, y="YoY Growth (%)") +
    theme_pub()

  p05 <- (p_cert + p_lts) / (p_cof + p_isg) +
    plot_annotation(
      title   = "FIGURE 05 — Deposit Channel Dynamics",
      caption = "Source: NCUA Form 5300 Call Report; FRB CCAR 2026 Baseline",
      theme   = theme(plot.title=element_text(face="bold",size=12))
    )
  save_plot(p05, "05_deposit_channel.png", w=13, h=10)
}

# =============================================================================
# CHART 06 — Asset Tier Response: dq_rate & netintmrg by tier
# =============================================================================
hdr("Chart 06: Asset tier response")

if ("asset_tier" %in% names(panel)) {
  tier_agg <- agg_group(panel,
                         intersect(c("dq_rate","netintmrg","costfds",
                                     "cert_share","insured_share_growth",
                                     "member_growth_yoy","pll_rate"),
                                   names(panel)),
                         "asset_tier")
  tier_agg <- merge(tier_agg,
                    mac_spine[,.(yyyyqq,pbrent,yoy_oil)],
                    by="yyyyqq", all.x=TRUE)

  make_tier_plot <- function(v, lab, y_fmt=waiver()) {
    if (!v %in% names(tier_agg)) return(NULL)
    ggplot(tier_agg[!is.na(get(v))],
           aes(x=cal_date, y=get(v), colour=asset_tier)) +
      ep_rects() +
      geom_line(linewidth=0.75) +
      scale_colour_manual(values=TIER_COLS, name="Asset Tier") +
      scale_x_date(date_breaks="3 years", date_labels="%Y") +
      scale_y_continuous(labels=y_fmt) +
      labs(title=lab, x=NULL, y=lab) +
      theme_pub() +
      theme(legend.position="right")
  }

  t_panels <- list(
    make_tier_plot("dq_rate",           "Delinquency Rate (%)"),
    make_tier_plot("pll_rate",          "PLL Rate (% of Avg Loans)",
                   number_format(accuracy=0.01)),
    make_tier_plot("netintmrg",         "Net Interest Margin (%)"),
    make_tier_plot("costfds",           "Cost of Funds (%)"),
    make_tier_plot("cert_share",        "Certificate Share",
                   percent_format(scale=100, accuracy=0.1)),
    make_tier_plot("member_growth_yoy", "Membership Growth (YoY%)",
                   number_format(accuracy=0.1))
  )
  t_panels <- Filter(Negate(is.null), t_panels)

  p06 <- wrap_plots(t_panels, ncol=2) +
    plot_annotation(
      title    = "FIGURE 06 — CU Outcomes by Asset Tier (2005–2025)",
      subtitle = "T1 < $10M | T2 $10-100M | T3 $100M-$1B | T4 > $1B",
      caption  = "Source: NCUA Form 5300 Call Report",
      theme    = theme(plot.title=element_text(face="bold",size=12),
                       plot.subtitle=element_text(size=9,colour="#555"))
    )
  save_plot(p06, "06_asset_tier_response.png", w=13, h=10)
}

# =============================================================================
# CHART 07 — Structural Break: pre vs post 2015Q1 (shale era)
# =============================================================================
hdr("Chart 07: Structural break 2015Q1")

agg_sb <- merge(agg_quarter(panel, cu_outcomes),
                mac_spine[,.(yyyyqq,yoy_oil)],
                by="yyyyqq", all.x=TRUE)
agg_sb[, era := fifelse(yyyyqq < 201501L, "Pre-Shale\n(2005–2014)",
                                           "Post-Shale\n(2015–2025)")]

sb_long <- melt(agg_sb[!is.na(yoy_oil)],
                id.vars     = c("yyyyqq","cal_date","era","yoy_oil"),
                measure.vars= intersect(cu_outcomes, names(agg_sb)),
                variable.name="outcome", value.name="value")
sb_long[, outcome_label := out_labels[as.character(outcome)]]
sb_long[is.na(outcome_label), outcome_label := as.character(outcome)]

# Scatter: yoy_oil vs each outcome, coloured by era
plot_sb_outcomes <- intersect(c("dq_rate","pll_rate","netintmrg","costfds",
                                 "insured_share_growth",
                                 "member_growth_yoy"), cu_outcomes)

sb_plots <- lapply(plot_sb_outcomes, function(v) {
  d <- sb_long[outcome==v & !is.na(value) & !is.na(yoy_oil)]
  if (nrow(d) < 10) return(NULL)
  ggplot(d, aes(x=yoy_oil, y=value, colour=era)) +
    geom_point(alpha=0.5, size=1.2) +
    geom_smooth(method="lm", se=TRUE, linewidth=0.8) +
    scale_colour_manual(values=c("Pre-Shale\n(2005–2014)"="#1a3a5c",
                                  "Post-Shale\n(2015–2025)"="#b5470a"),
                        name="Era") +
    labs(title  = unique(d$outcome_label),
         x      = "PBRENT YoY %",
         y      = unique(d$outcome_label)) +
    theme_pub() +
    theme(legend.position="right")
})
sb_plots <- Filter(Negate(is.null), sb_plots)

if (length(sb_plots) >= 2) {
  p07 <- wrap_plots(sb_plots, ncol=2) +
    plot_annotation(
      title    = "FIGURE 07 — Structural Break: Oil-CU Relationship Pre vs Post Shale Era (2015Q1)",
      subtitle = "Slope change between eras tests the shale revolution structural break hypothesis",
      caption  = "Source: NCUA Form 5300; FRB CCAR 2026 Baseline",
      theme    = theme(plot.title=element_text(face="bold",size=12),
                       plot.subtitle=element_text(size=9,colour="#555"))
    )
  save_plot(p07, "07_structural_break.png", w=13, h=9)
}

# =============================================================================
# CHART 08 — Spillover: non-oil CU response by spillover tercile
# =============================================================================
hdr("Chart 08: Spillover exposure terciles")

if ("spillover_exposure" %in% names(panel)) {
  non_oil <- panel[oil_exposure_bin == 0 & !is.na(spillover_exposure)]

  if (nrow(non_oil) > 100) {
    # Assign terciles — guard against non-unique breaks (low-variance spillover)
    breaks <- unique(quantile(non_oil$spillover_exposure,
                              probs=c(0,1/3,2/3,1), na.rm=TRUE))

    if (length(breaks) >= 3) {
      # Enough unique breaks for 3 groups
      n_labels <- length(breaks) - 1
      tier_labels <- c("Low","Medium","High")[1:n_labels]
      non_oil[, spill_tercile := cut(spillover_exposure, breaks=breaks,
                                      labels=tier_labels,
                                      include.lowest=TRUE)]
    } else {
      # Fallback: use median split when quantiles collapse
      med <- median(non_oil$spillover_exposure, na.rm=TRUE)
      non_oil[, spill_tercile := fifelse(spillover_exposure > med,
                                          "High", "Low")]
      non_oil[, spill_tercile := factor(spill_tercile,
                                         levels=c("Low","High"))]
      msg("  NOTE: spillover_exposure low variance — using median split")
    }

    spill_agg <- agg_group(non_oil,
                            intersect(c("dq_rate","pll_rate","netintmrg",
                                        "insured_share_growth",
                                        "costfds","member_growth_yoy"), names(non_oil)),
                            "spill_tercile")
    spill_agg <- merge(spill_agg,
                       mac_spine[,.(yyyyqq,yoy_oil)],
                       by="yyyyqq", all.x=TRUE)

    SPILL_COLS <- c("Low"="#a8d8a8","Medium"="#4a9a6a","High"="#1a5a3a",
                    "Low"="#a8d8a8","High"="#1a5a3a")  # handles 2-level fallback
    SPILL_COLS <- SPILL_COLS[!duplicated(names(SPILL_COLS))]

    sp_panels <- lapply(
      intersect(c("dq_rate","pll_rate","netintmrg","member_growth_yoy"),
                names(spill_agg)),
      function(v) {
        lab <- out_labels[v] %||% v
        ggplot(spill_agg[!is.na(get(v)) & !is.na(spill_tercile)],
               aes(x=cal_date, y=get(v), colour=spill_tercile)) +
          ep_rects() +
          geom_line(linewidth=0.8) +
          scale_colour_manual(values=SPILL_COLS,
                               name="Spillover\nTercile") +
          scale_x_date(date_breaks="3 years", date_labels="%Y") +
          labs(title=lab, x=NULL, y=lab) +
          theme_pub()
      })

    if (length(sp_panels) >= 1) {
      p08 <- wrap_plots(sp_panels, ncol=2) +
        plot_annotation(
          title    = "FIGURE 08 — Indirect Channel: Non-Oil CU Response by Spillover Exposure Tercile",
          subtitle = "High spillover = non-oil state with strong economic linkage to oil states",
          caption  = "Source: NCUA Form 5300; BLS QCEW adjacency-weighted spillover index",
          theme    = theme(plot.title=element_text(face="bold",size=12),
                           plot.subtitle=element_text(size=9,colour="#555"))
        )
      save_plot(p08, "08_spillover_terciles.png", w=13, h=6.5)
    }
  }
}

# =============================================================================
# CHART 09 — Episode Heatmap: oil cycle × CU outcome
# =============================================================================
hdr("Chart 09: Episode heatmap")

# Define oil price episodes
ep_def <- data.table(
  episode = c("Pre-GFC\n2005-07","GFC\n2008-09","Recovery\n2010-13",
              "Shale Bust\n2014-16","Rebound\n2017-19","COVID\n2020",
              "Surge\n2021-22","Post-Surge\n2023-25"),
  yr_from = c(2005, 2008, 2010, 2014, 2017, 2020, 2021, 2023),
  yr_to   = c(2007, 2009, 2013, 2016, 2019, 2020, 2022, 2025),
  q_from  = c(1, 1, 1, 1, 1, 1, 1, 1),
  q_to    = c(4, 4, 4, 4, 4, 4, 4, 4)
)
ep_def[, yyyyqq_from := yr_from * 100L + q_from]
ep_def[, yyyyqq_to   := yr_to   * 100L + q_to]

agg_ep <- agg_quarter(panel, cu_outcomes)
agg_ep <- merge(agg_ep, mac_spine[,.(yyyyqq,yoy_oil)], by="yyyyqq", all.x=TRUE)

# Compute mean per episode
heat_list <- lapply(1:nrow(ep_def), function(i) {
  ep  <- ep_def[i]
  sub <- agg_ep[yyyyqq >= ep$yyyyqq_from & yyyyqq <= ep$yyyyqq_to]
  if (nrow(sub) == 0) return(NULL)
  # Use sapply with explicit NA return when column missing or all-NA
  means <- sapply(cu_outcomes, function(v) {
    if (!v %in% names(sub)) return(NA_real_)
    val <- mean(sub[[v]], na.rm=TRUE)
    if (!is.finite(val)) NA_real_ else val  # NaN/Inf → NA so heatmap shows blank not crash
  })
  as.data.table(c(list(episode=ep$episode), as.list(means)))
})
heat_dt <- rbindlist(heat_list, fill=TRUE)

# ── Diagnostic: show which variables have NA in which episodes ────────────────
cat("\n  Heatmap NA coverage by outcome × episode:\n")
na_check <- heat_dt[, lapply(.SD, function(x) sum(is.na(x))),
                     .SDcols = intersect(cu_outcomes, names(heat_dt))]
if (any(unlist(na_check) > 0)) {
  problem_vars <- names(na_check)[unlist(na_check) > 0]
  cat(sprintf("  Variables with missing episode means: %s\n",
              paste(problem_vars, collapse=", ")))
  cat("  Check: are these variables populated in the call report for recent quarters?\n")
} else {
  cat("  All outcomes have complete episode coverage\n")
}

# Normalise each outcome 0-1 for heatmap colouring
heat_long <- melt(heat_dt, id.vars="episode",
                  variable.name="outcome", value.name="value")
heat_long[, outcome_label := out_labels[as.character(outcome)]]
heat_long[is.na(outcome_label), outcome_label := as.character(outcome)]
heat_long[, norm_val := {
  mn <- min(value, na.rm=TRUE); mx <- max(value, na.rm=TRUE)
  if (is.finite(mn) && is.finite(mx) && mx > mn)
    (value - mn)/(mx - mn)
  else
    rep(0.5, .N)   # constant series → neutral colour, not dropped
}, by=outcome]
heat_long[, episode := factor(episode, levels=ep_def$episode)]

# For outcomes where HIGH value = STRESS → invert so red = bad consistently
# pll_rate: high provisions = stress; dq_rate: high delinquency = stress
# costfds: high cost of funds = stress; chg_tot_lns_ratio: high charge-offs = stress
# pll_per_loan: high per-loan provision = stress
stress_outcomes <- c("dq_rate","chg_tot_lns_ratio","costfds",
                     "pll_rate","pll_per_loan")
heat_long[outcome %in% stress_outcomes, norm_val := 1 - norm_val]

# Display label for value: show 2 decimal places, handle NA gracefully
heat_long[, value_label := fifelse(
  is.na(value) | !is.finite(value),
  "",
  as.character(round(value, 2))
)]

p09 <- ggplot(heat_long,   # use ALL rows including NA norm_val
              aes(x=episode, y=outcome_label, fill=norm_val)) +
  geom_tile(colour="white", linewidth=0.5) +
  # Blank cells (NA norm_val) get light grey fill explicitly
  geom_tile(data=heat_long[is.na(norm_val)],
            fill="#eeeeee", colour="white", linewidth=0.5) +
  geom_text(aes(label=value_label), size=2.5, colour="#1a1a1a") +
  # Show "N/A" in blank cells so they are visibly missing not invisible
  geom_text(data=heat_long[is.na(norm_val)],
            aes(label="N/A"), size=2.2, colour="#aaaaaa", fontface="italic") +
  scale_fill_gradient2(
    low      = "#2166ac",    # blue  = favourable / low stress
    mid      = "#f7f7f7",    # white = neutral / average
    high     = "#d73027",    # red   = elevated stress / tight
    midpoint = 0.5,
    na.value = "#eeeeee",    # grey for NA cells
    name     = "Relative level\n(Blue = good\nRed = stress)"
  ) +
  scale_x_discrete(position="top") +
  labs(
    title    = "FIGURE 09 \u2014 CU Outcome Heatmap by Oil Price Episode",
    subtitle = paste(
      "Blue = favourable conditions | Red = elevated stress | White = average |",
      "Grey (N/A) = variable not available for that episode\n",
      "Colour direction: for stress variables (DQ rate, PLL rate, CoF)",
      "higher values shown as red; for performance variables higher = blue"
    ),
    caption  = "Source: NCUA Form 5300; FRB CCAR 2026 Baseline | Raw episode mean values shown",
    x=NULL, y=NULL
  ) +
  theme_pub() +
  theme(
    axis.text.x     = element_text(size=8.5, angle=0, face="bold"),
    axis.text.y     = element_text(size=8.5),
    legend.position = "right",
    legend.title    = element_text(size=8.5, face="bold"),
    legend.text     = element_text(size=8),
    plot.subtitle   = element_text(size=8, lineheight=1.3),
    panel.grid      = element_blank()
  )

save_plot(p09, "09_episode_heatmap.png", w=14, h=8)   # taller to fit all rows

# =============================================================================
# CHART 10 — Missingness Overview
# =============================================================================
hdr("Chart 10: Missingness")

check_miss <- c(
  "dq_rate","chg_tot_lns_ratio","netintmrg","pcanetworth",
  "networth","costfds","roa","insured_tot","dep_shrcert","acct_018",
  "insured_share_growth","cert_share","loan_to_share","nim_spread",
  "members","member_growth_yoy",
  "pll","lns_tot","lns_tot_n","pll_rate","pll_per_loan",
  "macro_base_pbrent","macro_base_lurc","macro_base_pcpi",
  "macro_base_rmtg","macro_base_phpi","macro_base_uypsav",
  "macro_base_yoy_oil","macro_base_yield_curve","macro_base_fomc_regime",
  "oil_exposure_cont","oil_exposure_bin","spillover_exposure",
  "oil_x_brent","fomc_x_brent","oil_bartik_iv"
)
fv  <- intersect(check_miss, names(panel))
pct <- sapply(panel[, ..fv],
              function(x) round(mean(is.na(x))*100, 1))
miss_tbl <- data.table(
  variable    = names(pct),
  pct_missing = pct,
  category    = fcase(
    names(pct) %like% "macro_", "CCAR Macro",
    names(pct) %like% "oil_|spill|brent|fomc_x|bartik", "Exposure/Interaction",
    default = "CU Outcome"
  )
)[order(-pct_missing)]

p10 <- ggplot(miss_tbl, aes(x=pct_missing,
                              y=reorder(variable, -pct_missing),
                              fill=category)) +
  geom_col(width=0.7) +
  geom_vline(xintercept=c(5,10,20), linetype="dashed",
             colour="#bbbbbb", linewidth=0.3) +
  scale_fill_manual(
    values=c("CU Outcome"="#1a3a5c","CCAR Macro"="#b5470a",
             "Exposure/Interaction"="#2d7a4a"),
    name="Category") +
  scale_x_continuous(labels=percent_format(scale=1),
                     limits=c(0, max(miss_tbl$pct_missing)*1.1)) +
  labs(title   = "FIGURE 10 — Variable Missingness Overview (panel_base)",
       subtitle = "Variables with >10% missing require investigation before modelling",
       caption  = "Source: NCUA Form 5300 + FRB CCAR 2026 Baseline merged panel",
       x="% Missing", y=NULL) +
  theme_pub() +
  theme(legend.position="right",
        panel.grid.major.y=element_blank())

save_plot(p10, "10_missingness.png", w=10, h=7)

# =============================================================================
# CHART NEW-A — Membership Growth: Time Series vs Oil Price
# =============================================================================
hdr("Chart NEW-A: Membership growth vs oil price")

if ("member_growth_yoy" %in% names(panel)) {

  mem_agg <- agg_quarter(panel, "member_growth_yoy")
  mem_agg <- merge(mem_agg,
                   mac_spine[, .(yyyyqq, pbrent, yoy_oil, fomc_regime)],
                   by="yyyyqq", all.x=TRUE)

  # Scale oil for dual axis
  mem_max <- max(abs(mem_agg$member_growth_yoy), na.rm=TRUE)
  pb_max  <- max(abs(mem_agg$pbrent), na.rm=TRUE)
  oil_sc  <- if (is.finite(mem_max) && pb_max > 0) mem_max / pb_max else 1

  # Panel 1: YoY membership growth bar + oil overlay
  pA1 <- ggplot(mem_agg[!is.na(member_growth_yoy)],
                aes(x=cal_date)) +
    ep_rects() +
    geom_col(aes(y=member_growth_yoy, fill=member_growth_yoy >= 0),
             width=70, show.legend=FALSE) +
    geom_line(aes(y=pbrent * oil_sc, colour="PBRENT (scaled)"),
              linewidth=0.7, linetype="dashed") +
    geom_hline(yintercept=0, linewidth=0.4) +
    scale_fill_manual(values=c("TRUE"=COL_POS, "FALSE"=COL_NEG)) +
    scale_colour_manual(values=c("PBRENT (scaled)"=COL_OIL), name=NULL) +
    scale_x_date(date_breaks="2 years", date_labels="%Y") +
    scale_y_continuous(labels=number_format(accuracy=0.1),
                       sec.axis=sec_axis(~./oil_sc,
                                         name="PBRENT ($/bbl)",
                                         labels=dollar_format(suffix="/bbl"))) +
    labs(title    = "Membership Growth (YoY%) vs Brent Oil Price",
         subtitle = "Green = membership expanding | Red = contracting | Dashed = oil price (right axis)",
         x=NULL, y="Member Growth YoY (%)") +
    theme_pub()

  # Panel 2: Direct vs Indirect membership growth comparison
  if ("plot_group" %in% names(panel04) && !is.null(panel04)) {
    mem_grp <- panel04[!is.na(member_growth_yoy) & !is.na(plot_group),
                        .(member_growth_yoy = mean(member_growth_yoy, na.rm=TRUE),
                          cal_date = first(cal_date)),
                        by=.(yyyyqq, plot_group)][order(yyyyqq)]

    pA2 <- ggplot(mem_grp, aes(x=cal_date, y=member_growth_yoy,
                                colour=plot_group)) +
      ep_rects() +
      geom_line(linewidth=0.85) +
      geom_hline(yintercept=0, linewidth=0.4, colour="#888") +
      scale_colour_manual(values=GRP_COLS04, name=NULL) +
      scale_x_date(date_breaks="2 years", date_labels="%Y") +
      scale_y_continuous(labels=number_format(accuracy=0.1)) +
      labs(title    = "Membership Growth: Oil-State vs Non-Oil CUs",
           subtitle = "Buffer hypothesis: oil-state CUs should gain members when oil rises",
           x=NULL, y="Member Growth YoY (%)") +
      theme_pub() +
      theme(legend.position="bottom")
  } else {
    pA2 <- NULL
  }

  if (!is.null(pA2)) {
    pNA <- pA1 / pA2
  } else {
    pNA <- pA1
  }

  pNA <- pNA + plot_annotation(
    title   = "FIGURE NEW-A — Membership Growth: Oil Price Transmission Channel",
    subtitle = "Membership growth is a leading indicator of financial stress — members leave before deposits do",
    caption = "Source: NCUA Form 5300 Call Report; FRB CCAR 2026 Baseline",
    theme   = theme(plot.title=element_text(face="bold",size=12),
                    plot.subtitle=element_text(size=9,colour="#555"))
  )
  save_plot(pNA, "NEW_A_member_growth_vs_oil.png", w=13, h=10)
  msg("  Saved: NEW_A_member_growth_vs_oil.png")
} else {
  msg("  SKIP NEW-A: member_growth_yoy not in panel")
}

# =============================================================================
# CHART NEW-B — Membership Growth Structural Break: Pre vs Post 2015
# =============================================================================
hdr("Chart NEW-B: Membership growth structural break")

if ("member_growth_yoy" %in% names(panel)) {

  # Use the same structural break data already built for Chart 07
  mem_sb <- merge(
    agg_quarter(panel, "member_growth_yoy"),
    mac_spine[, .(yyyyqq, yoy_oil)],
    by="yyyyqq", all.x=TRUE)
  mem_sb[, era := fifelse(yyyyqq < 201501L,
                           "Pre-Shale (2005\u20132014)",
                           "Post-Shale (2015\u20132025)")]

  pNB <- ggplot(mem_sb[!is.na(member_growth_yoy) & !is.na(yoy_oil)],
                aes(x=yoy_oil, y=member_growth_yoy, colour=era)) +
    geom_point(alpha=0.5, size=1.5) +
    geom_smooth(method="lm", se=TRUE, linewidth=1.0) +
    geom_hline(yintercept=0, linewidth=0.4, colour="#888") +
    geom_vline(xintercept=0, linewidth=0.4, colour="#888") +
    scale_colour_manual(
      values=c("Pre-Shale (2005\u20132014)"="#1a3a5c",
               "Post-Shale (2015\u20132025)"="#b5470a"),
      name="Era") +
    scale_x_continuous(labels=number_format(accuracy=1, suffix="%")) +
    scale_y_continuous(labels=number_format(accuracy=0.1, suffix="%")) +
    labs(
      title    = "FIGURE NEW-B — Membership Growth Structural Break: Pre vs Post Shale Era",
      subtitle = paste("Slope test: pre-shale oil rises should attract members (positive slope).",
                       "Post-shale: inflation erodes real income \u2192 slope may flatten or reverse."),
      caption  = "Source: NCUA Form 5300; FRB CCAR 2026 Baseline | Fitted lines = OLS with 95% CI",
      x        = "PBRENT YoY %",
      y        = "Membership Growth YoY %"
    ) +
    theme_pub() +
    theme(legend.position="right")

  save_plot(pNB, "NEW_B_member_growth_structural_break.png", w=11, h=7)
  msg("  Saved: NEW_B_member_growth_structural_break.png")
} else {
  msg("  SKIP NEW-B: member_growth_yoy not in panel")
}

# =============================================================================
# CHART NEW-C — Membership Growth by Asset Tier
# =============================================================================
hdr("Chart NEW-C: Membership growth by asset tier")

if ("member_growth_yoy" %in% names(panel) && "asset_tier" %in% names(panel)) {

  tier_mem <- agg_group(panel, "member_growth_yoy", "asset_tier")
  tier_mem <- merge(tier_mem,
                    mac_spine[,.(yyyyqq, yoy_oil)],
                    by="yyyyqq", all.x=TRUE)

  pNC1 <- ggplot(tier_mem[!is.na(member_growth_yoy)],
                 aes(x=cal_date, y=member_growth_yoy, colour=asset_tier)) +
    ep_rects() +
    geom_line(linewidth=0.8) +
    geom_hline(yintercept=0, linewidth=0.35, colour="#888") +
    scale_colour_manual(values=TIER_COLS, name="Asset Tier") +
    scale_x_date(date_breaks="2 years", date_labels="%Y") +
    scale_y_continuous(labels=number_format(accuracy=0.1)) +
    labs(title    = "Membership Growth by Asset Tier (Time Series)",
         subtitle = "Small CUs (T1/T2) more volatile — single employer/geography concentration",
         x=NULL, y="Member Growth YoY (%)") +
    theme_pub() +
    theme(legend.position="right")

  # Boxplot: distribution of membership growth by tier — shows spread not just mean
  pNC2 <- ggplot(
    panel[!is.na(member_growth_yoy) & !is.na(asset_tier)],
    aes(x=asset_tier, y=member_growth_yoy, fill=asset_tier)) +
    geom_boxplot(outlier.size=0.4, outlier.alpha=0.3,
                 linewidth=0.5, width=0.6) +
    geom_hline(yintercept=0, linewidth=0.4, linetype="dashed", colour="#888") +
    scale_fill_manual(values=TIER_COLS, guide="none") +
    scale_y_continuous(labels=number_format(accuracy=0.1)) +
    labs(title    = "Membership Growth Distribution by Asset Tier",
         subtitle = "Interquartile range reveals size-dependent volatility",
         x        = "Asset Tier",
         y        = "Member Growth YoY (%)") +
    theme_pub()

  pNC <- pNC1 / pNC2 + plot_annotation(
    title   = "FIGURE NEW-C — Membership Growth by Asset Tier (2005\u20132025)",
    subtitle = "T1 < $10M | T2 $10-100M | T3 $100M-$1B | T4 > $1B",
    caption = "Source: NCUA Form 5300 Call Report",
    theme   = theme(plot.title=element_text(face="bold",size=12),
                    plot.subtitle=element_text(size=9,colour="#555"))
  )
  save_plot(pNC, "NEW_C_member_growth_by_tier.png", w=13, h=10)
  msg("  Saved: NEW_C_member_growth_by_tier.png")
} else {
  msg("  SKIP NEW-C: member_growth_yoy or asset_tier not in panel")
}

# =============================================================================
# CHART NEW-D — Membership Growth Lead/Lag Cross-Correlation
# =============================================================================
hdr("Chart NEW-D: Membership growth cross-correlation")

if ("member_growth_yoy" %in% names(panel)) {

  mem_cc_agg <- merge(
    agg_quarter(panel, "member_growth_yoy"),
    mac_spine[, .(yyyyqq, yoy_oil)],
    by="yyyyqq", all.x=TRUE)

  # Cross-correlate member_growth_yoy with oil AND with insured_share_growth
  # to show whether membership leads or lags deposits

  # 1. Oil → membership timing
  cc_mem_oil <- {
    x   <- mem_cc_agg$yoy_oil
    y   <- mem_cc_agg$member_growth_yoy
    ok  <- !is.na(x) & !is.na(y)
    lags <- -4:8
    cors <- sapply(lags, function(k) {
      if (k >= 0) cor(x[ok][1:(sum(ok)-k)], y[ok][(k+1):sum(ok)],
                      use="complete.obs")
      else        cor(x[ok][(-k+1):sum(ok)], y[ok][1:(sum(ok)+k)],
                      use="complete.obs")
    })
    data.table(relationship="Oil YoY → Membership Growth", lag=lags, correlation=cors)
  }

  # 2. Membership → insured share growth (does membership predict deposit flows?)
  if ("insured_share_growth" %in% names(mem_cc_agg)) {
    cc_mem_dep <- {
      x  <- mem_cc_agg$member_growth_yoy
      y  <- mem_cc_agg$insured_share_growth
      ok <- !is.na(x) & !is.na(y)
      lags <- -4:8
      cors <- sapply(lags, function(k) {
        if (k >= 0) cor(x[ok][1:(sum(ok)-k)], y[ok][(k+1):sum(ok)],
                        use="complete.obs")
        else        cor(x[ok][(-k+1):sum(ok)], y[ok][1:(sum(ok)+k)],
                        use="complete.obs")
      })
      data.table(relationship="Membership Growth → Deposit Growth", lag=lags, correlation=cors)
    }
    cc_mem_all <- rbindlist(list(cc_mem_oil, cc_mem_dep))
  } else {
    cc_mem_all <- cc_mem_oil
  }

  REL_COLS <- c("Oil YoY → Membership Growth"         = COL_OIL,
                 "Membership Growth → Deposit Growth"   = COL_DIRECT)

  pND <- ggplot(cc_mem_all[!is.na(correlation)],
                aes(x=lag, y=correlation, colour=relationship)) +
    geom_hline(yintercept=0, linewidth=0.3, colour="#888") +
    geom_vline(xintercept=0, linewidth=0.4, linetype="dashed", colour=COL_OIL) +
    # Significance band: ±1.96/sqrt(N)
    geom_hline(yintercept= 1.96/sqrt(sum(!is.na(mem_cc_agg$member_growth_yoy))),
               linetype="dotted", colour="#aaa", linewidth=0.4) +
    geom_hline(yintercept=-1.96/sqrt(sum(!is.na(mem_cc_agg$member_growth_yoy))),
               linetype="dotted", colour="#aaa", linewidth=0.4) +
    geom_line(linewidth=0.9) +
    geom_point(size=2.2) +
    annotate("text", x=0.1, y=Inf,
             label="\u2190 Oil leads  |  Oil lags \u2192",
             vjust=1.5, hjust=0, size=2.8, colour="#888") +
    annotate("text",
             x=max(cc_mem_all$lag),
             y= 1.96/sqrt(sum(!is.na(mem_cc_agg$member_growth_yoy))) + 0.01,
             label="5% significance", hjust=1, size=2.5, colour="#aaa") +
    scale_colour_manual(values=REL_COLS, name=NULL) +
    scale_x_continuous(breaks=-4:8,
                        labels=c(paste0("Lag\n",4:1),"0",paste0("Lead\n",1:8))) +
    scale_y_continuous(labels=number_format(accuracy=0.01)) +
    labs(
      title    = "FIGURE NEW-D — Lead/Lag: Oil Price, Membership Growth & Deposit Flows",
      subtitle = paste("Does oil price predict membership changes? Does membership change lead deposit outflows?",
                       "Dotted lines = 5% significance threshold"),
      caption  = "Source: NCUA Form 5300; FRB CCAR 2026 Baseline",
      x        = "Quarter lag (positive = first variable leads)",
      y        = "Pearson Correlation"
    ) +
    theme_pub() +
    theme(legend.position="right",
          legend.key.height=unit(0.5,"cm"))

  save_plot(pND, "NEW_D_member_growth_crosscorr.png", w=11, h=6.5)
  msg("  Saved: NEW_D_member_growth_crosscorr.png")
} else {
  msg("  SKIP NEW-D: member_growth_yoy not in panel")
}

# =============================================================================
# CHART NEW-E — PLL Rate: Provision for Loan Losses Dynamics
# =============================================================================
hdr("Chart NEW-E: PLL Rate analysis")

if ("pll_rate" %in% names(panel)) {

  pll_agg <- agg_quarter(panel, c("pll_rate","dq_rate","pll_per_loan"))
  pll_agg <- merge(pll_agg,
                   mac_spine[, .(yyyyqq, pbrent, yoy_oil, fomc_regime)],
                   by="yyyyqq", all.x=TRUE)

  # ── Panel 1: PLL rate vs dq_rate vs oil — the lead/lag story ─────────────
  # PLL is forward-looking (management expectation); dq_rate is backward-looking
  # (loans already delinquent). PLL should LEAD dq_rate.
  pll_max <- max(abs(pll_agg$pll_rate), na.rm=TRUE)
  dq_max  <- max(abs(pll_agg$dq_rate),  na.rm=TRUE)
  pll_scale <- if (is.finite(dq_max) && pll_max > 0) dq_max / pll_max else 1

  pE1 <- ggplot(pll_agg[!is.na(pll_rate) & !is.na(dq_rate)],
                aes(x=cal_date)) +
    ep_rects() +
    geom_line(aes(y=dq_rate, colour="Delinquency Rate (actual)"),
              linewidth=0.85) +
    geom_line(aes(y=pll_rate * pll_scale,
                  colour="PLL Rate (scaled, forward-looking)"),
              linewidth=0.85, linetype="dashed") +
    scale_colour_manual(
      values=c("Delinquency Rate (actual)"              = "#c0392b",
               "PLL Rate (scaled, forward-looking)"     = "#1a3a5c"),
      name=NULL) +
    scale_x_date(date_breaks="2 years", date_labels="%Y") +
    scale_y_continuous(labels=number_format(accuracy=0.01)) +
    labs(title    = "PLL Rate vs Delinquency Rate — Forward vs Backward Looking",
         subtitle = "Dashed = PLL rate (management provision expectation) | Solid = actual delinquency rate",
         x=NULL, y="Rate (%%)") +
    theme_pub()

  # ── Panel 2: PLL rate by direct vs indirect CU group ─────────────────────
  if (!is.null(panel04) && "plot_group" %in% names(panel04)) {
    pll_grp <- panel04[!is.na(pll_rate) & !is.na(plot_group),
                        .(pll_rate = mean(pll_rate, na.rm=TRUE),
                          cal_date = first(cal_date)),
                        by=.(yyyyqq, plot_group)][order(yyyyqq)]

    pE2 <- ggplot(pll_grp, aes(x=cal_date, y=pll_rate, colour=plot_group)) +
      ep_rects() +
      geom_line(linewidth=0.85) +
      geom_hline(yintercept=0, linewidth=0.4, colour="#888") +
      scale_colour_manual(values=GRP_COLS04, name=NULL) +
      scale_x_date(date_breaks="2 years", date_labels="%Y") +
      scale_y_continuous(labels=number_format(accuracy=0.01)) +
      labs(title    = "PLL Rate: Oil-State vs Non-Oil CUs",
           subtitle = "Direct channel: oil bust → energy-sector stress → higher provisions in oil-state CUs",
           x=NULL, y="PLL Rate (%% of avg loans)") +
      theme_pub() + theme(legend.position="bottom")
  } else {
    pE2 <- NULL
  }

  # ── Panel 3: PLL rate structural break — pre/post 2015 ───────────────────
  pll_sb <- merge(
    agg_quarter(panel, "pll_rate"),
    mac_spine[, .(yyyyqq, yoy_oil)], by="yyyyqq", all.x=TRUE)
  pll_sb[, era := fifelse(yyyyqq < 201501L,
                           "Pre-Shale (2005\u20132014)",
                           "Post-Shale (2015\u20132025)")]

  pE3 <- ggplot(pll_sb[!is.na(pll_rate) & !is.na(yoy_oil)],
                aes(x=yoy_oil, y=pll_rate, colour=era)) +
    geom_point(alpha=0.5, size=1.5) +
    geom_smooth(method="lm", se=TRUE, linewidth=1.0) +
    geom_hline(yintercept=0, linewidth=0.4, colour="#888") +
    geom_vline(xintercept=0, linewidth=0.4, colour="#888") +
    scale_colour_manual(
      values=c("Pre-Shale (2005\u20132014)"="#1a3a5c",
               "Post-Shale (2015\u20132025)"="#b5470a"),
      name="Era") +
    labs(title    = "PLL Rate vs Oil Price: Structural Break at 2015Q1",
         subtitle = "Pre-shale: oil bust → PLL spikes. Post-shale: rate cycle dominates provisioning",
         x="PBRENT YoY %%", y="PLL Rate (%% of avg loans)") +
    theme_pub() + theme(legend.position="right")

  # Assemble
  pE_top <- if (!is.null(pE2)) pE1 | pE2 else pE1
  pNE    <- pE_top / pE3 + plot_annotation(
    title    = "FIGURE NEW-E — Provision for Loan Losses (PLL) Rate Dynamics",
    subtitle = paste("PLL rate = pll / avg(lns_tot) — forward-looking credit quality indicator.",
                     "Leads delinquency by 1\u20132 quarters; more sensitive to oil shocks than dq_rate."),
    caption  = "Source: NCUA Form 5300 | pll / ((lns_tot_t + lns_tot_{t-1})/2) \u00d7 100",
    theme    = theme(plot.title=element_text(face="bold",size=12),
                     plot.subtitle=element_text(size=9,colour="#555"))
  )
  save_plot(pNE, "NEW_E_pll_rate_dynamics.png", w=13, h=12)
  msg("  Saved: NEW_E_pll_rate_dynamics.png")
} else {
  msg("  SKIP NEW-E: pll_rate not in panel")
}


# =============================================================================
hdr("Descriptive statistics")

desc_vars <- intersect(c("dq_rate","chg_tot_lns_ratio","netintmrg",
                          "pcanetworth","costfds","roa",
                          "insured_share_growth","cert_share",
                          "loan_to_share","member_growth_yoy",
                          "pll_rate","pll_per_loan"), names(panel))

group_col <- if ("cu_group" %in% names(panel)) "cu_group" else "oil_exposure_bin"

desc_tbl <- panel[!is.na(get(group_col)),
  lapply(.SD, function(x)
    sprintf("%.3f (%.3f)", mean(x,na.rm=TRUE), sd(x,na.rm=TRUE))),
  by = group_col,
  .SDcols = desc_vars]

cat("\n  Descriptive statistics by CU group (mean (sd)):\n")
print(t(desc_tbl), quote=FALSE)


# =============================================================================
# SEVERELY ADVERSE SCENARIO CHARTS  (Charts 11-15)
# =============================================================================
hdr("SECTION: Severely Adverse Scenario Visuals")

# Load severely adverse macro if not already in memory
# Load macro_base if not in memory (EDA may run standalone)
if (!exists("macro_base") || !is.data.table(macro_base)) {
  if (file.exists("Data/macro_base.rds")) {
    macro_base <- readRDS("Data/macro_base.rds")
    setDT(macro_base)
    macro_base[, cal_date := as.Date(paste(year,
                                            Q_MONTH[as.character(quarter)],
                                            "01", sep="-"))]
    msg("  Loaded macro_base.rds")
  } else {
    msg("  WARNING: macro_base.rds not found — skipping severe scenario charts")
    macro_base  <- NULL
    macro_severe <- NULL
  }
}

# Load macro_severe
if (!exists("macro_severe") || !is.data.table(macro_severe)) {
  if (file.exists("Data/macro_severe.rds")) {
    macro_severe <- readRDS("Data/macro_severe.rds")
    setDT(macro_severe)
    macro_severe[, cal_date := as.Date(paste(year,
                                              Q_MONTH[as.character(quarter)],
                                              "01", sep="-"))]
    msg("  Loaded macro_severe.rds")
  } else {
    msg("  WARNING: macro_severe.rds not found — skipping severe scenario charts")
    macro_severe <- NULL
  }
}

# ── CRITICAL: ensure cal_date exists on both macro objects ────────────────────
# macro_base / macro_severe may arrive from 01_data_prep.R in the same session
# without cal_date (it is only added inside the load blocks above when the
# object didn't already exist). Always recompute if missing.
for (.nm in c("macro_base","macro_severe")) {
  .obj <- tryCatch(get(.nm), error=function(e) NULL)
  if (!is.null(.obj) && is.data.table(.obj) &&
      !"cal_date" %in% names(.obj)) {
    .obj[, cal_date := as.Date(paste(year,
                                      Q_MONTH[as.character(quarter)],
                                      "01", sep="-"))]
    assign(.nm, .obj)
    msg("  Added cal_date to %s", .nm)
  }
}
rm(.nm, .obj)

if (!is.null(macro_severe) && !is.null(macro_base) &&
    nrow(macro_severe) > 0 && nrow(macro_base) > 0) {

  # ── Build unified macro comparison spine ──────────────────────────────────
  # Historical: use macro_base up to last observed quarter
  # Projection: from first CCAR quarter onwards — both baseline and severe

  last_hist <- panel[, max(yyyyqq)]
  hist_date <- panel[yyyyqq == last_hist, first(cal_date)]

  # Key macro vars — build scenario comparison table
  base_vars   <- c("macro_base_pbrent","macro_base_lurc","macro_base_pcpi",
                   "macro_base_rmtg","macro_base_yield_curve",
                   "macro_base_fomc_regime","macro_base_yoy_oil",
                   "macro_base_real_rate","macro_base_hike_run")
  severe_vars <- gsub("macro_base_","macro_severe_", base_vars)

  base_spine   <- mac_spine[!is.na(pbrent)]
  severe_spine <- macro_severe[cal_date > hist_date]

  # Rename severe cols to common names for plotting
  sev_rename <- function(dt, pfx="macro_severe_") {
    nms <- names(dt)
    new <- gsub(pfx, "sev_", nms)
    setnames(dt, nms, new)
    dt
  }

  # Forward projection window
  proj_start <- min(macro_base[cal_date > hist_date, cal_date], na.rm=TRUE)
  proj_end   <- max(macro_base$cal_date, na.rm=TRUE)

  base_proj   <- macro_base[cal_date >= proj_start]
  severe_proj <- macro_severe[cal_date >= proj_start & !is.na(cal_date)]

  msg("  Historical through: %s | Projection: %s to %s",
      format(hist_date, "%Y-Q%q"),
      format(proj_start, "%Y-%m"),
      format(proj_end,   "%Y-%m"))

  # ── Shading for projection window ─────────────────────────────────────────
  proj_rect <- list(
    annotate("rect",
             xmin  = proj_start, xmax = proj_end,
             ymin  = -Inf,       ymax = Inf,
             fill  = "#f0f4ff",  alpha = 0.5),
    annotate("text",
             x     = proj_start + (proj_end - proj_start)/2,
             y     = Inf, vjust = 1.5, size = 2.8,
             label = "← CCAR 2026 Projection →",
             colour = "#6688aa", fontface = "italic")
  )

  vline_hist <- geom_vline(xintercept = as.numeric(proj_start),
                            linetype = "dashed",
                            colour = "#888888", linewidth = 0.5)

  # ==========================================================================
  # CHART 11 — PBRENT: Historical + Baseline vs Severely Adverse
  # ==========================================================================
  hdr("Chart 11: PBRENT scenario comparison")

  pb_base <- macro_base[!is.na(macro_base_pbrent),
                          .(cal_date, pbrent = macro_base_pbrent,
                            scenario = "Baseline")]
  pb_sev  <- macro_severe[!is.na(macro_severe_pbrent),
                            .(cal_date, pbrent = macro_severe_pbrent,
                              scenario = "Severely Adverse")]
  pb_hist <- mac_spine[!is.na(pbrent) & cal_date <= hist_date,
                        .(cal_date, pbrent, scenario = "Historical")]

  pb_all  <- rbindlist(list(pb_hist, pb_base, pb_sev), fill=TRUE)
  pb_all  <- pb_all[!duplicated(paste(cal_date, scenario))]

  SCEN_COLS <- c("Historical"        = "#1a1a1a",
                 "Baseline"          = COL_DIRECT,
                 "Severely Adverse"  = COL_NEG)
  SCEN_LT   <- c("Historical"        = "solid",
                 "Baseline"          = "dashed",
                 "Severely Adverse"  = "dashed")

  p11 <- ggplot(pb_all, aes(x=cal_date, y=pbrent,
                              colour=scenario, linetype=scenario)) +
    proj_rect +
    vline_hist +
    ep_rects() +
    geom_line(linewidth=0.85) +
    scale_colour_manual(values=SCEN_COLS, name="Scenario") +
    scale_linetype_manual(values=SCEN_LT, name="Scenario") +
    scale_x_date(date_breaks="2 years", date_labels="%Y") +
    scale_y_continuous(labels=dollar_format(prefix="$",suffix="/bbl")) +
    labs(title    = "FIGURE 11 — Brent Oil Price: Historical + CCAR 2026 Scenarios",
         subtitle = "Shaded blue = CCAR projection window | Severely Adverse: sharp PBRENT decline",
         x=NULL, y="$/barrel",
         caption  = "Source: FRB CCAR 2026 Baseline & Severely Adverse scenarios") +
    theme_pub()

  save_plot(p11, "11_pbrent_scenarios.png", w=11, h=6)

  # ==========================================================================
  # CHART 12 — Key Macro Paths: Baseline vs Severely Adverse (4-panel)
  # ==========================================================================
  hdr("Chart 12: Macro scenario comparison — 4 key variables")

  make_scen_plot <- function(base_col, sev_col, hist_col=NULL,
                              title, y_lab, y_fmt=waiver(),
                              hist_dt=NULL) {

    # Historical line (from macro_spine or passed dt)
    if (is.null(hist_dt)) {
      h_col <- intersect(hist_col, names(mac_spine))[1]
      hist_d <- if (!is.na(h_col) && !is.null(h_col))
        mac_spine[!is.na(get(h_col)) & cal_date <= hist_date,
                   .(cal_date, value=get(h_col), scenario="Historical")]
      else NULL
    } else {
      hist_d <- hist_dt
    }

    # Baseline projection
    b_col <- intersect(base_col, names(macro_base))[1]
    base_d <- if (!is.na(b_col))
      macro_base[!is.na(get(b_col)) & cal_date >= proj_start,
                  .(cal_date, value=get(b_col), scenario="Baseline")]
    else NULL

    # Severe projection
    s_col <- intersect(sev_col, names(macro_severe))[1]
    sev_d <- if (!is.na(s_col))
      macro_severe[!is.na(get(s_col)) & cal_date >= proj_start,
                    .(cal_date, value=get(s_col), scenario="Severely Adverse")]
    else NULL

    all_d <- rbindlist(Filter(Negate(is.null),
                               list(hist_d, base_d, sev_d)), fill=TRUE)
    if (nrow(all_d) == 0) return(NULL)

    ggplot(all_d, aes(x=cal_date, y=value,
                       colour=scenario, linetype=scenario)) +
      proj_rect +
      vline_hist +
      geom_line(linewidth=0.8) +
      scale_colour_manual(values=SCEN_COLS, name=NULL) +
      scale_linetype_manual(values=SCEN_LT,  name=NULL) +
      scale_x_date(date_breaks="3 years", date_labels="%Y") +
      scale_y_continuous(labels=y_fmt) +
      labs(title=title, x=NULL, y=y_lab) +
      theme_pub() +
      theme(legend.position="none")
  }

  p12_panels <- list(
    # PBRENT
    make_scen_plot("macro_base_pbrent",      "macro_severe_pbrent",
                   "pbrent",
                   "Brent Oil ($/bbl)",      "$/bbl",
                   dollar_format(prefix="$")),
    # Unemployment
    make_scen_plot("macro_base_lurc",        "macro_severe_lurc",
                   "lurc",
                   "Unemployment Rate (LURC)", "%",
                   number_format(accuracy=0.1, suffix="%")),
    # CPI
    make_scen_plot("macro_base_pcpi",        "macro_severe_pcpi",
                   "pcpi",
                   "Consumer Price Index (PCPI)", "Index",
                   number_format(accuracy=0.1)),
    # Mortgage rate
    make_scen_plot("macro_base_rmtg",        "macro_severe_rmtg",
                   "rmtg",
                   "30Y Fixed Mortgage Rate (RMTG)", "%",
                   number_format(accuracy=0.1, suffix="%"))
  )
  p12_panels <- Filter(Negate(is.null), p12_panels)

  # Shared legend
  leg_data <- data.frame(
    cal_date = Sys.Date(), value = 1,
    scenario = c("Historical","Baseline","Severely Adverse"))
  leg_p <- ggplot(leg_data, aes(cal_date, value,
                                  colour=scenario, linetype=scenario)) +
    geom_line() +
    scale_colour_manual(values=SCEN_COLS, name="Scenario") +
    scale_linetype_manual(values=SCEN_LT,  name="Scenario") +
    theme_pub() + theme(legend.position="bottom")

  p12 <- wrap_plots(p12_panels, ncol=2) +
    plot_annotation(
      title    = "FIGURE 12 — Key Macro Variables: Baseline vs Severely Adverse (CCAR 2026)",
      subtitle = "Solid = historical | Dashed = CCAR projection | Blue shading = projection window",
      caption  = "Source: FRB CCAR 2026 Baseline & Severely Adverse scenarios",
      theme    = theme(plot.title    = element_text(face="bold", size=12),
                       plot.subtitle = element_text(size=9, colour="#555"))
    )
  save_plot(p12, "12_macro_scenario_comparison.png", w=13, h=10)

  # ==========================================================================
  # CHART 13 — Yield Curve & Rate Environment: Baseline vs Severe
  # ==========================================================================
  hdr("Chart 13: Yield curve & rate scenario")

  p13_panels <- list(
    make_scen_plot("macro_base_yield_curve",  "macro_severe_yield_curve",
                   "yield_curve" %||% NULL,
                   "Yield Curve (10Y-3M, bps)", "bps",
                   number_format(accuracy=0.1)),
    make_scen_plot("macro_base_real_rate",    "macro_severe_real_rate",
                   "real_rate" %||% NULL,
                   "Real Fed Funds Rate (%)", "%",
                   number_format(accuracy=0.1, suffix="%")),
    make_scen_plot("macro_base_fomc_regime",  "macro_severe_fomc_regime",
                   NULL,
                   "FOMC Regime (+1 Hike / -1 Cut)", "",
                   number_format(accuracy=1))
  )
  p13_panels <- Filter(Negate(is.null), p13_panels)

  if (length(p13_panels) >= 2) {
    p13 <- wrap_plots(p13_panels, ncol=2) +
      plot_annotation(
        title    = "FIGURE 13 — Rate Environment & Yield Curve: CCAR 2026 Scenarios",
        subtitle = "NIM compression risk: flat/inverted curve + hiking regime = CoF squeeze",
        caption  = "Source: FRB CCAR 2026; Derived: yield_curve = RS10Y-RS3M; real_rate = RFF-PCPI_yoy",
        theme    = theme(plot.title    = element_text(face="bold", size=12),
                         plot.subtitle = element_text(size=9, colour="#555"))
      )
    save_plot(p13, "13_rate_scenario.png", w=13, h=8)
  }

  # ==========================================================================
  # CHART 14 — Scenario Divergence: How Different Are the Two Paths?
  # ==========================================================================
  hdr("Chart 14: Scenario divergence")

  # Compute gap: severe minus baseline for each key variable
  diverge_vars <- list(
    list(b="macro_base_pbrent",      s="macro_severe_pbrent",
         lab="PBRENT Gap (Severe-Base, $/bbl)"),
    list(b="macro_base_lurc",        s="macro_severe_lurc",
         lab="Unemployment Gap (pp)"),
    list(b="macro_base_pcpi",        s="macro_severe_pcpi",
         lab="CPI Gap (index points)"),
    list(b="macro_base_rmtg",        s="macro_severe_rmtg",
         lab="Mortgage Rate Gap (pp)")
  )

  div_list <- lapply(diverge_vars, function(v) {
    bc <- intersect(v$b, names(macro_base))[1]
    sc <- intersect(v$s, names(macro_severe))[1]
    if (is.na(bc) || is.na(sc)) return(NULL)

    base_d <- macro_base[cal_date >= proj_start & !is.na(get(bc)),
                          .(cal_date, base_val = get(bc))]
    sev_d  <- macro_severe[cal_date >= proj_start & !is.na(get(sc)),
                             .(cal_date, sev_val  = get(sc))]
    mrg    <- merge(base_d, sev_d, by="cal_date", all=FALSE)
    mrg[, `:=`(gap = sev_val - base_val, variable = v$lab)]
    mrg[, .(cal_date, gap, variable)]
  })
  div_dt <- rbindlist(Filter(Negate(is.null), div_list))

  if (nrow(div_dt) > 0) {
    p14 <- ggplot(div_dt, aes(x=cal_date, y=gap,
                               fill=gap < 0)) +
      geom_col(width=70, show.legend=FALSE) +
      geom_hline(yintercept=0, linewidth=0.4) +
      scale_fill_manual(values=c("TRUE"=COL_NEG, "FALSE"=COL_POS)) +
      scale_x_date(date_breaks="1 year", date_labels="%Y") +
      facet_wrap(~variable, scales="free_y", ncol=2) +
      labs(title    = "FIGURE 14 — Scenario Divergence: Severely Adverse minus Baseline",
           subtitle = "Red = severe worse than baseline | Green = severe better | Width of gap = stress magnitude",
           caption  = "Source: FRB CCAR 2026 Baseline & Severely Adverse",
           x=NULL, y="Difference (Severe - Baseline)") +
      theme_pub() +
      theme(strip.text=element_text(size=8, face="bold"))
    save_plot(p14, "14_scenario_divergence.png", w=13, h=8)
  }

  # ==========================================================================
  # CHART 15 — Implied CU Stress Under Severely Adverse
  # (Apply historical oil-CU relationship to severe PBRENT path)
  # ==========================================================================
  hdr("Chart 15: Implied CU stress under severely adverse")

  # Use historical cross-correlation coefficients to project implied stress
  # Simple approach: OLS slope from historical PBRENT YoY vs each CU outcome
  # Apply that slope to the CCAR severe scenario PBRENT path

  if ("macro_severe_yoy_oil" %in% names(macro_severe) &&
      nrow(agg) > 0) {

    # Historical slopes: regress each outcome on PBRENT YoY (aggregated)
    hist_agg <- merge(agg_quarter(panel, cu_outcomes),
                      mac_spine[, .(yyyyqq, yoy_oil)],
                      by="yyyyqq", all.x=TRUE)

    plot_outcomes_15 <- intersect(c("dq_rate","netintmrg",
                                    "insured_share_growth","costfds"),
                                   cu_outcomes)

    implied_list <- lapply(plot_outcomes_15, function(v) {
      d <- hist_agg[!is.na(get(v)) & !is.na(yoy_oil)]
      if (nrow(d) < 20) return(NULL)

      # Historical mean and OLS slope
      hist_mean  <- mean(d[[v]], na.rm=TRUE)
      fit        <- lm(as.formula(paste(v, "~ yoy_oil")), data=d)
      beta       <- coef(fit)["yoy_oil"]
      hist_sd    <- sd(d[[v]], na.rm=TRUE)

      # Apply slope to CCAR severe projected PBRENT YoY
      sev_proj <- macro_severe[!is.na(macro_severe_yoy_oil) &
                                  cal_date >= proj_start,
                                 .(cal_date,
                                   yoy_oil = macro_severe_yoy_oil,
                                   scenario = "Severely Adverse")]
      base_proj2 <- macro_base[!is.na(macro_base_yoy_oil) &
                                   cal_date >= proj_start,
                                  .(cal_date,
                                    yoy_oil = macro_base_yoy_oil,
                                    scenario = "Baseline")]

      proj_both <- rbindlist(list(sev_proj, base_proj2))
      proj_both[, implied := hist_mean + beta * yoy_oil]

      # Historical actual
      hist_line <- d[, .(cal_date, implied=get(v), scenario="Historical")]

      all_lines <- rbindlist(list(hist_line, proj_both[,.(cal_date,implied,scenario)]))
      all_lines[, outcome := v]
      all_lines[, beta_label := sprintf("β=%.3f", beta)]
      all_lines
    })

    implied_dt <- rbindlist(Filter(Negate(is.null), implied_list))

    if (nrow(implied_dt) > 0) {
      implied_dt[, outcome_label := out_labels[outcome]]
      implied_dt[is.na(outcome_label), outcome_label := outcome]

      p15 <- ggplot(implied_dt,
                    aes(x=cal_date, y=implied,
                        colour=scenario, linetype=scenario)) +
        proj_rect +
        vline_hist +
        geom_line(linewidth=0.8) +
        scale_colour_manual(values=SCEN_COLS, name="Scenario") +
        scale_linetype_manual(values=SCEN_LT,  name="Scenario") +
        scale_x_date(date_breaks="2 years", date_labels="%Y") +
        facet_wrap(~outcome_label, scales="free_y", ncol=2) +
        geom_text(data=implied_dt[scenario=="Severely Adverse",
                                    .SD[which.min(cal_date)],
                                    by=outcome_label],
                  aes(label=beta_label), hjust=0, vjust=-0.5,
                  size=2.5, colour=COL_NEG, show.legend=FALSE) +
        labs(title    = "FIGURE 15 — Implied CU Stress Under CCAR 2026 Severely Adverse Scenario",
             subtitle = "Projected using historical OLS slope (β) of CU outcome on PBRENT YoY % change",
             caption  = paste("Note: Simple single-factor projection for illustration.",
                               "β = historical OLS coefficient on PBRENT YoY.",
                               "Full VARX model projections in Script 05."),
             x=NULL, y="Outcome Value") +
        theme_pub() +
        theme(strip.text=element_text(size=8.5, face="bold"),
              legend.position="bottom")

      save_plot(p15, "15_implied_cu_stress_severe.png", w=13, h=10)
    }
  }

  msg("  ✓ Severely adverse charts 11-15 complete")

} else {
  msg("  Severely adverse scenario data not available — charts 11-15 skipped")
}

# =============================================================================
# COMPLETE
# =============================================================================
cat("\n=================================================================\n")
cat(" SCRIPT 02 COMPLETE — FIGURES SAVED TO Figures/\n")
cat("=================================================================\n")
cat("  NEW in this run (member_growth_yoy charts):\n")
cat("    NEW_A_member_growth_vs_oil.png          Time series + direct/indirect\n")
cat("    NEW_B_member_growth_structural_break.png Pre/post 2015 slope test\n")
cat("    NEW_C_member_growth_by_tier.png          By asset tier + boxplot\n")
cat("    NEW_D_member_growth_crosscorr.png        Lead/lag vs oil & deposits\n")
cat("=================================================================\n")
figs <- list.files("Figures", pattern="\\.png$", full.names=FALSE)
for (f in sort(figs)) cat(sprintf("  %s\n", f))
cat("=================================================================\n")
