# =============================================================================
# SCRIPT 04d — MEDIATION ANALYSIS: OIL → MACRO → CU BALANCE SHEET → OUTCOMES
# Version 5 | 2026-03-27
# Project: Oil Price Shock × Credit Union Research
# Working Dir: S:/Projects/Oil_Price_Shock_2026/
#
# PURPOSE:
#   Formally establish the three-link Baron-Kenny causal chain:
#     Link 1: Oil (PBRENT) → Macro variables (time-series OLS)
#     Link 2: Macro variables → CU balance sheet (panel FE, pooled fallback)
#     Link 3: CU balance sheet → CU outcomes (AR persistence, panel FE)
#   Sobel standard errors for indirect effects.
#   Causal chain diagram output in ggplot2.
#
# CRITICAL FIXEST RULE (enforced throughout):
#   NEVER use coef(summary(fit)) — returns empty matrix for feols objects.
#   ALWAYS use: coef(fit), se(fit), pvalue(fit)   [fixest native API]
#
# INPUTS:
#   Models/panel_model.rds        — CU panel with macro vars ALREADY merged
#   Results/04d_link1_oil_macro.csv  — if Link 1 was saved in a prior run
#
# OUTPUTS:
#   Results/04d_link1_oil_macro.csv
#   Results/04d_link2_macro_cu.csv
#   Results/04d_link3_ar_persistence.csv
#   Results/04d_sobel_indirect.csv
#   Figures/04d_link1_elasticities.png
#   Figures/04d_link2_macro_cu.png
#   Figures/04d_causal_chain_diagram.png
# =============================================================================

cat(sprintf("\n=== 04d_mediation.R | v5 | %s ===\n\n", Sys.time()))

# ── 0. Packages ───────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

# ── 0.1 Utility: extract coefficients from ANY model (fixest or lm) ──────────
# REASON THIS EXISTS: coef(summary(feols_object)) returns an empty named matrix.
# The correct fixest API is coef(), se(), pvalue().
# This helper unifies feols + lm into one interface.
get_cf <- function(fit) {
  if (is.null(fit)) return(NULL)
  if (inherits(fit, "fixest")) {
    # fixest native API — the ONLY safe way
    cf  <- coef(fit)
    se_ <- se(fit)
    pv  <- pvalue(fit)
    terms_ <- names(cf)
    data.table(
      term     = terms_,
      estimate = as.numeric(cf),
      std_err  = as.numeric(se_),
      p_value  = as.numeric(pv)
    )
  } else {
    # lm / standard model
    sm <- summary(fit)$coefficients
    data.table(
      term     = rownames(sm),
      estimate = sm[, 1L],
      std_err  = sm[, 2L],
      p_value  = sm[, 4L]
    )
  }
}

# Helper: find a coefficient row by partial name match (strips backticks)
find_coef <- function(cf_dt, pattern) {
  if (is.null(cf_dt) || nrow(cf_dt) == 0L) return(NULL)
  terms_clean <- gsub("`", "", cf_dt$term)
  idx <- grep(pattern, terms_clean, fixed = TRUE)
  if (length(idx) == 0L) idx <- grep(pattern, terms_clean)
  if (length(idx) == 0L) return(NULL)
  cf_dt[idx[1L]]
}

# msg helper
msg <- function(...) cat(sprintf(...), "\n")

# ── 1. Load panel ─────────────────────────────────────────────────────────────
msg("[04d] Loading panel_model.rds...")
panel_model <- readRDS("Models/panel_model.rds")
setDT(panel_model)

msg("  Rows: %s | CUs: %s | Quarters: %s",
    format(nrow(panel_model), big.mark = ","),
    format(uniqueN(panel_model$join_number), big.mark = ","),
    format(uniqueN(panel_model$yyyyqq), big.mark = ","))

# ── 1.1 Confirm macro vars are present in panel ───────────────────────────────
MACRO_VARS <- c("macro_base_lurc",    # unemployment
                "macro_base_pcpi",    # CPI
                "macro_base_yield_curve",
                "macro_base_rmtg",    # mortgage rate
                "hpi_yoy")            # HPI

present_macros <- intersect(MACRO_VARS, names(panel_model))
missing_macros <- setdiff(MACRO_VARS, names(panel_model))

if (length(missing_macros) > 0) {
  warning("MACRO_VARS missing from panel_model: ",
          paste(missing_macros, collapse = ", "),
          "\n  → Check that macro_base.rds merge happened in Script 02/03.")
}
msg("  Macro vars present in panel: %s / %s",
    length(present_macros), length(MACRO_VARS))

# Key outcome variables (CU balance sheet)
OUTCOME_VARS <- c("nim_ratio", "delinq_rate", "chgoff_rate",
                  "nw_ratio", "loan_to_share", "cost_of_funds")
present_outcomes <- intersect(OUTCOME_VARS, names(panel_model))

# Oil variable
OIL_VAR <- "macro_base_yoy_oil"   # YoY % change in PBRENT, from FRB CCAR
if (!OIL_VAR %in% names(panel_model))
  stop("Oil variable '", OIL_VAR, "' not found in panel_model. Check merge.")

# ── 1.2 Build a QUARTER-LEVEL time series for Links 1 & aggregate checks ─────
# Macro vars are constant within quarter → one row per quarter suffices for
# Link 1 (oil → macro). We also need it to confirm constant-within-quarter.
quarter_ts <- panel_model[,
  .(
    macro_base_yoy_oil    = first(macro_base_yoy_oil),
    macro_base_lurc       = first(get(intersect("macro_base_lurc",    names(panel_model))[1], .SD)),
    macro_base_pcpi       = first(get(intersect("macro_base_pcpi",    names(panel_model))[1], .SD)),
    macro_base_yield_curve = first(get(intersect("macro_base_yield_curve", names(panel_model))[1], .SD)),
    macro_base_rmtg       = first(get(intersect("macro_base_rmtg",    names(panel_model))[1], .SD)),
    hpi_yoy               = if ("hpi_yoy" %in% names(panel_model))
                               first(hpi_yoy) else NA_real_
  ),
  keyby = yyyyqq
]

# Safer version without using `first(get(...))` which can fail:
quarter_ts <- unique(
  panel_model[, c("yyyyqq", OIL_VAR, present_macros), with = FALSE],
  by = "yyyyqq"
)
setkey(quarter_ts, yyyyqq)

msg("  Quarter time-series: %d quarters", nrow(quarter_ts))

# =============================================================================
# LINK 1: OIL → MACRO VARIABLES (time-series OLS)
# y_t = α + β · oil_t + ε_t
# where oil_t = macro_base_yoy_oil (YoY % PBRENT from FRB CCAR)
# Dynamic spec: add lag(y, 1) to control for persistence
# =============================================================================
msg("\n[04d] LINK 1: Oil → Macro variables (time-series OLS)...")

link1_results <- rbindlist(lapply(present_macros, function(mv) {
  df1 <- quarter_ts[!is.na(get(OIL_VAR)) & !is.na(get(mv))]

  # Add lag of dependent variable for AR(1) correction
  df1[, y_lag := shift(get(mv), 1L)]
  df1 <- df1[!is.na(y_lag)]

  if (nrow(df1) < 20L) {
    msg("  %-35s : too few obs (%d), skipping", mv, nrow(df1))
    return(NULL)
  }

  fml <- as.formula(sprintf("`%s` ~ `%s` + y_lag", mv, OIL_VAR))
  fit <- tryCatch(lm(fml, data = df1), error = function(e) NULL)
  cf  <- get_cf(fit)
  row <- find_coef(cf, OIL_VAR)

  if (is.null(row)) {
    msg("  %-35s : OLS ran but oil coef not found", mv)
    return(data.table(macro_var = mv, beta = NA_real_,
                      se = NA_real_, pval = NA_real_, n = nrow(df1)))
  }

  sig <- dplyr::case_when(
    row$p_value < 0.001 ~ "***",
    row$p_value < 0.01  ~ "** ",
    row$p_value < 0.05  ~ "*  ",
    row$p_value < 0.10  ~ ".  ",
    TRUE                 ~ "   "
  )
  msg("  %-35s : β=%+.4f  SE=%.4f  p=%.3f %s",
      mv, row$estimate, row$std_err, row$p_value, sig)

  data.table(macro_var = mv,
             beta      = row$estimate,
             se        = row$std_err,
             pval      = row$p_value,
             stars     = trimws(sig),
             n         = nrow(df1))
}))

fwrite(link1_results, "Results/04d_link1_oil_macro.csv")
msg("  Link 1 → Results/04d_link1_oil_macro.csv")

# Link 1 plot
if (!is.null(link1_results) && nrow(link1_results) > 0 &&
    any(!is.na(link1_results$beta))) {

  p_link1 <- link1_results[!is.na(beta)] |>
    ggplot(aes(x = reorder(macro_var, beta),
               y = beta, fill = beta > 0)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_errorbar(aes(ymin = beta - 1.96 * se,
                      ymax = beta + 1.96 * se), width = 0.25) +
    geom_text(aes(label = paste0(round(beta, 3), " ", stars),
                  y = beta + sign(beta) * 0.002),
              size = 3.2, hjust = 0.5, vjust = -0.4) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#0a9396", "FALSE" = "#ae2012")) +
    labs(title    = "Link 1: Oil Price → Macro Variables",
         subtitle = "OLS with AR(1) correction | Dep var = FRB CCAR macro | *** p<0.001",
         x = NULL, y = "β (effect of YoY PBRENT change)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))

  ggsave("Figures/04d_link1_elasticities.png", p_link1,
         width = 9, height = 5, dpi = 150, bg = "white")
  msg("  Link 1 plot → Figures/04d_link1_elasticities.png")
}

# =============================================================================
# LINK 2: MACRO → CU BALANCE SHEET
#
# IDENTIFICATION PROBLEM:
#   Macro vars (unemployment, CPI, etc.) are CROSS-SECTIONALLY CONSTANT —
#   the same value for every CU in a given quarter.
#   CU fixed effects (| join_number) absorb all within-CU time variation,
#   but quarter FE (| yyyyqq) would absorb the macro vars entirely.
#
# SOLUTION:
#   Use TWO-WAY FE = CU FE + quarter FE is impossible with pure time-series
#   macro regressors. Instead:
#   (a) feols with CU FE only (| join_number), NO quarter FE
#       → macro vars vary over time, CU FE absorbs time-invariant CU traits
#   (b) Pooled OLS with clustered SE as robustness
#   The CU FE-only spec is the primary (it removes omitted CU heterogeneity
#   without absorbing the time-varying macro signal).
# =============================================================================
msg("\n[04d] LINK 2: Macro → CU balance sheet (CU FE, no quarter FE)...")

link2_results <- rbindlist(lapply(present_outcomes, function(ov) {
  rbindlist(lapply(present_macros, function(mv) {

    df2 <- panel_model[!is.na(get(ov)) & !is.na(get(mv)) & !is.na(join_number)]

    if (nrow(df2) < 500L) {
      msg("  %-25s ~ %-30s : too few obs", ov, mv)
      return(NULL)
    }

    # Primary: CU fixed effects only (macro var varies over time → identified)
    fml_fe <- as.formula(sprintf("`%s` ~ `%s` | join_number", ov, mv))

    fit_fe <- tryCatch(
      feols(fml_fe, data = df2, cluster = ~join_number, warn = FALSE),
      error = function(e) NULL
    )

    cf_fe  <- get_cf(fit_fe)
    row_fe <- find_coef(cf_fe, mv)

    # Fallback: pooled OLS with clustered SE (if FE fails or returns NULL)
    fit_pool <- NULL
    row_pool <- NULL
    if (is.null(row_fe)) {
      fml_pool <- as.formula(sprintf("`%s` ~ `%s`", ov, mv))
      fit_pool <- tryCatch(
        feols(fml_pool, data = df2, cluster = ~join_number, warn = FALSE),
        error = function(e) NULL
      )
      cf_pool  <- get_cf(fit_pool)
      row_pool <- find_coef(cf_pool, mv)
    }

    row  <- if (!is.null(row_fe)) row_fe else row_pool
    spec <- if (!is.null(row_fe)) "CU FE" else "Pooled OLS"

    if (is.null(row)) {
      msg("  %-25s ~ %-30s : BOTH specs returned NULL", ov, mv)
      return(data.table(outcome = ov, macro_var = mv,
                        beta = NA_real_, se = NA_real_,
                        pval = NA_real_, spec = "failed", n = nrow(df2)))
    }

    sig <- dplyr::case_when(
      row$p_value < 0.001 ~ "***",
      row$p_value < 0.01  ~ "** ",
      row$p_value < 0.05  ~ "*  ",
      row$p_value < 0.10  ~ ".  ",
      TRUE                 ~ "   "
    )
    msg("  %-25s ~ %-30s : β=%+.4f  p=%.3f %s [%s]",
        ov, mv, row$estimate, row$p_value, trimws(sig), spec)

    data.table(outcome   = ov,
               macro_var = mv,
               beta      = row$estimate,
               se        = row$std_err,
               pval      = row$p_value,
               stars     = trimws(sig),
               spec      = spec,
               n         = nrow(df2))
  }))
}))

fwrite(link2_results, "Results/04d_link2_macro_cu.csv")
msg("  Link 2 → Results/04d_link2_macro_cu.csv")
msg("  Rows with estimates: %d / %d",
    sum(!is.na(link2_results$beta)), nrow(link2_results))

# Link 2 heatmap
if (any(!is.na(link2_results$beta))) {

  p_link2 <- link2_results[!is.na(beta)] |>
    ggplot(aes(x = macro_var, y = outcome, fill = beta)) +
    geom_tile(color = "white") +
    geom_text(aes(label = paste0(round(beta, 3), "\n", stars)),
              size = 2.8) +
    scale_fill_gradient2(low = "#ae2012", mid = "white", high = "#0a9396",
                         midpoint = 0, name = "β") +
    labs(title    = "Link 2: Macro Variables → CU Balance Sheet",
         subtitle = "CU fixed effects | Cluster SE by CU",
         x = "Macro Variable", y = "CU Outcome") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          plot.title   = element_text(face = "bold"))

  ggsave("Figures/04d_link2_macro_cu.png", p_link2,
         width = 10, height = 6, dpi = 150, bg = "white")
  msg("  Link 2 plot → Figures/04d_link2_macro_cu.png")
}

# =============================================================================
# LINK 3: AR PERSISTENCE OF CU BALANCE SHEET OUTCOMES
# y_{i,t} = α_i + ρ · y_{i,t-1} + ε_{i,t}
# ρ is the persistence coefficient: how long a shock to the balance sheet lasts.
# Use CU FE only (same reasoning as Link 2 — no quarter FE to avoid
# eliminating the lagged-DV variation).
# =============================================================================
msg("\n[04d] LINK 3: AR(1) persistence of CU outcomes (CU FE)...")

link3_results <- rbindlist(lapply(present_outcomes, function(ov) {
  lag_col <- paste0(ov, "_lag1")

  df3 <- copy(panel_model[!is.na(get(ov)) & !is.na(join_number)])
  setkey(df3, join_number, yyyyqq)
  df3[, (lag_col) := shift(get(ov), 1L), by = join_number]
  df3 <- df3[!is.na(get(lag_col))]

  if (nrow(df3) < 500L) {
    msg("  %-30s : too few obs (%d)", ov, nrow(df3))
    return(NULL)
  }

  fml_ar <- as.formula(sprintf("`%s` ~ `%s` | join_number", ov, lag_col))

  fit_ar <- tryCatch(
    feols(fml_ar, data = df3, cluster = ~join_number, warn = FALSE),
    error = function(e) NULL
  )

  cf_ar  <- get_cf(fit_ar)
  row_ar <- find_coef(cf_ar, lag_col)

  # Within-R²
  r2_within <- tryCatch(
    as.numeric(r2(fit_ar, "wr2"))[[1L]],
    error = function(e) NA_real_
  )

  if (is.null(row_ar)) {
    msg("  %-30s : AR(1) coef not found in fixest output", ov)
    # Try fallback column name without backticks
    row_ar <- find_coef(cf_ar, gsub("`", "", lag_col))
  }

  if (is.null(row_ar)) {
    msg("  %-30s : AR(1) FAILED", ov)
    return(data.table(outcome = ov, ar1 = NA_real_,
                      se = NA_real_, pval = NA_real_,
                      r2_within = r2_within, n = nrow(df3)))
  }

  msg("  %-30s : AR1=%+.4f  p=%.3f  within-R²=%.3f  N=%s",
      ov, row_ar$estimate, row_ar$p_value,
      r2_within %||% NA_real_,
      format(nrow(df3), big.mark = ","))

  data.table(outcome   = ov,
             ar1       = row_ar$estimate,
             se        = row_ar$std_err,
             pval      = row_ar$p_value,
             r2_within = r2_within,
             n         = nrow(df3))
}))

fwrite(link3_results, "Results/04d_link3_ar_persistence.csv")
msg("  Link 3 → Results/04d_link3_ar_persistence.csv")

# =============================================================================
# SOBEL INDIRECT EFFECTS
# Indirect effect = β_link1 × β_link2
# SE (Sobel) = sqrt( β2² × se1² + β1² × se2² )
# This quantifies: how much of oil's effect on CU outcomes flows through
# each macro channel?
# =============================================================================
msg("\n[04d] Computing Sobel indirect effects...")

# For Sobel we need:
#   Link 1: oil → macro_var  (β1, se1)
#   Link 2: macro_var → CU outcome  (β2, se2)
# Indirect = β1 × β2 through each macro channel

sobel_results <- rbindlist(lapply(present_macros, function(mv) {
  l1 <- link1_results[macro_var == mv]
  if (nrow(l1) == 0 || is.na(l1$beta)) return(NULL)

  rbindlist(lapply(present_outcomes, function(ov) {
    l2 <- link2_results[macro_var == mv & outcome == ov]
    if (nrow(l2) == 0 || is.na(l2$beta)) return(NULL)

    b1 <- l1$beta;  se1 <- l1$se
    b2 <- l2$beta;  se2 <- l2$se

    indirect <- b1 * b2
    se_sobel <- sqrt(b2^2 * se1^2 + b1^2 * se2^2)
    z_sobel  <- indirect / se_sobel
    p_sobel  <- 2 * pnorm(-abs(z_sobel))

    data.table(
      macro_channel = mv,
      cu_outcome    = ov,
      beta_link1    = b1,
      beta_link2    = b2,
      indirect      = indirect,
      se_sobel      = se_sobel,
      z_sobel       = z_sobel,
      p_sobel       = p_sobel
    )
  }))
}))

if (!is.null(sobel_results) && nrow(sobel_results) > 0) {
  fwrite(sobel_results, "Results/04d_sobel_indirect.csv")
  msg("  Sobel results → Results/04d_sobel_indirect.csv")
  print(sobel_results[order(-abs(indirect))][1:min(10, .N)],
        digits = 4)
} else {
  msg("  Sobel: no complete Link 1 × Link 2 pairs found.")
}

# =============================================================================
# CAUSAL CHAIN DIAGRAM (ggplot2)
# Oil → Macro → CU Balance Sheet → CU Outcomes
# Shows β coefficients with significance stars at each arrow
# =============================================================================
msg("\n[04d] Building causal chain diagram...")

# Pull headline coefficients for diagram labels
# Link 1: oil → CPI (strongest expected)
b_oil_cpi <- {
  r <- link1_results[macro_var == "macro_base_pcpi"]
  if (nrow(r) > 0 && !is.na(r$beta))
    sprintf("β=%.3f%s", r$beta, r$stars) else "N/A"
}
b_oil_unemp <- {
  r <- link1_results[macro_var == "macro_base_lurc"]
  if (nrow(r) > 0 && !is.na(r$beta))
    sprintf("β=%.3f%s", r$beta, r$stars) else "N/A"
}
b_oil_mtg <- {
  r <- link1_results[macro_var == "macro_base_rmtg"]
  if (nrow(r) > 0 && !is.na(r$beta))
    sprintf("β=%.3f%s", r$beta, r$stars) else "N/A"
}

# Link 2: CPI → NIM (illustrative — pick most significant)
best_l2 <- if (nrow(link2_results[!is.na(beta)]) > 0)
  link2_results[!is.na(beta)][which.min(pval)]
else
  data.table(outcome = "?", macro_var = "?", beta = NA, stars = "")

b_cpi_nim <- sprintf("β=%.3f%s", best_l2$beta %||% 0, best_l2$stars %||% "")

# Link 3: NIM AR1
b_ar_nim <- {
  r <- link3_results[outcome == "nim_ratio"]
  if (!is.null(r) && nrow(r) > 0 && !is.na(r$ar1))
    sprintf("ρ=%.3f", r$ar1) else "ρ=N/A"
}

# Build the diagram as a ggplot annotation canvas
nodes <- data.frame(
  x     = c(1, 3, 3, 3, 5, 5, 5),
  y     = c(4, 6, 4, 2, 6, 4, 2),
  label = c("Oil\nPrice\n(PBRENT)",
            "Inflation\n(CPI)",
            "Mortgage\nRate",
            "Unemployment",
            "NIM",
            "Delinquency\nRate",
            "Net Worth\nRatio"),
  type  = c("oil", "macro", "macro", "macro",
            "cu", "cu", "cu"),
  stringsAsFactors = FALSE
)

col_map <- c(oil = "#005f73", macro = "#ee9b00", cu = "#0a9396")

arrows <- data.frame(
  x1 = c(1.3, 1.3, 1.3, 3.3, 3.3, 5.3),
  y1 = c(4.1, 4.0, 3.9, 6,   4,   4),
  x2 = c(2.7, 2.7, 2.7, 4.7, 4.7, 6.3),
  y2 = c(6,   4,   2,   6,   4,   4),
  label = c(b_oil_cpi, b_oil_mtg, b_oil_unemp,
            "", b_cpi_nim, b_ar_nim),
  stringsAsFactors = FALSE
)

p_chain <- ggplot() +
  # Arrows
  geom_segment(data = arrows,
    aes(x = x1, y = y1, xend = x2, yend = y2),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    color = "grey40", linewidth = 0.7) +
  geom_label(data = arrows[arrows$label != "", ],
    aes(x = (x1 + x2) / 2, y = (y1 + y2) / 2, label = label),
    size = 2.8, fill = "white", label.size = 0.2, color = "grey30") +
  # Nodes
  geom_point(data = nodes,
    aes(x = x, y = y, color = type), size = 18, alpha = 0.15) +
  geom_text(data = nodes,
    aes(x = x, y = y, label = label, color = type),
    size = 3, fontface = "bold", lineheight = 0.9) +
  scale_color_manual(values = col_map, guide = "none") +
  # Stage labels
  annotate("text", x = 1,   y = 7.5, label = "LINK 1\nOil Shock",
           fontface = "bold", size = 3.5, color = "#005f73") +
  annotate("text", x = 3,   y = 7.5, label = "LINK 2\nMacro Channel",
           fontface = "bold", size = 3.5, color = "#ee9b00") +
  annotate("text", x = 5,   y = 7.5, label = "LINK 3\nCU Persistence",
           fontface = "bold", size = 3.5, color = "#0a9396") +
  xlim(0, 7) + ylim(0.5, 8.5) +
  labs(title    = "Causal Chain: Oil Price → Macro → CU Balance Sheet",
       subtitle = paste("Baron-Kenny mediation | Sobel SEs |",
                        format(Sys.Date(), "%B %Y")),
       caption  = "Coefficients from Links 1–3. Stars: ***p<0.001 **p<0.01 *p<0.05") +
  theme_void(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "grey40"),
        plot.caption  = element_text(size = 8,  hjust = 0.5, color = "grey50"),
        plot.margin   = margin(10, 10, 10, 10))

ggsave("Figures/04d_causal_chain_diagram.png", p_chain,
       width = 12, height = 7, dpi = 150, bg = "white")
msg("  Causal chain diagram → Figures/04d_causal_chain_diagram.png")

# =============================================================================
# SUMMARY CONSOLE PRINT
# =============================================================================
cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║  04d MEDIATION SUMMARY                                   ║\n")
cat("╠══════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Link 1 (Oil → Macro)       : %2d estimates              ║\n",
    sum(!is.na(link1_results$beta))))
cat(sprintf("║  Link 2 (Macro → CU)        : %2d estimates              ║\n",
    sum(!is.na(link2_results$beta))))
cat(sprintf("║  Link 3 (AR persistence)    : %2d estimates              ║\n",
    sum(!is.na(link3_results$ar1))))
cat(sprintf("║  Sobel indirect effects     : %2d computed               ║\n",
    if (!is.null(sobel_results)) nrow(sobel_results) else 0))
cat("╚══════════════════════════════════════════════════════════╝\n")

cat("\n[04d] DONE —", format(Sys.time()), "\n\n")

# =============================================================================
# DEBUGGING CHECKLIST (run in console if Link 2 still returns 0)
# =============================================================================
# 1. Check macro vars are in panel:
#    intersect(MACRO_VARS, names(panel_model))
#
# 2. Check they are NOT all-NA:
#    panel_model[, lapply(.SD, function(x) sum(!is.na(x))), .SDcols=MACRO_VARS]
#
# 3. Check cross-sectional variation within a quarter (should be zero —
#    they are constant within quarter, which is expected):
#    panel_model[yyyyqq == "2010.2",
#                .(sd_lurc = sd(macro_base_lurc, na.rm=TRUE))]
#
# 4. Confirm get_cf() works on a test feols object:
#    test_fit <- feols(nim_ratio ~ macro_base_pcpi | join_number,
#                      data = panel_model[1:10000], warn=FALSE)
#    get_cf(test_fit)   # should return a data.table with pcpi row
# =============================================================================
