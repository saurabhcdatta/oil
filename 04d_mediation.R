# =============================================================================
# Script 04d — Mediation Analysis: Oil → Macro → CU Balance Sheet → Outcomes
# Baron-Kenny three-link causal chain
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(stringr)
})

hdr <- function(x) cat("\n", strrep("=",70), "\n##", x, "\n", strrep("=",70), "\n")
msg <- function(fmt, ...) cat(sprintf(paste0("  ", fmt, "\n"), ...))

dir.create("Results", showWarnings=FALSE)
dir.create("Figures", showWarnings=FALSE)

set.seed(20260101L)
t0 <- proc.time()

cat("\n")
cat("  ██████████████████████████████████████████████████████████\n")
cat("  ██  Script 04d — Mediation Analysis [v4 — 2026-03-27]  ██\n")
cat("  ██  If you do NOT see this banner — WRONG FILE          ██\n")
cat("  ██████████████████████████████████████████████████████████\n")
cat("\n")

# =============================================================================
# 0. CONFIGURATION — explicit, no auto-detection
# =============================================================================
hdr("SECTION 0: Configuration")

Y_VARS <- c("dq_rate","pll_rate","netintmrg","insured_share_growth",
            "member_growth_yoy","costfds","loan_to_share","pcanetworth")

OUTCOME_LABELS <- c(
  dq_rate="Delinquency Rate", pll_rate="PLL Rate",
  netintmrg="Net Interest Margin", insured_share_growth="Deposit Growth",
  member_growth_yoy="Membership Growth", costfds="Cost of Funds",
  loan_to_share="Loan-to-Share", pcanetworth="Net Worth Ratio"
)

# These will be resolved against actual data after loading
# Candidates listed in order of preference
OIL_CANDIDATES   <- c("macro_base_yoy_oil","yoy_oil","oil_yoy","pbrent_yoy")
MACRO_CANDIDATES <- list(
  unemp  = c("macro_base_lurc","lurc","unemployment"),
  cpi    = c("macro_base_pcpi","pcpi","cpi_level"),
  yield  = c("macro_base_yield_curve","yield_curve","term_spread"),
  mtg    = c("macro_base_rmtg","rmtg","mortgage_rate"),
  hpi    = c("hpi_yoy","hpi","house_price_yoy")
)
MACRO_NICE <- c(unemp="Unemployment",cpi="CPI Level",
                yield="Yield Curve",mtg="Mortgage Rate",hpi="HPI YoY")

P_LAG <- 2L

# =============================================================================
# 1. LOAD DATA
# =============================================================================
hdr("SECTION 1: Load Data")

panel      <- readRDS("Data/panel_model.rds")
macro_base <- readRDS("Data/macro_base.rds")
setDT(panel); setDT(macro_base)

# Print available columns for diagnosis
msg("panel columns (%d): %s ...", length(names(panel)),
    paste(head(names(panel), 20), collapse=", "))
msg("macro_base columns (%d): %s ...", length(names(macro_base)),
    paste(head(names(macro_base), 20), collapse=", "))

# ── Resolve column names against actual data ──────────────────────────────────
resolve_col <- function(candidates, data_names, label) {
  found <- intersect(candidates, data_names)
  if (length(found) > 0) {
    msg("  %-20s -> %s", label, found[1])
    return(found[1])
  }
  msg("  %-20s -> NOT FOUND in %s", label, deparse(substitute(data_names)))
  return(NULL)
}

# Search both panel AND macro_base for each variable
all_cols <- unique(c(names(panel), names(macro_base)))

msg("Resolving variable names:")
OIL_VAR <- resolve_col(OIL_CANDIDATES, all_cols, "Oil YoY")

MACRO_VARS   <- character(0)
MACRO_LABELS <- character(0)
for (nm in names(MACRO_CANDIDATES)) {
  col <- resolve_col(MACRO_CANDIDATES[[nm]], all_cols, MACRO_NICE[nm])
  if (!is.null(col)) {
    MACRO_VARS   <- c(MACRO_VARS, col)
    MACRO_LABELS <- c(MACRO_LABELS, MACRO_NICE[nm])
  }
}
names(MACRO_LABELS) <- MACRO_VARS

if (is.null(OIL_VAR)) stop("Oil variable not found — check OIL_CANDIDATES")
if (length(MACRO_VARS) == 0) stop("No macro variables found — check MACRO_CANDIDATES")
msg("Using oil var: %s | macro vars (%d): %s",
    OIL_VAR, length(MACRO_VARS), paste(MACRO_VARS, collapse=", "))

# ── Merge macro_base onto panel for any vars not already present ──────────────
vars_need_merge <- setdiff(c(OIL_VAR, MACRO_VARS), names(panel))
vars_need_merge <- intersect(vars_need_merge, names(macro_base))
if (length(vars_need_merge) > 0) {
  mc <- c("yyyyqq", vars_need_merge)
  panel <- merge(panel, macro_base[, ..mc], by="yyyyqq", all.x=TRUE)
  msg("Merged from macro_base: %s", paste(vars_need_merge, collapse=", "))
} else {
  msg("All variables already in panel")
}

# ── Build lag columns ─────────────────────────────────────────────────────────
setorder(panel, join_number, yyyyqq)
vars_to_lag <- intersect(c(Y_VARS, MACRO_VARS, OIL_VAR), names(panel))
for (v in vars_to_lag) {
  for (k in 1:P_LAG) {
    nm <- paste0(v, "_lag", k)
    if (!nm %in% names(panel))
      panel[, (nm) := shift(.SD[[v]], k, type="lag"),
             by=join_number, .SDcols=v]
  }
}
msg("Lags built for %d variables", length(vars_to_lag))

# ── Build panel_clean (complete cases on all required columns) ────────────────
req_cols <- unique(c("join_number","yyyyqq",
  Y_VARS, OIL_VAR, MACRO_VARS,
  paste0(OIL_VAR,  "_lag", 1:P_LAG),
  unlist(lapply(Y_VARS,      function(v) paste0(v,"_lag",1:P_LAG))),
  unlist(lapply(MACRO_VARS,  function(v) paste0(v,"_lag",1:P_LAG)))
))
req_cols   <- intersect(req_cols, names(panel))
panel_clean <- panel[complete.cases(panel[, ..req_cols])]
msg("panel_clean: %s obs | %s CUs",
    format(nrow(panel_clean),big.mark=","),
    format(uniqueN(panel_clean$join_number),big.mark=","))

# Verify key vars are present
missing_key <- setdiff(c(OIL_VAR, MACRO_VARS), names(panel_clean))
if (length(missing_key) > 0)
  warning("Missing from panel_clean: ", paste(missing_key, collapse=", "))

# ── Build quarterly aggregate for Link 1 (time-series) ───────────────────────
ts_cols <- intersect(c(OIL_VAR, MACRO_VARS,
  paste0(OIL_VAR, "_lag", 1:P_LAG),
  unlist(lapply(MACRO_VARS, function(v) paste0(v,"_lag",1:P_LAG)))),
  names(panel))
agg_ts  <- panel[, lapply(.SD, mean, na.rm=TRUE), by=yyyyqq, .SDcols=ts_cols]
setorder(agg_ts, yyyyqq)
msg("agg_ts: %d quarters | %d columns", nrow(agg_ts), ncol(agg_ts))

# =============================================================================
# 2. SAFE REGRESSION HELPERS
# =============================================================================
# safe_feols: run panel FE regression, return NULL on any error
safe_feols <- function(formula_str, data) {
  tryCatch(
    feols(as.formula(formula_str), data=data,
          se="cluster", cluster=~join_number, notes=FALSE),
    error=function(e) NULL
  )
}

# safe_r2: extract R² without crashing on degenerate models
safe_r2 <- function(fit) {
  tryCatch(r2(fit, type="within"),
    error=function(e) tryCatch(r2(fit, type="r2"),
      error=function(e2) NA_real_))
}

# safe_coef: extract named coefficient safely
safe_cf <- function(cf, name, col) {
  if (name %in% rownames(cf)) cf[name, col] else NA_real_
}

# =============================================================================
# 3. LINK 1 — OIL → MACRO (time-series OLS)
# =============================================================================
hdr("SECTION 3: Link 1 — Oil → Macro Variables")

link1 <- rbindlist(lapply(MACRO_VARS, function(m) {
  if (!m %in% names(agg_ts)) { msg("SKIP %s — not in agg_ts", m); return(NULL) }
  lag_oil <- intersect(paste0(OIL_VAR,"_lag",1:P_LAG), names(agg_ts))
  lag_m   <- intersect(paste0(m,"_lag",1:P_LAG), names(agg_ts))
  rhs     <- paste(c(OIL_VAR, lag_oil, lag_m), collapse=" + ")
  fit     <- tryCatch(
    lm(as.formula(paste0(m," ~ ",rhs)),
       data=agg_ts[complete.cases(agg_ts[, c(m,OIL_VAR,lag_oil,lag_m), with=FALSE])]),
    error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- summary(fit)$coefficients
  if (!OIL_VAR %in% rownames(cf)) return(NULL)
  ar1  <- if (paste0(m,"_lag1") %in% rownames(cf)) cf[paste0(m,"_lag1"),"Estimate"] else 0
  data.table(
    macro_var   = m,
    macro_label = MACRO_LABELS[m],
    a_path      = cf[OIL_VAR,"Estimate"],
    a_se        = cf[OIL_VAR,"Std. Error"],
    a_p         = cf[OIL_VAR,"Pr(>|t|)"],
    lr_mult     = if (abs(1-ar1)>0.01) cf[OIL_VAR,"Estimate"]/(1-ar1) else NA_real_,
    r2          = summary(fit)$r.squared,
    n_obs       = nobs(fit),
    sig         = fcase(cf[OIL_VAR,"Pr(>|t|)"]<0.01,"***",
                        cf[OIL_VAR,"Pr(>|t|)"]<0.05,"**",
                        cf[OIL_VAR,"Pr(>|t|)"]<0.10,"*", default="")
  )
}))

fwrite(link1, "Results/04d_link1_oil_macro.csv")
msg("Link 1: %d macro vars estimated", nrow(link1))
if (nrow(link1) > 0)
  print(link1[, .(macro_label, a_path=round(a_path,4), a_se=round(a_se,4),
                   a_p=round(a_p,4), sig, lr=round(lr_mult,3), r2=round(r2,3))])

# =============================================================================
# 4. LINK 2 — MACRO → CU OUTCOMES (panel FE)
# =============================================================================
hdr("SECTION 4: Link 2 — Macro Variables → CU Outcomes")

# DIAGNOSTIC: show exactly what's available
msg("panel_clean columns: %d total", length(names(panel_clean)))
msg("Y_VARS in panel_clean: %s",
    paste(intersect(Y_VARS, names(panel_clean)), collapse=", "))
msg("MACRO_VARS in panel_clean: %s",
    paste(intersect(MACRO_VARS, names(panel_clean)), collapse=", "))
msg("OIL_VAR in panel_clean: %s", as.character(OIL_VAR %in% names(panel_clean)))
macro_lags_present <- intersect(paste0(MACRO_VARS[1],"_lag1"), names(panel_clean))
msg("Sample lag check (%s_lag1): %s", MACRO_VARS[1],
    as.character(length(macro_lags_present) > 0))

# If MACRO_VARS not in panel_clean, try adding them now from panel
macro_missing_pc <- setdiff(MACRO_VARS, names(panel_clean))
if (length(macro_missing_pc) > 0) {
  msg("Adding missing macro vars to panel_clean from panel ...")
  avail <- intersect(macro_missing_pc, names(panel))
  if (length(avail) > 0) {
    merge_key <- c("join_number","yyyyqq",avail)
    panel_clean <- merge(panel_clean,
                          panel[, ..merge_key],
                          by=c("join_number","yyyyqq"), all.x=TRUE)
    msg("Added: %s", paste(avail, collapse=", "))
    # Build their lags too
    setorder(panel_clean, join_number, yyyyqq)
    for (v in avail) {
      for (k in 1:P_LAG) {
        nm <- paste0(v,"_lag",k)
        if (!nm %in% names(panel_clean))
          panel_clean[, (nm) := shift(.SD[[v]], k, type="lag"),
                       by=join_number, .SDcols=v]
      }
    }
  }
}

# NOTE: Macro vars vary only over time (not across CUs).
# CU fixed effects absorb them — fall back to pooled OLS if needed.

link2 <- rbindlist(lapply(Y_VARS, function(y) {
  rbindlist(lapply(MACRO_VARS, function(m) {
    if (!all(c(y,m) %in% names(panel_clean))) return(NULL)
    lag_y <- intersect(paste0(y,"_lag",1:P_LAG), names(panel_clean))
    rhs   <- paste(c(m, lag_y), collapse=" + ")

    # Try 1: CU fixed effects (ideal if macro var has residual variation)
    fit <- safe_feols(paste0(y," ~ ",rhs," | join_number"), panel_clean)

    # If collinear (fixest drops the macro var), fall back to pooled OLS
    if (!is.null(fit)) {
      cf <- coef(summary(fit))
      if (!m %in% rownames(cf)) fit <- NULL   # var was dropped — collinear
    }

    # Try 2: pooled OLS with clustered SE (no FE — macro var retained)
    if (is.null(fit)) {
      fit <- tryCatch(
        feols(as.formula(paste0(y," ~ ",rhs)),
              data=panel_clean, se="cluster",
              cluster=~join_number, notes=FALSE),
        error=function(e) NULL)
    }

    if (is.null(fit)) return(NULL)
    cf <- coef(summary(fit))
    if (!m %in% rownames(cf)) return(NULL)

    data.table(outcome=y, macro_var=m,
               b_path=cf[m,"Estimate"], b_se=cf[m,"Std. Error"],
               b_p=cf[m,"Pr(>|t|)"], r2=safe_r2(fit), n_obs=nobs(fit),
               sig=fcase(cf[m,"Pr(>|t|)"]<0.01,"***",
                         cf[m,"Pr(>|t|)"]<0.05,"**",
                         cf[m,"Pr(>|t|)"]<0.10,"*", default=""))
  }))
}))

fwrite(link2, "Results/04d_link2_macro_cu.csv")
msg("Link 2: %d regressions completed", nrow(link2))
if (nrow(link2) > 0) {
  sig2 <- link2[b_p < 0.10]
  if (nrow(sig2) > 0)
    print(sig2[, .(outcome=OUTCOME_LABELS[outcome], macro=MACRO_LABELS[macro_var],
                    b=round(b_path,5), p=round(b_p,4), sig)])
  else msg("No significant relationships at p < 0.10")
}

# =============================================================================
# 5. LINK 3 — LAGGED BALANCE SHEET → OUTCOMES (AR persistence)
# =============================================================================
hdr("SECTION 5: Link 3 — Lagged Balance Sheet Persistence")

# Ensure Y lags exist in panel_clean
y_lags_needed <- unlist(lapply(Y_VARS, function(v) paste0(v,"_lag",1:P_LAG)))
y_lags_missing <- setdiff(y_lags_needed, names(panel_clean))
if (length(y_lags_missing) > 0) {
  msg("Building %d missing Y lag columns in panel_clean ...", length(y_lags_missing))
  setorder(panel_clean, join_number, yyyyqq)
  for (v in Y_VARS) {
    for (k in 1:P_LAG) {
      nm <- paste0(v,"_lag",k)
      if (!nm %in% names(panel_clean) && v %in% names(panel_clean))
        panel_clean[, (nm) := shift(.SD[[v]], k, type="lag"),
                     by=join_number, .SDcols=v]
    }
  }
}

link3 <- rbindlist(lapply(Y_VARS, function(y) {
  lag_terms <- intersect(
    unlist(lapply(Y_VARS, function(v) paste0(v,"_lag",1:P_LAG))),
    names(panel_clean))
  if (length(lag_terms) == 0) return(NULL)
  rhs <- paste(lag_terms, collapse=" + ")
  fit <- safe_feols(paste0(y," ~ ",rhs," | join_number"), panel_clean)
  if (is.null(fit)) return(NULL)
  cf   <- coef(summary(fit))

  # Diagnostic for first outcome only
  if (y == Y_VARS[1])
    msg("  Link3 coef names sample: %s",
        paste(head(rownames(cf), 5), collapse=", "))

  # feols sometimes uses backtick-quoted names — search flexibly
  find_coef <- function(cf, name, col) {
    # exact match first
    if (name %in% rownames(cf)) return(cf[name, col])
    # try with backticks stripped
    clean_names <- gsub("`","", rownames(cf))
    idx <- which(clean_names == name)
    if (length(idx) > 0) return(cf[idx[1], col])
    return(NA_real_)
  }

  l1 <- paste0(y,"_lag1"); l2 <- paste0(y,"_lag2")
  data.table(
    outcome   = y,
    ar1_coef  = find_coef(cf, l1, "Estimate"),
    ar1_p     = find_coef(cf, l1, "Pr(>|t|)"),
    ar2_coef  = find_coef(cf, l2, "Estimate"),
    ar2_p     = find_coef(cf, l2, "Pr(>|t|)"),
    r2_within = safe_r2(fit),
    n_obs     = nobs(fit)
  )
}))

fwrite(link3, "Results/04d_link3_lag_cu.csv")
msg("Link 3: %d outcomes estimated", nrow(link3))
if (nrow(link3) > 0)
  print(link3[, .(outcome=OUTCOME_LABELS[outcome],
                   ar1=round(ar1_coef,3), ar1_p=round(ar1_p,4),
                   r2=round(r2_within,3))])

# =============================================================================
# 6. TOTAL EFFECTS — OIL → CU (path c)
# =============================================================================
hdr("SECTION 6: Total Effects (Path c)")
# Oil var is time-varying but cross-sectionally constant.
# CU FE absorbs it — fall back to pooled OLS with clustered SE if needed.

total_eff <- rbindlist(lapply(Y_VARS, function(y) {
  if (!OIL_VAR %in% names(panel_clean)) return(NULL)
  lag_y <- intersect(paste0(y,"_lag",1:P_LAG), names(panel_clean))
  rhs   <- paste(c(OIL_VAR, lag_y), collapse=" + ")
  fit   <- safe_feols(paste0(y," ~ ",rhs," | join_number"), panel_clean)
  if (!is.null(fit)) {
    cf_tmp <- coef(summary(fit))
    if (!OIL_VAR %in% rownames(cf_tmp)) fit <- NULL
  }
  if (is.null(fit))
    fit <- tryCatch(feols(as.formula(paste0(y," ~ ",rhs)),
                          data=panel_clean, se="cluster",
                          cluster=~join_number, notes=FALSE),
                    error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- coef(summary(fit))
  if (!OIL_VAR %in% rownames(cf)) return(NULL)
  data.table(outcome=y, c_total=cf[OIL_VAR,"Estimate"],
             c_se=cf[OIL_VAR,"Std. Error"], c_p=cf[OIL_VAR,"Pr(>|t|)"])
}))

msg("Total effects: %d outcomes", nrow(total_eff))
if (nrow(total_eff) > 0)
  print(total_eff[, .(outcome=OUTCOME_LABELS[outcome],
                       c=round(c_total,5), p=round(c_p,4))])

# =============================================================================
hdr("SECTION 7: Direct Effects (Path c')")
# Same fallback approach for direct effects

direct_eff <- rbindlist(lapply(Y_VARS, function(y) {
  rbindlist(lapply(MACRO_VARS, function(m) {
    if (!all(c(y,m,OIL_VAR) %in% names(panel_clean))) return(NULL)
    lag_y <- intersect(paste0(y,"_lag",1:P_LAG), names(panel_clean))
    rhs   <- paste(c(OIL_VAR, m, lag_y), collapse=" + ")
    fit   <- safe_feols(paste0(y," ~ ",rhs," | join_number"), panel_clean)
    if (!is.null(fit)) {
      cf_tmp <- coef(summary(fit))
      if (!OIL_VAR %in% rownames(cf_tmp) || !m %in% rownames(cf_tmp)) fit <- NULL
    }
    if (is.null(fit))
      fit <- tryCatch(feols(as.formula(paste0(y," ~ ",rhs)),
                            data=panel_clean, se="cluster",
                            cluster=~join_number, notes=FALSE),
                      error=function(e) NULL)
    if (is.null(fit)) return(NULL)
    cf <- coef(summary(fit))
    data.table(outcome=y, macro_var=m,
               c_prime=safe_cf(cf, OIL_VAR, "Estimate"),
               b_path =safe_cf(cf, m, "Estimate"),
               b_se   =safe_cf(cf, m, "Std. Error"))
  }))
}))

msg("Direct effects: %d pairs", nrow(direct_eff))


hdr("SECTION 8: Mediation Decomposition")

if (nrow(link1) == 0 || nrow(direct_eff) == 0 || nrow(total_eff) == 0) {
  msg("Skipping mediation — upstream steps incomplete")
  msg("  link1 rows: %d | direct_eff rows: %d | total_eff rows: %d",
      nrow(link1), nrow(direct_eff), nrow(total_eff))
  med_dt      <- data.table()
  prop_med    <- data.table()
} else {
  a_paths <- link1[, .(macro_var, a_path, a_se)]
  med_dt  <- merge(direct_eff, a_paths, by="macro_var")
  med_dt  <- merge(med_dt, total_eff, by="outcome")
  med_dt[, indirect_ab   := a_path * b_path]
  med_dt[, prop_mediated := indirect_ab / c_total]
  med_dt[, indirect_pct  := indirect_ab / c_total * 100]
  med_dt[, ab_se_sobel   := sqrt(b_path^2 * a_se^2 + a_path^2 * b_se^2)]
  med_dt[, ab_z          := indirect_ab / pmax(ab_se_sobel, 1e-10)]
  med_dt[, ab_p          := 2 * pnorm(-abs(ab_z))]
  med_dt[, outcome_label := OUTCOME_LABELS[outcome]]
  med_dt[, macro_label   := MACRO_LABELS[macro_var]]
  med_dt[, sig_ab        := fcase(ab_p<0.01,"***",ab_p<0.05,"**",
                                   ab_p<0.10,"*",default="")]

  prop_med <- med_dt[is.finite(prop_mediated),
    .(total_indirect_pct=sum(indirect_pct, na.rm=TRUE),
      n_sig=sum(ab_p < 0.10, na.rm=TRUE)),
    by=outcome]
  prop_med[, outcome_label := OUTCOME_LABELS[outcome]]

  fwrite(med_dt,   "Results/04d_mediation_summary.csv")
  fwrite(prop_med, "Results/04d_proportion_mediated.csv")
  msg("Mediation complete: %d rows", nrow(med_dt))
  cat("\n  Proportion of oil effect mediated through macro channels:\n\n")
  print(prop_med[order(-total_indirect_pct),
    .(outcome_label, pct_indirect=round(total_indirect_pct,1), n_sig)])
}

# =============================================================================
# 9. CHARTS
# =============================================================================
hdr("SECTION 9: Charts")

MACRO_COLS <- setNames(
  c("#2d7a4a","#c04828","#4a2080","#185FA5","#3B6D11"),
  names(MACRO_NICE)
)
# Map detected MACRO_VARS to colours via their code
macro_col_vec <- setNames(
  MACRO_COLS[match(names(MACRO_LABELS), names(MACRO_NICE))],
  MACRO_VARS
)
macro_col_vec[is.na(macro_col_vec)] <- "#888888"

ok <- function(dt) !is.null(dt) && nrow(dt) > 0

# ── Chart 1: Link 1 — oil → macro elasticities ───────────────────────────────
if (ok(link1) && "a_path" %in% names(link1)) {
  p1 <- ggplot(link1, aes(x=reorder(macro_label, abs(a_path)),
                            y=a_path, fill=macro_var)) +
    geom_col(width=0.65) +
    geom_errorbar(aes(ymin=a_path-1.96*a_se, ymax=a_path+1.96*a_se),
                  width=0.25, linewidth=0.6, colour="#333") +
    geom_text(aes(label=paste0(sig,"\nLR=",round(lr_mult,3))),
              hjust=ifelse(link1$a_path>=0,-0.1,1.1), size=3, fontface="bold") +
    geom_hline(yintercept=0, linewidth=0.4, colour="#555") +
    scale_fill_manual(values=macro_col_vec, guide="none") +
    scale_y_continuous(expand=expansion(mult=c(0.2,0.2))) +
    coord_flip() +
    labs(title="LINK 1 \u2014 Oil YoY \u2192 Macro Variables",
         subtitle="OLS quarterly time series | Error bars = 95% CI | LR = long-run multiplier",
         caption="*** p<0.01, ** p<0.05, * p<0.10",
         x=NULL, y="Coefficient per 1pp oil YoY") +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold",size=11),
          plot.subtitle=element_text(size=8.5,colour="#444"),
          plot.caption=element_text(size=7.5,colour="#888",hjust=0),
          panel.grid.major.y=element_blank())
  ggsave("Figures/04d_link1_elasticities.png", p1, width=11, height=6,
         dpi=300, bg="white")
  msg("Chart 1 saved")
}

# ── Chart 2: Link 2 heatmap — macro → CU ─────────────────────────────────────
if (ok(link2) && "b_path" %in% names(link2)) {
  link2[, outcome_label := factor(OUTCOME_LABELS[outcome], levels=OUTCOME_LABELS)]
  link2[, macro_label   := MACRO_LABELS[macro_var]]
  link2[, cell_lbl := fifelse(b_p < 0.10, sprintf("%.4f%s",b_path,sig), "")]
  # Winsorise for colour scale
  q99 <- quantile(abs(link2$b_path), 0.99, na.rm=TRUE)
  link2[, b_plot := pmax(pmin(b_path, q99), -q99)]

  p2 <- ggplot(link2[!is.na(outcome_label)],
               aes(x=outcome_label, y=macro_label, fill=b_plot)) +
    geom_tile(colour="white", linewidth=0.5) +
    geom_text(aes(label=cell_lbl), size=2.5, fontface="bold",
              colour=ifelse(abs(link2$b_plot[!is.na(link2$outcome_label)]) >
                              0.5*q99, "white","#333")) +
    scale_fill_gradient2(low="#185FA5", mid="white", high="#b5470a",
                         midpoint=0, name="Coeff\n(panel FE)") +
    scale_x_discrete(guide=guide_axis(angle=35)) +
    labs(title="LINK 2 \u2014 Macro Variables \u2192 CU Outcomes",
         subtitle="Panel FE | CU fixed effects + lagged Y | Significant cells (p<0.10) labelled",
         caption="Clustered SE by CU | feols()",
         x=NULL, y=NULL) +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold",size=11),
          plot.subtitle=element_text(size=8.5,colour="#444"),
          plot.caption=element_text(size=7.5,colour="#888",hjust=0),
          panel.grid=element_blank())
  ggsave("Figures/04d_link2_macro_heatmap.png", p2, width=13, height=7,
         dpi=300, bg="white")
  msg("Chart 2 saved")
}

# ── Chart 3: Mediation waterfall ──────────────────────────────────────────────
if (ok(med_dt) && ok(total_eff)) {
  med_agg <- med_dt[, .(direct=mean(c_prime,na.rm=TRUE),
                          indirect=sum(indirect_ab,na.rm=TRUE)),
                     by=outcome]
  med_agg[, outcome_label := factor(OUTCOME_LABELS[outcome], levels=OUTCOME_LABELS)]
  med_long <- melt(med_agg[!is.na(outcome_label)],
                    id.vars=c("outcome","outcome_label"),
                    variable.name="type", value.name="effect")
  med_long[, type_lbl := fifelse(type=="direct",
    "Direct (oil \u2192 outcome)",
    "Indirect (oil \u2192 macro \u2192 outcome)")]
  med_long[, type_lbl := factor(type_lbl,
    levels=c("Direct (oil \u2192 outcome)",
             "Indirect (oil \u2192 macro \u2192 outcome)"))]

  p3 <- ggplot(med_long, aes(x=outcome_label, y=effect, fill=type_lbl)) +
    geom_col(position="dodge", width=0.72, colour="white", linewidth=0.2) +
    geom_hline(yintercept=0, linewidth=0.4, colour="#444") +
    scale_fill_manual(values=c("Direct (oil \u2192 outcome)"="#b5470a",
                                "Indirect (oil \u2192 macro \u2192 outcome)"="#4a2080"),
                      name=NULL) +
    scale_x_discrete(guide=guide_axis(angle=35)) +
    scale_y_continuous(labels=function(x) sprintf("%+.5f",x)) +
    labs(title="MEDIATION \u2014 Direct vs Indirect Oil Effect",
         subtitle="Direct = oil \u2192 CU | Indirect = oil \u2192 macro \u2192 CU (summed across mediators)",
         caption="Baron-Kenny | Sobel SE | Panel FE",
         x=NULL, y="Effect size") +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold",size=11),
          plot.subtitle=element_text(size=8.5,colour="#444"),
          plot.caption=element_text(size=7.5,colour="#888",hjust=0),
          panel.grid.major.x=element_blank(), legend.position="top")
  ggsave("Figures/04d_mediation_waterfall.png", p3, width=13, height=7,
         dpi=300, bg="white")
  msg("Chart 3 saved")
}

# ── Chart 4: Proportion mediated ─────────────────────────────────────────────
if (ok(prop_med) && ok(total_eff)) {
  prop_plot <- merge(prop_med, total_eff[, .(outcome, c_p)], by="outcome")
  prop_plot[, outcome_label := factor(OUTCOME_LABELS[outcome], levels=OUTCOME_LABELS)]
  prop_plot[, sig_lbl := fcase(c_p<0.01,"***",c_p<0.05,"**",c_p<0.10,"*",default="(ns)")]
  prop_plot[, pct_capped := pmin(pmax(total_indirect_pct,-150),200)]

  p4 <- ggplot(prop_plot[!is.na(outcome_label)],
               aes(x=reorder(outcome_label,total_indirect_pct),
                   y=pct_capped, fill=total_indirect_pct>50)) +
    geom_col(width=0.72) +
    geom_hline(yintercept=c(0,50,100), linetype=c("solid","dashed","dashed"),
               colour=c("#444","#b5470a","#2d7a4a"), linewidth=c(0.4,0.5,0.5)) +
    geom_text(aes(label=paste0(round(total_indirect_pct,0),"% ",sig_lbl)),
              hjust=ifelse(prop_plot$pct_capped>=0,-0.1,1.1),
              size=3, fontface="bold") +
    scale_fill_manual(values=c("FALSE"="#888","TRUE"="#b5470a"), guide="none") +
    scale_y_continuous(labels=function(x) paste0(x,"%"),
                       expand=expansion(mult=c(0.05,0.25))) +
    coord_flip() +
    labs(title="MEDIATION \u2014 % of Oil Effect via Macro Channels",
         subtitle="% mediated = indirect(ab) / total(c) x 100 | Summed across macro mediators",
         caption="Baron-Kenny | Panel FE regressions",
         x=NULL, y="% of total oil effect mediated") +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(face="bold",size=11),
          plot.subtitle=element_text(size=8.5,colour="#444"),
          plot.caption=element_text(size=7.5,colour="#888",hjust=0),
          panel.grid.major.y=element_blank())
  ggsave("Figures/04d_proportion_mediated.png", p4, width=12, height=7,
         dpi=300, bg="white")
  msg("Chart 4 saved")
}

# ── Chart 5: Causal chain diagram ────────────────────────────────────────────
nodes <- data.table(
  label=c("Oil price\nYoY %","Macro vars\n(unemp, CPI,\nrates, HPI)",
          "CU balance\nsheet lags","CU outcomes\n(dq, pll, nim...)"),
  x=c(1,2.5,4,5.5), y=c(2,2,2,2),
  col=c("#b5470a","#4a2080","#185FA5","#2d7a4a")
)
mean_a  <- if (ok(link1)) mean(abs(link1$a_path),na.rm=TRUE) else NA
mean_b  <- if (ok(link2)) mean(abs(link2$b_path[link2$b_p<0.10]),na.rm=TRUE) else NA
mean_ar <- if (ok(link3)) mean(link3$ar1_coef,na.rm=TRUE) else NA

p5 <- ggplot() +
  geom_segment(aes(x=1.3,xend=2.2,y=2,yend=2),
               arrow=arrow(length=unit(0.3,"cm"),type="closed"),
               colour="#333",linewidth=1.2) +
  geom_segment(aes(x=2.8,xend=3.6,y=2,yend=2),
               arrow=arrow(length=unit(0.3,"cm"),type="closed"),
               colour="#333",linewidth=1.2) +
  geom_segment(aes(x=3.6,xend=4.3,y=2,yend=2),
               arrow=arrow(length=unit(0.3,"cm"),type="closed"),
               colour="#333",linewidth=1.2) +
  geom_segment(aes(x=4.3,xend=5.2,y=2,yend=2),
               arrow=arrow(length=unit(0.3,"cm"),type="closed"),
               colour="#333",linewidth=1.2) +
  annotate("curve",x=1.3,xend=5.2,y=1.5,yend=1.5,curvature=-0.2,
           colour="#b5470a",linewidth=0.9,
           arrow=arrow(length=unit(0.25,"cm"),type="closed")) +
  annotate("text",x=3.25,y=1.2,
           label=sprintf("Direct effect (path c')\nMean = %.5f",
                          if (!is.na(mean_b)) mean_b else 0),
           size=2.8,colour="#b5470a",fontface="bold") +
  annotate("label",x=1.75,y=2.35,
           label=sprintf("Link 1\na=%.4f",if(!is.na(mean_a)) mean_a else 0),
           size=2.8,fill="white",colour="#333",label.size=0.2,fontface="bold") +
  annotate("label",x=3.2,y=2.35,
           label=sprintf("Link 2\nb=%.4f",if(!is.na(mean_b)) mean_b else 0),
           size=2.8,fill="white",colour="#333",label.size=0.2,fontface="bold") +
  annotate("label",x=3.95,y=2.35,
           label=sprintf("Link 3\nAR=%.2f",if(!is.na(mean_ar)) mean_ar else 0),
           size=2.8,fill="white",colour="#333",label.size=0.2,fontface="bold") +
  geom_tile(data=nodes, aes(x=x,y=y,fill=col),
            width=0.9,height=0.55,colour="white",linewidth=1.2) +
  scale_fill_identity() +
  geom_text(data=nodes, aes(x=x,y=y,label=label),
            size=2.8,colour="white",fontface="bold",lineheight=1.2) +
  xlim(0.4,6.2) + ylim(0.8,2.9) +
  labs(title="CAUSAL CHAIN \u2014 Oil Price Shock Transmission to CU Outcomes",
       subtitle="Three-link mediation | All coefficients estimated from data",
       caption="Link 1: OLS time series | Link 2: Panel FE | Link 3: AR persistence") +
  theme_void(base_size=10) +
  theme(plot.title=element_text(face="bold",size=12,hjust=0.5),
        plot.subtitle=element_text(size=9,colour="#444",hjust=0.5),
        plot.caption=element_text(size=7.5,colour="#888",hjust=0.5),
        plot.margin=margin(10,10,10,10))
ggsave("Figures/04d_causal_chain.png", p5, width=14, height=7,
       dpi=300, bg="white")
msg("Chart 5 saved")

# =============================================================================
# 10. MANIFEST
# =============================================================================
hdr("SECTION 10: Output Manifest")

for (grp in c("Results","Figures")) {
  ff <- list.files(grp, pattern="04d_", full.names=TRUE)
  cat(sprintf("\n  %s:\n", grp))
  for (f in ff)
    cat(sprintf("    %-50s [%s bytes]\n", basename(f),
                format(file.size(f), big.mark=",")))
}

t_el <- (proc.time()-t0)["elapsed"]
cat(sprintf("\n  Total: %.1f sec | Script 04d complete\n", t_el))
