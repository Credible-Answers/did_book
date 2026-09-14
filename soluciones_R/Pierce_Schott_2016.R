# 0) Setup
pkgs <- c("haven", "data.table", "estimatr", "TwoWayFEWeights", "StuteTest")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install) > 0) install.packages(to_install)

library(haven)
library(data.table)
library(estimatr)
library(TwoWayFEWeights)
library(StuteTest)

# 0.1) Load data
paths <- c(
  "C:/Users/134476/C DE CHAISEMARTIN Dropbox/clément de chaisemartin/A Mini course DID/Applications/Data sets/Pierce and Schott 2016/pierce_schott_didtextbook.dta",
  "C:/Users/134476.SCIENCESPO/C DE CHAISEMARTIN Dropbox/clément de chaisemartin/A Mini course DID/Applications/Data sets/Pierce and Schott 2016/pierce_schott_didtextbook.dta"
)

dta_path <- paths[file.exists(paths)][1]
stopifnot(!is.na(dta_path))

df <- as.data.table(read_dta(dta_path))
# Helper: HC2 regression with df-adjusted inference
# Stata: reg Y X, vce(hc2, dfadjust)
run_hc2 <- function(y, x = "ntrgap", data = df) {
  f <- as.formula(paste0(y, " ~ ", x))
  lm_robust(f, data = data, se_type = "HC2")  # df-adjusted t by default
}

# 1) TWFE regressions (here: simple cross-sectional regressions)

ys_step1 <- c("delta2001", "delta2002", "delta2004", "delta2005")
ys_step1 <- ys_step1[ys_step1 %in% names(df)]

m_step1 <- setNames(lapply(ys_step1, run_hc2), ys_step1)

# Print summaries
for(nm in names(m_step1)) {
  cat("\n====================\nStep 1:", nm, "~ ntrgap\n====================\n")
  print(summary(m_step1[[nm]]))
}

# 2) Weights analysis
if(all(c("delta2001", "indusid", "cons", "ntrgap") %in% names(df))) {
  w_delta2001 <- twowayfeweights(
    data = df,
    Y    = "delta2001",
    G    = "indusid",
    T    = "cons",
    D    = "ntrgap",
    type = "fdTR",
    D0   = "ntrgap",
    summary_measures = TRUE
  )
  cat("\n====================\nStep 2: twowayfeweights (fdTR)\n====================\n")
  print(w_delta2001)
}

# 3) Test that NTR-gap is as good as randomly assigned
covs <- c("lemp1997", "lemp1998", "lemp1999", "lemp2000")
covs <- covs[covs %in% names(df)]

if(length(covs) > 0) {
  f3 <- as.formula(paste0("ntrgap ~ ", paste(covs, collapse = " + ")))
  m_random <- lm_robust(f3, data = df, se_type = "HC2")
  cat("\n====================\nStep 3: ntrgap ~ pre covariates (HC2)\n====================\n")
  print(summary(m_random))
}

# 4) Stute test (cross-section) + Joint test (panel mode after reshape)
run_stute_cs <- function(y, d = "ntrgap", seed = 1, order = 1, data = df) {
  stute_test(df = data, Y = y, D = d, seed = seed, order = order)
}

# Individual Stute tests for delta20xx
for(y in ys_step1) {
  cat("\n====================\nStep 4: Stute test (cross-section) for", y, "\n====================\n")
  print(run_stute_cs(y, seed = 1))
}

# Joint Stute test (reshape long for delta*)
delta_cols <- grep("^delta\\d{4}$", names(df), value = TRUE)

if(all(c("indusid", "ntrgap") %in% names(df)) && length(delta_cols) >= 2) {
  long_delta <- melt(
    df,
    id.vars = c("indusid", "ntrgap"),
    measure.vars = delta_cols,
    variable.name = "year",
    value.name = "delta"
  )
  long_delta[, year := as.integer(sub("^delta", "", year))]

  # Follow Stata condition: if year >= 2001
  long_delta_post <- long_delta[year >= 2001]

  # Enforce strong balance (StuteTest panel mode requires no gaps)
  yrs <- sort(unique(long_delta_post$year))
  balanced_ids <- long_delta_post[, .N, by = indusid][N == length(yrs), indusid]
  long_delta_post <- long_delta_post[indusid %in% balanced_ids]

  cat("\n====================\nStep 4: Joint Stute test for delta (year>=2001)\n====================\n")
  print(stute_test(
    df    = long_delta_post,
    Y     = "delta",
    D     = "ntrgap",
    group = "indusid",
    time  = "year",
    seed  = 1
  ))
}

# Quasi stayers statistic
if("ntrgap" %in% names(df)) {
  v <- sort(df$ntrgap)
  stat_test_qs <- v[1] / (v[2] - v[1])
  cat("\n====================\nStep 4: Quasi stayers statistic\n====================\n")
  cat("stat_test_qs =", stat_test_qs, "\n")
}

# 5) Pre-trends test: linear
ys_pre <- c("delta1999", "delta1998", "delta1997")
ys_pre <- ys_pre[ys_pre %in% names(df)]

for(y in ys_pre) {
  cat("\n====================\nStep 5:", y, "~ ntrgap (HC2)\n====================\n")
  print(summary(run_hc2(y)))
}

# 6) Pre-trends with industry-specific linear trends
ys_tr_pre <- c("deltalintrend1998", "deltalintrend1997")
ys_tr_pre <- ys_tr_pre[ys_tr_pre %in% names(df)]

for(y in ys_tr_pre) {
  cat("\n====================\nStep 6 (linear):", y, "~ ntrgap (HC2)\n====================\n")
  print(summary(run_hc2(y)))
}

for(y in ys_tr_pre) {
  cat("\n====================\nStep 6 (nonparam, order=0): Stute test for", y, "\n====================\n")
  print(run_stute_cs(y, seed = 1, order = 0))
}

# Joint Stute test for deltalintrend (year<=1998, order=0)
dlt_cols <- grep("^deltalintrend\\d{4}$", names(df), value = TRUE)

if(all(c("indusid", "ntrgap") %in% names(df)) && length(dlt_cols) >= 2) {
  long_dlt <- melt(
    df,
    id.vars = c("indusid", "ntrgap"),
    measure.vars = dlt_cols,
    variable.name = "year",
    value.name = "deltalintrend"
  )
  long_dlt[, year := as.integer(sub("^deltalintrend", "", year))]

  long_dlt_pre <- long_dlt[year <= 1998]

  yrs <- sort(unique(long_dlt_pre$year))
  balanced_ids <- long_dlt_pre[, .N, by = indusid][N == length(yrs), indusid]
  long_dlt_pre <- long_dlt_pre[indusid %in% balanced_ids]

  cat("\n====================\nStep 6 (joint, order=0): Stute test for deltalintrend (year<=1998)\n====================\n")
  print(stute_test(
    df    = long_dlt_pre,
    Y     = "deltalintrend",
    D     = "ntrgap",
    group = "indusid",
    time  = "year",
    order = 0,
    seed  = 1
  ))
}

# 7) Stute test, linear trends (post) + joint tests
ys_tr_post <- c("deltalintrend2001", "deltalintrend2002", "deltalintrend2004", "deltalintrend2005")
ys_tr_post <- ys_tr_post[ys_tr_post %in% names(df)]

for(y in ys_tr_post) {
  cat("\n====================\nStep 7: Stute test for", y, "(seed=1)\n====================\n")
  print(run_stute_cs(y, seed = 1))
}

# Joint post test (year>=2001)
if(exists("long_dlt")) {
  long_dlt_post <- long_dlt[year >= 2001]
  yrs <- sort(unique(long_dlt_post$year))
  balanced_ids <- long_dlt_post[, .N, by = indusid][N == length(yrs), indusid]
  long_dlt_post <- long_dlt_post[indusid %in% balanced_ids]

  cat("\n====================\nStep 7: Joint Stute test for deltalintrend (year>=2001)\n====================\n")
  print(stute_test(
    df    = long_dlt_post,
    Y     = "deltalintrend",
    D     = "ntrgap",
    group = "indusid",
    time  = "year",
    seed  = 1
  ))
}

# 9) Estimators with linear trends (HC2 regressions)
for(y in ys_tr_post) {
  cat("\n====================\nStep 9:", y, "~ ntrgap (HC2)\n====================\n")
  print(summary(run_hc2(y)))
}
