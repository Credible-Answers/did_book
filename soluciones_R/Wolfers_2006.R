# 0) Setup
* Install (run once)
pkgs <- c("haven", "fixest", "TwoWayFEWeights", "did", "DIDmultiplegtDYN", "didimputation", "data.table")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install) > 0) install.packages(to_install)

library(haven)
library(fixest)
library(TwoWayFEWeights)
library(did)
library(DIDmultiplegtDYN)
library(didimputation)
library(data.table)

# 0.1) Load data
paths <- c(
  "C:/Users/134476.SCIENCESPO/C DE CHAISEMARTIN Dropbox/clément de chaisemartin/A Mini course DID/Applications/Data sets/Wolfers 2006/wolfers2006_didtextbook.dta",
  "C:/Users/134476/C DE CHAISEMARTIN Dropbox/clément de chaisemartin/A Mini course DID/Applications/Data sets/Wolfers 2006/wolfers2006_didtextbook.dta"
)

dta_path <- paths[file.exists(paths)][1]
stopifnot(!is.na(dta_path))

df <- read_dta(dta_path)
df <- as.data.table(df)

* Helpful: ensure identifiers are clean
# (fixest can handle numeric or character ids; you can keep as-is)
 df[, state := as.factor(state)]
 df[, year  := as.integer(year)]

# 1) Static TWFE regression 
m_twfe <- feols(
  div_rate ~ udl | state + year,
  data    = df,
  weights = ~ stpop,
  cluster = ~ state
)

etable(m_twfe)

# 2) Decomposing the static TWFE regression (TwoWayFEWeights)
w_twfe <- twowayfeweights(
  data               = df,
  Y                  = "div_rate",
  G                  = "state",
  T                  = "year",
  D                  = "udl",
  type               = "feTR",
  test_random_weights= "exposurelength",
  weights            = "stpop",
  summary_measures   = TRUE
)

print(w_twfe)
# w_twfe$weights_data (name may vary by version) usually contains the cell weights.

# 3) Testing randomized treatment timing
df_sub <- df[cohort != 1956 & year <= 1968]

m_timing <- feols(
  div_rate ~ i(early_late_never) | 0,
  data    = df_sub,
  weights = ~ stpop,
  cluster = ~ state
)

etable(m_timing)

# 4) Event-study TWFE regression
rel_vars <- grep("^rel_time", names(df), value = TRUE)
stopifnot(length(rel_vars) > 0)

f_es <- as.formula(
  paste0("div_rate ~ ", paste(rel_vars, collapse = " + "), " | state + year")
)

m_es <- feols(
  f_es,
  data    = df,
  weights = ~ stpop,
  cluster = ~ state
)

etable(m_es)

# Joint test of pre-trends (like Stata "test rel_timeminus1 ... rel_timeminus9")
# fixest::wald tests joint nullity of coefficients that match a regex :contentReference[oaicite:5]{index=5}
wald(m_es, keep = "^rel_timeminus")

# Optional: event-study plot (only if your rel_time dummies are named cleanly)
# coefplot(m_es, keep = "^rel_time")

# 5) Decomposing the first estimated effect in the event-study TWFE regression
other_treats <- paste0("rel_time", 2:16)
controls_pre <- paste0("rel_timeminus", 1:9)

w_es1 <- twowayfeweights(
  data            = df,
  Y               = "div_rate",
  G               = "state",
  T               = "year",
  D               = "rel_time1",
  type            = "feTR",
  test_random_weights = "year",
  weights         = "stpop",
  other_treatments= other_treats,
  controls        = controls_pre
)

print(w_es1)

# 6) Sun and Abraham event-study estimators
df[, cohort0 := cohort]
df[is.na(cohort0), cohort0 := 0]

never_val <- max(df$year, na.rm = TRUE) + 1000
df[, cohort_sa := fifelse(cohort0 == 0, never_val, cohort0)]

m_sa <- feols(
  div_rate ~ sunab(cohort_sa, year, ref.p = -1) | state + year,
  data    = df,
  weights = ~ stpop,
  cluster = ~ state
)

etable(m_sa)
iplot(m_sa)      # dynamic ATT by event time
summary(m_sa, agg = "att")  # overall ATT

# 7) Callaway and Sant'Anna estimators (did package)
df[, cohort_cs := cohort0]  # 0 for never-treated

cs_attgt <- att_gt(
  yname         = "div_rate",
  tname         = "year",
  idname        = "state",
  gname         = "cohort_cs",
  data          = df,
  panel         = TRUE,
  control_group = "notyettreated",
  weightsname   = "stpop",
  clustervars   = "state",
  bstrap        = TRUE
)

# Event-study aggregation (similar to agg(event))
cs_event <- aggte(cs_attgt, type = "dynamic")
summary(cs_event)

# Plot (optional)
# ggdid(cs_event)

# 8) de Chaisemartin and D'Haultfoeuille (did_multiplegt_dyn in R)
dc_dh <- did_multiplegt_dyn(
  df        = df,
  outcome   = "div_rate",
  group     = "state",
  time      = "year",
  treatment = "udl",
  effects   = 16,
  placebo   = 9,
  weight    = "stpop",
  cluster   = "state"
)

# Results + plot object are stored in the returned list :contentReference[oaicite:8]{index=8}
print(dc_dh)
# dc_dh$plot  (if you want to display the ggplot)

# 9) Borusyak et al. imputation estimator (didimputation::did_imputation)
# Option A: never-treated coded as 0 (fine per docs)
df[, cohort_imp := cohort_cs]  # 0 for never treated

imp <- did_imputation(
  data        = df,
  yname       = "div_rate",
  gname       = "cohort_imp",
  tname       = "year",
  idname      = "state",
  wname       = "stpop",
  horizon     = 0:15,
  pretrends   = -1:-9,
  cluster_var = "state"
)

print(imp)
