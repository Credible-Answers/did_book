# Gentzkow et al. (2011) — DID Textbook replication in R
# --- Packages ---
pkgs <- c("haven", "fixest", "TwoWayFEWeights", "DIDmultiplegtDYN", "dplyr")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(to_install)) install.packages(to_install)

library(haven)
library(fixest)
library(TwoWayFEWeights)
library(DIDmultiplegtDYN)
library(dplyr)

# Optional: did_multiplegt_stat is on GitHub
if(!requireNamespace("didmultiplegtstat", quietly = TRUE) &&
   !("didmultiplegtstat" %in% installed.packages()[,"Package"])) {
  install.packages("remotes")
  # GitHub repo contains an R package in subfolder /R
  remotes::install_github("chaisemartinPackages/did_multiplegt_stat", subdir = "R")
}
# If installed, load it:
suppressWarnings(suppressMessages(try(library(didmultiplegtstat), silent = TRUE)))

# Load data 
paths <- c(
  "C:/Users/134476.SCIENCESPO/C DE CHAISEMARTIN Dropbox/clément de chaisemartin/A Mini course DID/Applications/Data sets/Gentzkow et al 2011/gentzkowetal_didtextbook.dta",
  "C:/Users/134476/C DE CHAISEMARTIN Dropbox/clément de chaisemartin/A Mini course DID/Applications/Data sets/Gentzkow et al 2011/gentzkowetal_didtextbook.dta"
)
path_ok <- paths[file.exists(paths)][1]
if(is.na(path_ok)) stop("Could not find the .dta file. Update `paths`.")
df <- read_dta(path_ok) |> as.data.frame()

# Make sure key variables exist
needed <- c("prestout","cnty90","year","numdailies","styr",
            "changeprestout","changedailies","lag_numdailies","lag_ishare_urb")
missing <- setdiff(needed, names(df))
if(length(missing)) message("WARNING: missing columns: ", paste(missing, collapse=", "))

# Chapter 5

# 1) TWFE Regression (areg prestout i.year numdailies, absorb(cnty90) cluster(cnty90))
m_twfe <- feols(prestout ~ numdailies | cnty90 + year, data = df, cluster = ~cnty90)
summary(m_twfe)

# Decomposition (twowayfeweights prestout cnty90 year numdailies, type(feTR))
w_twfe <- twowayfeweights(
  data = df, Y = "prestout", G = "cnty90", T = "year", D = "numdailies",
  type = "feTR"
)
print(w_twfe)

# 2) TWFE with state-specific trends (your Stata creates dummies for styr and includes them)
# Estimation in R: easiest is to absorb styr as an extra FE (equivalent to adding all styr dummies)
m_twfe_styr <- feols(prestout ~ numdailies | cnty90 + year + styr, data = df, cluster = ~cnty90)
summary(m_twfe_styr)

# If you want to print coefficient + SE like Stata's di _b[numdailies], _se[numdailies]
cbind(beta = coef(m_twfe_styr)["numdailies"], se = se(m_twfe_styr)["numdailies"])

# Decomposition with styr dummies as "controls" (like controls(styr1-styr683) in Stata):
styr_dum <- model.matrix(~ 0 + factor(styr), data = df)
colnames(styr_dum) <- paste0("styr", seq_len(ncol(styr_dum)))
df2 <- cbind(df, styr_dum)
styr_controls <- colnames(styr_dum)

w_twfe_styr <- twowayfeweights(
  data = df2, Y = "prestout", G = "cnty90", T = "year", D = "numdailies",
  type = "feTR", controls = styr_controls
)
print(w_twfe_styr)

# 3) FD Regression with state-specific trends
m_fd <- feols(changeprestout ~ changedailies | styr, data = df, cluster = ~cnty90)
summary(m_fd)

# Decomposition:
w_fd <- twowayfeweights(
  data = df2, Y = "changeprestout", G = "cnty90", T = "year",
  D = "changedailies", D0 = "numdailies",
  type = "fdTR", controls = styr_controls
)
print(w_fd)

# Assessing if weights correlated with year (test_random_weights(year))
w_fd_testyr <- twowayfeweights(
  data = df2, Y = "changeprestout", G = "cnty90", T = "year",
  D = "changedailies", D0 = "numdailies",
  type = "fdTR", controls = styr_controls,
  test_random_weights = "year"
)
print(w_fd_testyr)

# Chapter 8

# 1) Testing whether change in daily newspapers as good as random
m_rand1 <- feols(changedailies ~ lag_numdailies, data = df, cluster = ~cnty90)
m_rand2 <- feols(changedailies ~ lag_ishare_urb, data = df, cluster = ~cnty90)
summary(m_rand1)
summary(m_rand2)

# 2) Distributed lag TWFE Regression
m_dlag <- feols(prestout ~ numdailies + lag_numdailies | cnty90 + year, data = df, cluster = ~cnty90)
summary(m_dlag)

# Decomposition with other_treatments (two runs, like Stata)
w_dlag_num <- twowayfeweights(
  data = df, Y = "prestout", G = "cnty90", T = "year", D = "numdailies",
  type = "feTR", other_treatments = "lag_numdailies"
)
w_dlag_lag <- twowayfeweights(
  data = df, Y = "prestout", G = "cnty90", T = "year", D = "lag_numdailies",
  type = "feTR", other_treatments = "numdailies"
)
print(w_dlag_num)
print(w_dlag_lag)

# 3) Non-normalized event-study effects (did_multiplegt_dyn ... effects(4) placebo(4) effects_equal(all))
res_dyn <- did_multiplegt_dyn(
  data = df,
  outcome = "prestout", group = "cnty90", time = "year", treatment = "numdailies",
  effects = 4, placebo = 4,
  effects_equal = TRUE,
  graph_off = TRUE
)
print(res_dyn)

# 4) Analyzing the paths averaged in the non-normalized effects (design(0.8,console))
res_dyn_d1 <- did_multiplegt_dyn(
  data = df,
  outcome = "prestout", group = "cnty90", time = "year", treatment = "numdailies",
  effects = 1, design = c(0.8, "console"),
  graph_off = TRUE
)
res_dyn_d2 <- did_multiplegt_dyn(
  data = df,
  outcome = "prestout", group = "cnty90", time = "year", treatment = "numdailies",
  effects = 2, design = c(0.8, "console"),
  graph_off = TRUE
)
res_dyn_d4 <- did_multiplegt_dyn(
  data = df,
  outcome = "prestout", group = "cnty90", time = "year", treatment = "numdailies",
  effects = 4, design = c(0.8, "console"),
  graph_off = TRUE
)

# 5) Normalized event-study effects
res_dyn_norm <- did_multiplegt_dyn(
  data = df,
  outcome = "prestout", group = "cnty90", time = "year", treatment = "numdailies",
  effects = 4, placebo = 4,
  normalized = TRUE, normalized_weights = TRUE,
  effects_equal = TRUE,
  graph_off = TRUE
)
print(res_dyn_norm)
#With graph
p_norm <- res_dyn_norm$plot +
  labs(
    x = "Relative time to change in newspapers",
    y = "Effect",
    title = "Normalized DID estimates"
  ) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(-4, 4, 1)) +
  scale_y_continuous(breaks = seq(-0.01, 0.02, 0.01)) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-0.01, 0.02))

print(p_norm)

# 6) Testing if the lagged number of newspapers affects turnout (subset + same_switchers)
if(all(c("first_change","same_treat_after_first_change") %in% names(df))) {
  df_sub <- df |> filter(year <= first_change | same_treat_after_first_change == 1)

  res_lag_test <- did_multiplegt_dyn(
    data = df_sub,
    outcome = "prestout", group = "cnty90", time = "year", treatment = "numdailies",
    effects = 2,
    effects_equal = TRUE,
    same_switchers = TRUE,
    graph_off = TRUE
  )
  print(res_lag_test)
} else {
  message("Skipping step 6: first_change / same_treat_after_first_change not found in data.")
}

# 7) Estimators assuming away effects of lagged treatments on the outcome (did_multiplegt_stat)
# Create election_number = group(year) like egen group(year)
df <- df |> mutate(election_number = as.integer(factor(year)))

# Run did_multiplegt_stat with placebo(1).
# NOTE: Stata has option exact_match; R package may or may not expose it.
if("did_multiplegt_stat" %in% ls(getNamespace("didmultiplegtstat"))) {
  out_stat <- tryCatch(
    did_multiplegt_stat(df, Y = "prestout", ID = "cnty90", Time = "election_number",
                       D = "numdailies", placebo = 1, exact_match = TRUE),
    error = function(e) {
      message("Exact_match not accepted in your R version; running without it.")
      did_multiplegt_stat(df, Y = "prestout", ID = "cnty90", Time = "election_number",
                         D = "numdailies", placebo = 1)
    }
  )
  print(out_stat)
} else {
  message("did_multiplegt_stat() not available. Check GitHub install step above.")
}

# Reproduce the two Stata tabulations:
if("lag_numdailies" %in% names(df)) {
  if(all(c("first_change","year") %in% names(df))) {
    print(table(df$lag_numdailies[df$year == df$first_change], useNA = "ifany"))
  }
  if(all(c("changedailies","year") %in% names(df))) {
    idx <- !is.na(df$changedailies) & df$changedailies != 0 & df$year != 1868
    print(table(df$lag_numdailies[idx], useNA = "ifany"))
  }
}
