# Solution_R_moser_voena_didtextbook.dta
# Moser & Voena (2012) – DID textbook replication in R
# install.packages(c(
#   "haven","data.table","fixest","ggplot2","stringr","broom",
#   "didimputation","DIDmultiplegtDYN","HonestDiD","synthdid","remotes","fect"
# ))
# remotes::install_github("jonathandroth/pretrends")  # if not installed

library(haven)
library(data.table)
library(fixest)
library(ggplot2)
library(stringr)
library(broom)
library(didimputation)
library(DIDmultiplegtDYN)
library(pretrends)
library(HonestDiD)
library(synthdid)
library(fect)

* 0) Load data

* Choose ONE path (as in your Stata code)
dta_path <- "C:/Users/134476.SCIENCESPO/C DE CHAISEMARTIN Dropbox/clément de chaisemartin/A Mini course DID/Applications/Data sets/Moser and Voena 2012/moser_voena_didtextbook.dta"
# dta_path <- "C:/Users/134476/C DE CHAISEMARTIN Dropbox/clément de chaisemartin/A Mini course DID/Applications/Data sets/Moser and Voena 2012/moser_voena_didtextbook.dta"

dt <- as.data.table(read_dta(dta_path))

* Detect event-study dummy variables (assumed to exist in the dataset)
pre_vars  <- grep("^reltimeminus[0-9]+$", names(dt), value = TRUE)
post_vars <- grep("^reltimeplus[0-9]+$", names(dt), value = TRUE)

* Helper: build ES dataframe (point estimates + CI) from a fixest model
es_df_from_fixest <- function(mod, pre_vars, post_vars){
  ct <- as.data.table(broom::tidy(mod, conf.int = TRUE))
  ct <- ct[term %in% c(pre_vars, post_vars)]

  ct[, rel_time := fifelse(term %in% pre_vars,
                           -as.integer(str_extract(term, "[0-9]+")),
                           as.integer(str_extract(term, "[0-9]+")))]

  * Add the reference period (t = 0) with effect 0 
  zero <- data.table(term="t0", estimate=0, conf.low=0, conf.high=0, rel_time=0L)
  out <- rbindlist(list(ct[, .(rel_time, estimate, conf.low, conf.high)], zero), fill=TRUE)
  setorder(out, rel_time)
  out[]
}

plot_es <- function(esdf, title, xlab, ylab){
  ggplot(esdf, aes(x=rel_time, y=estimate)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_line() +
    geom_point(size=2) +
    geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0) +
    labs(title=title, x=xlab, y=ylab) +
    theme_minimal(base_size=14)
}


* CHAPTER 3
* 1) Static TWFE regression
m_twfe_static <- feols(patents ~ twea | subclass + year, data=dt, cluster=~subclass)
etable(m_twfe_static)

* 2) Equivalence between static TWFE and DID
m_did_static <- feols(patents ~ treatmentgroup + post + twea, data=dt, cluster=~subclass)
etable(m_did_static)

* 3) Testing randomized treatment
m_rand_test <- feols(patents ~ treatmentgroup, data=dt[year <= 1918], cluster=~subclass)
etable(m_rand_test)

* 4) Event-study TWFE regression (with pre and post)
f_es_full <- as.formula(paste(
  "patents ~ treatmentgroup +",
  paste(c(pre_vars, post_vars), collapse=" + "),
  "| year"
))
m_es_full <- feols(f_es_full, data=dt, cluster=~subclass)

# Joint test of pre-trends (Wald)
pre_test_h0 <- paste0(pre_vars, " = 0", collapse=" & ")
wald(m_es_full, pre_test_h0)

# ES plot data
es_full_df <- es_df_from_fixest(m_es_full, pre_vars, post_vars)
res_post_df <- es_full_df[rel_time >= 0 & rel_time <= 21]

p1 <- plot_es(
  es_full_df,
  title="TWFE Event-study estimates",
  xlab="Relative time to year before TWEA",
  ylab="Effect"
)
p1
# ggsave("C:/.../graphES_moser1.pdf", p1, width=8, height=5)

* Verify “by hand” that the 1919 coefficient is a simple DID
m1 <- dt[year==1919 & treatmentgroup==1, mean(patents, na.rm=TRUE)]
m2 <- dt[year==1918 & treatmentgroup==1, mean(patents, na.rm=TRUE)]
m3 <- dt[year==1919 & treatmentgroup==0, mean(patents, na.rm=TRUE)]
m4 <- dt[year==1918 & treatmentgroup==0, mean(patents, na.rm=TRUE)]
did_by_hand <- (m1 - m2) - (m3 - m4)
did_by_hand

* 5) Event-study TWFE regression WITHOUT pre-period estimates (only post)
* Create a "yearpost" factor: collapse all pre-years into a single category
dt[, yearpost := ifelse(year <= 1918, 0L, as.integer(year))]
dt[, yearpost := factor(yearpost)]

f_es_postonly <- as.formula(paste(
  "patents ~ treatmentgroup +",
  paste(post_vars, collapse=" + "),
  "| yearpost"
))
m_es_postonly <- feols(f_es_postonly, data=dt, cluster=~subclass)

* Post-only ES dataframe
ct_postonly <- as.data.table(broom::tidy(m_es_postonly, conf.int=TRUE))
ct_postonly <- ct_postonly[term %in% post_vars]
ct_postonly[, rel_time := as.integer(str_extract(term, "[0-9]+"))]
es_postonly_df <- rbind(
  ct_postonly[, .(rel_time, estimate, conf.low, conf.high)],
  data.table(rel_time=0L, estimate=0, conf.low=0, conf.high=0)
)
setorder(es_postonly_df, rel_time)

* Overlay plot: without pre (blue) vs with pre (red)
df_overlay <- rbind(
  es_postonly_df[, .(rel_time, estimate, conf.low, conf.high, group="Without Pre-Periods")],
  res_post_df[,     .(rel_time, estimate, conf.low, conf.high, group="With Pre-Periods")]
)

p2 <- ggplot(df_overlay, aes(x=rel_time, y=estimate, color=group)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_line() + geom_point(size=2) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0) +
  labs(title="TWFE Event-study estimates", x="Years after TWEA", y="Effect", color=NULL) +
  theme_minimal(base_size=14)
p2
# ggsave("C:/.../graphES_moser2.pdf", p2, width=8, height=5)

* Borusyak et al. DID imputation (numerically close to Stata did_imputation)
dt[, cohort := ifelse(treatmentgroup==1, 1919L, 0L)]  # 0 = never treated

imp_es <- did_imputation(
  data   = dt,
  yname  = "patents",
  gname  = "cohort",
  tname  = "year",
  idname = "subclass",
  horizon = 0:20
)
imp_es

* Equivalent to including all year FE in regression
f_es_all_yearFE_post <- as.formula(paste(
  "patents ~ treatmentgroup +",
  paste(post_vars, collapse=" + "),
  "| year"
))
m_es_all_yearFE_post <- feols(f_es_all_yearFE_post, data=dt, cluster=~subclass)
etable(m_es_all_yearFE_post)

* 6) Linear pre-trends we could fail to detect (power analysis)
pre6_vars <- pre_vars[pre_vars %in% paste0("reltimeminus", 1:6)]
dt_win <- dt[year >= 1912 & year <= 1939]

f_es_pre6 <- as.formula(paste(
  "patents ~ treatmentgroup +",
  paste(c(pre6_vars, post_vars), collapse=" + "),
  "| year"
))
m_es_pre6 <- feols(f_es_pre6, data=dt_win, cluster=~subclass)

* Extract betahat and sigma for ES coefficients
coef_names <- names(coef(m_es_pre6))
es_terms <- c(pre6_vars, post_vars)
keep <- intersect(es_terms, coef_names)

betahat <- coef(m_es_pre6)[keep]
sigma   <- vcov(m_es_pre6)[keep, keep, drop=FALSE]

* Build tVec: negative pre, positive post (assumes reltimeminus1..6 and reltimeplus1..K)
tVec <- c(-1:-length(pre6_vars), 1:length(post_vars))
referencePeriod <- 0

slope50 <- slope_for_power(
  sigma = sigma,
  targetPower = 0.5,
  tVec = tVec,
  referencePeriod = referencePeriod
)
slope50

* Plot with dashed line y = x * slope50
es6_df <- es_df_from_fixest(m_es_pre6, pre6_vars, post_vars)

p3 <- ggplot(es6_df, aes(x=rel_time, y=estimate)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_line() + geom_point(size=2) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0) +
  stat_function(fun = function(x) x * slope50, linetype="dashed") +
  labs(title="TWFE Event-study estimates",
       x="Relative time to year before TWEA", y="Effect") +
  theme_minimal(base_size=14)
p3
# ggsave("C:/.../graphES_moser3.pdf", p3, width=8, height=5)

* 7) Variance of the effect of being exposed for 14 years (year==1932)
dt_1932 <- dt[year == 1932]
sd_t1 <- dt_1932[treatmentgroup==1, sd(diffpatentswrt1918, na.rm=TRUE)]
sd_t0 <- dt_1932[treatmentgroup==0, sd(diffpatentswrt1918, na.rm=TRUE)]
sd_effects <- sd_t1 - sd_t0

m_var <- feols(diffpatentswrt1918 ~ treatmentgroup, data=dt_1932)
b <- coef(m_var)["treatmentgroup"]
ci_manual <- c(b - 1.96*sd_effects, b + 1.96*sd_effects)
list(b=b, sd_effects=sd_effects, ci_manual=ci_manual)

* 8) Placebo test for the variance assumptions
placebo_sd <- function(y){
  dty <- dt[year==y]
  data.table(
    year=y,
    sd_treated = dty[treatmentgroup==1, sd(diffpatentswrt1918, na.rm=TRUE)],
    sd_control = dty[treatmentgroup==0, sd(diffpatentswrt1918, na.rm=TRUE)]
  )
}
placebo_1904 <- placebo_sd(1904)
placebo_1904

placebo_all <- rbindlist(lapply(1900:1939, placebo_sd))
placebo_all[, sd_diff := sd_treated - sd_control]
placebo_all

* CHAPTER 4

* 1) Estimators with controls
f_controls <- as.formula(paste(
  "patents ~",
  paste(c(pre_vars, post_vars), collapse=" + "),
  "+ i(year, patents1900) | treatmentgroup"
))
m_es_controls <- feols(f_controls, data=dt, cluster=~subclass)

* Test pretrends
pre_test_h0_all <- paste0(pre_vars, " = 0", collapse=" & ")
wald(m_es_controls, pre_test_h0_all)

* Plot
es_controls_df <- es_df_from_fixest(m_es_controls, pre_vars, post_vars)
p_controls <- plot_es(es_controls_df,
                      title="TWFE Event-study estimates (controls)",
                      xlab="Relative time to year before TWEA",
                      ylab="Effect")
p_controls
# ggsave("C:/.../graphES_moser_controls.pdf", p_controls, width=8, height=5)

* Test that treatment indicator and baseline covariate are correlated
m_corr <- feols(patents1900 ~ treatmentgroup, data=dt[year==1900])
etable(m_corr)

* DIDmultiplegtDYN with nonparametric trends in patents1900
res_mgt <- DIDmultiplegtDYN::did_multiplegt_dyn(
  data = dt,
  Y = "patents",
  G = "subclass",
  T = "year",
  D = "twea",
  effects = 21,
  placebo = 18,
  trends_nonparam = "patents1900",
  cluster = "subclass"
)
res_mgt

* 2) Interactive fixed effects (IFE) using fect
fit_cv <- fect(
  formula = patents ~ twea,
  data = dt,
  index = c("subclass","year"),
  method = "ife",
  r = 0:4,
  CV = TRUE
)
fit_cv

set.seed(1)
fit_ife <- fect(
  formula = patents ~ twea,
  data = dt,
  index = c("subclass","year"),
  method = "ife",
  r = 2,
  se = TRUE
)

* Inspect outputs (structure depends on package version)
names(fit_ife)
summary(fit_ife)

* 3) Synthetic control (SC) and 4) Synthetic DID (SDID) using synthdid
* In R, we build a balanced outcome matrix Y (units x time), order controls first,
* define N0 (# controls) and T0 (# pre periods), then bootstrap the effect curve.

* Build wide matrix Y: subclasses x years
years <- sort(unique(dt$year))
wide <- dcast(dt, subclass ~ year, value.var="patents")
Y <- as.matrix(wide[, -1])
rownames(Y) <- wide$subclass

* Treated indicator per unit: max(twea)
unit_treat <- dt[, .(treated = max(twea, na.rm=TRUE)), by=subclass]
setkey(unit_treat, subclass)
treated_flag <- unit_treat[.(rownames(Y)), treated]

* Order: controls first, treated last
ord <- order(treated_flag)
Y <- Y[ord, , drop=FALSE]
treated_flag <- treated_flag[ord]

N0 <- sum(treated_flag == 0)
T0 <- sum(years <= 1918)  # pre period ends in 1918

boot_effect_curve <- function(est_fun, Y, N0, T0, B=200, alpha=0.05, seed=1){
  set.seed(seed)
  N <- nrow(Y); TT <- ncol(Y)
  H <- TT - T0
  curves <- matrix(NA_real_, nrow=B, ncol=H)

  for(b in 1:B){
    idx_c <- sample(1:N0, N0, replace=TRUE)
    idx_t <- sample((N0+1):N, N-N0, replace=TRUE)
    Yb <- rbind(Y[idx_c, , drop=FALSE], Y[idx_t, , drop=FALSE])
    estb <- est_fun(Yb, N0, T0)
    curves[b, ] <- synthdid_effect_curve(estb)
  }

  est0 <- est_fun(Y, N0, T0)
  point <- synthdid_effect_curve(est0)

  lo <- apply(curves, 2, quantile, probs=alpha/2, na.rm=TRUE)
  hi <- apply(curves, 2, quantile, probs=1-alpha/2, na.rm=TRUE)

  data.table(h=1:H, estimate=point, conf.low=lo, conf.high=hi)
}

* SC (standard)
sc_df <- boot_effect_curve(synthdid::sc_estimate, Y, N0, T0, B=200, seed=1)

* SC (demeaned)
pre_mean <- rowMeans(Y[, 1:T0, drop=FALSE], na.rm=TRUE)
Y_dm <- Y - pre_mean
sc_dm_df <- boot_effect_curve(synthdid::sc_estimate, Y_dm, N0, T0, B=200, seed=1)

* SDID
sdid_df <- boot_effect_curve(synthdid::synthdid_estimate, Y, N0, T0, B=200, seed=1)

* Compare with TWFE post (h = 1..21)
twfe_post_plotdf <- res_post_df[rel_time >= 1, .(h=rel_time, estimate, conf.low, conf.high)]
twfe_post_plotdf[, group := "TWFE"]

plot_compare <- function(dfA, title){
  df <- rbind(dfA, twfe_post_plotdf, fill=TRUE)
  ggplot(df, aes(x=h, y=estimate, color=group)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_line() + geom_point(size=2) +
    geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0) +
    labs(title=title, x="Years after TWEA", y="Effect", color=NULL) +
    theme_minimal(base_size=14)
}

p_sc   <- plot_compare(sc_df[, .(h, estimate, conf.low, conf.high, group="SC")],
                       "Event-study estimates: SC vs TWFE")
p_scdm <- plot_compare(sc_dm_df[, .(h, estimate, conf.low, conf.high, group="Demeaned SC")],
                       "Event-study estimates: Demeaned SC vs TWFE")
p_sdid <- plot_compare(sdid_df[, .(h, estimate, conf.low, conf.high, group="SDID")],
                       "Event-study estimates: SDID vs TWFE")
p_sc
p_scdm
p_sdid
# ggsave("C:/.../graphES_moser_sc_did.pdf", p_sc, width=8, height=5)
# ggsave("C:/.../graphES_moser_demeaned_sc_did.pdf", p_scdm, width=8, height=5)
# ggsave("C:/.../graphES_moser_sd_did.pdf", p_sdid, width=8, height=5)

* 5) Sensitivity analysis (Rambachan & Roth) using HonestDiD
m_honest <- feols(
  as.formula(paste(
    "patents ~", paste(c(pre_vars, post_vars), collapse=" + "),
    "| year + treatmentgroup"
  )),
  data=dt,
  cluster=~subclass
)

betahat <- coef(m_honest)[c(pre_vars, post_vars)]
sigma   <- vcov(m_honest)[c(pre_vars, post_vars), c(pre_vars, post_vars)]

numPre  <- length(pre_vars)
numPost <- length(post_vars)

* l_vec selects horizon 14 (same as Stata: l_vec[14]=1)
l_vec <- rep(0, numPost)
l_vec[14] <- 1

sens <- HonestDiD::createSensitivityResults_relativeMagnitudes(
  betahat = betahat,
  sigma = sigma,
  numPrePeriods  = numPre,
  numPostPeriods = numPost,
  relativeMagnitudes = seq(0.5, 2, by=0.5),
  l_vec = l_vec
)

p_honest <- HonestDiD::createSensitivityPlot_relativeMagnitudes(sens)
p_honest
# ggsave("C:/.../coefplot.pdf", p_honest, width=7, height=5)

* Subset version (year==1918 or 1904 or 1932)
dt_sub <- dt[year %in% c(1918, 1904, 1932)]
m_honest2 <- feols(
  patents ~ reltimeminus14 + reltimeplus14 | year + treatmentgroup,
  data=dt_sub,
  cluster=~subclass
)

betahat2 <- coef(m_honest2)[c("reltimeminus14","reltimeplus14")]
sigma2   <- vcov(m_honest2)[c("reltimeminus14","reltimeplus14"),
                            c("reltimeminus14","reltimeplus14")]

sens2 <- HonestDiD::createSensitivityResults_relativeMagnitudes(
  betahat = betahat2,
  sigma = sigma2,
  numPrePeriods  = 1,
  numPostPeriods = 1,
  relativeMagnitudes = seq(2, 10, by=2)
)

p_honest2 <- HonestDiD::createSensitivityPlot_relativeMagnitudes(sens2)
p_honest2
# ggsave("C:/.../coefplot_subset.pdf", p_honest2, width=7, height=5)
