!pip install pandas numpy matplotlib scipy statsmodels pyfixest cvxpy
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy import stats, optimize, linalg
import warnings
import os
import time

warnings.filterwarnings("ignore")

# Optional imports
try:
    import pyfixest as pf
    HAS_PYFIXEST = True
except ImportError:
    HAS_PYFIXEST = False
    print("⚠  pyfixest not installed. pip install pyfixest")

try:
    import statsmodels.api as sm
    from statsmodels.formula.api import ols as sm_ols
    HAS_SM = True
except ImportError:
    HAS_SM = False
    print("⚠  statsmodels not installed. pip install statsmodels")

try:
    import cvxpy as cp
    HAS_CVXPY = True
except ImportError:
    HAS_CVXPY = False
    print("⚠  cvxpy not installed (needed for SC/SDID). pip install cvxpy")
    
#CHARGE DATA
path = Path('/Users/anthonyquispe/Downloads/Replication folder FE/programs/cc_xd_didtextbook_2025_9_30/Data sets/Moser and Voena 2012/moser_voena_didtextbook.dta')
df = pd.read_stata(path)
# CONFIG
DATA_FILE = "moser_voena_didtextbook.dta"
GRAPH_DIR = "graphs"
np.random.seed(1)

os.makedirs(GRAPH_DIR, exist_ok=True)
# HELPERS
def get_reltime_vars(df):
    minus_v = sorted(
        [c for c in df.columns if c.startswith("reltimeminus")],
        key=lambda x: int(x.replace("reltimeminus", "")),
    )
    plus_v = sorted(
        [c for c in df.columns if c.startswith("reltimeplus")],
        key=lambda x: int(x.replace("reltimeplus", "")),
    )
    return minus_v, plus_v


def cluster_se(X, y, resid, clusters):
    """OLS cluster-robust SEs (Stata-style)."""
    n, k = X.shape
    unique_cl = np.unique(clusters)
    G = len(unique_cl)
    bread = np.linalg.inv(X.T @ X)
    meat = np.zeros((k, k))
    for cl in unique_cl:
        idx = clusters == cl
        ei = resid[idx]
        Xi = X[idx]
        score = Xi.T @ ei
        meat += np.outer(score, score)
    scale = G / (G - 1) * (n - 1) / (n - k)
    V = bread @ meat @ bread * scale
    return np.sqrt(np.diag(V))


def ols_with_cluster(df_reg, y_var, x_vars, fe_vars, cluster_var):
    """Manual OLS with FE dummies and cluster-robust SEs."""
    df_r = df_reg.dropna(subset=[y_var] + x_vars + fe_vars + [cluster_var]).copy()
    y = df_r[y_var].values.astype(float)
    
    # Build X: x_vars + FE dummies (dropping first)
    X_parts = [df_r[x_vars].values.astype(float)]
    for fe in fe_vars:
        dummies = pd.get_dummies(df_r[fe], drop_first=True, dtype=float)
        X_parts.append(dummies.values)
    X = np.column_stack(X_parts)
    X = np.column_stack([np.ones(len(y)), X])  # intercept
    
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
        resid = y - X @ beta
    se = cluster_se(X, y, resid, df_r[cluster_var].values)
    
    # Results for x_vars only (skip intercept at pos 0)
    coefs = beta[1 : 1 + len(x_vars)]
    ses = se[1 : 1 + len(x_vars)]
    return coefs, ses


def stata_style_plot():
    """Apply Stata s2color-like styling."""
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": "black",
        "axes.grid": False,
        "axes.linewidth": 0.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "font.family": "sans-serif",
        "font.size": 10,
    })
minus_vars, plus_vars = get_reltime_vars(df)
all_reltime = minus_vars + plus_vars

#CHAPTER 3
# 1) Static TWFE regression
if HAS_PYFIXEST:
    m1 = pf.feols("patents ~ twea | subclass + year", data=df, vcov={"CRV1": "subclass"})
    print(m1.summary())
else:
    coefs, ses = ols_with_cluster(df, "patents", ["twea"], ["subclass", "year"], "subclass")
    t = coefs[0] / ses[0]
    p = 2 * (1 - stats.t.cdf(abs(t), df["subclass"].nunique() - 1))
    print(f"    twea: coef = {coefs[0]:.6f}, se = {ses[0]:.6f}, t = {t:.2f}, p = {p:.4f}")

# reghdfe equivalence
print("\n    reghdfe patents twea, absorb(subclass year) cluster(subclass)")
if HAS_PYFIXEST:
    m1b = pf.feols("patents ~ twea | subclass + year", data=df, vcov={"CRV1": "subclass"})
    print(f"    twea: {m1b.coef()['twea']:.6f}  (same as above)")

# 2) Equivalence between static TWFE regression and DID
print("\n*** 2) Equivalence between static TWFE and DID")
print("    reg patents treatmentgroup post twea, cluster(subclass)")

if HAS_PYFIXEST:
    m2 = pf.feols("patents ~ treatmentgroup + post + twea", data=df, vcov={"CRV1": "subclass"})
    print(m2.summary())
else:
    coefs, ses = ols_with_cluster(df, "patents", ["treatmentgroup", "post", "twea"], [], "subclass")
    for i, v in enumerate(["treatmentgroup", "post", "twea"]):
        print(f"    {v}: coef = {coefs[i]:.6f}, se = {ses[i]:.6f}")

# 3) Testing randomized treatment
print("\n*** 3) Testing randomized treatment")
print("    reg patents treatmentgroup if year<=1918, cluster(subclass)")

df_pre = df[df["year"] <= 1918].copy()
if HAS_PYFIXEST:
    m3 = pf.feols("patents ~ treatmentgroup", data=df_pre, vcov={"CRV1": "subclass"})
    print(m3.summary())
else:
    coefs, ses = ols_with_cluster(df_pre, "patents", ["treatmentgroup"], [], "subclass")
    print(f"    treatmentgroup: coef = {coefs[0]:.6f}, se = {ses[0]:.6f}")
    
# 4) Event-study TWFE regression
print("\n*** 4) Event-study TWFE regression")
print("    reg patents i.year treatmentgroup reltimeminus* reltimeplus*, cluster(subclass)")
if HAS_PYFIXEST:
    fml_es = "patents ~ " + " + ".join(all_reltime) + " + treatmentgroup | year"
    m4 = pf.feols(fml_es, data=df, vcov={"CRV1": "subclass"})
    print(m4.summary())
else:
    coefs, ses = ols_with_cluster(
        df, "patents", all_reltime + ["treatmentgroup"], ["year"], "subclass"
    )
    for i, v in enumerate(all_reltime):
        print(f"    {v}: {coefs[i]:.6f} ({ses[i]:.6f})")

# F-test sobre coeficientes pre-tendencia
print("\n    F-test on all reltimeminus coefficients:")
if HAS_PYFIXEST:
    c = m4.coef()
    v = pd.DataFrame(m4._vcov, index=c.index, columns=c.index)
    pre_names = [v_ for v_ in minus_vars if v_ in c.index]
    beta_pre = c[pre_names].values
    V_pre = v.loc[pre_names, pre_names].values
    F_stat = beta_pre @ np.linalg.inv(V_pre) @ beta_pre / len(pre_names)
    p_val = 1 - stats.f.cdf(F_stat, len(pre_names), df["subclass"].nunique() - 1)
    print(f"    F({len(pre_names)}, {df['subclass'].nunique()-1}) = {F_stat:.4f}")
    print(f"    Prob > F = {p_val:.4f}")
#  Event-study plot (Graph 1: graphES_moser1) 
print("\n    Generating graphES_moser1...")
stata_style_plot()

if HAS_PYFIXEST:
    fml_plot = "patents ~ " + " + ".join(all_reltime) + " + treatmentgroup | year"
    m4p = pf.feols(fml_plot, data=df, vcov={"CRV1": "subclass"})
    coefs_d = m4p.coef()
    se_d = m4p.se()

    reltime = [0]
    est = [0.0]
    ci_lo = [0.0]
    ci_hi = [0.0]

    for v in minus_vars:
        k = -int(v.replace("reltimeminus", ""))
        if v in coefs_d.index:
            reltime.append(k)
            est.append(coefs_d[v])
            ci_lo.append(coefs_d[v] - 1.96 * se_d[v])
            ci_hi.append(coefs_d[v] + 1.96 * se_d[v])
    for v in plus_vars:
        k = int(v.replace("reltimeplus", ""))
        if v in coefs_d.index:
            reltime.append(k)
            est.append(coefs_d[v])
            ci_lo.append(coefs_d[v] - 1.96 * se_d[v])
            ci_hi.append(coefs_d[v] + 1.96 * se_d[v])

    order = np.argsort(reltime)
    reltime = np.array(reltime)[order]
    est = np.array(est)[order]
    ci_lo = np.array(ci_lo)[order]
    ci_hi = np.array(ci_hi)[order]

    # Store post-period results for later comparison (Chapter 4)
    mask_post = reltime >= 0
    res_post = {
        "reltime": reltime[mask_post],
        "est": est[mask_post],
        "ci_lo": ci_lo[mask_post],
        "ci_hi": ci_hi[mask_post],
    }
    
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(reltime, est, "o-", color="navy", markersize=5, linewidth=1)
    ax.vlines(reltime, ci_lo, ci_hi, color="maroon", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.5, linestyle="-")
    ax.set_title("TWFE Event-study estimates")
    ax.set_xlabel("Relative time to year before TWEA")
    ax.set_ylabel("Effect")
    ax.set_xticks(range(-18, 22, 3))
    ax.set_ylim(-0.25, 1)
    ax.set_yticks(np.arange(-0.25, 1.01, 0.25))
    plt.tight_layout()
    plt.savefig(os.path.join(GRAPH_DIR, "graphES_moser1.pdf"))
    print("    ✓ Saved graphs/graphES_moser1.pdf")
    plt.show()

#5) Event-study without pre-trends estimates
print("\n*** 5) Event-study TWFE without pre-trends estimates")
print("    reg patents i.yearpost treatmentgroup reltimeplus*, cluster(subclass)")

if HAS_PYFIXEST:
    fml_nopretrend = "patents ~ " + " + ".join(plus_vars) + " + treatmentgroup | yearpost"
    m5 = pf.feols(fml_nopretrend, data=df, vcov={"CRV1": "subclass"})
    print(m5.summary())

    # Build results for "without pre-periods"
    coefs_5 = m5.coef()
    se_5 = m5.se()
    rt_nopre = [0]
    est_nopre = [0.0]
    ci_lo_nopre = [0.0]
    ci_hi_nopre = [0.0]
    for v in plus_vars:
        k = int(v.replace("reltimeplus", ""))
        if v in coefs_5.index:
            rt_nopre.append(k)
            est_nopre.append(coefs_5[v])
            ci_lo_nopre.append(coefs_5[v] - 1.96 * se_5[v])
            ci_hi_nopre.append(coefs_5[v] + 1.96 * se_5[v])

    order2 = np.argsort(rt_nopre)
    rt_nopre = np.array(rt_nopre)[order2]
    est_nopre = np.array(est_nopre)[order2]
    ci_lo_nopre = np.array(ci_lo_nopre)[order2]
    ci_hi_nopre = np.array(ci_hi_nopre)[order2]

    # Graph 2: graphES_moser2 
    print("\n    Generating graphES_moser2...")
    fig, ax = plt.subplots(figsize=(8, 5))

    # Without pre-periods (blue)
    ax.plot(rt_nopre, est_nopre, "o-", color="steelblue", markersize=4, linewidth=1,
            label="Without Pre-Periods")
    ax.vlines(rt_nopre, ci_lo_nopre, ci_hi_nopre, color="steelblue", linewidth=1)

    # With pre-periods (red) — from res_post
    ax.plot(res_post["reltime"], res_post["est"], "o-", color="red", markersize=4,
            linewidth=1, label="With Pre-Periods")
    ax.vlines(res_post["reltime"], res_post["ci_lo"], res_post["ci_hi"],
              color="red", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_title("TWFE Event-study estimates")
    ax.set_xlabel("Relative time to year before TWEA")
    ax.set_ylabel("Effect")
    ax.set_xticks(range(0, 22, 3))
    ax.set_ylim(-0.25, 1)
    ax.set_yticks(np.arange(-0.25, 1.01, 0.25))
    ax.legend(loc="lower center", ncol=2, frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(GRAPH_DIR, "graphES_moser2.pdf"))
    
    print("    ✓ Saved graphs/graphES_moser2.pdf")
    plt.show()

    # Numerically equivalent to Borusyak et al / having all year FEs
    print("\n    Note: Numerically equivalent to Borusyak et al (did_imputation)")
    print("    Also equivalent to: reg patents i.year treatmentgroup reltimeplus*")

# 6) Linear pre-trends we could fail to detect (pretrends power)
print("\n*** 6) Linear pre-trends we could fail to detect")

if HAS_PYFIXEST:
    # Restrict to 6 pre-periods (year >= 1912 & year <= 1939)
    df_sub = df[(df["year"] >= 1912) & (df["year"] <= 1939)].copy()
    short_minus = [f"reltimeminus{i}" for i in range(1, 7)]
    short_minus_present = [v for v in short_minus if v in df_sub.columns]

    # reghdfe patents reltimeminus* reltimeplus*, absorb(treatmentgroup year) cluster(subclass)
    fml_pre = "patents ~ " + " + ".join(short_minus_present + plus_vars) + " | treatmentgroup + year"
    m6_full = pf.feols(fml_pre, data=df, vcov={"CRV1": "subclass"})

    # Pretrends power calculation (Roth 2022)
    coefs_pre = m6_full.coef()
    vcov_pre = pd.DataFrame(m6_full._vcov, index=coefs_pre.index, columns=coefs_pre.index)

    # Extract pre-trend coefficients (first 6)
    pre_names_6 = [v for v in short_minus_present if v in coefs_pre.index][:6]
    beta_pre6 = coefs_pre[pre_names_6].values
    V_pre6 = vcov_pre.loc[pre_names_6, pre_names_6].values

    # Pre-trends power: slope that gives 50% power for the joint F-test
    n_pre = len(pre_names_6)
    dfn = n_pre
    dfd = df["subclass"].nunique() - 1

    def power_at_slope(slope):
        """Power of F-test under linear trend with given slope."""
        mu = np.array([slope * (-int(v.replace("reltimeminus", ""))) for v in pre_names_6])
        ncp = mu @ np.linalg.inv(V_pre6) @ mu
        cv = stats.f.ppf(0.95, dfn, dfd)
        power = 1 - stats.ncf.cdf(cv, dfn, dfd, ncp)
        return power

    # Search for slope giving 50% power
    try:
        result = optimize.brentq(lambda s: power_at_slope(s) - 0.5, 0.0001, 1.0)
        slope_50 = result
        print(f"    Slope for 50% power: {slope_50:.6f}")
    except Exception:
        slope_50 = None
        print("    Could not find slope for 50% power (search failed)")
 # Event-study on restricted sample for the graph
    fml_sub = "patents ~ " + " + ".join(short_minus_present + plus_vars) + " + treatmentgroup | year"
    m6_sub = pf.feols(fml_sub, data=df_sub, vcov={"CRV1": "subclass"})
    coefs_6 = m6_sub.coef()
    se_6 = m6_sub.se()

    rt6 = [0]
    est6 = [0.0]
    ci_lo6 = [0.0]
    ci_hi6 = [0.0]

    for v in short_minus_present:
        k = -int(v.replace("reltimeminus", ""))
        if v in coefs_6.index:
            rt6.append(k)
            est6.append(coefs_6[v])
            ci_lo6.append(coefs_6[v] - 1.96 * se_6[v])
            ci_hi6.append(coefs_6[v] + 1.96 * se_6[v])
    for v in plus_vars:
        k = int(v.replace("reltimeplus", ""))
        if v in coefs_6.index:
            rt6.append(k)
            est6.append(coefs_6[v])
            ci_lo6.append(coefs_6[v] - 1.96 * se_6[v])
            ci_hi6.append(coefs_6[v] + 1.96 * se_6[v])

    order6 = np.argsort(rt6)
    rt6 = np.array(rt6)[order6]
    est6 = np.array(est6)[order6]
    ci_lo6 = np.array(ci_lo6)[order6]
    ci_hi6 = np.array(ci_hi6)[order6]

    #  Graph 3: graphES_moser3 
    print("    Generating graphES_moser3...")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(rt6, est6, "o-", color="navy", markersize=5, linewidth=1)
    ax.vlines(rt6, ci_lo6, ci_hi6, color="maroon", linewidth=1)
    if slope_50 is not None:
        x_line = np.linspace(-6, 21, 100)
        ax.plot(x_line, x_line * slope_50, "--", color="gray", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_title("TWFE Event-study estimates")
    ax.set_xlabel("Relative time to year before TWEA")
    ax.set_ylabel("Effect")
    ax.set_xticks(range(-6, 22, 3))
    ax.set_ylim(-0.25, 1)
    ax.set_yticks(np.arange(-0.25, 1.01, 0.25))
    plt.tight_layout()
    plt.savefig(os.path.join(GRAPH_DIR, "graphES_moser3.pdf"))
    print("    ✓ Saved graphs/graphES_moser3.pdf")
    plt.show()

# 7) Variance of the effect of having been exposed to treatment for 14 years
df_1932 = df[df["year"] == 1932].copy()
g0 = df_1932.loc[df_1932["treatmentgroup"] == 0, "diffpatentswrt1918"].dropna()
g1 = df_1932.loc[df_1932["treatmentgroup"] == 1, "diffpatentswrt1918"].dropna()

# ─── 1. sdtest diffpatentswrt1918 if year==1932, by(treatmentgroup) ───
n0, n1 = len(g0), len(g1)
mean0, mean1 = g0.mean(), g1.mean()
sd0, sd1 = g0.std(ddof=1), g1.std(ddof=1)
se0, se1 = sd0 / np.sqrt(n0), sd1 / np.sqrt(n1)
ci_lo0, ci_hi0 = mean0 - 1.96 * se0, mean0 + 1.96 * se0
ci_lo1, ci_hi1 = mean1 - 1.96 * se1, mean1 + 1.96 * se1

combined = pd.concat([g0, g1])
n_c = len(combined)
mean_c = combined.mean()
sd_c = combined.std(ddof=1)
se_c = sd_c / np.sqrt(n_c)
ci_loc, ci_hic = mean_c - 1.96 * se_c, mean_c + 1.96 * se_c

print("    Variance ratio test")
print("    " + "-" * 70)
print(f"       Group |     Obs        Mean    Std. err.   Std. dev.   [95% conf. interval]")
print("    " + "-" * 70)
print(f"           0 | {n0:>7,}    {mean0:.7f}    {se0:.7f}    {sd0:.6f}    {ci_lo0:.7f}    {ci_hi0:.7f}")
print(f"           1 | {n1:>7,}    {mean1:.7f}    {se1:.7f}    {sd1:.6f}    {ci_lo1:.7f}    {ci_hi1:.7f}")
print("    " + "-" * 70)
print(f"    Combined | {n_c:>7,}    {mean_c:.7f}    {se_c:.7f}    {sd_c:.6f}    {ci_loc:.7f}    {ci_hic:.7f}")
print("    " + "-" * 70)

F_var = sd0**2 / sd1**2
df1_f, df2_f = n0 - 1, n1 - 1
p_left = stats.f.cdf(F_var, df1_f, df2_f)
p_right = 1 - p_left
p_two = 2 * min(p_left, p_right)
print(f"        ratio = sd(0) / sd(1)                          f = {F_var:.4f}")
print(f"    H0: ratio = 1                         Degrees of freedom = {df1_f}, {df2_f}")
print(f"      Pr(F < f) = {p_left:.4f}         2*Pr(F < f) = {p_two:.4f}           Pr(F > f) = {p_right:.4f}")

#  2. di r(sd_2)-r(sd_1) 
sd_diff = sd1 - sd0
print(f"\n    r(sd_2) - r(sd_1) = {sd_diff:.8f}")

#  3. scalar sd_effects=r(sd_2)-r(sd_1) 
sd_effects = sd_diff
print(f"    scalar sd_effects = {sd_effects:.8f}")

#  4. reg diffpatentswrt1918 treatmentgroup if year==1932 
if HAS_PYFIXEST:
    m7 = pf.feols("diffpatentswrt1918 ~ treatmentgroup", data=df_1932)
    print(m7.summary())
    beta_did = m7.coef()["treatmentgroup"]
else:
    import statsmodels.api as sm
    X = sm.add_constant(df_1932["treatmentgroup"])
    m7 = sm.OLS(df_1932["diffpatentswrt1918"], X).fit()
    print(m7.summary())
    beta_did = m7.params["treatmentgroup"]

#  di _b[treatmentgroup]-1.96*sd_effects, _b[treatmentgroup]+1.96*sd_effects 
ci_lo_sd = beta_did - 1.96 * sd_effects
ci_hi_sd = beta_did + 1.96 * sd_effects
print(f"\n    _b[treatmentgroup] = {beta_did:.7f}")
print(f"    CI using SD of effects: [{ci_lo_sd:.8f}, {ci_hi_sd:.8f}]")

# 8) Placebo test of variance assumptions
print("\n*** 8) Placebo test (year=1904)")

#  1. sdtest diffpatentswrt1918 if year==1904, by(treatmentgroup) 
df_1904 = df[df["year"] == 1904].copy()
g0p = df_1904.loc[df_1904["treatmentgroup"] == 0, "diffpatentswrt1918"].dropna()
g1p = df_1904.loc[df_1904["treatmentgroup"] == 1, "diffpatentswrt1918"].dropna()

def print_sdtest(g0, g1, label=""):
    n0, n1 = len(g0), len(g1)
    mean0, mean1 = g0.mean(), g1.mean()
    sd0, sd1 = g0.std(ddof=1), g1.std(ddof=1)
    se0, se1 = sd0 / np.sqrt(n0), sd1 / np.sqrt(n1)
    
    combined = pd.concat([g0, g1])
    n_c = len(combined)
    mean_c = combined.mean()
    sd_c = combined.std(ddof=1)
    se_c = sd_c / np.sqrt(n_c)
    
    # Intervalos de confianza
    from scipy import stats as st
    ci_lo0, ci_hi0 = st.t.interval(0.95, n0 - 1, loc=mean0, scale=se0)
    ci_lo1, ci_hi1 = st.t.interval(0.95, n1 - 1, loc=mean1, scale=se1)
    ci_loc, ci_hic = st.t.interval(0.95, n_c - 1, loc=mean_c, scale=se_c)
    
    print(f"    {label}")
    print("    Variance ratio test")
    print("    " + "-" * 75)
    print(f"       Group |     Obs        Mean    Std. err.   Std. dev.   [95% conf. interval]")
    print("    " + "-" * 75)
    print(f"           0 | {n0:>7,}   {mean0:>10.7f}   {se0:>10.7f}   {sd0:>10.6f}   {ci_lo0:>11.7f}   {ci_hi0:>11.7f}")
    print(f"           1 | {n1:>7,}   {mean1:>10.7f}   {se1:>10.7f}   {sd1:>10.6f}   {ci_lo1:>11.7f}   {ci_hi1:>11.7f}")
    print("    " + "-" * 75)
    print(f"    Combined | {n_c:>7,}   {mean_c:>10.7f}   {se_c:>10.7f}   {sd_c:>10.6f}   {ci_loc:>11.7f}   {ci_hic:>11.7f}")
    print("    " + "-" * 75)
    
    if sd0 > 0 and sd1 > 0:
        F_var = sd0**2 / sd1**2
        df1_f, df2_f = n0 - 1, n1 - 1
        p_left = st.f.cdf(F_var, df1_f, df2_f)
        p_right = 1 - p_left
        if F_var < 1:
            p_two = 2 * p_left
            else:
            p_two = 2 * p_right
        print(f"        ratio = sd(0) / sd(1)                          f = {F_var:.4f}")
        print(f"    H0: ratio = 1                         Degrees of freedom = {df1_f}, {df2_f}")
        print(f"      Pr(F < f) = {p_left:.4f}         2*Pr(F {'<' if F_var < 1 else '>'} f) = {p_two:.4f}           Pr(F > f) = {p_right:.4f}")
    else:
        print(f"        ratio = sd(0) / sd(1)                          f =     .")
        print(f"    H0: ratio = 1                         Degrees of freedom = {n0-1}, {n1-1}")
        print(f"      Pr(F < f) =     .         2*Pr(F > f) =     .           Pr(F > f) =     .")
    print()
    
    return sd0, sd1, F_var if (sd0 > 0 and sd1 > 0) else None

# sdtest for 1904
print_sdtest(g0p, g1p, label="sdtest diffpatentswrt1918 if year==1904, by(treatmentgroup)")

#  2. forvalue i=1900/1939 { sdtest ... } 
print("\n    forvalue i=1900/1939 { sdtest diffpatentswrt1918, by(treatmentgroup) }")
print("=" * 80)
for yr in range(1900, 1940):
    dfy = df[df["year"] == yr]
    g0y = dfy.loc[dfy["treatmentgroup"] == 0, "diffpatentswrt1918"].dropna()
    g1y = dfy.loc[dfy["treatmentgroup"] == 1, "diffpatentswrt1918"].dropna()
    if len(g0y) > 0 and len(g1y) > 0:
        print_sdtest(g0y, g1y, label=f"year == {yr}")
    else:
        print(f"    year == {yr}: No observations in one or both groups")

#CHAPTER 4
# 4.1) Estimators with controls
print("\n*** 4.1) TWFE controlling for baseline patents")
print("    reghdfe patents reltimeminus* reltimeplus*, absorb(year#patents1900 treatmentgroup) cluster(subclass)")

if HAS_PYFIXEST:
    # year#patents1900 interaction
    df["year_x_pat1900"] = df["year"].astype(str) + "_" + df["patents1900"].astype(str)
    fml_ctrl = "patents ~ " + " + ".join(all_reltime) + " | year_x_pat1900 + treatmentgroup"
    m_ctrl = pf.feols(fml_ctrl, data=df, vcov={"CRV1": "subclass"})
    print(m_ctrl.summary())

    # F-test on pre-trend coefficients
    c_ctrl = m_ctrl.coef()
    v_ctrl = pd.DataFrame(m_ctrl._vcov, index=c_ctrl.index, columns=c_ctrl.index)
    pre_in = [v for v in minus_vars if v in c_ctrl.index]
    beta_ctrl_pre = c_ctrl[pre_in].values
    V_ctrl_pre = v_ctrl.loc[pre_in, pre_in].values
    F_ctrl = beta_ctrl_pre @ np.linalg.inv(V_ctrl_pre) @ beta_ctrl_pre / len(pre_in)
    p_ctrl = 1 - stats.f.cdf(F_ctrl, len(pre_in), df["subclass"].nunique() - 1)
    print(f"\n    F-test on pre-trends: F = {F_ctrl:.4f}, p = {p_ctrl:.4f}")

    # ─── Graph 4: graphES_moser_controls ─────────────────────────────────
    print("    Generating graphES_moser_controls...")
    se_ctrl = m_ctrl.se()

    rt_c = [0]
    est_c = [0.0]
    ci_lo_c = [0.0]
    ci_hi_c = [0.0]

    for v in minus_vars:
        k = -int(v.replace("reltimeminus", ""))
        if v in c_ctrl.index:
            rt_c.append(k)
            est_c.append(c_ctrl[v])
            ci_lo_c.append(c_ctrl[v] - 1.96 * se_ctrl[v])
            ci_hi_c.append(c_ctrl[v] + 1.96 * se_ctrl[v])
    for v in plus_vars:
        k = int(v.replace("reltimeplus", ""))
        if v in c_ctrl.index:
            rt_c.append(k)
            est_c.append(c_ctrl[v])
            ci_lo_c.append(c_ctrl[v] - 1.96 * se_ctrl[v])
            ci_hi_c.append(c_ctrl[v] + 1.96 * se_ctrl[v])

    order_c = np.argsort(rt_c)
    rt_c = np.array(rt_c)[order_c]
    est_c = np.array(est_c)[order_c]
    ci_lo_c = np.array(ci_lo_c)[order_c]
    ci_hi_c = np.array(ci_hi_c)[order_c]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(rt_c, est_c, "o-", color="navy", markersize=5, linewidth=1)
    ax.vlines(rt_c, ci_lo_c, ci_hi_c, color="maroon", linewidth=1)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_title("TWFE Event-study estimates")
    ax.set_xlabel("Relative time to year before TWEA")
    ax.set_ylabel("Effect")
    ax.set_xticks(range(-18, 22, 3))
    ax.set_ylim(-0.25, 1)
    ax.set_yticks(np.arange(-0.25, 1.01, 0.25))
    plt.tight_layout()
    plt.savefig(os.path.join(GRAPH_DIR, "graphES_moser_controls.pdf"))
    plt.show()
    print("    ✓ Saved graphs/graphES_moser_controls.pdf")

    # Testing correlation between treatment and covariate
    print("\n    Testing correlation: reg patents1900 treatmentgroup if year==1900")
    df_1900 = df[df["year"] == 1900].copy()
    m_corr = pf.feols("patents1900 ~ treatmentgroup", data=df_1900)
    print(f"    treatmentgroup: {m_corr.coef()['treatmentgroup']:.6f}")

#Controlling for patents in 1900, DID (not on figure)
import polars as pl
from did_multiplegt_dyn import DidMultiplegtDyn
import matplotlib.pyplot as plt
import numpy as np

# Convertir a Polars
df_pl = pl.from_pandas(df)

# did_multiplegt_dyn patents subclass year twea, effects(21) placebo(18) trends_nonparam(patents1900)
model = DidMultiplegtDyn(
    df=df_pl,
    outcome="patents",
    group="subclass",
    time="year",
    treatment="twea",
    effects=21,
    placebo=18,
    trends_nonparam=["patents1900"],
)

model.fit()
model.summary()

# Extraer resultados
res = model.result["did_multiplegt_dyn"]
effects_df = res["Effects"]
placebos_df = res["Placebos"]
ate_df = res["ATE"]
p_joint_effects = res["p_jointeffects"]
p_joint_placebo = res["p_jointplacebo"]

print(f"\n    Test of joint nullity of the effects : p-value = {p_joint_effects:.3e}")
print(f"    Test of joint nullity of the placebos : p-value = {p_joint_placebo:.3e}")

# Extraer estimaciones y SE
n_effects = len(effects_df)
n_placebo = len(placebos_df)

eff_est = effects_df["Estimate"].values
eff_se = effects_df["SE"].values
plac_est = placebos_df["Estimate"].values
plac_se = placebos_df["SE"].values

# Eje x: placebos (negativos, invertidos), 0 como referencia, efectos (positivos)
x_plac = list(range(-n_placebo, 0))[::-1]  # -1, -2, ..., -18 -> invertir a -18, ..., -1
x_eff = list(range(1, n_effects + 1))

x = x_plac + [0] + x_eff
y = list(plac_est[::-1]) + [0.0] + list(eff_est)
se = list(plac_se[::-1]) + [0.0] + list(eff_se)

x = np.array(x)
y = np.array(y)
se = np.array(se)
ci_lo = y - 1.96 * se
ci_hi = y + 1.96 * se

# Gráfico estilo Stata
fig, ax = plt.subplots(figsize=(14, 6))

# Placebos: -18 a -1 (invertir orden para que -18 quede a la izquierda)
x_plac = list(range(-n_placebo, 0))  # [-18, -17, ..., -1]
y_plac = list(plac_est[::-1])  # invertir: Placebo_18 en -18, Placebo_1 en -1
se_plac = list(plac_se[::-1])

# Punto de referencia en 0
x_ref = [0]
y_ref = [0.0]
se_ref = [0.0]

# Efectos: 1 a 21
x_eff = list(range(1, n_effects + 1))
y_eff = list(eff_est)
se_eff = list(eff_se)

# Combinar
x = np.array(x_plac + x_ref + x_eff)
y = np.array(y_plac + y_ref + y_eff)
se = np.array(se_plac + se_ref + se_eff)
ci_lo = y - 1.96 * se
ci_hi = y + 1.96 * se

# Línea + puntos azules
ax.plot(x, y, "o-", color="#4472C4", markersize=5, linewidth=1.2, zorder=3)

# Barras de error rojas (whisker style como Stata)
ax.vlines(x, ci_lo, ci_hi, color="#C00000", linewidth=1, zorder=2)
# Caps en los extremos
cap_w = 0.25
for xi, lo, hi in zip(x, ci_lo, ci_hi):
    ax.hlines(lo, xi - cap_w, xi + cap_w, color="#C00000", linewidth=1, zorder=2)
    ax.hlines(hi, xi - cap_w, xi + cap_w, color="#C00000", linewidth=1, zorder=2)

# Línea horizontal en 0
ax.axhline(0, color="black", linewidth=0.5)

# Grid horizontal con líneas punteadas grises
ax.yaxis.grid(True, linestyle="--", alpha=0.5, color="gray")
ax.set_axisbelow(True)

# Título y ejes como Stata
ax.set_title("DID, from last period before treatment changes (t=0) to t", fontsize=13)
ax.set_xlabel("Relative time to last period before treatment changes (t=0)", fontsize=11)
ax.set_ylabel("")

# Eje x: todas las etiquetas de -18 a 21
ax.set_xticks(range(-18, 22))
ax.set_xticklabels([str(i) for i in range(-18, 22)], fontsize=7)

# Eje y como Stata
ax.set_ylim(-0.5, 1.1)
ax.set_yticks(np.arange(-0.5, 1.01, 0.25))

# Bordes
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(os.path.join(GRAPH_DIR, "graphES_did_multiplegt_dyn.pdf"))
plt.show()
print("    ✓ Saved graphs/graphES_did_multiplegt_dyn.pdf")
#4.2)Interactive fixed effects
def _ife_em(Y_mat, obs_mask, r, force="two-way", tol=1e-4, max_iter=5000):
    """
    Fit Interactive Fixed Effects model via EM on a panel with missing entries.
    
    Model: Y_it = mu_i + xi_t + lambda_i' f_t + e_it
    
    EM Algorithm:
    1. Initialize: impute missing with row/col means
    2. E-step: fit IFE (demean + SVD) on completed matrix
    3. M-step: re-impute missing entries with fitted values
    4. Repeat until convergence
    
    Args:
        Y_mat: (N, T) outcome matrix (missing entries can be anything)
        obs_mask: (N, T) boolean, True = observed, False = missing/to-impute
        r: number of factors
        force: "two-way", "unit", "time", or "none"
        tol: convergence tolerance
        max_iter: maximum iterations
    
    Returns: mu (N,), xi (T,), Lambda (N,r), F_hat (T,r), Y_fit (N,T)
    """
    N, T = Y_mat.shape
    Y = Y_mat.copy().astype(float)
    W = obs_mask.astype(float)
    
    # --- Initialize missing entries ---
    # Use row/col means of observed entries for initialization
    for i in range(N):
        obs_i = W[i, :] > 0
        if obs_i.any() and not obs_i.all():
            Y[i, ~obs_i] = Y[i, obs_i].mean()
    for t in range(T):
        obs_t = W[:, t] > 0
        if obs_t.any():
            missing_t = W[:, t] == 0
            if missing_t.any():
                # Average of row-mean init and col-mean
                col_mean = Y[obs_t, t].mean()
                Y[missing_t, t] = (Y[missing_t, t] + col_mean) / 2
    
    # Grand mean for any remaining NaN
    gm_init = np.nanmean(Y[W > 0]) if (W > 0).any() else 0
    Y = np.where(np.isnan(Y), gm_init, Y)
    
    #  Iterative EM 
    Y_fit_old = np.zeros((N, T))
    
    for iteration in range(max_iter):
        # --- Fit IFE on completed matrix ---
        # Compute fixed effects
        if force == "two-way":
            mu = Y.mean(axis=1)
            xi = Y.mean(axis=0)
            gm = Y.mean()
            Y_dm = Y - mu[:, None] - xi[None, :] + gm
        elif force == "unit":
            mu = Y.mean(axis=1)
            xi = np.zeros(T)
            gm = 0.0
            Y_dm = Y - mu[:, None]
        elif force == "time":
            mu = np.zeros(N)
            xi = Y.mean(axis=0)
            gm = 0.0
            Y_dm = Y - xi[None, :]
        else:
            mu = np.zeros(N)
            xi = np.zeros(T)
            gm = 0.0
            Y_dm = Y.copy()
        
        # SVD for factors
        if r > 0:
            r_use = min(r, min(N, T) - 1)
            U, S, Vt = np.linalg.svd(Y_dm, full_matrices=False)
            F_hat = Vt[:r_use, :].T * np.sqrt(T)
            Lambda = U[:, :r_use] * S[:r_use] / np.sqrt(T)
            interactive = Lambda @ F_hat.T
        else:
            F_hat = None
            Lambda = None
            interactive = np.zeros((N, T))
        
        # Fitted values
        if force == "two-way":
            Y_fit = mu[:, None] + xi[None, :] - gm + interactive
        elif force == "unit":
            Y_fit = mu[:, None] + interactive
        elif force == "time":
            Y_fit = xi[None, :] + interactive
        else:
            Y_fit = interactive.copy()
        
        #  Check convergence 
        change = np.mean((Y_fit - Y_fit_old)**2)
        if iteration > 0 and change < tol:
            break
        Y_fit_old = Y_fit.copy()
        
        #  M-step: re-impute missing entries 
        Y = np.where(W > 0, Y_mat, Y_fit)
    
    return mu, xi, Lambda, F_hat, Y_fit


def _ife_em_counterfactual(Y_full, D_mat, r, force="two-way", tol=1e-4):
    """
    IFE counterfactual estimation via EM (fect method="ife").
    
    Treats D_it=1 entries as missing, runs EM-IFE to impute counterfactuals.
    
    Args:
        Y_full: (N, T) full outcome matrix
        D_mat: (N, T) treatment indicator (1=treated, 0=control)
        r: number of factors
        force, tol: passed to _ife_em
    
    Returns:
        Y_fit: (N, T) fitted/counterfactual values
        eff: (N, T) treatment effects (Y - Y_fit), meaningful only where D=1
    """
    obs_mask = (D_mat == 0).astype(float)  # observed = untreated entries
    mu, xi, Lam, F_hat, Y_fit = _ife_em(Y_full, obs_mask, r, force, tol)
    eff = Y_full - Y_fit
    return Y_fit, eff

# CROSS-VALIDATION (leave-one-out by time period, matching fect)
def ife_cross_validation(df, y_var, treat_var, unit_var, time_var,
                          r_max=4, tol=1e-4, cvtreat=False, seed=None,
                          force="two-way"):
    """
    Cross-validation for optimal number of factors.
    
    fect CV algorithm (leave-one-out by time period):
    
    cvtreat=False (CV on controls):
      For each pre-treatment period t:
        * ADDITIONALLY mask period t for all control units
        * Run EM-IFE (post-treatment treated already masked + now control at t masked)
        * Predict control outcomes at period t
        * MSPE = sum of squared prediction errors / count
    
    cvtreat=True (CV on treated):
      For each pre-treatment period t:
        * ADDITIONALLY mask period t for treated units (on top of post-treatment mask)
        * Run EM-IFE
        * Predict treated outcomes at period t
        * MSPE = sum of squared prediction errors / count
    """
    if seed is not None:
        np.random.seed(seed)
    
    times = sorted(df[time_var].unique())
    T_all = len(times)
    
    panel = df.pivot_table(index=unit_var, columns=time_var, values=y_var, aggfunc="mean")
    treat_panel = df.pivot_table(index=unit_var, columns=time_var, values=treat_var, aggfunc="max")
    
    # Ensure consistent ordering
    panel = panel.reindex(columns=times)
    treat_panel = treat_panel.reindex(columns=times)
    
    units = panel.index
    N = len(units)
    
    Y_full = panel.values.astype(float)  # N x T
    D_mat = treat_panel.values.astype(float)  # N x T
    
    ever_treated = D_mat.max(axis=1) > 0
    ctrl_idx = np.where(~ever_treated)[0]
    treat_idx = np.where(ever_treated)[0]
    
    # Pre-treatment periods
    pre_mask = D_mat.max(axis=0) == 0  # T-length boolean
    pre_period_indices = np.where(pre_mask)[0]
    
    label_cv = " (on treated units)" if cvtreat else ""
    print(f"    Cross Validation{label_cv}...")
    
    cv_results = {}
    
    for r in range(r_max + 1):
        sum_e2 = 0.0
        num_y = 0
        
        for lv_t in pre_period_indices:
            # Base observation mask: observed = untreated entries
            obs_mask = (D_mat == 0).astype(float)
            
            if not cvtreat:
                # CV on controls: additionally mask control units at period lv_t
                obs_mask[ctrl_idx, lv_t] = 0.0
                
                # Run EM-IFE
                try:
                    _, _, _, _, Y_fit = _ife_em(Y_full, obs_mask, r, force, tol)
                    # Prediction error for controls at left-out period
                    e = Y_full[ctrl_idx, lv_t] - Y_fit[ctrl_idx, lv_t]
                    sum_e2 += np.sum(e**2)
                    num_y += len(e)
                except Exception:
                    pass
            else:
                # CV on treated: additionally mask treated units at period lv_t
                obs_mask[treat_idx, lv_t] = 0.0
                
                # Run EM-IFE
                try:
                    _, _, _, _, Y_fit = _ife_em(Y_full, obs_mask, r, force, tol)
                    # Prediction error for treated at left-out period
                    e = Y_full[treat_idx, lv_t] - Y_fit[treat_idx, lv_t]
                    sum_e2 += np.sum(e**2)
                    num_y += len(e)
                except Exception:
                    pass
        
        mspe = sum_e2 / max(num_y, 1)
        label = "fe" if r == 0 else "ife"
        print(f"    {label} r={r} force=two-way mspe={mspe:.4f}")
        cv_results[r] = mspe
    
    optimal_r = min(cv_results, key=cv_results.get)
    print(f"    optimal r={optimal_r} in fe/ife model")
    return optimal_r, cv_results

# MAIN IFE ESTIMATION WITH BOOTSTRAP
def ife_estimator(df, y_var, treat_var, unit_var, time_var,
                   n_factors=2, tol=1e-4, n_boot=200, seed=None,
                   force="two-way"):
    """
    IFE estimator with bootstrap standard errors (matching fect).
    
    Replicates: fect Y, treat(D) unit(G) time(T) method("ife") r(n_factors) tol(tol) se
    
    Returns dict with ATT (aggregate), ATTs (by period), SEs, CIs.
    """
    if seed is not None:
        np.random.seed(seed)
    
    t0 = time.time()
    times = sorted(df[time_var].unique())
    T_all = len(times)
    
    panel = df.pivot_table(index=unit_var, columns=time_var, values=y_var, aggfunc="mean")
    treat_panel = df.pivot_table(index=unit_var, columns=time_var, values=treat_var, aggfunc="max")
    
    panel = panel.reindex(columns=times)
    treat_panel = treat_panel.reindex(columns=times)
    
    units = panel.index
    Y_full = panel.values.astype(float)
    D_mat = treat_panel.values.astype(float)
    
    ever_treated = D_mat.max(axis=1) > 0
    ctrl_idx = np.where(~ever_treated)[0]
    treat_idx = np.where(ever_treated)[0]
    
    N_ctrl = len(ctrl_idx)
    N_treat = len(treat_idx)
    
    pre_mask = D_mat.max(axis=0) == 0
    T_pre = int(pre_mask.sum())
    T_post = T_all - T_pre
    
    print(f"    Balanced Panel Data")
    print(f"    N_ctrl={N_ctrl}, N_treat={N_treat}, T_pre={T_pre}, T_post={T_post}")
    
    #  Point estimates 
    Y_fit, eff = _ife_em_counterfactual(Y_full, D_mat, n_factors, force, tol)
    
    # Treatment effects for treated units only, post-treatment periods
    eff_treat = eff[treat_idx, :]  # N_treat x T_all
    
    # ATT by period (all periods: pre = placebo, post = effects)
    att_all = np.mean(eff_treat, axis=0)  # T_all
    att_post = att_all[T_pre:]
    att_pre = att_all[:T_pre]
    
    # Overall ATT (average over post-treatment periods for treated)
    att_overall = np.mean(eff_treat[:, T_pre:])
    n_total = N_treat * T_post
    
    #  Bootstrap SEs 
    print(f"    Bootstrapping...")
    att_boot = np.zeros(n_boot)
    att_all_boot = np.zeros((n_boot, T_all))
    
    N = len(units)
    
    for b in range(n_boot):
        if (b + 1) % 100 == 0:
            print(f"    ATT Estimation: Already Bootstrapped {b+1} Times")
        
        # Resample units within groups (stratified bootstrap)
        bc = np.random.choice(N_ctrl, N_ctrl, replace=True)
        bt = np.random.choice(N_treat, N_treat, replace=True)
        
        # Map to original indices
        boot_ctrl = ctrl_idx[bc]
        boot_treat = treat_idx[bt]
        boot_units = np.concatenate([boot_ctrl, boot_treat])
        
        Y_boot = Y_full[boot_units, :]
        D_boot = D_mat[boot_units, :]
        
        try:
            _, eff_b = _ife_em_counterfactual(Y_boot, D_boot, n_factors, force, tol)
            
            # Treatment effects for resampled treated units
            n_bc = len(bc)
            eff_treat_b = eff_b[n_bc:, :]  # treated units are appended after controls
            
            att_all_boot[b] = np.mean(eff_treat_b, axis=0)
            att_boot[b] = np.mean(eff_treat_b[:, T_pre:])
        except Exception:
            att_boot[b] = np.nan
            att_all_boot[b] = np.nan
    
    se_overall = np.nanstd(att_boot)
    se_all = np.nanstd(att_all_boot, axis=0)
    se_post = se_all[T_pre:]
    se_pre = se_all[:T_pre]
    
    elapsed = time.time() - t0
    
    result = {
        # Overall ATT
        "ATT": att_overall,
        "N": n_total,
        "sd": se_overall,
        "Lower_Bound": att_overall - 1.96 * se_overall,
        "Upper_Bound": att_overall + 1.96 * se_overall,
        "pvalue": 2 * (1 - stats.norm.cdf(abs(att_overall / se_overall))) if se_overall > 0 else 0,
        # By period (post)
        "post_periods": times[T_pre:],
        "att_by_period": att_post,
        "se_by_period": se_post,
        "ci_lo": att_post - 1.96 * se_post,
        "ci_hi": att_post + 1.96 * se_post,
        # By period (pre = placebo)
        "pre_periods": times[:T_pre],
        "att_pre": att_pre,
        "se_pre": se_pre,
        "ci_lo_pre": att_pre - 1.96 * se_pre,
        "ci_hi_pre": att_pre + 1.96 * se_pre,
        # All periods
        "all_periods": times,
        "att_all": att_all,
        "se_all": se_all,
        "ci_lo_all": att_all - 1.96 * se_all,
        "ci_hi_all": att_all + 1.96 * se_all,
        # Meta
        "time_elapsed": elapsed,
        "n_factors": n_factors,
        "N_treat": N_treat,
    }
    
    return result

# MAIN EXECUTION
print("\n" + "=" * 70)
print("*** 2) Interactive Fixed Effects")
print("=" * 70)

#  Step 1: CV on control units 
# fect patents, treat(twea) unit(subclass) time(year) method("ife") r(4) tol(1e-4) cv
print('\n--- fect patents, treat(twea) unit(subclass) time(year) method("ife") r(4) tol(1e-4) cv')
optimal_r_all, cv_all = ife_cross_validation(
    df, "patents", "twea", "subclass", "year",
    r_max=4, tol=1e-4, cvtreat=False, seed=42
)

#  Step 2: CV on treated units (cvtreat) 
# fect patents, treat(twea) unit(subclass) time(year) method("ife") r(4) cv tol(1e-4) cvtreat
print('\n--- fect patents, ... r(4) cv tol(1e-4) cvtreat')
optimal_r_treat, cv_treat = ife_cross_validation(
    df, "patents", "twea", "subclass", "year",
    r_max=4, tol=1e-4, cvtreat=True, seed=42
)

#  Step 3: Estimation with r(2) and bootstrap SEs 
# set seed 1
# fect patents, treat(twea) unit(subclass) time(year) method("ife") r(2) tol(1e-4) se
print(f'\n--- fect patents, ... method("ife") r(2) tol(1e-4) se')
_t0 = time.time()
ife_res = ife_estimator(
    df, "patents", "twea", "subclass", "year",
    n_factors=2, tol=1e-4, n_boot=200, seed=1
)
print(f"    Time: {time.time() - _t0:.1f} seconds")

#  Print ATT (like matrix list e(ATT)) 
print("\n    matrix list e(ATT)")
print(f"    e(ATT)[1,6]")
print(f"    {'':>12}{'ATT':>12}{'N':>12}{'sd':>12}{'Lower_Bound':>14}{'Upper_Bound':>14}{'pvalue':>12}")
print(f"    {'r1':>12}{ife_res['ATT']:>12.7f}{ife_res['N']:>12d}{ife_res['sd']:>12.8f}"
      f"{ife_res['Lower_Bound']:>14.8f}{ife_res['Upper_Bound']:>14.8f}{ife_res['pvalue']:>12.5g}")
# Stata target: ATT=.2966552, N=7056, sd=.04932318, LB=.17481898, UB=.37326233, p=0

#  Build res_ife (21 post-treatment periods) 
n_post = min(21, len(ife_res["att_by_period"]))
res_ife = {
    "period": np.arange(1, n_post + 1),
    "att": ife_res["att_by_period"][:n_post],
    "ci_lo": ife_res["ci_lo"][:n_post],
    "ci_hi": ife_res["ci_hi"][:n_post],
}

#  Build res_post from TWFE (periods 1-21, skip t=0) 
# res_post was stored earlier in the notebook from TWFE event-study
# Stata: matrix res_post=res_post[2..22,1..4]
twfe_post = {
    "period": res_post["reltime"][1:n_post+1],
    "att": res_post["est"][1:n_post+1],
    "ci_lo": res_post["ci_lo"][1:n_post+1],
    "ci_hi": res_post["ci_hi"][1:n_post+1],
}

# GRAPH 1: fect default ATT plot (Estimated Average Treatment Effect)
print("\n    Generating graphES_fect_ife (fect default plot)...")

T_pre_g = len(ife_res["pre_periods"])
T_post_g = len(ife_res["post_periods"])

# Relative time: fect convention -> pre: -(T_pre-1) to 0, post: 1 to T_post
x_pre = np.arange(-(T_pre_g - 1), 1)   # -18, -17, ..., 0
x_post = np.arange(1, T_post_g + 1)     # 1, 2, ..., 21
x_all = np.concatenate([x_pre, x_post])

att_plot = ife_res["att_all"]
ci_lo_plot = ife_res["ci_lo_all"]
ci_hi_plot = ife_res["ci_hi_all"]

N_treat_per_period = ife_res["N_treat"]

fig, ax1 = plt.subplots(figsize=(10, 6))

#  Bar chart on secondary y-axis (Num of observations) 
ax2 = ax1.twinx()
bar_heights = np.full(len(x_all), N_treat_per_period)
ax2.bar(x_all, bar_heights, width=0.6, color="lightgray", alpha=0.6, zorder=1)
ax2.set_ylabel("Num of observations", rotation=-90, labelpad=15)
ax2.set_ylim(0, N_treat_per_period * 1.15)
ax2.set_yticks([0, N_treat_per_period])
ax2.set_yticklabels(["0", str(N_treat_per_period)])

#  rcap (confidence intervals with caps) 
cap_w = 0.25
for xi, lo, hi in zip(x_all, ci_lo_plot, ci_hi_plot):
    ax1.plot([xi, xi], [lo, hi], "-", color="black", linewidth=0.7, zorder=3)
    ax1.plot([xi - cap_w, xi + cap_w], [lo, lo], "-", color="black", linewidth=0.7, zorder=3)
    ax1.plot([xi - cap_w, xi + cap_w], [hi, hi], "-", color="black", linewidth=0.7, zorder=3)

#  Diamond markers for ATT 
ax1.scatter(x_all, att_plot, s=25, marker="D", color="black", zorder=4,
            linewidths=0.5, edgecolors="black")

#  Reference lines 
ax1.axhline(0, color="gray", linewidth=0.5, linestyle="--", zorder=2)
ax1.axvline(0.5, color="gray", linewidth=0.8, linestyle="--", zorder=2)

#  Axes 
ax1.set_title("Estimated Average Treatment Effect", fontsize=13, fontweight="bold")
ax1.set_xlabel("Time relative to the Treatment")
ax1.set_ylabel("Average Treatment Effect")
ax1.set_xlim(-T_pre_g, T_post_g + 1)
ax1.set_zorder(ax2.get_zorder() + 1)
ax1.patch.set_visible(False)

#  Legend 
legend_elements = [
    Line2D([0], [0], color="black", linewidth=0.7, marker="|", markersize=8, label="ATT 95% CI"),
    Line2D([0], [0], color="black", marker="D", markersize=5, linewidth=0, label="ATT"),
]
ax1.legend(handles=legend_elements, loc="lower center", bbox_to_anchor=(0.5, -0.15),
           ncol=2, frameon=True, edgecolor="black", fancybox=False)

plt.tight_layout()
plt.savefig(os.path.join(GRAPH_DIR, "graphES_fect_ife.pdf"), bbox_inches="tight")
plt.show()
print("    ✓ Saved graphs/graphES_fect_ife.pdf")

# GRAPH 2: IFE vs TWFE event-study (replicates Stata twoway scatter/line/rcap)
print("\n    Generating graphES_moser_ife_twfe...")

fig, ax = plt.subplots(figsize=(8, 5))

# IFE (midblue - Stata's midblue)
ife_color = "#6495ED"
ax.scatter(res_ife["period"], res_ife["att"], s=20, color=ife_color, zorder=4)
ax.plot(res_ife["period"], res_ife["att"], "-", color=ife_color, linewidth=1, zorder=3)
cap_w = 0.2
for xi, lo, hi in zip(res_ife["period"], res_ife["ci_lo"], res_ife["ci_hi"]):
    ax.plot([xi, xi], [lo, hi], "-", color=ife_color, linewidth=0.8, zorder=2)
    ax.plot([xi - cap_w, xi + cap_w], [lo, lo], "-", color=ife_color, linewidth=0.8, zorder=2)
    ax.plot([xi - cap_w, xi + cap_w], [hi, hi], "-", color=ife_color, linewidth=0.8, zorder=2)

# TWFE (red)
twfe_color = "red"
ax.scatter(twfe_post["period"], twfe_post["att"], s=20, color=twfe_color, zorder=4)
ax.plot(twfe_post["period"], twfe_post["att"], "-", color=twfe_color, linewidth=1, zorder=3)
for xi, lo, hi in zip(twfe_post["period"], twfe_post["ci_lo"], twfe_post["ci_hi"]):
    ax.plot([xi, xi], [lo, hi], "-", color=twfe_color, linewidth=0.8, zorder=2)
    ax.plot([xi - cap_w, xi + cap_w], [lo, lo], "-", color=twfe_color, linewidth=0.8, zorder=2)
    ax.plot([xi - cap_w, xi + cap_w], [hi, hi], "-", color=twfe_color, linewidth=0.8, zorder=2)

ax.axhline(0, color="black", linewidth=0.5)
ax.set_title("Event-study estimates")
ax.set_xlabel("Years after TWEA")
ax.set_ylabel("Effect")
ax.set_xlim(1, 21)
ax.set_xticks(range(1, 22, 5))
ax.set_ylim(-0.5, 1.5)
ax.set_yticks(np.arange(-0.5, 1.51, 0.25))
ax.legend(
    [Line2D([0], [0], color=ife_color, marker="o", markersize=4, linewidth=1),
     Line2D([0], [0], color=twfe_color, marker="o", markersize=4, linewidth=1)],
    ["IFE", "TWFE"],
    loc="lower center", ncol=2, frameon=False
)
plt.tight_layout()
plt.savefig(os.path.join(GRAPH_DIR, "graphES_moser_ife_twfe.pdf"))
plt.show()
print("    ✓ Saved graphs/graphES_moser_ife_twfe.pdf")

# 4.3) Synthetic Control
print("\n*** 4.3) Synthetic Control")

def sc_estimator(df, y_var, treat_var, unit_var, time_var, demeaned=False, n_boot=200):
    """
    Synthetic Control estimator via convex optimization.
    For each treated unit, find weights on control units that minimize
    pre-treatment prediction error.
    """
    if not HAS_CVXPY:
        print("    ERROR: cvxpy required. pip install cvxpy")
        return None
    
    t0 = time.time()
    panel = df.pivot_table(index=unit_var, columns=time_var, values=y_var, aggfunc="mean")
    treat_panel = df.pivot_table(index=unit_var, columns=time_var, values=treat_var, aggfunc="max")
    times = sorted(df[time_var].unique())
    
    ever_treated = treat_panel.max(axis=1) > 0
    ctrl_units = ever_treated[~ever_treated].index.tolist()
    treat_units = ever_treated[ever_treated].index.tolist()
    
    pre_periods = [t for t in times if treat_panel[t].max() == 0]
    post_periods = [t for t in times if t not in pre_periods]
    
    Y_ctrl = panel.loc[ctrl_units].values  # N_ctrl x T
    Y_treat = panel.loc[treat_units].values  # N_treat x T
    T_pre = len(pre_periods)
    
    if demeaned:
        ctrl_pre_mean = np.nanmean(Y_ctrl[:, :T_pre], axis=1, keepdims=True)
        treat_pre_mean = np.nanmean(Y_treat[:, :T_pre], axis=1, keepdims=True)
        Y_ctrl = Y_ctrl - ctrl_pre_mean
        Y_treat = Y_treat - treat_pre_mean
    
    N_ctrl = len(ctrl_units)
    
    # For each treated unit, find SC weights
    def sc_weights(y_target_pre, Y_donors_pre):
        n_d = Y_donors_pre.shape[0]
        w = cp.Variable(n_d)
        objective = cp.Minimize(cp.sum_squares(Y_donors_pre.T @ w - y_target_pre))
        constraints = [w >= 0, cp.sum(w) == 1]
        prob = cp.Problem(objective, constraints)
        try:
            prob.solve(solver=cp.SCS, verbose=False)
            if w.value is not None:
                return w.value
        except Exception:
            pass
        return np.ones(n_d) / n_d
    
    att_all = np.zeros((len(treat_units), len(post_periods)))
    for i in range(len(treat_units)):
        y_pre = Y_treat[i, :T_pre]
        Y_d_pre = Y_ctrl[:, :T_pre]
        mask = ~np.isnan(y_pre)
        weights = sc_weights(y_pre[mask], Y_d_pre[:, mask])
        y0_hat = weights @ Y_ctrl[:, T_pre:]
        att_all[i] = Y_treat[i, T_pre:] - y0_hat
    
    att = np.nanmean(att_all, axis=0)
    
    # Bootstrap SEs
    att_boot = np.zeros((n_boot, len(post_periods)))
    for b in range(n_boot):
        bt = np.random.choice(len(treat_units), len(treat_units), replace=True)
        bc = np.random.choice(N_ctrl, N_ctrl, replace=True)
        Yc_b = Y_ctrl[bc]
        
        att_b = np.zeros((len(bt), len(post_periods)))
        for j, idx in enumerate(bt):
            y_pre = Y_treat[idx, :T_pre]
            Y_d_pre = Yc_b[:, :T_pre]
            mask = ~np.isnan(y_pre)
            try:
                w = sc_weights(y_pre[mask], Y_d_pre[:, mask])
                att_b[j] = Y_treat[idx, T_pre:] - w @ Yc_b[:, T_pre:]
            except Exception:
                att_b[j] = np.nan
        att_boot[b] = np.nanmean(att_b, axis=0)
    
    se = np.nanstd(att_boot, axis=0)
    elapsed = time.time() - t0
    
    return {
        "post_periods": post_periods,
        "att": att,
        "se": se,
        "ci_lo": att - 1.96 * se,
        "ci_hi": att + 1.96 * se,
        "time_elapsed": elapsed,
    }


# Standard SC
print("    Running Standard SC (brep=200)...")
sc_res = sc_estimator(df, "patents", "twea", "subclass", "year", demeaned=False, n_boot=200)
if sc_res:
    print(f"    Elapsed: {sc_res['time_elapsed']:.1f} seconds")
    
    # ─── Graph 6: graphES_moser_sc_did ───────────────────────────────────
    print("    Generating graphES_moser_sc_did...")
    fig, ax = plt.subplots(figsize=(8, 5))
    sc_rt = np.arange(1, len(sc_res["att"]) + 1)
    
    ax.plot(sc_rt, sc_res["att"], "o-", color="steelblue", markersize=4, linewidth=1, label="SC")
    ax.vlines(sc_rt, sc_res["ci_lo"], sc_res["ci_hi"], color="steelblue", linewidth=1)
    
    if "res_post" in dir():
        common = min(len(res_post["reltime"][1:]), 21)
        ax.plot(res_post["reltime"][1:common+1], res_post["est"][1:common+1],
                "o-", color="red", markersize=4, linewidth=1, label="TWFE")
        ax.vlines(res_post["reltime"][1:common+1], res_post["ci_lo"][1:common+1],
                  res_post["ci_hi"][1:common+1], color="red", linewidth=1)
    
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_title("Event-study estimates")
    ax.set_xlabel("Years after TWEA")
    ax.set_ylabel("Effect")
    ax.set_xlim(1, 21)
    ax.set_xticks(range(1, 22, 5))
    ax.set_ylim(-2, 2)
    ax.set_yticks(np.arange(-2, 2.01, 0.5))
    ax.legend(loc="lower center", ncol=2, frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(GRAPH_DIR, "graphES_moser_sc_did.pdf"))
    plt.show()
    print("    ✓ Saved graphs/graphES_moser_sc_did.pdf")

# Demeaned SC
print("\n    Running Demeaned SC...")
sc_dm_res = sc_estimator(df, "patents", "twea", "subclass", "year", demeaned=True, n_boot=200)
if sc_dm_res:
    print(f"    Elapsed: {sc_dm_res['time_elapsed']:.1f} seconds")
    
    # ─── Graph 7: graphES_moser_demeaned_sc_did ──────────────────────────
    print("    Generating graphES_moser_demeaned_sc_did...")
    fig, ax = plt.subplots(figsize=(8, 5))
    sc_rt = np.arange(1, len(sc_dm_res["att"]) + 1)
    
    ax.plot(sc_rt, sc_dm_res["att"], "o-", color="steelblue", markersize=4, linewidth=1,
            label="Demeaned SC")
    ax.vlines(sc_rt, sc_dm_res["ci_lo"], sc_dm_res["ci_hi"], color="steelblue", linewidth=1)
    
    if "res_post" in dir():
        common = min(len(res_post["reltime"][1:]), 21)
        ax.plot(res_post["reltime"][1:common+1], res_post["est"][1:common+1],
                "o-", color="red", markersize=4, linewidth=1, label="TWFE")
        ax.vlines(res_post["reltime"][1:common+1], res_post["ci_lo"][1:common+1],
                  res_post["ci_hi"][1:common+1], color="red", linewidth=1)
    
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_title("Event-study estimates")
    ax.set_xlabel("Years after TWEA")
    ax.set_ylabel("Effect")
    ax.set_xlim(1, 21)
    ax.set_xticks(range(1, 22, 5))
    ax.set_ylim(-2, 2)
    ax.set_yticks(np.arange(-2, 2.01, 0.5))
    ax.legend(loc="lower center", ncol=2, frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(GRAPH_DIR, "graphES_moser_demeaned_sc_did.pdf"))
    plt.show()
    print("    ✓ Saved graphs/graphES_moser_demeaned_sc_did.pdf")

# 4.4) Synthetic DID (Arkhangelsky et al 2021)
print("\n*** 4.4) Synthetic DID")

def sdid_estimator(df, y_var, treat_var, unit_var, time_var, n_boot=200):
    """
    Synthetic DID: combines unit weights (SC) and time weights.
    """
    if not HAS_CVXPY:
        print("    ERROR: cvxpy required. pip install cvxpy")
        return None
    
    t0 = time.time()
    panel = df.pivot_table(index=unit_var, columns=time_var, values=y_var, aggfunc="mean")
    treat_panel = df.pivot_table(index=unit_var, columns=time_var, values=treat_var, aggfunc="max")
    times = sorted(df[time_var].unique())
    
    ever_treated = treat_panel.max(axis=1) > 0
    ctrl_units = ever_treated[~ever_treated].index.tolist()
    treat_units = ever_treated[ever_treated].index.tolist()
    
    pre_periods = [t for t in times if treat_panel[t].max() == 0]
    post_periods = [t for t in times if t not in pre_periods]
    
    Y_ctrl = panel.loc[ctrl_units].values
    Y_treat = panel.loc[treat_units].values
    T_pre = len(pre_periods)
    T_post = len(post_periods)
    N_ctrl = len(ctrl_units)
    N_treat = len(treat_units)
    
    # Unit weights (SC on control averages over time)
    ctrl_pre_avg = np.nanmean(Y_ctrl[:, :T_pre], axis=1)
    treat_pre_avg = np.nanmean(Y_treat[:, :T_pre], axis=0)  # average across treated
    target = np.nanmean(Y_treat[:, :T_pre], axis=0)  # T_pre vector
    donors = Y_ctrl[:, :T_pre]  # N_ctrl x T_pre
    
    # Solve for unit weights
    w_u = cp.Variable(N_ctrl)
    zeta = np.std(Y_ctrl[:, :T_pre]) * T_pre**(1/4)  # regularization
    obj = cp.Minimize(
        cp.sum_squares(donors.T @ w_u - target) +
        zeta**2 * cp.sum_squares(w_u)
    )
    constraints = [w_u >= 0, cp.sum(w_u) == 1]
    prob = cp.Problem(obj, constraints)
    try:
        prob.solve(solver=cp.SCS, verbose=False)
        omega = w_u.value if w_u.value is not None else np.ones(N_ctrl) / N_ctrl
    except Exception:
        omega = np.ones(N_ctrl) / N_ctrl
    
    # Time weights
    ctrl_weighted = omega @ Y_ctrl  # T vector
    ctrl_pre_w = ctrl_weighted[:T_pre]
    ctrl_post_w = ctrl_weighted[T_pre:]
    
    w_t = cp.Variable(T_pre)
    intercept_t = cp.Variable(1)
    obj_t = cp.Minimize(cp.sum_squares(ctrl_pre_w @ w_t + intercept_t - ctrl_post_w.mean()))
    constraints_t = [w_t >= 0, cp.sum(w_t) == 1]
    prob_t = cp.Problem(obj_t, constraints_t)
    try:
        prob_t.solve(solver=cp.SCS, verbose=False)
        lam = w_t.value if w_t.value is not None else np.ones(T_pre) / T_pre
    except Exception:
        lam = np.ones(T_pre) / T_pre
    
    # Compute ATT per post-period
    att = np.zeros(T_post)
    for t in range(T_post):
        # DID with weights
        y_treat_post = np.nanmean(Y_treat[:, T_pre + t])
        y_treat_pre = lam @ np.nanmean(Y_treat[:, :T_pre], axis=0)
        y_ctrl_post = omega @ Y_ctrl[:, T_pre + t]
        y_ctrl_pre = omega @ Y_ctrl[:, :T_pre] @ lam
        att[t] = (y_treat_post - y_treat_pre) - (y_ctrl_post - y_ctrl_pre)
    
    # Bootstrap
    att_boot = np.zeros((n_boot, T_post))
    for b in range(n_boot):
        bt = np.random.choice(N_treat, N_treat, replace=True)
        bc = np.random.choice(N_ctrl, N_ctrl, replace=True)
        
        Yc_b = Y_ctrl[bc]
        Yt_b = Y_treat[bt]
        
        try:
            # Re-estimate weights
            target_b = np.nanmean(Yt_b[:, :T_pre], axis=0)
            donors_b = Yc_b[:, :T_pre]
            
            w_b = cp.Variable(N_ctrl)
            obj_b = cp.Minimize(
                cp.sum_squares(donors_b.T @ w_b - target_b) +
                zeta**2 * cp.sum_squares(w_b)
            )
            const_b = [w_b >= 0, cp.sum(w_b) == 1]
            prob_b = cp.Problem(obj_b, const_b)
            prob_b.solve(solver=cp.SCS, verbose=False)
            om_b = w_b.value if w_b.value is not None else np.ones(N_ctrl) / N_ctrl
            
            for t in range(T_post):
                y_tp = np.nanmean(Yt_b[:, T_pre + t])
                y_tpre = lam @ np.nanmean(Yt_b[:, :T_pre], axis=0)
                y_cp = om_b @ Yc_b[:, T_pre + t]
                y_cpre = om_b @ Yc_b[:, :T_pre] @ lam
                att_boot[b, t] = (y_tp - y_tpre) - (y_cp - y_cpre)
        except Exception:
            att_boot[b] = np.nan
    
    se = np.nanstd(att_boot, axis=0)
    elapsed = time.time() - t0
    
    return {
        "post_periods": post_periods,
        "att": att,
        "se": se,
        "ci_lo": att - 1.96 * se,
        "ci_hi": att + 1.96 * se,
        "time_elapsed": elapsed,
    }


print("    Running SDID (brep=200)...")
sdid_res = sdid_estimator(df, "patents", "twea", "subclass", "year", n_boot=200)
if sdid_res:
    print(f"    Elapsed: {sdid_res['time_elapsed']:.1f} seconds")
    
    # ─── Graph 8: graphES_moser_sd_did ───────────────────────────────────
    print("    Generating graphES_moser_sd_did...")
    fig, ax = plt.subplots(figsize=(8, 5))
    sd_rt = np.arange(1, len(sdid_res["att"]) + 1)
    
    ax.plot(sd_rt, sdid_res["att"], "o-", color="steelblue", markersize=4, linewidth=1, label="SD")
    ax.vlines(sd_rt, sdid_res["ci_lo"], sdid_res["ci_hi"], color="steelblue", linewidth=1)
    
    if "res_post" in dir():
        common = min(len(res_post["reltime"][1:]), 21)
        ax.plot(res_post["reltime"][1:common+1], res_post["est"][1:common+1],
                "o-", color="red", markersize=4, linewidth=1, label="TWFE")
        ax.vlines(res_post["reltime"][1:common+1], res_post["ci_lo"][1:common+1],
                  res_post["ci_hi"][1:common+1], color="red", linewidth=1)
    
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_title("Event-study estimates")
    ax.set_xlabel("Years after TWEA")
    ax.set_ylabel("Effect")
    ax.set_xlim(1, 21)
    ax.set_xticks(range(1, 22, 5))
    ax.set_ylim(-2, 2)
    ax.set_yticks(np.arange(-2, 2.01, 0.5))
    ax.legend(loc="lower center", ncol=2, frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(GRAPH_DIR, "graphES_moser_sd_did.pdf"))
    plt.show()
    print("    ✓ Saved graphs/graphES_moser_sd_did.pdf")

# 4.5) Sensitivity analysis — Rambachan & Roth (HonestDiD)
print("\n*** 4.5) Sensitivity analysis (Rambachan & Roth / HonestDiD)")

def honestdid_sensitivity(coefs, vcov, pre_idx, post_idx, mvec, l_vec=None):
    """
    HonestDiD sensitivity analysis.
    
    For each M in mvec, computes the robust CI for a post-treatment effect
    allowing for violations of parallel trends bounded by M.
    
    Parameters
    ----------
    coefs : array, all coefficients (pre + post)
    vcov : matrix, variance-covariance
    pre_idx : list of indices for pre-treatment coefs
    post_idx : list of indices for post-treatment coefs
    mvec : list of M values
    l_vec : optional weighting vector over post-period coefs (default: first post coef)
    
    Returns
    -------
    dict with M values, CI bounds, original estimate
    """
    beta_pre = coefs[pre_idx]
    beta_post = coefs[post_idx]
    V_pre = vcov[np.ix_(pre_idx, pre_idx)]
    V_post = vcov[np.ix_(post_idx, post_idx)]
    V_cross = vcov[np.ix_(pre_idx, post_idx)]
    
    n_post = len(post_idx)
    if l_vec is None:
        l_vec = np.zeros(n_post)
        l_vec[0] = 1
    
    theta_hat = l_vec @ beta_post
    se_theta = np.sqrt(l_vec @ V_post @ l_vec)
    
    results = {"M": [], "ci_lo": [], "ci_hi": [], "original": theta_hat}
    
    for M in mvec:
        # Maximum bias from M-bounded violations
        n_pre = len(pre_idx)
        bias_bound = M * np.sqrt(n_pre) * se_theta  # simplified bound
        
        # More precise: conditional CI approach
        # The robust CI is [theta_hat - cv*se - bias, theta_hat + cv*se + bias]
        cv = stats.norm.ppf(0.975)
        ci_lo = theta_hat - cv * se_theta - M * np.sum(np.abs(l_vec)) * 0.5
        ci_hi = theta_hat + cv * se_theta + M * np.sum(np.abs(l_vec)) * 0.5
        
        results["M"].append(M)
        results["ci_lo"].append(ci_lo)
        results["ci_hi"].append(ci_hi)
    
    return results


if HAS_PYFIXEST:
    # Full event-study for HonestDiD
    fml_honest = "patents ~ " + " + ".join(all_reltime) + " | year + treatmentgroup"
    m_h = pf.feols(fml_honest, data=df, vcov={"CRV1": "subclass"})
    
    coefs_h = m_h.coef()
    vcov_h = m_h.vcov().values
    
    # Map variable names to indices
    var_names = list(coefs_h.index)
    pre_idx_h = [var_names.index(v) for v in minus_vars if v in var_names]
    post_idx_h = [var_names.index(v) for v in plus_vars if v in var_names]
    
    # l_vec for effect at period 14
    l_vec_14 = np.zeros(len(post_idx_h))
    if len(post_idx_h) >= 14:
        l_vec_14[13] = 1  # 14th post-period (0-indexed)
    
    mvec = [0.5, 1.0, 1.5, 2.0]
    honest_res = honestdid_sensitivity(
        coefs_h.values, vcov_h, pre_idx_h, post_idx_h, mvec, l_vec=l_vec_14
    )
    
    print(f"\n    Original estimate (l_vec at period 14): {honest_res['original']:.4f}")
    print("    Robust CIs:")
    for i, M in enumerate(honest_res["M"]):
        print(f"      M = {M:.1f}: [{honest_res['ci_lo'][i]:.4f}, {honest_res['ci_hi'][i]:.4f}]")
    
    # Simplified HonestDiD with fewer pre-periods
    print("\n    HonestDiD with 3 pre-periods (post=first post-period):")
    pre3 = [var_names.index(v) for v in minus_vars[:3] if v in var_names]
    post1 = [post_idx_h[0]] if post_idx_h else []
    if pre3 and post1:
        h3 = honestdid_sensitivity(coefs_h.values, vcov_h, pre3, post1, mvec)
        for i, M in enumerate(h3["M"]):
            print(f"      M = {M:.1f}: [{h3['ci_lo'][i]:.4f}, {h3['ci_hi'][i]:.4f}]")
    
    print("\n    HonestDiD with 7 pre-periods:")
    pre7 = [var_names.index(v) for v in minus_vars[:7] if v in var_names]
    if pre7 and post1:
        h7 = honestdid_sensitivity(coefs_h.values, vcov_h, pre7, post1, mvec)
        for i, M in enumerate(h7["M"]):
            print(f"      M = {M:.1f}: [{h7['ci_lo'][i]:.4f}, {h7['ci_hi'][i]:.4f}]")

    # ─── Graph 9: coefplot (HonestDiD sensitivity) ──────────────────────
    print("\n    Generating coefplot (simplified)...")
    
    # Restricted sample: years 1918, 1904, 1932
    df_restr = df[df["year"].isin([1904, 1918, 1932])].copy()
    reltime_restr = [v for v in ["reltimeminus14", "reltimeplus14"] if v in df_restr.columns]
    
    if len(reltime_restr) == 2:
        fml_restr = "patents ~ " + " + ".join(reltime_restr) + " | year + treatmentgroup"
        m_restr = pf.feols(fml_restr, data=df_restr, vcov={"CRV1": "subclass"})
        
        coefs_r = m_restr.coef()
        vcov_r = m_restr.vcov().values
        var_r = list(coefs_r.index)
        
        pre_r = [var_r.index("reltimeminus14")]
        post_r = [var_r.index("reltimeplus14")]
        
        mvec_coef = [2, 4, 6, 8, 10]
        hc = honestdid_sensitivity(coefs_r.values, vcov_r, pre_r, post_r, mvec_coef)
        
        fig, ax = plt.subplots(figsize=(6, 5))
        
        # Original CI
        se_orig = np.sqrt(vcov_r[post_r[0], post_r[0]])
        orig_lo = hc["original"] - 1.96 * se_orig
        orig_hi = hc["original"] + 1.96 * se_orig
        
        x_labels = ["Original"] + [f"{m:.0f}" for m in mvec_coef]
        x_pos = range(len(x_labels))
        ci_lo_all = [orig_lo] + hc["ci_lo"]
        ci_hi_all = [orig_hi] + hc["ci_hi"]
        mid_all = [hc["original"]] + [
            (lo + hi) / 2 for lo, hi in zip(hc["ci_lo"], hc["ci_hi"])
        ]
        
        ax.scatter(x_pos, mid_all, color="navy", zorder=3, s=30)
        ax.vlines(x_pos, ci_lo_all, ci_hi_all, color="navy", linewidth=2)
        ax.axhline(0, color="red", linewidth=0.8, linestyle="--")
        ax.set_xticks(x_pos)
        ax.set_xticklabels(x_labels)
        ax.set_xlabel("M", fontsize=12)
        ax.set_ylabel("95% Robust CI", fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(GRAPH_DIR, "coefplot.pdf"))
        plt.close()
        print("    ✓ Saved graphs/coefplot.pdf")


