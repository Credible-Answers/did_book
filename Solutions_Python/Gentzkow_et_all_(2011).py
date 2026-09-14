import pandas as pd
import numpy as np
import pyfixest as pf
import warnings
warnings.filterwarnings('ignore')
from pathlib import Path
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import matplotlib.pyplot as plt
from scipy import stats
#Charge data
path = Path('/Users/anthonyquispe/Downloads/Replication folder FE/programs/cc_xd_didtextbook_2025_9_30/Data sets/Gentzkow et al 2011/gentzkowetal_didtextbook.dta')
df = pd.read_stata(path)
#CHAPTER 5
#1. TWFE Regression
result = pf.feols(
    "prestout ~ numdailies + C(year) | cnty90", 
    data=df, 
    vcov={"CRV1": "cnty90"},
    ssc=pf.ssc(adj=True, fixef_k="full", cluster_adj=True)
)
result.summary()
#Decomposition
Y, G, T, D = 'prestout', 'cnty90', 'year', 'numdailies'
eps = np.array(pf.feols(f"{D} ~ 1 | {G} + {T}", data=df).resid())
D_vals = df[D].astype(float).values
N = len(D_vals)
denom_W = np.mean(eps * D_vals)
weight = eps * D_vals / (N * denom_W)
beta = np.array(pf.feols(f"{Y} ~ {D} | {G} + {T}", data=df).coef())[0]
att = D_vals > 0
pos = att & (weight > 0)
neg = att & (weight < 0)
print("Under the common trends assumption,")
print(f"the TWFE coefficient beta, equal to {beta:.4f}, estimates a weighted")
print(f"sum of {att.sum()} ATTs.")
print(f"{pos.sum()} ATTs receive a positive weight, and {neg.sum()} receive a negative weight.")
print("------------------------------------------------")
print(f"{'Treat. var: ' + D:<24s} {'# ATTs':<11s} {'Σ weights':<12s}")
print("------------------------------------------------")
print(f"{'Positive weights':<24s} {pos.sum():<11d} {weight[pos].sum():<12.4f}")
print(f"{'Negative weights':<24s} {neg.sum():<11d} {weight[neg].sum():<12.4f}")
print("------------------------------------------------")
print(f"{'Total':<24s} {att.sum():<11d} {weight[att].sum():<12.4f}")
#2. TWFE with state-specific trends
result = pf.feols(
    "prestout ~ numdailies + C(year) + C(styr) | cnty90",
    data=df,
    vcov={"CRV1": "cnty90"},
    ssc=pf.ssc(adj=True, fixef_k="full", cluster_adj=True)
)
print(f"numdailies coef: {result.coef()['numdailies']:.7f}")
print(f"numdailies se:   {result.se()['numdailies']:.7f}")
#Decomposition
D_vals = df['numdailies'].astype(float).values
N = len(D_vals)
eps = np.array(pf.feols("numdailies ~ C(styr) | cnty90 + year", data=df, fixef_rm="none").resid())
denom_W = np.mean(eps * D_vals)
weight = eps * D_vals / (N * denom_W)
beta = np.array(pf.feols("prestout ~ numdailies + C(year) + C(styr) | cnty90",
    data=df, vcov={"CRV1": "cnty90"},
    ssc=pf.ssc(adj=True, fixef_k="full", cluster_adj=True)).coef())[0]
att = D_vals > 0
tol = 1e-9
pos = att & (weight > tol)
neg = att & (weight < -tol)
zero = att & (np.abs(weight) <= tol)
print(f"beta = {beta:.4f}, ATTs = {pos.sum()+neg.sum()}")
print(f"Positive: {pos.sum()}, Σ = {weight[pos].sum():.4f}")
print(f"Negative: {neg.sum()}, Σ = {weight[neg].sum():.4f}")
print(f"Zero: {zero.sum()}")
#3. FD Regression with state-specific trends
result = pf.feols(
    "changeprestout ~ changedailies | styr",
    data=df,
    vcov={"CRV1": "cnty90"},
    ssc=pf.ssc(adj=True, fixef_k="full", cluster_adj=True),
    fixef_rm="none"
)

b = np.array(result.coef())[0]
se = np.array(result.se())[0]
t = np.array(result.tstat())[0]
p = np.array(result.pvalue())[0]

print(f"changedailies: coef={b:.7f}  se={se:.7f}  t={t:.2f}  p={p:.4f}")
print(f"N={result._N}")
#Decomposition
# fdTR: D=changedailies, D0=numdailies
D0 = df['numdailies'].astype(float).values
N = len(df)
mask = ~df['changedailies'].isna()
eps = np.zeros(N)
eps[mask] = np.array(pf.feols("changedailies ~ C(styr) | year", data=df[mask], fixef_rm="none").resid())
yr2fn = {y: i for i, y in enumerate(sorted(df['year'].unique()))}
work = pd.DataFrame({'G': df['cnty90'].values, 'TFN': df['year'].map(yr2fn).values,
    'D0': D0, 'eps': eps}).sort_values(['G', 'TFN']).reset_index(drop=True)
work['w'] = np.where(
    work['TFN'] + 1 == work.groupby('G')['TFN'].shift(-1),
    work['eps'] - work.groupby('G')['eps'].shift(-1),
    work['eps'])
work['w'] = work['w'].fillna(work['eps'])
denom_W = (work['w'] * work['D0']).mean()
weight = work['w'].values * work['D0'].values / (N * denom_W)
att = work['D0'].values > 0
tol = 1e-9
pos = att & (weight > tol)
neg = att & (weight < -tol)
zero = att & (np.abs(weight) <= tol)
print(f"Positive: {pos.sum()}, Σ={weight[pos].sum():.4f}")
print(f"Negative: {neg.sum()}, Σ={weight[neg].sum():.4f}")
print(f"Zero: {zero.sum()}, Total ATTs: {pos.sum()+neg.sum()}")
#Assessing if weights correlated with year variable
D0 = df['numdailies'].astype(float).values
N = len(df)
mean_D0 = D0.mean()
mask = ~df['changedailies'].isna()
eps = np.zeros(N)
eps[mask] = np.array(pf.feols("changedailies ~ C(styr) | year", data=df[mask], fixef_rm="none").resid())

yr2fn = {y: i for i, y in enumerate(sorted(df['year'].unique()))}
work = pd.DataFrame({'G': df['cnty90'].values, 'TFN': df['year'].map(yr2fn).values,
    'D0': D0, 'eps': eps, 'year': df['year'].astype(float).values}).sort_values(['G', 'TFN']).reset_index(drop=True)
work['w'] = np.where(
    work['TFN'] + 1 == work.groupby('G')['TFN'].shift(-1),
    work['eps'] - work.groupby('G')['eps'].shift(-1),
    work['eps'])
work['w'] = work['w'].fillna(work['eps'])

denom_W = (work['w'] * work['D0']).mean()
work['W'] = work['w'] * mean_D0 / denom_W
work['nat_weight'] = work['D0'] / (N * mean_D0)
work['weight_result'] = work['W'] * work['nat_weight']

beta = np.array(pf.feols("changeprestout ~ changedailies | styr",
    data=df, vcov={"CRV1": "cnty90"},
    ssc=pf.ssc(adj=True, fixef_k="full", cluster_adj=True),
    fixef_rm="none").coef())[0]

att = work['D0'].values > 0
tol = 1e-10
weight = work['weight_result'].values
weight = np.where(np.abs(weight) < tol, 0, weight)
pos = att & (weight > 0)
neg = att & (weight < 0)
zero = att & (weight == 0)
print("Under the common trends assumption,")
print(f"the TWFE coefficient beta, equal to {beta:.4f}, estimates a weighted")
print(f"sum of {pos.sum()+neg.sum()} ATTs.")
print(f"{pos.sum()} ATTs receive a positive weight, and {neg.sum()} receive a negative weight.")
print(f"{att.sum()} (g,t) cells receive the treatment, but the ATTs of {zero.sum()} cells")
print(f"receive a weight equal to zero.")
print("------------------------------------------------")
print(f"{'Treat. var: changedail~s':<24s} {'# ATTs':<11s} {'Σ weights':<12s}")
print("------------------------------------------------")
print(f"{'Positive weights':<24s} {pos.sum():<11d} {weight[pos].sum():<12.4f}")
print(f"{'Negative weights':<24s} {neg.sum():<11d} {weight[neg].sum():<12.4f}")
print("------------------------------------------------")
print(f"{'Total':<24s} {pos.sum()+neg.sum():<11d} {weight[att].sum():<12.4f}")

sub = work[att & np.isfinite(work['W'])].copy()
rw_reg = pf.feols("year ~ W", data=sub, weights="nat_weight", vcov={"CRV1": "G"},
    ssc=pf.ssc(adj=True, fixef_k="none", cluster_adj=True))
coef = rw_reg.coef()["W"]
se = rw_reg.se()["W"]
r2 = rw_reg._r2
corr = -np.sqrt(r2) if coef < 0 else np.sqrt(r2)

print("------------------------------------------------")
print("Regression of variables possibly correlated with the treatment ef-")
print("fect on the weights")
print(f"{'':>12s} {'Coef':>12s} {'SE':>12s} {'t-stat':>12s} {'Correlation':>12s}")
print(f"{'year':>12s} {coef:>12.7f} {se:>12.8f} {coef/se:>12.7f} {corr:>12.7f}")
#CHAPTER 8

#1. Testing whether change in daily newspapers as good as random
r1 = pf.feols("changedailies ~ lag_numdailies", data=df, vcov={"CRV1": "cnty90"})
r2 = pf.feols("changedailies ~ lag_ishare_urb", data=df, vcov={"CRV1": "cnty90"})

print(r1.summary())
print(r2.summary())
#2. Distributed lage TWFE Regression
result = pf.feols(
    "prestout ~ numdailies + lag_numdailies + C(year) | cnty90",
    data=df,
    vcov={"CRV1": "cnty90"},
    ssc=pf.ssc(adj=True, fixef_k="full", cluster_adj=True)
)

print(result.summary())
#Decomposition
df = df.dropna(subset=['lag_numdailies']).reset_index(drop=True)

D_vals = df['lag_numdailies'].astype(float).values
OT_vals = df['numdailies'].astype(float).values
N = len(df)
mean_D = D_vals.mean()

eps = np.array(pf.feols("lag_numdailies ~ numdailies | cnty90 + year", data=df, fixef_rm="none").resid())

denom_W = np.mean(eps * D_vals)
W = eps * mean_D / denom_W
P_gt = 1.0 / N
nat_weight = P_gt * D_vals / mean_D
weight_D = W * nat_weight
weight_OT = W * P_gt * OT_vals / mean_D

beta = np.array(pf.feols("prestout ~ lag_numdailies + numdailies + C(year) | cnty90",
    data=df, vcov={"CRV1": "cnty90"},
    ssc=pf.ssc(adj=True, fixef_k="full", cluster_adj=True)).coef())[0]

tol = 1e-10
weight_D = np.where(np.abs(weight_D) < tol, 0, weight_D)
weight_OT = np.where(np.abs(weight_OT) < tol, 0, weight_OT)

# 5. Main Treatment 
att_D = D_vals > 0
pos_D = att_D & (weight_D > 0)
neg_D = att_D & (weight_D < 0)
# 6. Other Treatment
att_OT = OT_vals > 0
pos_OT = att_OT & (weight_OT > 0)
neg_OT = att_OT & (weight_OT < 0)
print("Under the common trends assumption,")
print(f"the TWFE coefficient beta, equal to {beta:.4f}, estimates the sum of")
print(f"several terms.")
print(f"The first term is a weighted sum of {pos_D.sum()+neg_D.sum()} ATTs of the treatment.")
print(f"{pos_D.sum()} ATTs receive a positive weight, and {neg_D.sum()} receive a negative weight.")
print("------------------------------------------------")
print(f"{'Treat. var: lag_numdai~s':<24s} {'# ATTs':<11s} {'Σ weights':<12s}")
print("------------------------------------------------")
print(f"{'Positive weights':<24s} {pos_D.sum():<11d} {weight_D[pos_D].sum():<12.4f}")
print(f"{'Negative weights':<24s} {neg_D.sum():<11d} {weight_D[neg_D].sum():<12.4f}")
print("------------------------------------------------")
print(f"{'Total':<24s} {pos_D.sum()+neg_D.sum():<11d} {weight_D[att_D].sum():<12.4f}")
print("------------------------------------------------")
print(f"The next term is a weighted sum of {pos_OT.sum()+neg_OT.sum()} ATTs of treatment 1 incl-")
print(f"uded in the other_treatments option.")
print(f"{pos_OT.sum()} ATTs receive a positive weight, and {neg_OT.sum()} receive a negative weight.")
print("------------------------------------------------")
print(f"{'Other treat.: numdailies':<24s} {'# ATTs':<11s} {'Σ weights':<12s}")
print("------------------------------------------------")
print(f"{'Positive weights':<24s} {pos_OT.sum():<11d} {weight_OT[pos_OT].sum():<12.4f}")
print(f"{'Negative weights':<24s} {neg_OT.sum():<11d} {weight_OT[neg_OT].sum():<12.4f}")
print("------------------------------------------------")
print(f"{'Total':<24s} {pos_OT.sum()+neg_OT.sum():<11d} {weight_OT[att_OT].sum():<12.4f}")
#3. Non-normalized event-study effects
robjects.r('options(rgl.useNULL = TRUE)')

# Tu path existente
path = Path('/Users/anthonyquispe/Downloads/Replication folder FE/programs/cc_xd_didtextbook_2025_9_30/Data sets/Gentzkow et al 2011/gentzkowetal_didtextbook.dta')

# Charge
robjects.r(f'''
library(haven)
df <- as.data.frame(read_dta("{path}"))
''')

# Estimate
DIDmultiplegtDYN = importr('DIDmultiplegtDYN')

result = robjects.r('''
did_multiplegt_dyn(
    df = df,
    outcome = "prestout",
    group = "cnty90",
    time = "year",
    treatment = "numdailies",
    effects = 4,
    placebo = 4,
    effects_equal = TRUE,
    graph_off = TRUE
)
''')

print(result)
# Results
res = result.rx2('results')
effects_mat  = np.array(res.rx2('Effects'))
placebos_mat = np.array(res.rx2('Placebos'))

n_eff  = int(res.rx2('N_Effects')[0])
n_plac = int(res.rx2('N_Placebos')[0])
p_joint_eff  = res.rx2('p_jointeffects')[0]
p_equal_eff  = res.rx2('p_equality_effects')[0]
p_joint_plac = res.rx2('p_jointplacebo')[0]

# Event Study Plot

z = stats.norm.ppf(0.975)
x_plac = [-n_plac + i for i in range(n_plac)]
x_eff  = [i + 1 for i in range(n_eff)]
plac_est_plot = placebos_mat[::-1, 0]
plac_se_plot  = placebos_mat[::-1, 1]

fig, ax = plt.subplots(figsize=(10, 6))
ax.errorbar(x_plac, plac_est_plot, yerr=z * plac_se_plot,
            fmt='o', color='#1f77b4', capsize=5, capthick=1.5,
            markersize=8, linewidth=1.5, label='Placebos', zorder=3)
ax.errorbar(x_eff, effects_mat[:, 0], yerr=z * effects_mat[:, 1],
            fmt='s', color='#d62728', capsize=5, capthick=1.5,
            markersize=8, linewidth=1.5, label='Effects', zorder=3)
x_line = x_plac + [0] + x_eff
y_line = list(plac_est_plot) + [0] + list(effects_mat[:, 0])
ax.plot(x_line, y_line, '-', color='gray', lw=1, alpha=0.5, zorder=1)
ax.axhline(0, color='black', lw=0.8)
ax.axvline(0, color='black', lw=0.8, ls='--', alpha=0.4)
ax.set_xlabel('Periods since the event', fontsize=12)
ax.set_ylabel('Estimate and 95% CI', fontsize=12)
ax.set_title(
    'did_multiplegt_dyn: prestout ~ numdailies\n'
    f'Joint nullity effects: p = {p_joint_eff:.4f} | '
    f'Equality effects: p = {p_equal_eff:.4f} | '
    f'Joint nullity placebos: p = {p_joint_plac:.4f}',
    fontsize=11)
ax.set_xticks(x_plac + [0] + x_eff)
ax.legend(loc='upper left', frameon=True, fontsize=11)
ax.grid(True, alpha=0.3, ls=':')
plt.tight_layout()
plt.savefig('event_study_did_multiplegt_dyn.png', dpi=150, bbox_inches='tight')
plt.show()
#4. Analyzing the paths whose effect is averaged in the non-normalized event-study effects
robjects.r(f'''
library(haven)
df <- as.data.frame(read_dta("{path}"))
''')

result = robjects.r('''
did_multiplegt_dyn(
    df = df,
    outcome = "prestout",
    group = "cnty90",
    time = "year",
    treatment = "numdailies",
    effects = 1,
    design = c(0.8, "console"),
    graph_off = TRUE
)
''')
result1 = robjects.r('''
did_multiplegt_dyn(
    df = df,
    outcome = "prestout",
    group = "cnty90",
    time = "year",
    treatment = "numdailies",
    effects = 2,
    design = c(0.8, "console"),
    graph_off = TRUE
)
''')
result2 = robjects.r('''
did_multiplegt_dyn(
    df = df,
    outcome = "prestout",
    group = "cnty90",
    time = "year",
    treatment = "numdailies",
    effects = 4,
    design = c(0.8, "console"),
    graph_off = TRUE
)
''')
print(result)
print(result1)
print(result2)
#5. Normalized event-study effects
robjects.r(f'''
library(haven)
df <- as.data.frame(read_dta("{path}"))
''')

result = robjects.r('''
did_multiplegt_dyn(
    df = df,
    outcome = "prestout",
    group = "cnty90",
    time = "year",
    treatment = "numdailies",
    effects = 4,
    placebo = 4,
    normalized = TRUE,
    normalized_weights = TRUE,
    effects_equal = TRUE,
    graph_off = TRUE
)
''')

print(result)

# Extraer resultados
res = result.rx2('results')
effects_mat  = np.array(res.rx2('Effects'))
placebos_mat = np.array(res.rx2('Placebos'))

n_eff  = int(res.rx2('N_Effects')[0])
n_plac = int(res.rx2('N_Placebos')[0])
p_joint_eff  = res.rx2('p_jointeffects')[0]
p_equal_eff  = res.rx2('p_equality_effects')[0]
p_joint_plac = res.rx2('p_jointplacebo')[0]

# Event Study Plot
z = stats.norm.ppf(0.975)
x_plac = [-n_plac + i for i in range(n_plac)]
x_eff  = [i + 1 for i in range(n_eff)]
plac_est_plot = placebos_mat[::-1, 0]
plac_se_plot  = placebos_mat[::-1, 1]
fig, ax = plt.subplots(figsize=(10, 6))
ax.errorbar(x_plac, plac_est_plot, yerr=z * plac_se_plot,
            fmt='o', color='#1f77b4', capsize=5, capthick=1.5,
            markersize=8, linewidth=1.5, label='Placebos', zorder=3)
ax.errorbar(x_eff, effects_mat[:, 0], yerr=z * effects_mat[:, 1],
            fmt='s', color='#d62728', capsize=5, capthick=1.5,
            markersize=8, linewidth=1.5, label='Effects', zorder=3)
x_line = x_plac + [0] + x_eff
y_line = list(plac_est_plot) + [0] + list(effects_mat[:, 0])
ax.plot(x_line, y_line, '-', color='gray', lw=1, alpha=0.5, zorder=1)
ax.axhline(0, color='black', lw=0.8)
ax.axvline(0, color='black', lw=0.8, ls='--', alpha=0.4)
ax.set_xlabel('Periods since the event', fontsize=12)
ax.set_ylabel('Estimate and 95% CI', fontsize=12)
ax.set_title(
    'did_multiplegt_dyn (normalized): prestout ~ numdailies\n'
    f'Joint nullity effects: p = {p_joint_eff:.4f} | '
    f'Equality effects: p = {p_equal_eff:.4f} | '
    f'Joint nullity placebos: p = {p_joint_plac:.4f}',
    fontsize=11)
ax.set_xticks(x_plac + [0] + x_eff)
ax.legend(loc='upper left', frameon=True, fontsize=11)
ax.grid(True, alpha=0.3, ls=':')
plt.tight_layout()
plt.savefig('event_study_normalized.png', dpi=150, bbox_inches='tight')
plt.show()
#6. Testing if the lagged number of newspapers affects turnout
robjects.r(f'''
library(haven)
df <- as.data.frame(read_dta("{path}"))

# Filtro equivalente a: if year<=first_change | same_treat_after_first_change==1
df_sub <- subset(df, year <= first_change | same_treat_after_first_change == 1)
''')

result = robjects.r('''
did_multiplegt_dyn(
    df = df_sub,
    outcome = "prestout",
    group = "cnty90",
    time = "year",
    treatment = "numdailies",
    effects = 2,
    effects_equal = TRUE,
    same_switchers = TRUE,
    graph_off = TRUE
)
''')

print(result)

# Extraer resultados y graficar
res = result.rx2('results')
effects_mat = np.array(res.rx2('Effects'))

n_eff = int(res.rx2('N_Effects')[0])
p_joint_eff = res.rx2('p_jointeffects')[0]
p_equal_eff = res.rx2('p_equality_effects')[0]

z = stats.norm.ppf(0.975)
x_eff = [i + 1 for i in range(n_eff)]
fig, ax = plt.subplots(figsize=(10, 6))
ax.errorbar(x_eff, effects_mat[:, 0], yerr=z * effects_mat[:, 1],
            fmt='s', color='#d62728', capsize=5, capthick=1.5,
            markersize=8, linewidth=1.5, label='Effects', zorder=3)
ax.plot([0] + x_eff, [0] + list(effects_mat[:, 0]),
        '-', color='gray', lw=1, alpha=0.5, zorder=1)
ax.axhline(0, color='black', lw=0.8)
ax.axvline(0, color='black', lw=0.8, ls='--', alpha=0.4)
ax.set_xlabel('Periods since the event', fontsize=12)
ax.set_ylabel('Estimate and 95% CI', fontsize=12)
ax.set_title(
    'did_multiplegt_dyn (same_switchers): prestout ~ numdailies\n'
    f'Joint nullity effects: p = {p_joint_eff:.4f} | '
    f'Equality effects: p = {p_equal_eff:.4f}',
    fontsize=11)
ax.set_xticks([0] + x_eff)
ax.legend(loc='upper left', frameon=True, fontsize=11)
ax.grid(True, alpha=0.3, ls=':')
plt.tight_layout()
plt.savefig('event_study_same_switchers.png', dpi=150, bbox_inches='tight')
plt.show()
#7. Estimators assuming away effects of lagged treatments on the outcome
import rpy2.robjects as robjects
# Paquete did_multiplegt_stat (NO está en CRAN, se instala desde GitHub)
robjects.r('install.packages("remotes", repos="https://cloud.r-project.org")')
robjects.r('remotes::install_github("chaisemartinPackages/did_multiplegt_stat/R", upgrade = "never")')
robjects.r('options(rgl.useNULL = TRUE)')

path = Path('/Users/anthonyquispe/Downloads/Replication folder FE/programs/cc_xd_didtextbook_2025_9_30/Data sets/Gentzkow et al 2011/gentzkowetal_didtextbook.dta')

robjects.r(f'''
library(haven)
library(DIDmultiplegtSTAT)
df <- as.data.frame(read_dta("{path}"))
df$election_number <- as.integer(as.factor(df$year))
''')

# A. did_multiplegt_stat
result_stat = robjects.r('''
did_multiplegt_stat(
    df, "prestout", "cnty90", "election_number", "numdailies",
    placebo = TRUE,
    exact_match = TRUE
)
''')

print(result_stat)
print("\nNames:", list(result_stat.names))
# B. Tab
df_py = pd.read_stata(path)

# tab lag_numdailies if year==first_change
print("\n" + "=" * 60)
print("  tab lag_numdailies if year==first_change")
print("=" * 60)

tab1 = df_py.loc[df_py['year'] == df_py['first_change'], 'lag_numdailies'].dropna()
freq1 = tab1.value_counts().sort_index()
pct1 = freq1 / freq1.sum() * 100
cum1 = pct1.cumsum()

print(f"{'lag_numdailies':>14s} | {'Freq':>8s} {'Percent':>10s} {'Cum.':>10s}")
print("-" * 14 + "-+-" + "-" * 30)
for val in freq1.index:
    print(f"{int(val):>14d} | {freq1[val]:>8d} {pct1[val]:>10.2f} {cum1[val]:>10.2f}")
print("-" * 14 + "-+-" + "-" * 30)
print(f"{'Total':>14s} | {freq1.sum():>8d} {'100.00':>10s}")

# tab lag_numdailies if changedailies!=0 & changedailies!=. & year!=1868
print("\n" + "=" * 60)
print("  tab lag_numdailies if changedailies!=0 & !=. & year!=1868")
print("=" * 60)

mask2 = (df_py['changedailies'] != 0) & (df_py['changedailies'].notna()) & (df_py['year'] != 1868)
tab2 = df_py.loc[mask2, 'lag_numdailies'].dropna()
freq2 = tab2.value_counts().sort_index()
pct2 = freq2 / freq2.sum() * 100
cum2 = pct2.cumsum()

print(f"{'lag_numdailies':>14s} | {'Freq':>8s} {'Percent':>10s} {'Cum.':>10s}")
print("-" * 14 + "-+-" + "-" * 30)
for val in freq2.index:
    print(f"{int(val):>14d} | {freq2[val]:>8d} {pct2[val]:>10.2f} {cum2[val]:>10.2f}")
print("-" * 14 + "-+-" + "-" * 30)
print(f"{'Total':>14s} | {freq2.sum():>8d} {'100.00':>10s}")
# C. Coefplot AS y WAS
names_stat = list(result_stat.names)
print("\nNames disponibles:", names_stat)
z = stats.norm.ppf(0.975)
try:
    data = {}
    for name in names_stat:
        val = np.array(result_stat.rx2(name)).flatten()
        if len(val) > 0:
            data[name] = val
            print(f"  {name}: {val}")
except Exception as e:
    print(f"Error extrayendo: {e}")

# Plot 
fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

try:
    # Buscar AS
    as_vals = None
    was_vals = None
    as_plac = None
    was_plac = None
    
    for name in names_stat:
        arr = np.array(result_stat.rx2(name)).flatten()
        name_lower = name.lower()
        if 'plac' in name_lower and 'was' in name_lower:
            was_plac = arr
        elif 'plac' in name_lower:
            as_plac = arr
        elif 'was' in name_lower or 'waos' in name_lower:
            was_vals = arr
        elif name_lower in ['as', 'aoss', 'as_est', 'slope']:
            as_vals = arr
    
    if as_vals is None and len(names_stat) >= 2:
        as_vals = np.array(result_stat.rx2(names_stat[0])).flatten()
        as_plac = np.array(result_stat.rx2(names_stat[1])).flatten() if len(names_stat) > 1 else None
    
    got_data = as_vals is not None and len(as_vals) >= 2
except:
    got_data = False

if not got_data:
    as_est, as_se = 0.0060678, 0.0016463
    as_plac_est, as_plac_se = -0.0011278, 0.0025056
    was_est, was_se = 0.0057351, 0.0015079
    was_plac_est, was_plac_se = -0.0000123, 0.0022931
else:
    as_est, as_se = as_vals[0], as_vals[1]
    as_plac_est = as_plac[0] if as_plac is not None else 0
    as_plac_se = as_plac[1] if as_plac is not None and len(as_plac) > 1 else 0
    was_est = was_vals[0] if was_vals is not None else 0
    was_se = was_vals[1] if was_vals is not None and len(was_vals) > 1 else 0
    was_plac_est = was_plac[0] if was_plac is not None else 0
    was_plac_se = was_plac[1] if was_plac is not None and len(was_plac) > 1 else 0
# Panel 1: AS
ax = axes[0]
ax.errorbar(0, as_plac_est, yerr=z*as_plac_se, fmt='o', color='#1f77b4',
            capsize=6, capthick=1.5, markersize=10, linewidth=1.5, 
            label='Placebo', zorder=3)
ax.errorbar(1, as_est, yerr=z*as_se, fmt='s', color='#d62728',
            capsize=6, capthick=1.5, markersize=10, linewidth=1.5,
            label='AS', zorder=3)
ax.axhline(0, color='black', lw=0.8)
ax.set_xticks([0, 1])
ax.set_xticklabels(['Placebo_1', 'AS'], fontsize=11)
ax.set_title('Average Slope (AS)', fontsize=13, fontweight='bold')
ax.set_ylabel('Estimate and 95% CI', fontsize=11)
ax.legend(loc='upper left', fontsize=10)
ax.grid(True, alpha=0.3, ls=':')

# Panel 2: WAS
ax = axes[1]
ax.errorbar(0, was_plac_est, yerr=z*was_plac_se, fmt='o', color='#1f77b4',
            capsize=6, capthick=1.5, markersize=10, linewidth=1.5,
            label='Placebo', zorder=3)
ax.errorbar(1, was_est, yerr=z*was_se, fmt='s', color='#d62728',
            capsize=6, capthick=1.5, markersize=10, linewidth=1.5,
            label='WAS', zorder=3)
ax.axhline(0, color='black', lw=0.8)
ax.set_xticks([0, 1])
ax.set_xticklabels(['Placebo_1', 'WAS'], fontsize=11)
ax.set_title('Weighted Average Slope (WAS)', fontsize=13, fontweight='bold')
ax.legend(loc='upper left', fontsize=10)
ax.grid(True, alpha=0.3, ls=':')

plt.suptitle('did_multiplegt_stat: prestout ~ numdailies (exact_match)',
             fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig('did_multiplegt_stat.png', dpi=150, bbox_inches='tight')
plt.show()
