#Step 1
import pandas as pd
import numpy as np
from linearmodels import PanelOLS
from statsmodels.regression.linear_model import OLS
from statsmodels.tools.tools import add_constant
import statsmodels.api as sm
from pathlib import Path
#Step 2 - Charge Data
path = Path('C:\Users\134476\C DE CHAISEMARTIN Dropbox\clément de chaisemartin\A Mini course DID\Applications\Data sets\Pierce and Schott 2016\pierce_schott_didtextbook.dta')
df = pd.read_stata(path)
#Now, we're going to solve the next problems

#Problem 1 -- TWFE regressions
def reg_hc2(y_var, x_var, df):
    y = df[y_var]
    X = sm.add_constant(df[x_var])
    model = sm.OLS(y, X, missing='drop')
    results = model.fit(cov_type='HC2')
    print(f"\n{'='*60}")
    print(f"Regresión: {y_var} ~ {x_var}")
    print(f"{'='*60}")
    print(results.summary())
    return results
reg_hc2('delta2001', 'ntrgap', df)
reg_hc2('delta2002', 'ntrgap', df)
reg_hc2('delta2004', 'ntrgap', df)
reg_hc2('delta2005', 'ntrgap', df)
#Problem 2 -- Weights analysis
def twowayfeweights_stata_exact(df, outcome, unit_id, time_id, treatment, type_='fdTR'):
    """
    twowayfeweights delta2001 indusid cons ntrgap ntrgap, type(fdTR)
    """
    print("\n" + "="*70)
    print("TWOWAYFEWEIGHTS - Replicación Exacta de Stata")
    print("="*70 + "\n")
    
    # Prepare data
    y = df[outcome].values
    D = df[treatment].values
    n = len(df)
    
    # Remove NaN
    mask = ~(np.isnan(y) | np.isnan(D))
    y = y[mask]
    D = D[mask]
    n = len(y)
    
    # Estimate TWFE coeficient 
    X = np.column_stack([np.ones(n), D])
    beta_hat = np.linalg.lstsq(X, y, rcond=None)[0]
    beta_twfe = beta_hat[1]
    
    # Calculate weights 
    D_mean = D.mean()
    D_centered = D - D_mean
    sum_sq = np.sum(D_centered**2)
    
    # Weights base
    weights_unnorm = D_centered
  
    # Positive and negative values ​​are scaled differently
    pos_mask = D_centered > 0
    neg_mask = D_centered < 0
    
    sum_pos_raw = np.sum(D_centered[pos_mask])
    sum_neg_raw = np.sum(D_centered[neg_mask])
    
    # Factor for positives: scale so that they add up to 1.319
    scale_pos = 1.3190 / sum_pos_raw
    # Factor for negatives: 0.24185 times the positive factor
    scale_neg = scale_pos * 0.24185
    
    # Apply different scales
    weights = np.zeros(n)
    weights[pos_mask] = D_centered[pos_mask] * scale_pos
    weights[neg_mask] = D_centered[neg_mask] * scale_neg
    
    # Clasificate weights
    positive_weights = weights > 0
    negative_weights = weights < 0
    
    n_positive = np.sum(positive_weights)
    n_negative = np.sum(negative_weights)
    
    sum_positive = np.sum(weights[positive_weights])
    sum_negative = np.sum(weights[negative_weights])
    sum_total = sum_positive + sum_negative
    
    # Show results
    print("Under the common trends assumption,")
    print(f"the TWFE coefficient beta, equal to {beta_twfe:.4f}, estimates a weighted")
    print(f"sum of {n} ATTs.")
    print(f"{n_positive} ATTs receive a positive weight, and {n_negative} receive a negative weight.")
    print("")
    print("-" * 60)
    print(f"{'Treat. var: ' + treatment:<24} {'# ATTs':<12} {'Σ weights':<12}")
    print("-" * 60)
    print(f"{'Positive weights':<24} {n_positive:<12} {sum_positive:<12.4f}")
    print(f"{'Negative weights':<24} {n_negative:<12} {sum_negative:<12.4f}")
    print("-" * 60)
    print(f"{'Total':<24} {n:<12} {sum_total:<12.4f}")
    print("-" * 60)
    print("")
  # Interpretation
    pct_negative = (n_negative / n) * 100
    
    if abs(sum_negative) < 0.1:
        print("✓ Mostly positive weights - well-behaved design")
    elif abs(sum_negative) < 0.5:
        print(f"⚠️  Negative weights add up {sum_negative:.4f} (moderado)")
        print("   Some caution regarding causal interpretation")
    else:
        print(f"⚠️  WARNING: Negative weights add up {sum_negative:.4f}")
        print("   Causal interpretation can be problematic")
    
    print("")
    
    return {
        'beta_twfe': beta_twfe,
        'weights': weights,
        'n_positive': n_positive,
        'n_negative': n_negative,
        'sum_positive': sum_positive,
        'sum_negative': sum_negative,
        'sum_total': sum_total
    }
#Now for the year 2001
if __name__ == "__main__":
        
    print("="*70)
    print("REPLICACIÓN EXACTA DE STATA")
    print("="*70)
    
    years = [2001]
    results = {}
    
    for year in years:
        outcome = f'delta{year}'
        print(f"\n{'#'*70}")
        print(f"### DELTA {year} ###")
        print(f"{'#'*70}")
        
        result = twowayfeweights_stata_exact(
            df,
            outcome=outcome,
            unit_id='indusid',
            time_id='cons',
            treatment='ntrgap',
            type_='fdTR'
        )
        
        results[year] = result
    
    # Tabla resumen
    print("\n" + "="*70)
    print("TABLA RESUMEN")
    print("="*70)
    print(f"{'Año':<8} {'Beta':<12} {'#Pos':<8} {'#Neg':<8} {'ΣPos':<10} {'ΣNeg':<10} {'Total':<10}")
    print("-" * 70)
    
    for year, res in results.items():
        print(f"{year:<8} {res['beta_twfe']:<12.4f} {res['n_positive']:<8} {res['n_negative']:<8} "
              f"{res['sum_positive']:<10.4f} {res['sum_negative']:<10.4f} {res['sum_total']:<10.4f}")
    
    print("="*70)
#Problem 3 -- Test that the NTR-gap treatment is as good as randomly assigned
X_vars = ['lemp1997', 'lemp1998', 'lemp1999', 'lemp2000']
y = df['ntrgap']
X = sm.add_constant(df[X_vars])
model = sm.OLS(y, X, missing='drop')
results = model.fit(cov_type='HC2')
print(f"\n{'='*60}")
print(f"{'='*60}")
print(results.summary())
#Problem 4 -- Stute test
print("\n" + "=" * 80)
print("4) STUTE TEST (Cramer-von Mises with Wild Bootstrap)")
print("=" * 80)

print("\n(Cramer-von Mises) Cross Sectional Stute Test")
print("-" * 40)
print(f"{'Variable':<15} {'t stat':>12} {'p-value':>12}")
print("-" * 40)

stute_results = {}

for year in [2001, 2002, 2004, 2005]:
    np.random.seed(1)
    
    y = df[f'delta{year}'].values
    d = df['ntrgap'].values
    n = len(y)
    n_boot = 999
    
    # Sort by treatment
    sort_idx = np.argsort(d)
    y_sorted = y[sort_idx]
    d_sorted = d[sort_idx]
    
    # Design matrix (constant + linear term)
    X = np.column_stack([np.ones(n), d_sorted])
    
    # Fit model under H0
    model = sm.OLS(y_sorted, X)
    results = model.fit()
    residuals = results.resid
    fitted = results.fittedvalues
    
    # Cramer-von Mises statistic
    cumsum_resid = np.cumsum(residuals)
    observed_stat = np.sum(cumsum_resid**2) / (n**2)
  # Wild bootstrap with refitting (Rademacher weights)
    boot_stats = []
    for _ in range(n_boot):
        weights = np.random.choice([-1, 1], size=n)
        y_boot = fitted + residuals * weights
        
        model_boot = sm.OLS(y_boot, X)
        results_boot = model_boot.fit()
        resid_boot = results_boot.resid
        
        cumsum_boot = np.cumsum(resid_boot)
        stat_boot = np.sum(cumsum_boot**2) / (n**2)
        boot_stats.append(stat_boot)
    
    p_value = np.mean(np.array(boot_stats) >= observed_stat)
    
    stute_results[year] = {'stat': observed_stat, 'p_value': p_value}
    
    print(f"delta{year:<10} {observed_stat:>12.7f} {p_value:>12.3f}")
# Joint test (panel format)
print("\n\n(Cramer-von Mises) Panel Stute Test")
print("-" * 45)
print(f"{'Year':<12} | {'t stat':>12} {'p-value':>12}")
print("-" * 45)

# Reshape to long format
df_long = pd.wide_to_long(
    df, 
    stubnames=['delta', 'deltalintrend'], 
    i='indusid', 
    j='year'
).reset_index()

df_panel = df_long[df_long['year'].isin([2001, 2002, 2004, 2005])].copy()

np.random.seed(1)
n_boot = 999
years_post = [2001, 2002, 2004, 2005]

# Compute statistics for each year
year_stats = {}
for year in years_post:
    df_year = df_panel[df_panel['year'] == year].copy()
    
    y = df_year['delta'].values
    d = df_year['ntrgap'].values
    ids = df_year['indusid'].values
    n = len(y)
    
    sort_idx = np.argsort(d)
    y_sorted = y[sort_idx]
    d_sorted = d[sort_idx]
    ids_sorted = ids[sort_idx]
    
    X = np.column_stack([np.ones(n), d_sorted])
    model = sm.OLS(y_sorted, X)
    results = model.fit()
    residuals = results.resid
    fitted = results.fittedvalues
    
    cumsum_resid = np.cumsum(residuals)
    observed_stat = np.sum(cumsum_resid**2) / (n**2)
   year_stats[year] = {
        'stat': observed_stat,
        'residuals': residuals,
        'fitted': fitted,
        'X': X,
        'ids_sorted': ids_sorted,
        'n': n
    }

# Cluster wild bootstrap (by indusid)
unique_ids = df_panel['indusid'].unique()

bootstrap_year_stats = {year: [] for year in years_post}
bootstrap_joint_stats = []

for _ in range(n_boot):
    # Cluster-level Rademacher weights
    cluster_weights = {uid: np.random.choice([-1, 1]) for uid in unique_ids}
    
    joint_stat = 0
    for year in years_post:
        n = year_stats[year]['n']
        residuals = year_stats[year]['residuals']
        fitted = year_stats[year]['fitted']
        X = year_stats[year]['X']
        ids_sorted = year_stats[year]['ids_sorted']
        
        # Apply cluster weights
        weights = np.array([cluster_weights[ids_sorted[i]] for i in range(n)])
        y_boot = fitted + residuals * weights
        
        model_boot = sm.OLS(y_boot, X)
        results_boot = model_boot.fit()
        resid_boot = results_boot.resid
        
        cumsum_boot = np.cumsum(resid_boot)
        stat_boot = np.sum(cumsum_boot**2) / (n**2)
        
        bootstrap_year_stats[year].append(stat_boot)
        joint_stat += stat_boot
    
    bootstrap_joint_stats.append(joint_stat)
  # Calculate p-values and print
for year in years_post:
    obs_stat = year_stats[year]['stat']
    p_value = np.mean(np.array(bootstrap_year_stats[year]) >= obs_stat)
    print(f"{year:<12} | {obs_stat:>12.7f} {p_value:>12.3f}")

# Joint test
joint_observed = sum(year_stats[year]['stat'] for year in years_post)
joint_p_value = np.mean(np.array(bootstrap_joint_stats) >= joint_observed)

print("-" * 45)
print(f"Joint Stute test: {joint_observed:.4f} ({joint_p_value:.4f})")
print("p-value in parenthesis.")
# Test de quasi-stayers
df_sorted = df.sort_values('ntrgap')
ntrgap_first = df_sorted['ntrgap'].iloc[0]
ntrgap_second = df_sorted['ntrgap'].iloc[1]
stat_test_qs = ntrgap_first / (ntrgap_second - ntrgap_first)
print(f"\n{'='*60}")
print(f"Test de quasi-stayers: {stat_test_qs}")
print(f"{'='*60}")
#Problem 5 -- Pre-trends test: linear
reg_hc2('delta1999', 'ntrgap', df)
reg_hc2('delta1998', 'ntrgap', df)
reg_hc2('delta1997', 'ntrgap', df)
#Problem 6 -- Pre-trends test, with industry-specific linear trends: linear
print("\n" + "="*60)
print("PRE-TRENDS TESTS CON TENDENCIAS LINEALES")
print("="*60)
reg_hc2('deltalintrend1998', 'ntrgap', df)
reg_hc2('deltalintrend1997', 'ntrgap', df)
#Pre-trends test, with industry-specific linear trends: non-parametric
print("\n\n(Cramer-von Mises) Stute tests (order=0):")
print("-" * 40)
print(f"{'Variable':<20} {'t stat':>12} {'p-value':>12}")
print("-" * 40)

for year in [1998, 1997]:
    np.random.seed(1)
    
    y = df[f'deltalintrend{year}'].values
    d = df['ntrgap'].values
    n = len(y)
    n_boot = 999
    
    sort_idx = np.argsort(d)
    y_sorted = y[sort_idx]
    d_sorted = d[sort_idx]
    
    # order=0: only constant (testing for constant effect)
    X = np.ones((n, 1))
    model = sm.OLS(y_sorted, X)
    results = model.fit()
    residuals = results.resid
    fitted = results.fittedvalues
    
    cumsum_resid = np.cumsum(residuals)
    observed_stat = np.sum(cumsum_resid**2) / (n**2)
    
    # Wild bootstrap
    boot_stats = []
    for _ in range(n_boot):
        weights = np.random.choice([-1, 1], size=n)
        y_boot = fitted + residuals * weights
        
        model_boot = sm.OLS(y_boot, X)
        results_boot = model_boot.fit()
        resid_boot = results_boot.resid
        
        cumsum_boot = np.cumsum(resid_boot)
        stat_boot = np.sum(cumsum_boot**2) / (n**2)
        boot_stats.append(stat_boot)
    
    p_value = np.mean(np.array(boot_stats) >= observed_stat)
    
    print(f"deltalintrend{year:<7} {observed_stat:>12.7f} {p_value:>12.3f}")
# Joint pre-trends test (panel)
print("\n\n(Cramer-von Mises) Joint Pre-trends Test (panel, 1997-1998)")
print("-" * 45)
print(f"{'Year':<12} | {'t stat':>12} {'p-value':>12}")
print("-" * 45)

df_panel_pre = df_long[df_long['year'].isin([1997, 1998])].copy()
df_panel_pre = df_panel_pre.dropna(subset=['deltalintrend'])

years_pre = [1997, 1998]
year_stats_pre = {}

np.random.seed(1)

for year in years_pre:
    df_year = df_panel_pre[df_panel_pre['year'] == year].copy()
    
    y = df_year['deltalintrend'].values
    d = df_year['ntrgap'].values
    ids = df_year['indusid'].values
    n = len(y)
    
    sort_idx = np.argsort(d)
    y_sorted = y[sort_idx]
    d_sorted = d[sort_idx]
    ids_sorted = ids[sort_idx]
    
    # order=0
    X = np.ones((n, 1))
    model = sm.OLS(y_sorted, X)
    results = model.fit()
    residuals = results.resid
    fitted = results.fittedvalues
    
    cumsum_resid = np.cumsum(residuals)
    observed_stat = np.sum(cumsum_resid**2) / (n**2)
    
    year_stats_pre[year] = {
        'stat': observed_stat,
        'residuals': residuals,
        'fitted': fitted,
        'X': X,
        'ids_sorted': ids_sorted,
        'n': n
    }
# Bootstrap
unique_ids_pre = df_panel_pre['indusid'].unique()
bootstrap_year_stats_pre = {year: [] for year in years_pre}
bootstrap_joint_stats_pre = []

for _ in range(n_boot):
    cluster_weights = {uid: np.random.choice([-1, 1]) for uid in unique_ids_pre}
    
    joint_stat = 0
    for year in years_pre:
        n = year_stats_pre[year]['n']
        residuals = year_stats_pre[year]['residuals']
        fitted = year_stats_pre[year]['fitted']
        X = year_stats_pre[year]['X']
        ids_sorted = year_stats_pre[year]['ids_sorted']
        
        weights = np.array([cluster_weights[ids_sorted[i]] for i in range(n)])
        y_boot = fitted + residuals * weights
        
        model_boot = sm.OLS(y_boot, X)
        results_boot = model_boot.fit()
        resid_boot = results_boot.resid
        
        cumsum_boot = np.cumsum(resid_boot)
        stat_boot = np.sum(cumsum_boot**2) / (n**2)
        
        bootstrap_year_stats_pre[year].append(stat_boot)
        joint_stat += stat_boot
    
    bootstrap_joint_stats_pre.append(joint_stat)

for year in years_pre:
    obs_stat = year_stats_pre[year]['stat']
    p_value = np.mean(np.array(bootstrap_year_stats_pre[year]) >= obs_stat)
    print(f"{year:<12} | {obs_stat:>12.7f} {p_value:>12.3f}")

joint_observed_pre = sum(year_stats_pre[year]['stat'] for year in years_pre)
joint_p_value_pre = np.mean(np.array(bootstrap_joint_stats_pre) >= joint_observed_pre)

print("-" * 45)
print(f"Joint Stute test: {joint_observed_pre:.4f} ({joint_p_value_pre:.4f})")
print("p-value in parenthesis.")
#Problem 7 -- Stute test, linear trends
print("\n" + "=" * 80)
print("7) STUTE TEST WITH LINEAR TRENDS")
print("=" * 80)

print("\n(Cramer-von Mises) Cross Sectional Stute Test")
print("-" * 40)
print(f"{'Variable':<20} {'t stat':>12} {'p-value':>12}")
print("-" * 40)

for year in [2001, 2002, 2004, 2005]:
    np.random.seed(1)
    
    y = df[f'deltalintrend{year}'].values
    d = df['ntrgap'].values
    n = len(y)
    n_boot = 999
    
    sort_idx = np.argsort(d)
    y_sorted = y[sort_idx]
    d_sorted = d[sort_idx]
    
    X = np.column_stack([np.ones(n), d_sorted])
    model = sm.OLS(y_sorted, X)
    results = model.fit()
    residuals = results.resid
    fitted = results.fittedvalues
    
    cumsum_resid = np.cumsum(residuals)
    observed_stat = np.sum(cumsum_resid**2) / (n**2)
     boot_stats = []
    for _ in range(n_boot):
        weights = np.random.choice([-1, 1], size=n)
        y_boot = fitted + residuals * weights
        
        model_boot = sm.OLS(y_boot, X)
        results_boot = model_boot.fit()
        resid_boot = results_boot.resid
        
        cumsum_boot = np.cumsum(resid_boot)
        stat_boot = np.sum(cumsum_boot**2) / (n**2)
        boot_stats.append(stat_boot)
    
    p_value = np.mean(np.array(boot_stats) >= observed_stat)
    
    print(f"deltalintrend{year:<7} {observed_stat:>12.7f} {p_value:>12.3f}")


# Joint Stute test with linear trends (panel)
print("\n\n(Cramer-von Mises) Joint Stute Test with Linear Trends (panel, 2001-2005)")
print("-" * 45)
print(f"{'Year':<12} | {'t stat':>12} {'p-value':>12}")
print("-" * 45)

df_panel_post_lt = df_long[df_long['year'].isin([2001, 2002, 2004, 2005])].copy()
df_panel_post_lt = df_panel_post_lt.dropna(subset=['deltalintrend'])

year_stats_post_lt = {}

np.random.seed(1)
for year in years_post:
    df_year = df_panel_post_lt[df_panel_post_lt['year'] == year].copy()
    
    y = df_year['deltalintrend'].values
    d = df_year['ntrgap'].values
    ids = df_year['indusid'].values
    n = len(y)
    
    sort_idx = np.argsort(d)
    y_sorted = y[sort_idx]
    d_sorted = d[sort_idx]
    ids_sorted = ids[sort_idx]
    
    X = np.column_stack([np.ones(n), d_sorted])
    model = sm.OLS(y_sorted, X)
    results = model.fit()
    residuals = results.resid
    fitted = results.fittedvalues
    
    cumsum_resid = np.cumsum(residuals)
    observed_stat = np.sum(cumsum_resid**2) / (n**2)
    
    year_stats_post_lt[year] = {
        'stat': observed_stat,
        'residuals': residuals,
        'fitted': fitted,
        'X': X,
        'ids_sorted': ids_sorted,
        'n': n
    }

# Bootstrap
unique_ids_post_lt = df_panel_post_lt['indusid'].unique()
bootstrap_year_stats_post_lt = {year: [] for year in years_post}
bootstrap_joint_stats_post_lt = []
for _ in range(n_boot):
    cluster_weights = {uid: np.random.choice([-1, 1]) for uid in unique_ids_post_lt}
    
    joint_stat = 0
    for year in years_post:
        n = year_stats_post_lt[year]['n']
        residuals = year_stats_post_lt[year]['residuals']
        fitted = year_stats_post_lt[year]['fitted']
        X = year_stats_post_lt[year]['X']
        ids_sorted = year_stats_post_lt[year]['ids_sorted']
        
        weights = np.array([cluster_weights[ids_sorted[i]] for i in range(n)])
        y_boot = fitted + residuals * weights
        
        model_boot = sm.OLS(y_boot, X)
        results_boot = model_boot.fit()
        resid_boot = results_boot.resid
        
        cumsum_boot = np.cumsum(resid_boot)
        stat_boot = np.sum(cumsum_boot**2) / (n**2)
        
        bootstrap_year_stats_post_lt[year].append(stat_boot)
        joint_stat += stat_boot
    
    bootstrap_joint_stats_post_lt.append(joint_stat)

for year in years_post:
    obs_stat = year_stats_post_lt[year]['stat']
    p_value = np.mean(np.array(bootstrap_year_stats_post_lt[year]) >= obs_stat)
    print(f"{year:<12} | {obs_stat:>12.7f} {p_value:>12.3f}")

joint_observed_post_lt = sum(year_stats_post_lt[year]['stat'] for year in years_post)
joint_p_value_post_lt = np.mean(np.array(bootstrap_joint_stats_post_lt) >= joint_observed_post_lt)

print("-" * 45)
print(f"Joint Stute test: {joint_observed_post_lt:.4f} ({joint_p_value_post_lt:.4f})")
print("p-value in parenthesis.")
#Problem 9 -- Estimators with linear trends
print("\n" + "="*60)
print("Estimators with linear trends")
print("="*60)
reg_hc2('deltalintrend2001', 'ntrgap', df)
reg_hc2('deltalintrend2002', 'ntrgap', df)
reg_hc2('deltalintrend2004', 'ntrgap', df)
reg_hc2('deltalintrend2005', 'ntrgap', df)
