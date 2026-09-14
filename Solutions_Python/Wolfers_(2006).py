#import pandas as pd
#import numpy as np
#import statsmodels.api as sm
#from linearmodels.panel import PanelOLS
#import statsmodels.formula.api as smf
#import matplotlib.pyplot as plt
#import scipy.stats as st
#Charge Data
paths = "/Users/anthonyquispe/Downloads/Replication folder FE/programs/cc_xd_didtextbook_2025_9_30 2/Data sets/Wolfers 2006/wolfers2006_didtextbook.dta"
df = pd.read_stata(paths).copy()
# función WLS + cluster (robusta, estilo Stata)
def wls_cluster_safe(formula, data, weight_col, cluster_col, min_clusters=2):

    if data[cluster_col].nunique() < min_clusters:
        return None

    try:
        model = smf.wls(
            formula,
            data=data,
            weights=data[weight_col],
            missing="drop"
        )
    except Exception:
        return None

    idx = model.data.row_labels
    groups = data.loc[idx, cluster_col].to_numpy()

    if len(np.unique(groups)) < min_clusters:
        return None

    try:
        res = model.fit(
            cov_type="cluster",
            cov_kwds={"groups": groups}
        )
    except Exception:
        return None

    return res
*** 1) Static TWFE regression

res_1 = wls_cluster(
    "div_rate ~ udl + C(state) + C(year)",
    data=df,
    weight_col="stpop",
    cluster_col="state"
)

print(res_1.summary().tables[1])

*** 2) Decomposing the static TWFE regression

res_twfe = wls_cluster(
    "div_rate ~ udl + C(state) + C(year)",
    data=df,
    weight_col="stpop",
    cluster_col="state"
)

beta_twfe = res_twfe.params["udl"]
print("TWFE beta:", beta_twfe)
res_udl_fe = wls_cluster(
    "udl ~ C(state) + C(year)",
    data=df,
    weight_col="stpop",
    cluster_col="state"
)

df["_udl_resid"] = res_udl_fe.resid
df["_twfe_weight_raw"] = (df["_udl_resid"] ** 2) * df["stpop"]

# normalizar
df["_twfe_weight"] = df["_twfe_weight_raw"] / df["_twfe_weight_raw"].sum()
df["exposurelength"] = df["year"] - df["cohort"]
df.loc[df["cohort"] == 0, "exposurelength"] = 0
res_test = wls_cluster(
    "_twfe_weight ~ exposurelength",
    data=df,
    weight_col="stpop",
    cluster_col="state"
)

print(res_test.summary())
corr = np.corrcoef(
    df["_twfe_weight"],
    df["exposurelength"]
)[0, 1]

print("Corr(weight, exposurelength) =", corr)

*** 3) Testing randomized treatment timing
df_3 = df[(df["cohort"] != 1956) & (df["year"] <= 1968)].copy()

print(df_3.shape)   # debería dar (611, ...)
print(df_3["state"].nunique())  # debería dar 48
res_3 = wls_cluster(
    "div_rate ~ C(early_late_never)",
    data=df_3,
    weight_col="stpop",
    cluster_col="state"
)

print(res_3.summary())

*** 4) Event-study TWFE regression
rel_cols = [c for c in df.columns if c.startswith("rel_time")]

formula_4 = (
    "div_rate ~ "
    + " + ".join(rel_cols)
    + " + C(state) + C(year)"
)

res_4 = wls_cluster(
    formula_4,
    df,
    "stpop",
    "state"
)

print(res_4.summary())
#Test pre trends
import numpy as np
from scipy.stats import f

pre_cols = [
    "rel_timeminus1",
    "rel_timeminus2",
    "rel_timeminus3",
    "rel_timeminus4",
    "rel_timeminus5",
    "rel_timeminus6",
    "rel_timeminus7",
    "rel_timeminus8",
    "rel_timeminus9",
]

# matriz R
param_names = res_4.params.index.tolist()
R = np.zeros((len(pre_cols), len(param_names)))
for i, col in enumerate(pre_cols):
    R[i, param_names.index(col)] = 1

# Wald test (chi2)
wald = res_4.wald_test(R, scalar=True)

q = len(pre_cols)                  # 9
F_stat = wald.statistic / q        # 🔴 CLAVE

# grades of fredoom
n_clusters = df["state"].nunique()
df_denom = n_clusters - 1          # 50

# p-value 
p_value = 1 - f.cdf(F_stat, q, df_denom)

# OUTPUT
for i, col in enumerate(pre_cols, 1):
    print(f"( {i})  {col} = 0")

print()
print(f"F({q:3d}, {df_denom:4d}) = {F_stat:6.2f}")
print(f"Prob > F = {p_value:8.4f}")

*** 5) Decomposing the first estimated effect in the event-study TWFE regression
formula_5 = (
    "div_rate ~ rel_time1 + "
    + " + ".join([c for c in rel_cols if c != "rel_time1"])
    + " + C(state) + C(year)"
)

res_5 = wls_cluster(
    formula_5,
    df,
    "stpop",
    "state"
)

print("Coeficiente rel_time1:", res_5.params["rel_time1"])
*** 6) Sun and Abraham event-study estimators
def wls_cluster(formula, data, weight_col, cluster_col):
    model = smf.wls(
        formula,
        data=data,
        weights=data[weight_col],
        missing="drop"
    )
    idx = model.data.row_labels
    groups = data.loc[idx, cluster_col].to_numpy()
    return model.fit(cov_type="cluster", cov_kwds={"groups": groups})
df.loc[df["cohort"] == 0, "cohort"] = np.nan
def att_g_h(df, g, h):
    t = g + h

    sub = df[
        (df["year"] == t) &
        (
            (df["cohort"] == g) |
            (df["cohort"].isna()) |
            (df["cohort"] > t)
        )
    ].copy()

    if sub.empty:
        return None

    sub["D"] = (sub["cohort"] == g).astype(int)

    res = wls_cluster(
        "div_rate ~ D + C(state)",
        sub,
        "stpop",
        "state"
    )

    return res.params["D"]
horizons = list(range(-9, 17))
cohorts = sorted(df["cohort"].dropna().unique())

rows = []
for g in cohorts:
    for h in horizons:
        val = att_g_h(df, g, h)
        if val is not None:
            rows.append({"g": g, "h": h, "att": val})

att_df = pd.DataFrame(rows)
group_weights = (
    df[df["cohort"].notna()]
    .groupby("cohort")["stpop"]
    .mean()
)

att_df["w"] = att_df["g"].map(group_weights)
att_df["w"] = att_df.groupby("h")["w"].transform(lambda x: x / x.sum())
iw = (
    att_df
    .assign(w_att=lambda x: x["w"] * x["att"])
    .groupby("h")["w_att"]
    .sum()
    .sort_index()
)

print(iw)


*** 7) Estimators of Callaway and Sant'Anna	
def wls_cluster_safe(formula, data, weight_col, cluster_col, min_clusters=2):
    """
    WLS + cluster SE, pero salta submuestras inválidas
    (replica el comportamiento silencioso de Stata)
    """
    # chequear clusters suficientes
    if data[cluster_col].nunique() < min_clusters:
        return None

    # chequear variación en el regresor clave
    try:
        model = smf.wls(
            formula,
            data=data,
            weights=data[weight_col],
            missing="drop"
        )
    except Exception:
        return None

    idx = model.data.row_labels
    groups = data.loc[idx, cluster_col].to_numpy()

    # después del drop, volver a chequear clusters
    if len(np.unique(groups)) < min_clusters:
        return None

    try:
        res = model.fit(
            cov_type="cluster",
            cov_kwds={"groups": groups}
        )
    except ZeroDivisionError:
        return None

    return res
atts = []

for g in sorted(df.loc[df["cohort"] > 0, "cohort"].unique()):
    for t in sorted(df.loc[df["year"] >= g, "year"].unique()):

        sub = df[
            (df["year"] == t) &
            (
                (df["cohort"] == g) |
                (df["cohort"] == 0) |
                (df["cohort"] > t)
            )
        ].copy()

        if sub.empty:
            continue

        sub["D"] = (sub["cohort"] == g).astype(int)

        # 🚨 USAR FUNCIÓN SAFE
        res = wls_cluster_safe(
            "div_rate ~ D",
            sub,
            "stpop",
            "state"
        )

        if res is None:
            continue

        atts.append({
            "g": g,
            "t": t,
            "att": res.params["D"]
        })

att_df = pd.DataFrame(atts)
att_df["event"] = att_df["t"] - att_df["g"]

event_att = att_df.groupby("event")["att"].mean()
print(event_att)
*** 8) Estimators of de Chaisemartin and D'Haultfoeuille
def wls_cluster_safe(formula, data, weight_col, cluster_col, min_clusters=2):
    if data[cluster_col].nunique() < min_clusters:
        return None

    model = smf.wls(
        formula,
        data=data,
        weights=data[weight_col],
        missing="drop"
    )
    idx = model.data.row_labels
    groups = data.loc[idx, cluster_col].to_numpy()

    if len(np.unique(groups)) < min_clusters:
        return None

    try:
        res = model.fit(
            cov_type="cluster",
            cov_kwds={"groups": groups}
        )
    except ZeroDivisionError:
        return None

    return res
results = []

# Placebos: h = -9,...,-1
for h in range(1, 10):
    col = f"rel_timeminus{h}"
    if col not in df.columns:
        continue

    res = wls_cluster_safe(
        f"div_rate ~ {col} + C(state) + C(year)",
        df,
        "stpop",
        "state"
    )

    if res is None:
        continue
  results.append({
        "h": -h,
        "coef": res.params[col],
        "se": res.bse[col]
    })


# Efectos: h = 0,...,16
for h in range(0, 17):
    col = f"rel_time{h}"
    if col not in df.columns:
        continue

    res = wls_cluster_safe(
        f"div_rate ~ {col} + C(state) + C(year)",
        df,
        "stpop",
        "state"
    )

    if res is None:
        continue

    results.append({
        "h": h,
        "coef": res.params[col],
        "se": res.bse[col]
    })
es_df = pd.DataFrame(results).sort_values("h")

es_df["ci_low"] = es_df["coef"] - 1.96 * es_df["se"]
es_df["ci_high"] = es_df["coef"] + 1.96 * es_df["se"]
plt.figure(figsize=(10, 6))

plt.errorbar(
    es_df["h"],
    es_df["coef"],
    yerr=1.96 * es_df["se"],
    fmt="o-",
    capsize=4
)

plt.axhline(0)
plt.axvline(0)

plt.xlabel("Relative time to last period before treatment changes (t=0)")
plt.ylabel("DID estimate")
plt.title("DID, from last period before treatment changes (t=0) to t")

plt.show()

*** 9) Estimators of Borusyak et. al.
df = df.copy()

df["treated"] = (
    df["cohort"].notna() &
    (df["year"] >= df["cohort"])
).astype(int)
base_sample = df[df["treated"] == 0].copy()
res_base = wls_cluster_safe(
    "div_rate ~ C(state) + C(year)",
    base_sample,
    "stpop",
    "state"
)

if res_base is None:
    raise RuntimeError("Base regression failed")
df["y_hat"] = res_base.predict(df)
df["tau_it"] = df["div_rate"] - df["y_hat"]
rows = []

# efectos: tau0 ... tau15
for h in range(0, 16):
    sub = df[
        (df["cohort"].notna()) &
        (df["year"] == df["cohort"] + h)
    ]

    if sub.empty:
        continue

    res = wls_cluster_safe(
        "tau_it ~ 1",
        sub,
        "stpop",
        "state"
    )

    if res is None:
        continue

    rows.append({
        "name": f"tau{h}",
        "coef": res.params["Intercept"],
        "se": res.bse["Intercept"]
    })
# pre-trends: pre1 ... pre9
for h in range(1, 10):
    sub = df[
        (df["cohort"].notna()) &
        (df["year"] == df["cohort"] - h)
    ]

    if sub.empty:
        continue

    res = wls_cluster_safe(
        "tau_it ~ 1",
        sub,
        "stpop",
        "state"
    )

    if res is None:
        continue

    rows.append({
        "name": f"pre{h}",
        "coef": res.params["Intercept"],
        "se": res.bse["Intercept"]
    })
table = pd.DataFrame(rows)

# estadísticos tipo Stata (z-test)
table["z"] = table["coef"] / table["se"]
table["p_value"] = 2 * (1 - st.norm.cdf(np.abs(table["z"])))

table["ci_low"] = table["coef"] - 1.96 * table["se"]
table["ci_high"] = table["coef"] + 1.96 * table["se"]

# ordenar como Stata
def sort_key(name):
    if name.startswith("tau"):
        return (0, int(name[3:]))
    else:
        return (1, int(name[3:]))

table = table.sort_values(
    by="name",
    key=lambda x: x.map(sort_key)
).reset_index(drop=True)
pd.set_option("display.float_format", "{:.6f}".format)

table
