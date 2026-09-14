import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Charge data
datos = pd.read_stata("/Users/anthonyquispe/Downloads/Replication folder FE/programs/cc_xd_didtextbook_2025_9_30/Data sets/Gentzkow et al 2011/gentzkowetal_didtextbook.dta")

# Create treat_timing (captures staggered adoption)
datos['treat_timing'] = np.where(
    datos['first_change'] >= 1e+06, 0,
    np.where(datos['year'] >= datos['first_change'], 1, 0)
)

# Identify groups
ever_treated = datos.groupby('cnty90', observed=True)['treat_timing'].max().reset_index()
treated_ids = ever_treated[ever_treated['treat_timing'] == 1]['cnty90'].values
control_ids = ever_treated[ever_treated['treat_timing'] == 0]['cnty90'].values

# First year of treatment per county
first_treat = datos[datos['treat_timing'] == 1].groupby('cnty90', observed=True)['year'].min().reset_index()
first_treat.columns = ['cnty90', 'first_treat']

# Create status: 0=Control, 1=Treated Pre, 2=Treated Post
datos['status'] = np.where(~datos['cnty90'].isin(treated_ids), 0,
                  np.where(datos['treat_timing'] == 1, 2, 1))

# PLOT 1a: TREATMENT STATUS - ALL COUNTIES (1195)

def plot_treatment_status(datos, subset_ids, title_suffix="", filename_suffix="", axis_y_gap=100):
    order_treated_df = first_treat[first_treat['cnty90'].isin(subset_ids)].sort_values('first_treat')
    order_treated = order_treated_df['cnty90'].values
    order_control = [c for c in control_ids if c in subset_ids]
    order_all = np.concatenate([order_treated, order_control])

    pivot = datos[datos['cnty90'].isin(subset_ids)].pivot_table(
        index='cnty90', columns='year', values='status', aggfunc='first', observed=True
    )
    pivot = pivot.reindex(order_all)

    years = pivot.columns.values
    n_units = len(order_all)

    fig, ax = plt.subplots(figsize=(12, 10))
    cmap = plt.cm.colors.ListedColormap(['#D5D5D5', '#7BB3E0', '#2A6496'])
    ax.imshow(pivot.values, aspect='auto', cmap=cmap, interpolation='none', vmin=0, vmax=2)

    # Axis X
    xticks = np.arange(0, len(years), 4)
    ax.set_xticks(xticks)
    ax.set_xticklabels(years[xticks].astype(int), rotation=45, fontsize=8)
    ax.set_xlabel('Year', fontsize=12)

    # Axis Y
    ytick_step = axis_y_gap
    yticks = np.arange(0, n_units, ytick_step)
    ax.set_yticks(yticks)
    ax.set_yticklabels([str(order_all[i]) for i in yticks], fontsize=5)
    ax.set_ylabel('County', fontsize=12)

    # Gridlines
    for y in np.arange(-0.5, n_units, 1):
        ax.axhline(y, color='white', linewidth=0.15)
    for x in np.arange(-0.5, len(years), 1):
        ax.axvline(x, color='white', linewidth=0.3)

    # Legend
    patches = [mpatches.Patch(color='#D5D5D5', label='Control'),
               mpatches.Patch(color='#7BB3E0', label='Treated (Pre)'),
               mpatches.Patch(color='#2A6496', label='Treated (Post)')]
    ax.legend(handles=patches, loc='lower center', bbox_to_anchor=(0.5, -0.1),
              ncol=3, fontsize=10)

    ax.set_title(f'Treatment Status: Gentzkow et al. (2011){title_suffix}', fontsize=16)
    plt.tight_layout()
    plt.savefig(f'/Users/anthonyquispe/Documents/Libro de Science Po/Solution_Python_(PanelView)/Gentzkow_et_all/PanelView_gentzkow{filename_suffix}.png',
                dpi=150, bbox_inches='tight')
    plt.show()


# All counties
all_ids = np.concatenate([treated_ids, control_ids])
plot_treatment_status(datos, all_ids, title_suffix="", filename_suffix="_all", axis_y_gap=100)

# With 100 units (60 treated + 40 control)
subs = np.concatenate([treated_ids[:60], control_ids[:40]])
plot_treatment_status(datos, subs, title_suffix=" (100 units)", filename_suffix="_100", axis_y_gap=10)

# PLOT 2: OUTCOMES PER COHORT

first_treat_dict = first_treat.set_index('cnty90')['first_treat'].to_dict()
datos['cohort_year'] = datos['cnty90'].map(first_treat_dict)

fig, ax = plt.subplots(figsize=(10, 6))

# Control (never-treated)
avg_control = datos[datos['cnty90'].isin(control_ids)].groupby('year', observed=True)['prestout'].mean()
ax.plot(avg_control.index, avg_control.values, color='grey', linewidth=1.5, label='Controls')

# One line per cohort
cohort_years = sorted(first_treat['first_treat'].unique())
colors_pre = plt.cm.Oranges(np.linspace(0.3, 0.7, len(cohort_years)))
colors_post = plt.cm.Reds(np.linspace(0.4, 0.9, len(cohort_years)))

for i, cy in enumerate(cohort_years):
    cohort_counties = first_treat[first_treat['first_treat'] == cy]['cnty90'].values
    cohort_data = datos[datos['cnty90'].isin(cohort_counties)].groupby('year', observed=True)['prestout'].mean()

    pre = cohort_data[cohort_data.index < cy]
    post = cohort_data[cohort_data.index >= cy]

    ax.plot(pre.index, pre.values, color=colors_pre[i], linewidth=1.2)
    ax.plot(post.index, post.values, color=colors_post[i], linewidth=1.2)

# Manual legend
patches = [mpatches.Patch(color='grey', label='Controls'),
           mpatches.Patch(color='#F4A460', label='Treated (Pre)'),
           mpatches.Patch(color='#CC0000', label='Treated (Post)')]
ax.legend(handles=patches, fontsize=10)

ax.set_xlabel('Year', fontsize=12)
ax.set_ylabel('Presidential Turnout', fontsize=12)
ax.set_title('Average Outcomes by Cohort: Gentzkow et al. (2011)', fontsize=14)
plt.tight_layout()
plt.savefig('/Users/anthonyquispe/Documents/Libro de Science Po/Solution_Python_(PanelView)/Gentzkow_et_all/Cohort_Gentzkow(PanelView).png',
            dpi=150, bbox_inches='tight')
plt.show()
