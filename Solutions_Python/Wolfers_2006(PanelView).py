import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Charge data
datos = pd.read_stata("/Users/anthonyquispe/Downloads/Replication folder FE/programs/cc_xd_didtextbook_2025_9_30/Data sets/Wolfers 2006/wolfers2006_didtextbook.dta")

# 1. TREATMENT STATUS PLOT
# Identify groups
ever_treated = datos.groupby('state', observed=True)['udl'].max().reset_index()
treated_ids = ever_treated[ever_treated['udl'] == 1]['state'].values
control_ids = ever_treated[ever_treated['udl'] == 0]['state'].values

# Create status: 0=Control, 1=Treated Pre, 2=Treated Post
datos['status'] = np.where(~datos['state'].isin(treated_ids), 0,
                  np.where(datos['udl'] == 1, 2, 1))

# First year of treatment per state
first_treat = datos[datos['udl'] == 1].groupby('state', observed=True)['year'].min().reset_index()
first_treat.columns = ['state', 'first_treat']

# To order: treated above by timing, controls below
order_treated = first_treat.sort_values('first_treat')['state'].values
order_all = np.concatenate([order_treated, control_ids])

# Pivot for heatmap
pivot = datos.pivot_table(index='state', columns='year', values='status', aggfunc='first', observed=True)
pivot = pivot.reindex(order_all)

fig, ax = plt.subplots(figsize=(12, 10))
cmap = plt.cm.colors.ListedColormap(['#D5D5D5', '#7BB3E0', '#2A6496'])
ax.imshow(pivot.values, aspect='auto', cmap=cmap, interpolation='none', vmin=0, vmax=2)

# Axis X
years = pivot.columns.values
xticks = np.arange(0, len(years), 5)
ax.set_xticks(xticks)
ax.set_xticklabels(years[xticks].astype(int), rotation=45, fontsize=8)
ax.set_xlabel('Year', fontsize=12)

# Axis Y
yticks = np.arange(0, len(order_all))
ax.set_yticks(yticks)
ax.set_yticklabels(order_all, fontsize=6)
ax.set_ylabel('State', fontsize=12)

# Gridlines para separar celdas
for y in np.arange(-0.5, len(order_all), 1):
    ax.axhline(y, color='white', linewidth=0.3)
for x in np.arange(-0.5, len(years), 1):
    ax.axvline(x, color='white', linewidth=0.3)

# Legend
patches = [mpatches.Patch(color='#D5D5D5', label='Control'),
           mpatches.Patch(color='#7BB3E0', label='Treated (Pre)'),
           mpatches.Patch(color='#2A6496', label='Treated (Post)')]
ax.legend(handles=patches, loc='lower center', bbox_to_anchor=(0.5, -0.1),
          ncol=3, fontsize=10)

ax.set_title('Treatment Status: Wolfers (2006)', fontsize=16)
plt.tight_layout()
plt.savefig('/Users/anthonyquispe/Documents/Libro de Science Po/Solution_Python_(PanelView)/Wolfers/PanelView_wolfers.png', dpi=150, bbox_inches='tight')
plt.show()

# 2. OUTCOMES per COHORT
# Cohorts per year of adoption
first_treat_dict = first_treat.set_index('state')['first_treat'].to_dict()
datos['cohort_year'] = datos['state'].map(first_treat_dict)

fig, ax = plt.subplots(figsize=(10, 6))

# Control (never-treated)
avg_control = datos[datos['state'].isin(control_ids)].groupby('year', observed=True)['div_rate'].mean()
ax.plot(avg_control.index, avg_control.values, color='grey', linewidth=1.5, label='Controls')

# One line per cohort
cohort_years = sorted(first_treat['first_treat'].unique())
colors_pre = plt.cm.Oranges(np.linspace(0.3, 0.7, len(cohort_years)))
colors_post = plt.cm.Reds(np.linspace(0.4, 0.9, len(cohort_years)))

for i, cy in enumerate(cohort_years):
    cohort_states = first_treat[first_treat['first_treat'] == cy]['state'].values
    cohort_data = datos[datos['state'].isin(cohort_states)].groupby('year', observed=True)['div_rate'].mean()

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
ax.set_ylabel('Divorce Rate', fontsize=12)
ax.set_title('Average Outcomes by Cohort: Wolfers (2006)', fontsize=14)
plt.tight_layout()
plt.savefig('/Users/anthonyquispe/Documents/Libro de Science Po/Solution_Python_(PanelView)/Wolfers/Cohort_Wolfers(PanelView).png', dpi=150, bbox_inches='tight')
plt.show()
