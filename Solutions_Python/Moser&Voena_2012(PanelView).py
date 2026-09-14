import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Charge data
datos = pd.read_stata("/Users/anthonyquispe/Downloads/Replication folder FE/programs/cc_xd_didtextbook_2025_9_30/Data sets/Moser and Voena 2012/moser_voena_didtextbook.dta")
# We have to create a treatment's variable
datos['D'] = datos['treatmentgroup'] * datos['post']

# 1. TREATMENT STATUS PLOT
# Create status
datos['status'] = np.where(datos['D'] == 1, 'Treatment', 'Control')

# To order: treated above by timing, controls below
order = datos.groupby('subclass')['D'].max().reset_index()
order.columns = ['subclass', 'ever_treated']
order = order.sort_values(['ever_treated', 'subclass'], ascending=[False, True])

# Top 400 units
top400 = order.head(400)['subclass'].values
df = datos[datos['subclass'].isin(top400)].copy()

# Pivot for heatmap
pivot = df.pivot_table(index='subclass', columns='year', values='D', aggfunc='first')
pivot = pivot.reindex(top400)

fig, ax = plt.subplots(figsize=(12, 10))
cmap = plt.cm.colors.ListedColormap(['#D5D5D5', '#4A90D9'])
ax.imshow(pivot.values, aspect='auto', cmap=cmap, interpolation='none')

# Axis X
years = pivot.columns.values
xticks = np.arange(0, len(years), 5)
ax.set_xticks(xticks)
ax.set_xticklabels(years[xticks].astype(int), rotation=45, fontsize=8)
ax.set_xlabel('Year', fontsize=12)

# Axis Y
yticks = np.arange(0, len(top400), 10)
ax.set_yticks(yticks)
ax.set_yticklabels(top400[yticks].astype(int), fontsize=4)
ax.set_ylabel('Subclass', fontsize=12)

# Leyend
patches = [mpatches.Patch(color='#D5D5D5', label='Control'),
           mpatches.Patch(color='#4A90D9', label='Treatment')]
ax.legend(handles=patches, loc='lower center', bbox_to_anchor=(0.5, -0.12),
          ncol=2, fontsize=10)

ax.set_title('Treatment Status (Top 400 units, sorted)', fontsize=16)
plt.tight_layout()
plt.savefig('/Users/anthonyquispe/Documents/Libro de Science Po/Solution_Python_(PanelView)/Moser&Voena/treatment_status_moser.png', dpi=150, bbox_inches='tight')
plt.show()

# 2. OUTCOMES POR COHORTE
# Identify cohorts
treated = datos[datos['treatmentgroup'] == 1]
control = datos[datos['treatmentgroup'] == 0]

# Average per year
avg_treated_post = treated[treated['post'] == 1].groupby('year')['patents'].mean()
avg_treated_pre = treated[treated['post'] == 0].groupby('year')['patents'].mean()
avg_control = control.groupby('year')['patents'].mean()

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(avg_control.index, avg_control.values, color='grey', linewidth=1.5, label='Controls')
ax.plot(avg_treated_pre.index, avg_treated_pre.values, color='#F4A460', linewidth=1.5, label='Treated (Pre)')
ax.plot(avg_treated_post.index, avg_treated_post.values, color='#CC0000', linewidth=1.5, label='Treated (Post)')

ax.axvline(x=1919, color='black', linestyle='--', alpha=0.3, label='TWEA (1919)')
ax.set_xlabel('Year', fontsize=12)
ax.set_ylabel('Patents', fontsize=12)
ax.set_title('Average Outcomes by Cohort: Moser & Voena (2012)', fontsize=14)
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig('/Users/anthonyquispe/Documents/Libro de Science Po/Solution_Python_(PanelView)/Moser&Voena/cohort_moser.png', dpi=150, bbox_inches='tight')
plt.show()
