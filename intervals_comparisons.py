# This script evaluate the coverage of several methods
# and generates figure 4 and 5.

import numpy as np
import pandas as pd
from scipy.stats import norm
import seaborn as sns
import matplotlib.pyplot as plt
from calculate_intervals import bonferroni_interval, sidak_interval
from calculate_intervals import fuentes_intervals, fuentes_symmetric_interval
from calculate_intervals import SoS_shortest_interval, SoS_intervals, SoS_all_interval

alpha = 0.05
m = 100
k_values = np.arange(1,m)
sns.set(style="darkgrid", rc={'figure.figsize':(11,8)}, font_scale=1.75)

# %%

SoS_symmetric = []
SoS_shortest = []
bonferroni = []
sidak = []
fuentes_shortest = []
fuentes_symmetric = []

for k in np.arange(1,100):
    SoS_symmetric.append(SoS_intervals(k,m, delta = float(m)/(k+m)))
    SoS_shortest.append(SoS_shortest_interval(k,m))
    fuentes_shortest.append(fuentes_intervals(k,m))
    fuentes_symmetric.append(fuentes_symmetric_interval(k,m))
    bonferroni.append(bonferroni_interval(m))
    sidak.append(sidak_interval(m))

SoS_symmetric_int = np.array([u[0] for u in SoS_symmetric])
SoS_shortest_int = np.array([u[0] for u in SoS_shortest])
fuentes_shortest_int = np.array([u[0] for u in fuentes_shortest])
fuentes_symmetric_int = np.array([u[0] for u in fuentes_symmetric])
bonferroni_shortest_int = np.array([u[0] for u in bonferroni])
sidak_intervals = np.array([u[0] for u in sidak])

# %%
d = pd.DataFrame({'SoS symmetric':SoS_symmetric_int,
              'SoS shortest':SoS_shortest_int,
              'FCW symmetric':fuentes_symmetric_int,
              'FCW shortest':fuentes_shortest_int,
              'Bonferroni m':bonferroni_shortest_int,
              'Sidak m': sidak_intervals,
              'k':k_values,
                })
    
df = pd.melt(d, id_vars='k')
df.columns = ['k', 'Method', 'value']

plt.figure()  
sns_plot = sns.lineplot(x="k", y="value",
             hue="Method", style="Method",linewidth=2,
             data=df)
sns_plot.set(xlabel='$k$', ylabel='Length', title="Selection of $k$ out of $m=100$")
handles, labels = sns_plot.get_legend_handles_labels()
sns_plot.legend(handles=handles[1:], labels=labels[1:])
plt.savefig("figures/top_out_of_m_intervals.png", dpi=300)

# %%
k_20_m_100, _, deltas = SoS_all_interval(20,100)
k_10_m_1000, _, deltas = SoS_all_interval(10,1000)
k_10_m_100, _, deltas = SoS_all_interval(10,100)
k_20_m_1000, _, deltas = SoS_all_interval(20,1000)

min_values = [(deltas[u.argmin()],u.min()) for u in [k_10_m_100, k_10_m_1000, k_20_m_100, k_20_m_1000]]

tmp_df = pd.DataFrame({'k = 20, m = 100':k_20_m_100, 
                       'k = 10, m = 1000':k_10_m_1000,
                       'k = 10, m = 100':k_10_m_100,
                       'k = 20, m = 1000':k_20_m_1000,
                       'Delta':deltas,
                       })

df = pd.melt(tmp_df,id_vars='Delta')
df.columns = ['Delta', 'Selection', 'value']

plt.figure()
sns_plot = sns.lineplot(x="Delta", y="value",
             hue="Selection", style="Selection",linewidth=2,
             data=df)
sns_plot.set(xlabel=r'$\delta$', ylabel='Length', title="Interval length as a function of $\delta$")
handles, labels = sns_plot.get_legend_handles_labels()
sns_plot.legend(handles=handles[1:], labels=labels[1:])
plt.scatter(np.array(min_values)[:,0], np.array(min_values)[:,1], c=sns.color_palette()[:4], marker='o', s=100)

plt.savefig("figures/delta_interval_length.png", dpi=300)
