# This script calculate dependency simulations for several methods. 
# Final results is Figure 11, but the script is also calculating the 
# coverage other methods. Might take several minutes. 

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
k = 10
sns.set(style="darkgrid", rc={'figure.figsize':(11,7)}, font_scale=1.75)

# %%

def create_AR_cov_mat(m, rho):
    cov_mat = np.zeros((m,m))
    for i in np.arange(m):
        for j in np.arange(m):
            cov_mat[i,j] = rho**(np.abs(i-j))
    return cov_mat

def create_model_7_cov_mat(m):
    sigma = np.zeros((m,m))
    for i in np.arange(m):
        for j in np.arange(m):
            if i==j:
                sigma[i,j]  = 1
            else:
                sigma[i,j] = (0.5)/((np.abs(i-j))**(5))
                
    d = np.diag(np.random.uniform(1,3,m))
    d_square = np.sqrt(d)
    cov_mat = np.matmul(np.matmul(d_square,sigma),d_square)
    return cov_mat

def create_block_cov_mat(m, block_values = np.array([0,0.2,0.5,0.75,0.9])):
    sigma = np.zeros((m,m))
    num_corr_values = len(block_values)
    steps = int(m/num_corr_values)
    for k in np.arange(num_corr_values):
        sigma[k*steps:(k+1)*(steps),k*steps:(k+1)*(steps)] = block_values[k]
    np.fill_diagonal(sigma,1)
    return sigma

def simulate_diff(mu, cov_mat, num_samples = 150000):
    y = np.random.multivariate_normal(mu, cov_mat, size=num_samples)
    y_max_ix = y.argmax(1)
    y_max = y.max(1)
    
    max_diff = (y_max - mu[y_max_ix])
    return max_diff

def simulate_top_k_diff(mu, cov_mat, k, num_samples = 50000):
    y = np.random.multivariate_normal(mu, cov_mat, size=num_samples)
    sort_ix = y.argsort(1)
    sort_ix.shape
    y[sort_ix][0]
    
    y_sort = np.take_along_axis(y, sort_ix, axis=1)
    mu_sort = mu[sort_ix]
    y_diff = y_sort - mu_sort
    
    top_k_diff = y_diff[:,-k:]
    return top_k_diff

def top_k_coverage(top_k_diff, interval):
    coverage = np.mean(np.logical_and(top_k_diff.min(1) > -interval[1], top_k_diff.max(1) < interval[0]))
    return coverage

def cov_probability(max_diff, interval):
    coverage = np.mean(np.logical_and(max_diff > -interval[1], max_diff < interval[0]))
    return coverage

def interval_coverage(mu, cov_matrix, k):
    m = mu.shape[0]
    
    _,SoS_interval = SoS_shortest_interval(k,m)
    _,SoS_sym_interval = SoS_intervals(k,m, delta = float(m)/(k+m))
    _,fuentes_int = fuentes_intervals(k,m)
    _,fuentes_sym_int = fuentes_symmetric_interval(k,m)
    _,bonferroni_int = bonferroni_interval(m)
    _,sidak_int = sidak_interval(m)
    
    top_k_diff = simulate_top_k_diff(mu, cov_matrix, k)
    SoS_coverage = top_k_coverage(top_k_diff, SoS_interval)
    SoS_sym_coverage = top_k_coverage(top_k_diff, SoS_sym_interval)
    fuentes_coverage = top_k_coverage(top_k_diff, fuentes_int)
    fuentes_sym_coverage = top_k_coverage(top_k_diff, fuentes_sym_int)
    bonferroni_coverage = top_k_coverage(top_k_diff, bonferroni_int)
    sidak_coverage = top_k_coverage(top_k_diff, sidak_int)
    
    return([SoS_coverage,SoS_sym_coverage,fuentes_coverage,fuentes_sym_coverage, bonferroni_coverage, sidak_coverage])

# %%
eta_values = np.linspace(0,40,21)
rho_values = np.linspace(0,1,11)
mu = np.random.uniform(-1,1, m)
rho = 0.3

df_all = []

# %%
rho_values = [0, 0.3,0.7]
for rho in rho_values: 
    interval_results = []
    for eta in eta_values:
        cov_mat = create_AR_cov_mat(len(mu),rho)
        interval_results.append(interval_coverage(mu*eta,cov_mat,k))
        #print(eta)
    
    d = None
    d = pd.DataFrame(np.array(interval_results))
    d.columns = ["SoS shortest","SoS symmetric","FCW shortest","FCW symmetric", "Bonferroni m", "Sidak m"]
    d['Eta'] = eta_values
    d['Rho'] = np.repeat(rho,d.shape[0])
    d['Correlation'] = np.repeat("AR rho=%.1f"%rho,d.shape[0])
    
    df = pd.melt(d, id_vars=['Eta','Rho', 'Correlation'])
    df_all.append(df)
    #print('---------------')

# %%
cov_mat_model_seven=create_model_7_cov_mat(m)
corr_mat_model_seven = np.corrcoef(cov_mat_model_seven)

interval_results = []
for eta in eta_values:   
    interval_results.append(interval_coverage(mu*eta,corr_mat_model_seven,k))
    #print(eta)


d = None
d = pd.DataFrame(np.array(interval_results))
d.columns = ["SoS shortest","SoS symmetric","FCW shortest","FCW symmetric", "Bonferroni m", "Sidak m"]
d['Eta'] = eta_values
d['Rho'] = np.repeat('NA',d.shape[0])
d['Correlation'] = np.repeat("TD",d.shape[0])

df = pd.melt(d, id_vars=['Eta', 'Correlation','Rho'])
df_all.append(df)

# %%
cov_block_mat=create_block_cov_mat(m)
corr_block_mat = np.corrcoef(cov_block_mat)

interval_results = []
for eta in eta_values:   
    interval_results.append(interval_coverage(mu*eta,corr_block_mat,k))
    #print(eta)

d = None
d = pd.DataFrame(np.array(interval_results))
d.columns = ["SoS shortest","SoS symmetric","FCW shortest","FCW symmetric", "Bonferroni m", "Sidak m"]
d['Eta'] = eta_values
d['Rho'] = np.repeat('NA',d.shape[0])
d['Correlation'] = np.repeat("Block",d.shape[0])

df = pd.melt(d, id_vars=['Eta', 'Correlation','Rho'])
df_all.append(df)
# %%
df_full = pd.concat(df_all, sort=True)
df_full.columns = [u'Correlation', u'Eta', u'Rho', u'value', u'Method']
df_full = df_full.replace('AR rho=0.0', 'Uncorrelated')


# %%
df_sliced = pd.DataFrame(df_full[df_full['Method'].isin(['SoS shortest', 'SoS symmetric'])])

g = sns.FacetGrid(df_sliced, col="Method", hue="Correlation", height=5, aspect=1, col_wrap=2)
g.map(sns.lineplot, "Eta", "value", linewidth=2, alpha=.7)
for axlist in g.axes:
    axlist.set_ylim(.95,1)

g.set_axis_labels(r'$\eta$', "Coverage");
g.add_legend(title='');
g.set_titles("{col_name}")

plt.savefig("figures/Dependency_simulations_SoS.png", dpi=300)
