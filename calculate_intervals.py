import numpy as np
import pandas as pd
from scipy.stats import norm

def bonferroni_interval(m, alpha = 0.05):
	#Calculate the bonferroni m interval

    c_under_bar = norm.ppf(1-alpha/(2*m))
    c_upper_bar = -1*norm.ppf(alpha/(2*m))
    
    interval = [c_under_bar,c_upper_bar]
    interval_length = c_under_bar + c_upper_bar
    
    return [interval_length, interval]

def sidak_interval(m, alpha = 0.05):
	#Calculate the sidak m interval

    c_sidak = norm.ppf(1-(1-(1-alpha)**(1/float(m)))/2)
    c_upper_bar = c_under_bar = c_sidak
    
    interval = [c_under_bar,c_upper_bar]
    interval_length = c_under_bar + c_upper_bar
    
    return [interval_length, interval]

def fuentes_coverage(c,d,k,m):
	#Calculate the FCW coverage for a given c and d for the selection of k out of m

    return ((norm.cdf(c) - norm.cdf(-d))**(k-1))*(norm.cdf(c)**(m-k+1)-norm.cdf(-d)**(m-k+1))

def fuentes_intervals(k,m,alpha = 0.05, x_grid_values = 3001, y_grid_values = 3001):
	#Find the shortest FCW interval for the selection of k out of m by grid search

    xaxis = np.linspace(1.88, 5, x_grid_values)
    yaxis = np.linspace(1.88, 5, y_grid_values)
    result = fuentes_coverage(xaxis[:,None], yaxis[None,:],k,m)
    
    x_ix, y_ix = np.where(result > (1-alpha))
    
    min_ix = (xaxis[x_ix] + yaxis[y_ix]).argmin()
    x_ix_min = xaxis[x_ix][min_ix]
    y_ix_min = yaxis[y_ix][min_ix]
    
    interval_length = y_ix_min + x_ix_min
    
    return([interval_length,[x_ix_min,y_ix_min]])

def fuentes_symmetric_interval(k,m,alpha = 0.05, x_grid_values = 7001):
	#Find the symmetric FCW interval 
    c_values = np.linspace(1, 5, x_grid_values)
    result = fuentes_coverage(c_values,c_values,k,m)
    
    x_ix = np.where(result > (1-alpha))[0]
    c_min = c_values[x_ix.min()]
    result[x_ix.min()]    

    interval_length = c_min + c_min
    
    return([interval_length,[c_min,c_min]])

def SoS_intervals(k,m,delta, alpha = 0.05):
	#Find the SoS interval for a given delta for the selection of k out of m 

    c_under_bar = norm.ppf(1-(delta*alpha)/m)
    c_upper_bar = -1*norm.ppf((1-delta)*alpha/k)
    
    interval = [c_under_bar,c_upper_bar]
    interval_length = c_under_bar + c_upper_bar
    
    return [interval_length, interval]

def SoS_shortest_interval(k,m, alpha=0.05):
	#Find the shortest SoS interval for the selection of k out of m by grid search

    delta = np.linspace(0.001,.999,5001)
    SoS_results = SoS_intervals(k,m,delta, alpha=alpha)

    shortest_interval_ix = np.argmin(SoS_results[0])
    interval_length = SoS_results[0][shortest_interval_ix]
    shortest_c_under_bar = SoS_results[1][0][shortest_interval_ix]
    shortest_c_upper_bar = SoS_results[1][1][shortest_interval_ix]

    return([interval_length, [shortest_c_under_bar,shortest_c_upper_bar]])

def SoS_all_interval(k,m, alpha=0.05):
	#Return all possible SoS intervals for the selection of k out of m 

    delta = np.linspace(0.001,.999,5001)
    SoS_results = SoS_intervals(k,m,delta, alpha=alpha)

    interval_lengths = SoS_results[0]
    c_under_bar_vec = SoS_results[1][0]
    c_upper_bar_vec = SoS_results[1][1]

    return([interval_lengths, [c_under_bar_vec,c_upper_bar_vec], delta])
