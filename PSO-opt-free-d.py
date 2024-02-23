#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibration code for the paper:
    
A. Cinfrignini, D. Petturiti, B. Vantaggi (2023). 
Market consistent bid-ask option pricing under Dempster-Shafer uncertainty.
"""

"""
EXPLANATION OF THE CODE:
Market calibration of parameters u, d and b_d relying on the PSO technique.
"""

from pyswarm import pso
import pandas as pd
import math as m
import numpy as np

np.random.seed(12345)

# Maturity in trading days
T = 20
  
# Risk-free return
R = (1 + 0.0555)**(1/250)

# Names of the dataset files
file_calls = 'META_init_date_2023-09-29/META_calls_2023-10-27_s_1_0_LOWER.csv'
file_puts = 'META_init_date_2023-09-29/META_puts_2023-10-27_s_1_0_LOWER.csv'
    
# Load STOCK calls
STOCK_calls = pd.read_csv('./datasets/' + file_calls)[['strike','bid','ask']]
    
# Insert the stock as a degenerate call with strike 0
S_0_bid = 300.11
S_0_ask = 300.30
STOCK_calls.loc[len(STOCK_calls)] = [0, S_0_bid, S_0_ask]
    
# Load STOCK puts
STOCK_puts = pd.read_csv('./datasets/' + file_puts)[['strike','bid','ask']]
    
# Get indices of S_T values
I_S_T_vals = list(range(T+1))
    
    
# Get indices of calls and puts
I_calls = list(range(len(STOCK_calls)))

print('Call and put options dataset:')

print('# calls:', len(I_calls))

I_puts = list(range(len(STOCK_puts)))

print('# puts:', len(I_puts))
    
C_0_bid = STOCK_calls['bid']
C_0_ask = STOCK_calls['ask']
K_C_T = STOCK_calls['strike']
P_0_bid = STOCK_puts['bid']
P_0_ask = STOCK_puts['ask']
K_P_T = STOCK_puts['strike']



# Error function with the tree parameters u, d and b_d
def E(x):
    u = x[0]
    d = x[1]
    b_d = x[2]
    
    b_u = (R - d) / (u - d)
    
    E = 0
    
    for i in I_calls:
        E += ((1 / (R**T)) * sum(m.comb(T, k) * b_u**k * (1 - b_u)**(T - k) * max(u**k * d**(T - k) * S_0_bid - K_C_T[i], 0) for k in I_S_T_vals) - C_0_bid[i])**2
        E += ((1 / (R**T)) * sum(m.comb(T, k) * (1 - b_d)**k * b_d**(T - k) * max(u**k * d**(T - k) * S_0_bid - K_C_T[i], 0) for k in I_S_T_vals) - C_0_ask[i])**2
        
    for i in I_puts:
        E += ((1 / (R**T)) * sum(m.comb(T, k) * (1 - b_d)**k * b_d**(T - k) * max(K_P_T[i] - u**k * d**(T - k) * S_0_bid, 0) for k in I_S_T_vals) - P_0_bid[i])**2
        E += ((1 / (R**T)) * sum(m.comb(T, k) * b_u**k * (1 - b_u)**(T - k) * max(K_P_T[i] - u**k * d**(T - k) * S_0_bid, 0) for k in I_S_T_vals) - P_0_ask[i])**2

    return E


# Constraint function to assure b_d in (0, 1 - b_u]
def con(x):
    u = x[0]
    d = x[1]
    b_d = x[2]
    return [1 - (R - d) / (u - d) - b_d]

# Bounds for the oprimization variables u, d and b_d
lb = [R + 0.00001, 0, 0.3]
ub = [1.056, R - 0.00001, 0.6]


# OPTIMIZATION WITH PSO

xopt, fopt = pso(E, lb, ub, f_ieqcons=con, maxiter=2000)


print('[u, d, b_d] = ', xopt)
print('Minimum squared error:', fopt)

(u, d, b_d) = xopt
b_u = (R - d) / (u - d)
print('b_u = ', b_u)
print('b_d = ', b_d)
print('b_u + b_d = ', b_u + b_d)
