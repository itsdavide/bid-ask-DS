#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibration code for the paper:
    
A. Cinfrignini, D. Petturiti, B. Vantaggi (2023). 
Market consistent bid-ask option pricing under Dempster-Shafer uncertainty.
"""

"""
EXPLANATION OF THE CODE:
This code plots the theoretical bid-ask price curves and the observed market 
prices of European call and put options on META stock with maturity T = 20 
trading days as a function of K.
"""

from pyswarm import pso
import pandas as pd
import math as m
import numpy as np
import matplotlib.pyplot as plt

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


# Error function with the tree parameters u and b_d (assuming d = 1/ u)
def E(x):
    u = x[0]
    b_d = x[1]
    
    d = 1 / u
    
    b_u = (R - d) / (u - d)
    
    
    E = 0
    
    for i in I_calls:
        E += ((1 / (R**T)) * sum(m.comb(T, k) * b_u**k * (1 - b_u)**(T - k) * max(u**k * d**(T - k) * S_0_bid - K_C_T[i], 0) for k in I_S_T_vals) - C_0_bid[i])**2
        E += ((1 / (R**T)) * sum(m.comb(T, k) * (1 - b_d)**k * b_d**(T - k) * max(u**k * d**(T - k) * S_0_bid - K_C_T[i], 0) for k in I_S_T_vals) - C_0_ask[i])**2
    
    for i in I_puts:
        E += ((1 / (R**T)) * sum(m.comb(T, k) * (1 - b_d)**k * b_d**(T - k) * max(K_P_T[i] - u**k * d**(T - k) * S_0_bid, 0) for k in I_S_T_vals) - P_0_bid[i])**2
        E += ((1 / (R**T)) * sum(m.comb(T, k) * b_u**k * (1 - b_u)**(T - k) * max(K_P_T[i] - u**k * d**(T - k) * S_0_bid, 0) for k in I_S_T_vals) - P_0_ask[i])**2
    
    return E

# Border of the regionf of the admissible values of b_d
def phi(u):
    return 1 - (R - (1 / u)) / (u - (1 / u))

# Constraint function to assure b_d in (0, 1 - b_u]
def con(x):
    u = x[0]
    b_d = x[1]
    return [1 - (R - (1 / u)) / (u - (1 / u)) - b_d]

lb = [R + 0.00005, 0.3]
ub = [1.056, 0.6]


# OPTIMIZATION WITH PSO

xopt, fopt = pso(E, lb, ub, f_ieqcons=con, maxiter=2000)

print('[u, b_d] = ', xopt)
print('Minimum squared error:', fopt)

(u, b_d) = xopt
d = 1 / u
b_u = (R - d) / (u - d)
print('b_u = ', b_u)
print('b_d = ', b_d)
print('b_u + b_d = ', b_u + b_d)




############################################################
#
# Graph of the error surface and of the contour lines
#
############################################################


def F(u, b_d):
    return E([u, b_d])

# Vectorize the function
vE = np.vectorize(F)
    

# Make data
x = np.linspace(1.022, 1.036, 100)
y = np.linspace(0.47, 0.53, 100)


X, Y = np.meshgrid(x, y)
zs = np.array(vE(np.ravel(X), np.ravel(Y)))
Z = zs.reshape(X.shape)


levels = np.arange(0, 500, 10)


# Fontsize of the tick labels
plt.rc('xtick', labelsize=7)    
plt.rc('ytick', labelsize=7)

# Surface plot
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(X, Y, Z, rstride=10, cstride=10, cmap='viridis', edgecolors='k', lw=0.6)

plt.title('Graph of $E$')
ax.set_xlabel('$u$')
ax.set_ylabel('$\widehat{b_d}$', rotation=0)

ax.zaxis.set_rotate_label(False)
ax.set_zlabel('$E(u,\widehat{b_d})$', rotation=90)

ax.view_init(30, 250)
plt.savefig('Error.png', dpi=300)


# Contour lines plot
plt.clf()
plt.rc('xtick', labelsize=8)    
plt.rc('ytick', labelsize=8)
fig = plt.figure(figsize=(5,5))
plt.title('Contour lines of $E$')
plt.xlabel('$u$')
plt.ylabel('$\widehat{b_d}$')
plt.contour(X, Y, Z, levels)

y = phi(x)
plt.fill_between(x, y, 0.53, color='red', alpha=0.8)

x_s = [u]
y_s = [b_d]
plt.plot(x_s, y_s, marker="o", markersize=3, markeredgecolor='blue', markerfacecolor='blue')

plt.savefig('Error-cl.png', dpi=300)
