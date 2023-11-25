#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibration code for the paper:
    
A. Cinfrignini, D. Petturiti, B. Vantaggi (2023). 
Market consistent bid-ask option pricing under Dempster-Shafer uncertainty.
"""

import numpy as np
from BinomialTree import plotPutTree

S0 = 300.11
K = 350
T = 20
u = 1.02354592
d = 0.96425885
R = (1 + 0.0555)**(1/250)

b_u = (R - d) / (u - d)

b_d = 0.39260853

print('b_u = ', b_u)
print('b_d = ', b_d)
print('b_u + b_d = ', b_u + b_d)


S = np.zeros((T + 1, T + 1))
S[0, 0] = S0
z = 1
for t in range(1, T+1):
    for j in range(0, z):
        S[j, t] = S[j, t - 1] * u
        S[j + 1, t] = S[j, t - 1] * d
    z += 1
    
print('Stock:\n', S)

PA_bid = np.zeros_like(S)
PA_ask = np.zeros_like(S)
for j in range(0, T + 1):
    PA_bid[j, T] = max(K - S[j, T], 0)
    PA_ask[j, T] = max(K - S[j, T], 0)
    

EE_bid = np.zeros_like(S)
z = T
for t in range(T - 1, -1, -1):
    for j in range(0, z):
        P_hat = (1 / R) * ((1 - b_d) * PA_bid[j, t + 1] + b_d * PA_bid[j + 1, t + 1])
        e_P = max(K - S[j, t], 0)
        PA_bid[j, t] = max(P_hat, e_P)
        if e_P > P_hat:
            EE_bid[j, t] = 1
    z -= 1
    
EE_ask = np.zeros_like(S)
z = T
for t in range(T - 1, -1, -1):
    for j in range(0, z):
        P_hat = (1 / R) * (b_u * PA_ask[j, t + 1] + (1 - b_u) * PA_ask[j + 1, t + 1])
        e_P = max(K - S[j, t], 0)
        PA_ask[j, t] = max(P_hat, e_P)
        if e_P > P_hat:
            EE_ask[j, t] = 1
    z -= 1
    
print('EE_bid:\n', EE_bid)
print('PA_bid:\n', np.round(PA_bid, 2))

print('Bid-ask prices:')
print('PA_bid(0) = ', round(PA_bid[0 ,0],2))
print('PA_ask(0) = ', round(PA_ask[0 ,0],2))


# Tree representations
plotPutTree(EE_bid, PA_bid, False, T, 'Put-tree-bid.png')
plotPutTree(EE_ask, PA_ask, True, T, 'Put-tree-ask.png')



