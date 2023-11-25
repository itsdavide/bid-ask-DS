#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibration code for the paper:
    
A. Cinfrignini, D. Petturiti, B. Vantaggi (2023). 
Market consistent bid-ask option pricing under Dempster-Shafer uncertainty.
"""

import numpy as np
import matplotlib.pyplot as plt
import math as m

S0 = 300.11
T = 20
u = 1.02354592
d = 0.96425885
R = (1 + 0.0555)**(1/250)

b_u = (R - d) / (u - d)

b_d = 0.39260853

print('b_u = ', b_u)
print('b_d = ', b_d)
print('b_u + b_d = ', b_u + b_d)


def bidaskAmericanPut(u, d, b_u, b_d, S0, K, T, R):
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
    return (PA_bid[0, 0], PA_ask[0, 0])
    
Ks = np.arange(200, 400.25, 0.25)
PA_bids = []
PA_asks = []

for K in Ks:
   bid, ask = bidaskAmericanPut(u, d, b_u, b_d, S0, K, T, R)
   PA_bids.append(bid)
   PA_asks.append(ask)
   
   
plt.figure(figsize=(6, 4))
plt.xlabel(r'$K$')
plt.ylabel('Bid and ask prices')
plt.title(r'Theoretical bid and ask prices of an American put as a function of $K$')

plt.plot(Ks, PA_bids, color='blue')
plt.plot(Ks, PA_asks, color='red')
plt.savefig('Theor-PA.png', dpi=300)

spreads = np.array(PA_asks) - np.array(PA_bids)

plt.figure(figsize=(6, 4))
plt.xlabel(r'$K$')
plt.ylabel('Bid-ask spread')
plt.title(r'Theoretical bid-ask spread of an American put as a function of $K$')

plt.plot(Ks, spreads, color='green')
plt.savefig('Theor-spread-A.png', dpi=300)


def bidaskEuropeanPut(u, d, b_u, b_d, S0, K, T, R):
    PE_bid = (1 / (R**T)) * sum(m.comb(T, k) * (1 - b_d)**k * b_d**(T - k) * max(K - u**k * d**(T - k) * S0, 0) for k in range(0, T + 1))
    PE_ask = (1 / (R**T)) * sum(m.comb(T, k) * b_u**k * (1 - b_u)**(T - k) * max(K - u**k * d**(T - k) * S0, 0) for k in range(0, T + 1))
    return (PE_bid, PE_ask)
    
Ks = np.arange(200, 400.25, 0.25)
PE_bids = []
PE_asks = []

for K in Ks:
   bid, ask = bidaskEuropeanPut(u, d, b_u, b_d, S0, K, T, R)
   PE_bids.append(bid)
   PE_asks.append(ask)
   

plt.figure(figsize=(6, 4))
plt.xlabel(r'$K$')
plt.ylabel('Bid and ask prices')
plt.title(r'Theoretical bid and ask prices of a European put as a function of $K$')

plt.plot(Ks, PE_bids, color='blue')
plt.plot(Ks, PE_asks, color='red')
plt.savefig('Theor-PE.png', dpi=300)

spreads = np.array(PE_asks) - np.array(PE_bids)

plt.figure(figsize=(6, 4))
plt.xlabel(r'$K$')
plt.ylabel('Bid-ask spread')
plt.title(r'Theoretical bid-ask spread of a European put as a function of $K$')

plt.plot(Ks, spreads, color='green')
plt.savefig('Theor-spread-E.png', dpi=300)

