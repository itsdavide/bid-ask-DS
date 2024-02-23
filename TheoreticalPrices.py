#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibration code for the paper:
    
A. Cinfrignini, D. Petturiti, B. Vantaggi (2023). 
Market consistent bid-ask option pricing under Dempster-Shafer uncertainty.
"""

"""
EXPLANATION OF THE CODE:
Market calibration of parameters u and b_d (assuming that d = 1 / u) relying on
the PSO technique.
"""

import math as m
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file_calls = 'META_init_date_2023-09-29/META_calls_2023-10-27_s_1_0_LOWER.csv'
file_puts = 'META_init_date_2023-09-29/META_puts_2023-10-27_s_1_0_LOWER.csv'
    
# Load STOCK calls
STOCK_calls = pd.read_csv('./datasets/' + file_calls)[['strike','bid','ask']]
# Load STOCK puts
STOCK_puts = pd.read_csv('./datasets/' + file_puts)[['strike','bid','ask']]
    
print('K calls: [', min(STOCK_calls['strike']), ',', max(STOCK_calls['strike']), ']')
print('K puts: [', min(STOCK_puts['strike']), ',', max(STOCK_puts['strike']), ']')



T = 20
  
R = (1 + 0.0555)**(1/250)


S_0_bid = 300.11
S_0_ask = 300.30


u = 1.02354592
d = 0.96425885

b_u = (R - d) / (u - d)

b_d = 0.39260853

S_0_bid_teor = (1 / (R**T)) * sum(m.comb(T, k) * b_u**k * (1 - b_u)**(T - k) * max(u**k * d**(T - k) * S_0_bid - 0, 0) for k in range(T + 1))
S_0_ask_teor = (1 / (R**T)) * sum(m.comb(T, k) * (1 - b_d)**k * b_d**(T - k) * max(u**k * d**(T - k) * S_0_bid - 0, 0) for k in range(T + 1))

print('S_0_bid_teor =', S_0_bid_teor)
print('S_0_ask_teor =', S_0_ask_teor)

print('b_u:', round(b_u,4))

print('b_d:', round(b_d,4))

print('b_u + b_d:', b_u + b_d)


Ks = np.arange(min(STOCK_puts['strike']), max(STOCK_puts['strike']), 0.5)

P_bids = []
P_asks = []

for K in Ks:
    P_0_bid_teor = (1 / (R**T)) * sum(m.comb(T, k) * (1 - b_d)**k * b_d**(T - k) * max(K - u**k * d**(T - k) * S_0_bid, 0) for k in range(T + 1))
    P_0_ask_teor = (1 / (R**T)) * sum(m.comb(T, k) * b_u**k * (1 - b_u)**(T - k) * max(K - u**k * d**(T - k) * S_0_bid, 0) for k in range(T + 1))
    P_bids.append(P_0_bid_teor)
    P_asks.append(P_0_ask_teor)


fig = plt.figure(figsize=(6,4))
plt.title('Theoretical bid prices European puts')
plt.xlabel('$K$')
plt.ylabel('Bid prices') 
plt.plot(Ks, P_bids, linewidth=2, color='blue')
plt.scatter(STOCK_puts['strike'], STOCK_puts['bid'], c='green',s=10)
plt.savefig('ThPrice_bid_puts.png', dpi=300)

fig = plt.figure(figsize=(6,4))
plt.title('Theoretical ask prices European puts')
plt.xlabel('$K$')
plt.ylabel('Ask prices') 
plt.plot(Ks, P_asks, linewidth=2, color='red')
plt.scatter(STOCK_puts['strike'], STOCK_puts['ask'], c='orange',s=10)
plt.savefig('ThPrice_ask_puts.png', dpi=300)



Ks = np.arange(min(STOCK_calls['strike']), max(STOCK_calls['strike']), 0.5)

C_bids = []
C_asks = []

for K in Ks:
    C_0_bid_teor = (1 / (R**T)) * sum(m.comb(T, k) * b_u**k * (1 - b_u)**(T - k) * max(u**k * d**(T - k) * S_0_bid - K, 0) for k in range(T + 1))
    C_0_ask_teor = (1 / (R**T)) * sum(m.comb(T, k) * (1 - b_d)**k * b_d**(T - k) * max(u**k * d**(T - k) * S_0_bid - K, 0) for k in range(T + 1))
    C_bids.append(C_0_bid_teor)
    C_asks.append(C_0_ask_teor)

fig = plt.figure(figsize=(6,4))
plt.title('Theoretical bid prices European calls')
plt.xlabel('$K$')
plt.ylabel('Bid prices') 
plt.plot(Ks, C_bids, linewidth=2, color='blue')
plt.scatter(STOCK_calls['strike'], STOCK_calls['bid'], c='green',s=10)
plt.savefig('ThPrice_bid_calls.png', dpi=300)

fig = plt.figure(figsize=(6,4))
plt.title('Theoretical ask prices European calls')
plt.xlabel('$K$')
plt.ylabel('Ask prices') 
plt.plot(Ks, C_asks, linewidth=2, color='red')
plt.scatter(STOCK_calls['strike'], STOCK_calls['ask'], c='orange',s=10)
plt.savefig('ThPrice_ask_calls.png', dpi=300)













