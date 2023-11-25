#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibration code for the paper:
    
A. Cinfrignini, D. Petturiti, B. Vantaggi (2023). 
Market consistent bid-ask option pricing under Dempster-Shafer uncertainty.
"""

from matplotlib import pyplot as plt
import numpy as np

def plotPutTree(EE, PA, ask, T, filename):
    
    xlabels = []
    for i in range(0,T + 1):
        if ask:
            xlabels.append('$\overline{P}_{' + str(i) + '}$' )
        else:
            xlabels.append('$P_{' + str(i) + '}$' )


    fig = plt.figure(figsize=[14, 12])
    ax = fig.gca()
    ax.tick_params(bottom=False, top=True, left=False, right=False)
    ax.tick_params(labelbottom=False, labeltop=True, labelleft=False, labelright=False)
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xticklabels(xlabels, fontsize=18)

    for i in range(T):
        print('Checking i:', i)
        x = [1, 0, 1]
        row = 0
        for j in range(i):
            print('Checking j:', j)
            x.append(0)
            x.append(1)
        print('Orig. x:', x)
        x_orig = x
        x = np.array(x) + i
        y = np.arange(-(i+1), i+2)[::-1]
        print('x:', x)
        print('y:', y)
        plt.plot(x, y, 'o-', color='lightskyblue')
        # Draws EE
        x_EE = []
        y_EE = []
        row_EE = 0
        print(EE)
        for k in range(len(x)):
            if x_orig[k] == 0:
                if EE[row_EE, i] == 1:
                    x_EE.append(x[k - 1])
                    y_EE.append(y[k - 1])
                    x_EE.append(x[k])
                    y_EE.append(y[k])
                    x_EE.append(x[k + 1])
                    y_EE.append(y[k + 1])
                row_EE += 1
        plt.plot(x_EE, y_EE, 'o--', color='orange')       
        strings = []
        # Plot the strings
        if i < T - 1:
            for k in range(len(x)):
                if x_orig[k] == 0:
                    strings.append(str(round(PA[row, i], 2)))
                    row += 1
                else:
                    strings.append('')
            for k in range(len(x)):
                ax.annotate(strings[k], (x[k], y[k]), fontsize=12)
        else:
            strings = ['' for k in range(len(x))]
            x_last = []
            y_last = []
            for k in range(len(x)):
                if x_orig[k] == 0:
                    x_last.append(x[k-1])
                    y_last.append(y[k-1])
                    x_last.append(x[k+1])
                    y_last.append(y[k+1])
                    strings[k - 1] = str(round(PA[row, i + 1], 2))
                    strings[k] = str(round(PA[row, i], 2))
                    strings[k + 1] = str(round(PA[row + 1, i + 1], 2))
                    row += 1
            plt.plot(x_last, y_last, 'o', color='lightskyblue') 
            for k in range(len(x)):
                ax.annotate(strings[k], (x[k], y[k]), fontsize=12)
            
    
    plt.xticks(list(range(0,T+1)), xlabels)
    
    plt.savefig(filename, dpi=300)