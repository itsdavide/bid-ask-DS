# bid-ask-DS
Calibration code for the paper:
    
A. Cinfrignini, D. Petturiti, and B. Vantaggi. _Market consistent bid-ask option pricing under Dempster-Shafer uncertainty_. **Quantitative Finance**, DOI: 10.1080/14697688.2024.2353318, 2024. 


# Requirements
The code requires the PySwarm library available at: https://pypi.org/project/pyswarm/.

# File inventory
PSO-opt.py: Executes the 2D calibration on market data.

PSO-opt-free-d.py: Executes the 3D calibration on market data.

BinomialTree.py: Contains a utility function to plot bid-ask binomial trees for American put options.

AmericanPut.py: Computes the price of American options according to the calibration.

AmericanSpreads.py: Computes the theoretical bid-ask spreads as functions of K for American and European put options.

TheoreticalPrices.py: Plot the theoretical bid-ask price curves and the observed market prices for European options.
