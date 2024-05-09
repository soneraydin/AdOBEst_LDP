This repository contains the codes for reproducing the numerical results in the paper entitled

**Adaptive Online Bayesian Estimation of Frequency Distributions with Local Differential Privacy**

by Soner Aydın, Sinan Yıldırım

The main file is "experiments_for_SGDL.m" by which one can run the experiments and obtain the accuracy results (in terms of total variation distance) and the numbers $k$ (average cardinality of the subsets chosen by each algorithm). 

Each result in the paper is generated by a different combination of K $\in$ (10, 20), eps1_coeff $\in$ (0.8, 0.9).

After obtaining these results, one can also run the file "plot_results.m" to visualize them as barplots and heatmaps.

Brief information about the other files:
- "AdOBEst_LDP.m" contains the main algorithm.
- "calculate_utility.m" computes the estimation error for 6 different utility functions
- "calculate_FIM.m" computes the Fisher information matrix
- "make_G.m" computes the stochastic matrix G
- "make_p_yx_vec.m" computes the likelihood vector P(Y|X)
- "SGLD_LDP.m" implements the Stochastic Gradient Langevin Dynamics method for approximate MCMC moves.

More information about each of these codes can be found in the comments inside them.
