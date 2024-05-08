This repository contains the codes for reproducing the numerical results in "Adaptive Online Bayesian Estimation of Frequency Distributions
with Local Differential Privacy" paper.

The main file is "experiments_for_SGDL.m" by which one can run the experiments and obtain the accuracy results (in terms of total variation distance) and the numbers k (average cardinality of the subsets chosen by each algorithm). Here we advise the user to set the values for K {10, 20} and eps1_coeff {0.8, 0.9} manually, since these values were used in the paper.
After obtaining these results, one can also run the file "plot_results.m" to visualize them as barplots and heatmaps.
Brief information about the other files:
- "adaptive_density_est_DP_SGLD.m" contains the main algorithm.
- "calculate_est_error.m" computes the estimation error for 6 different loss functions
- "calculate_FIM.m" computes the Fisher information matrix
- "make_G.m" computes the stochastic matrix G
- "make_p_yx_vec.m" computes the likelihood vector P(Y|X)
- "SGLD_DP.m" implements the Stochastic Gradient Langevin Dynamics method for approximate MCMC moves
More information about each of these codes can be found in the comments inside them.
