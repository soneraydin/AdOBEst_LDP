% The main file to conduct the numerical experiments
clear; clc; close all;

rng(1); % Setting the seed to obtain reproducible results

% DP parameters
eps_DP_vec = [0.5, 1, 5]; L_e = length(eps_DP_vec);

% Loss types (the first five values are zero, since they 
% correspond to semi-adaptive approach that doesn't use any loss function)
loss_type_vec = [0 0 0 0 0 0 1 2 3 4 5 6]; L_m = length(loss_type_vec);

% Alpha parameters for semi-adaptive approach
% The last five values are zero, since they correspond to adaptive 
% approach which doesn't use the alpha parameter
alpha_vec = [1 0.95 0.9 0.8 0.6 0.2 0 0 0 0 0 0];

K_vec = 20; L_K = length(K_vec); % Number of categories

M = 20; % Number of MCMC runs

% Dirichlet parameters
rho_coeff_vec = [0.01 0.1 1]; L_r = length(rho_coeff_vec); 
MC_run = 50; % Number of Monte Carlo runs

S = 50; % Number of SGLD steps 

% Initialization of cells for TV results 
% and cardinalities of the selected subsets
TV = repmat({zeros(L_e, L_m, MC_run)}, L_K, L_r);
K_selected = cell(L_K, L_r, L_e, L_m, MC_run);

% Coefficient that relates eps1 to eps_DP
% eps1 = eps1_coeff*eps_DP
eps1_coeff = 0.80; 

% Name of the file for saving the results
datanametosave = ['K_vec_' sprintf('%d_', K_vec) sprintf('eps1coeff_%02d_', 100*eps1_coeff) ...
    'methods' sprintf('_%d', loss_type_vec)]; 

for mc = 1:MC_run
    for i1 = 1:L_K
        K = K_vec(i1);
        T = 500*K; 
        a = 0.5/T;
        rho0 = ones(1, K);
        for i2 = 1:L_r
            rho_coeff = rho_coeff_vec(i2);

            rho_x = rho_coeff*ones(1, K);
            theta_true = gamrnd(rho_x, 1);
            theta_true = theta_true/sum(theta_true);
            X = randsample(1:K, T, 'true', theta_true);

            for i3 = 1:L_e
                eps_DP = eps_DP_vec(i3);
                for i4 = 1:L_m
                    loss_type = loss_type_vec(i4);
                    alpha = alpha_vec(i4);
                    fprintf('Monte Carlo run no %d for \n epsilon = %.2f, rho = %.2f, K=%d, method=%d... \n', mc, eps_DP, rho_coeff, K, loss_type);
                    [theta_est, thetas, Y, k_selected] = adaptive_density_est_DP_SGLD(X, eps_DP, eps1_coeff, rho0, M, alpha, loss_type, S, a);
                    TV_current = 0.5*sum(abs(theta_est' - theta_true));
                    TV{i1, i2}(i3, i4, mc) = TV_current;
                    K_selected{i1, i2}(i3, i4, mc) = mean(k_selected);
                    disp(TV_current);
                end
            end
            save(datanametosave);
        end
    end
end


%% results
% Computing the means and standard deviations of the TV results
TV_mean = cell(L_K, L_r); 
TV_stdev = cell(L_K, L_r);
for i1 = 1:L_K
    for i2 = 1:L_r
        TV_mean{i1, i2} = mean(TV{i1, i2}, 3);
        TV_stdev{i1, i2} = std(TV{i1, i2}, [], 3);
    end
end

save(datanametosave);
