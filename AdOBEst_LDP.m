function [theta_est, thetas, Y, k_selected] = AdOBEst_LDP(X, eps_DP, eps1_coeff, rho0, M, alpha, loss_type, S, a)

% [theta_est, thetas, Y, k_selected] = AdOBEst_LDP(X, 
% eps_DP, eps1_coeff, rho0, M, alpha, loss_type, S, a)
%
% This function implements the main algorithm for adaptive, non-adaptive,
% and semi-adaptive approaches.
% X: real data
% eps_DP: DP parameter
% eps1_coeff: coefficient that relates eps1 to eps_DP 
% (eps1 = eps1_coeff*eps_DP)
% rho0: prior hyperparameter vector for Gamma distribution
% M: number of MCMC runs
% alpha: parameter for the semi-adaptive approach
% loss_type: a value between 1 and 6 to choose one of loss functions below:
% 1:FIM, 2:Entropy, 3:TV1, 4:TV2, 5:Expected MSE, 6. Prob(Y = x)
% S: subsample size
% a: step length for gradient update


T = length(X);
K = length(rho0);
Y = zeros(1, T);

% determine eps1 and eps2
if loss_type == 0
    eps1 = eps_DP;
    eps2_vec = zeros(1, K);
else
    eps1 = eps1_coeff*eps_DP;
    Sc_min = ceil(exp(eps_DP - eps1));
    k_vec = Sc_min:K;
    eps2_UB = eps_DP*ones(1, K);
    eps2_UB(k_vec) = (k_vec - 1)./(exp(eps1-eps_DP).*k_vec - 1);    
    eps2_vec = min(eps_DP, eps2_UB);
end

GG = cell(1, K);
for k = 1:K
    GG{k} = make_G(K, k, eps1, eps2_vec(k));
end

M_final = M*100;
phi_samp = ones(K, 1)/K;
theta_samp = phi_samp/sum(phi_samp);

P_yx = zeros(K, T);
ord_ind_inv = zeros(1, K);

k_selected = zeros(1, T);

for t = 1:T
    % Get x
    x = X(t);
    
    % A: Construct star set
    [theta_ord, ord_ind] = sort(theta_samp, 'descend');
    ord_ind_inv(ord_ind) = 1:K;

    if loss_type == 0 % semi-adaptive version (if alpha = 0, it implements the non-adaptive version)
        k_best = min(K, sum((cumsum(theta_ord)) < alpha) + 1); % for precision reasons
    else
        L = zeros(1, K);
        for k = 1:K
            L(k) = calculate_est_error(theta_ord, GG{k}, loss_type);
        end
        [~, k_best] = min(L);
    end
    star_set = ord_ind(1:k_best);
    k_selected(t) = k_best;
    
    % B: Sample y
    x_ord = ord_ind_inv(x);
    y_ord = find(rand <= cumsum(GG{k_best}(:, x_ord)), 1);
    Y(t) = ord_ind(y_ord);

    % C. Sample theta
    % The likelihood vector
    P_yx(:, t) = make_p_yx_vec(Y(t), K, star_set, k_best, eps1, eps2_vec(k_best)); 

    M_current = (t < T)*M + (t == T)*M_final;
    [phis] = SGLD_DP(phi_samp, P_yx(:, 1:t), rho0, M_current, S, a);
    phi_samp = phis(:, end);
    theta_samp = phi_samp/sum(phi_samp);
    
end

m_burn = M_final/2;
thetas = phis./sum(phis, 1);
theta_est = mean(thetas(:, m_burn+1:end), 2);
