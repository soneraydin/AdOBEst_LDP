function [phis] = SGLD_DP(phi, P_yx, rho0, M, S, a)

% [phis] = SGLD_DP(phi, P_yx, rho0, M, S, a)
% 
% This function implements the Stochastic Gradient Langevin Dynamics moves 
% to sample thetas given P_yx.
% phi: surrogate variable for theta
% P_yx: collection of likelihood vectors, P(y|x)
% rho0: prior hyperparameter vector for Gamma distribution
% M: number of MCMC runs
% S: subsample size
% a: step length for gradient update


[K, T] = size(P_yx);
phis = zeros(K, M);
sub_size = min(T, S);

for m = 1:M
    sub_ind = randi([1, T], 1, sub_size);
	P_yx_s = P_yx(:, sub_ind);

	theta = phi/sum(phi);

    J = eye(K, K-1)/sum(phi) - phi(1:K-1)'/sum(phi)^2;
    dfdtheta = (P_yx_s(1:K-1, :) - P_yx_s(K, :))./(theta'*P_yx_s);
    
    grad_phi = mean(J*dfdtheta, 2)*T + (rho0'-1)./phi - 1;

    phi = phi + (a/2)*grad_phi + randn(K, 1)*sqrt(a);
    phi = abs(phi);

    phis(:, m) = phi;
end