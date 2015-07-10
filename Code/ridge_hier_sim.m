function [lambda,estim_risk,pred_risk,lambda_th,estim_risk_th,pred_risk_th] = ridge_hier_sim(gamma,alpha,depth,n_lambda)
%run an BinaryTree simulation with parameters gamma, alpha, depth
%
% Inputs:
% depth - depth of tree
% gamma - aspect ration p/n
% alpha - signal strength
% n  - (optional) sample size of matrix used in simulation
% n_lambda - (optional) parameter determining number of lambdas

% Outputs
% lambda, estim_risk,pred_risk - a grid of regularization parameters, and empirical risk
%       estimates
% lambda_th, estim_risk_th,pred_risk_th - a larger grid of regularization parameters, and
%       theoretical risk calculations
%% calculate AR(1) parameters
if ~exist('n_lambda','var')
    n_lambda = 1e2;
end 
if ~exist('ub_only','var')
    ub_only = 1;
end

p = 2^depth - 1;
n = floor(p / gamma);

eigs_arr = {};
for idx = 1:depth
    eigs_arr{idx} = repmat(2^idx, [2^(depth - idx) 1]);
end
eigs = cat(1, eigs_arr{:});
eigs = eigs ./ (sum(eigs) / p);
Sigma = diag(eigs);

lambda = linspace(0.1,4,n_lambda).^2';

%% simulation
[estim_risk,pred_risk] = simulate_ridge_risk_cov(Sigma,n,p,alpha,lambda);
%% theoretical risk formula
t = eig(Sigma);
w = ones(p,1)/p;
[lambda_th,estim_risk_th,pred_risk_th] = compute_ridge_risk(w,t,gamma,alpha);
%% truncate the lambda-s to a common interval common interval

if ub_only
    %% common interval - upper bound only
    maxi = min(max(lambda), max(lambda_th));
    
    ind1 = ( lambda <=maxi);
    lambda = lambda(ind1);
    estim_risk = estim_risk(ind1);
    pred_risk = pred_risk(ind1);
    
    ind2 = ( lambda_th <=maxi) ;
    lambda_th = lambda_th(ind2);
    estim_risk_th = estim_risk_th(ind2);
    pred_risk_th = pred_risk_th(ind2);
else
    maxi = min(max(lambda), max(lambda_th));
    mini = max(min(lambda), min(lambda_th));
    
    ind1 = ( lambda <=maxi) & (lambda >=mini);
    lambda = lambda(ind1);
    estim_risk = estim_risk(ind1);
    pred_risk = pred_risk(ind1);
     
    ind2 = ( lambda_th <=maxi) & (lambda_th >=mini);
    lambda_th = lambda_th(ind2);
    estim_risk_th = estim_risk_th(ind2);
    pred_risk_th = pred_risk_th(ind2);
end

