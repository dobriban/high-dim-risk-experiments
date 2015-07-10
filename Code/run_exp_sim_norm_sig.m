function [lambda,risk,lambda_th,risk_th] = run_exp_sim_norm_sig(exp_lambda,gamma,oracle_error_rate,n_lambda,n,ub_only)
%run an EXP simulation
%
% Inputs:
% exp_lambda - rate parameter of exponential distribution
% gamma - aspect ration p/n
% percent_error - desired mis-classification error
% n  - (optional) sample size of matrix used in simulation
% ub_only - (optional) parameter determining how lambda=grids should be
% truncated (default = 1), only upper grids

% Outputs
% lambda, risk - a grid of regularization parameters, and empirical risk
%       estimates
% lambda_th, risk_th - a larger grid of regularization parameters, and
%       theoretical risk calculations
%% Set defaults
if ~exist('n','var')
    n = 1e3;
end
if ~exist('n_lambda','var')
  n_lambda = 1e2;
end
if ~exist('ub_only','var')
    ub_only = 1;
end

%% calculate EXP parameters
p = floor(n*gamma);
g = linspace(1/(2*p),1-1/(2*p),p);
g = 1/exp_lambda*log(1./g);
Sigma = diag(g);

%% grids, calibration of alpha level
Maxi = 10;
power = 2;
lambda = linspace(0,Maxi^1/power,n_lambda).^power';
c = sqrt(p/ trace(Sigma^(-1)));
alpha = c * abs(norminv(oracle_error_rate));

%% simulation
risk = simulate_rda_risk_cov(Sigma,n,p,alpha, lambda);
%% theoretical risk formula
t = eig(Sigma);
w = ones(p,1)/p;
[lambda_th,risk_th] = compute_rda_risk(w,t,gamma,alpha);

%% truncate the lambda-s to a common interval common interval

if ub_only
    %% common interval - upper bound only
    maxi = min(max(lambda), max(lambda_th));
    
    ind1 = ( lambda <=maxi);
    lambda = lambda(ind1);
    risk = risk(ind1);
    
    ind2 = ( lambda_th <=maxi) ;
    lambda_th = lambda_th(ind2);
    risk_th = risk_th(ind2);
else
    maxi = min(max(lambda), max(lambda_th));
    mini = max(min(lambda), min(lambda_th));
    
    ind1 = ( lambda <=maxi) & (lambda >=mini);
    lambda = lambda(ind1);
    risk = risk(ind1);
    
    ind2 = ( lambda_th <=maxi) & (lambda_th >=mini);
    lambda_th = lambda_th(ind2);
    risk_th = risk_th(ind2);
end

