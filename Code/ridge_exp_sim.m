function [lambda,estim_risk,pred_risk,lambda_th,estim_risk_th,pred_risk_th] = ridge_exp_sim(gamma,alpha,rate,n_lambda,n)
%run an Exponential simulation with parameters gamma, alpha, depth
%
% Inputs:
% rate - rate of exponential
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
if ~exist('n','var')
    n = 1e3;
end
if ~exist('n_lambda','var')
    n_lambda = 1e2;
end 
if ~exist('ub_only','var')
    ub_only = 1;
end

p = floor(n*gamma);
g = linspace(1/(2*p),1-1/(2*p),p);
g = 1/rate*log(1./g);
Sigma = diag(g);
 
lambda = linspace(0.1,2.5,n_lambda).^2';

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

