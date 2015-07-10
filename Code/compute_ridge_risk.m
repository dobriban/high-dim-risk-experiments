function [lambda,estim_risk,pred_risk] = compute_ridge_risk(w,t,gamma,alpha)
%Compute the prediction error of Ridge Regression
%
%Inputs
%w,t - input spectral distribution is a mixture of point masses H = sum
%delta_{t_i} * w_i
%gamma - aspect ratio
%alpha - signal strength
%
%Output
%lambda - grid of regularization parameters
%estim_risk - estimation error of ridge on the lambda grid
%pred_risk - estimation error of ridge on the lambda grid

[lambda,m,v,m_prime,v_prime] = compute_ST(w,t,gamma);

%% compute risk
estim_risk = gamma*m + gamma*lambda.*(lambda*alpha^2/gamma-1).*abs(m_prime);
pred_risk = (1 + (lambda*alpha^2/gamma-1).*(1-lambda.*v_prime./v))./(lambda.*v);


 
