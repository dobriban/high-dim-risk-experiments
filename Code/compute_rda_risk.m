function [lambda,risk] = compute_rda_risk(w,t,gamma,alpha)
%Compute the prediction error of Regularized Discriminant Analysis
%
%Inputs
%w,t - input spectral distribution is a mixture of point masses H = sum
%delta_{t_i} * w_i
%gamma - aspect ratio
%alpha - signal strength
%
%Output
%lambda - grid of regularization parameters
%risk - prediction error of RDA on the lambda grid

[lambda,m,v,~,v_prime] = compute_ST(w,t,gamma);

%% compute risk
Theta = alpha^2*m.*v.*lambda./sqrt(alpha^2*(v-lambda.*v_prime)/gamma+(v_prime-v.^2)./v.^2);
risk = normcdf(-Theta);

 
