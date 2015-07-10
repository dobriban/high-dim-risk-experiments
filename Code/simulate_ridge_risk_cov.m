function [estim_risk,pred_risk] = simulate_ridge_risk_cov(Sigma,n,p,alpha,lambda_arr)
%Simulate the Risk of Ridge Regression

%% Monte Carlo Evaluation
num_monte_in = 1e2;
num_monte = 5*1e2;
rng(0)
pred_risk = zeros(length(lambda_arr),1);
estim_risk = zeros(length(lambda_arr),1);
S = Sigma^(1/2);
tic
for k=1:length(lambda_arr)
    timer = toc;
    fprintf('Lambda: %d/%d; Time: %f\n',k,length(lambda_arr),timer);
    lambda = lambda_arr(k);
    pred_err = zeros(num_monte,num_monte_in);
    estim_err = zeros(num_monte,1);
    for i=1:num_monte
        X = randn(n,p)*S;
        w = alpha*randn(p,1)/sqrt(p);
        epsi = randn(n,1);
        y = X*w + epsi;
        
        ip = X'* y;
        w_hat = (X'*X +  n*lambda*eye(p)) \ ip;
        estim_err(i) = sum((w-w_hat).^2);
        %inner loop, generate random test data: x_0, y_0
        for j=1:num_monte_in
            x_0 = S*randn(p,1);
            y_0 = x_0'*w + randn;
            y_hat = x_0'*w_hat;
            pred_err(i,j) = (y_0 - y_hat)^2;
        end
    end
    pred_risk(k) = sum(sum(pred_err))/(num_monte_in*num_monte);
    estim_risk(k)  = mean(estim_err);
end


