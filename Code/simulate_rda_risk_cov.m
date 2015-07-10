function risk = simulate_rda_risk_cov(Sigma,n,p,alpha,lambda_arr,timing)
%Simulate the prediction error of Regularized Discriminant Analysis

if ~exist('timing','var')
    timing = 0;
end
%% Monte Carlo Evaluation
num_monte = 1e4;
rng(0)
risk = zeros(length(lambda_arr),1);
S = Sigma^(1/2);
if timing==1
    tic
end
for k=1:length(lambda_arr)
    if timing==1
        timer = toc;
        fprintf('Lambda: %d/%d; Time: %f\n',k,length(lambda_arr),timer);
    end
    lambda = lambda_arr(k);
    pred_err = zeros(num_monte, 1);
    mu = alpha*randn(p,1)/sqrt(p);
    y = sign(randn(n,1));
    X = randn(n,p)*S + y*mu';
    mu_hat = 1/n*X'*y;
    Sigma_hat = 1/n*(X'*X) - mu_hat*mu_hat';
    if lambda>0
        if lambda<Inf %lambda is a finite positive number
            beta_hat = (Sigma_hat+ lambda*eye(p)) \ mu_hat;
        else %lambda is +Infty - IR
            beta_hat = mu_hat;
        end
    else %lambda is 0 - LDA
        beta_hat = pinv(Sigma_hat+ lambda*eye(p)) * mu_hat;
    end
    
    %inner loop, generate random test data: x_0, y_0
    for j=1:num_monte
        y_0 = sign(randn);
        x_0 = y_0*mu+S*randn(p,1);
        y_hat = sign(x_0'*beta_hat);
        pred_err(j) = 1/2*abs(y_0 - y_hat);
    end
    risk(k) = mean(pred_err);
    
end



