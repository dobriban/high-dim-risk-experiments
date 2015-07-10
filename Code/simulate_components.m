function [new_lambda_arr,m_emp,m_prime_emp,K_emp,K_prime_emp,Psi_emp] = simulate_components(w,t,gamma,lambda_array)

%% downscale lambda
l = 20;
new_lambda_arr  = linspace(min(lambda_array),max(lambda_array),l);

%% Monte Carlo Evaluation
%set number of iterations in inner and outer MC loop
num_monte = 100;
n = 200;
p = floor(gamma*n);
rng(0)
m_emp = zeros(l,1);
m_prime_emp = zeros(l,1);
K_emp = zeros(l,1);
K_prime_emp = zeros(l,1);
Psi_emp = zeros(l,1);

tic
for k=1:l
  timer = toc;
  fprintf('Lambda: %d/%d; Time: %f\n',k,l,timer);
  lambda = new_lambda_arr(k);
  m_MC = zeros(num_monte,1);
  m_prime_MC = zeros(num_monte,1);
  K_MC = zeros(num_monte,1);
  K_prime_MC = zeros(num_monte,1);
  Psi_MC = zeros(num_monte,1);
  s = sqrt(random_discrete_draw(t,w,p));
  
  for i=1:num_monte
    mu = randn(p,1)/sqrt(p);
    y = sign(randn(n,1));
    
    %X = randn(n,p)*diag(s);
    X = randn(n,p)*diag(s) + y*mu';
    mu_hat = 1/n*X'*y;  % sigma^{1/2}
    %Sigma_hat = 1/n*(X'*X);
    Sigma_hat = 1/n*(X'*X) - mu_hat*mu_hat';
    beta_hat = (Sigma_hat+ lambda*eye(p)) \ mu_hat; % sigma^{1/2}, sigma_hat^{-1}
    a = (Sigma_hat+ lambda*eye(p)) \ mu;   %  sigma_hat^{-1}
    c = (Sigma_hat+ lambda*eye(p)) \ (s.*mu);  % sigma^{1}, sigma_hat^{-1}
    
    m_MC(i) = mu'*a; %sigma_hat^{-1}
    m_prime_MC(i) = a'*a; %sigma_hat^{-2}
    K_MC(i) = (s.*mu)'*(s.*a); % sigma^{1}, sigma_hat^{-1}
    K_prime_MC(i) = -a'*diag(s.^2)*a; %sigma^{1}, sigma_hat^{-2}
    Psi_MC(i) = c'*diag(s.^2)*c; %sigma^{1}, sigma_hat^{-2}
  end
  
  m_emp(k) = mean(m_MC);
  m_prime_emp(k) = mean(m_prime_MC);
  K_emp(k) = mean(K_MC);
  K_prime_emp(k) = mean(K_prime_MC);
  Psi_emp(k) = mean(Psi_MC);
  
end



