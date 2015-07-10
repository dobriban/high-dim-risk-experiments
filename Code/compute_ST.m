function [lambda,m,v,m_prime,v_prime] = compute_ST(w,t,gamma,grid_size)
%Compute the  Stieltjes transform of the ESD of sample covariance matrix
%
%Inputs
%w,t - input spectral distribution is a mixture of point masses H = sum
%delta_{t_i} * w_i
%gamma - aspect ratio
%
%Output
%lambda - grid
% m - Stieltjes transform of the ESD of sample covariance matrix
% with PSD equal to H
% m_prime - derivative of m
% v - companion Stieltjes transform
% v_prime - derivative of v

%% compute lambda, v, m
if ~exist('grid_size','var')
grid_size = 1e5;
end
%the v-grid
v = linspace(1/grid_size,1e3,grid_size)';
%v = 1./v;
%the inverse of v;
z = zeros(grid_size,1);
for i=1:grid_size
  z(i) = -1./v(i) + gamma* sum(w.*t./(1 + v(i)*t));
end
%find the region where z<0, this corresponds to lambda>0
v = v(z<0);
z = z(z<0);
lambda = -z;

ind = (lambda<10)&(lambda>1e-2);
lambda = lambda(ind);
v = v(ind);
z = z(ind);
m = v/gamma + (1/gamma-1)./z;

%% compute m',v'
L = length(lambda);

v_prime = zeros(L,1);
for i=1:L
  v_prime(i) = 1/[1/v(i)^2 - gamma* sum( w.*t.^2./(1 + t.*v(i)).^2) ];
end

m_prime = v_prime/gamma - (1/gamma-1)./z.^2;
