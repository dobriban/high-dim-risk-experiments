function [X, TrueSpectrum] = generateCorrelatedData(n,p,testType,rho,diff)
% Generate Data with specified population covariance

% Inputs
% n - sample size
% p - dimension
% other inputs -parameters of covariance matrix
% testType - the name of the spectrum, see list specified below
% rho - covariance of AR(1) model
% diff - separation between spikes in two-cluster model
% 
%
% Outputs
% X - n x p data matrix, normalized by 1/sqrt(n)
% TrueSectrum - spectrum of population covariance matrix



switch testType
    
    case 'ones'
        TrueSpectrum = ones(p,1);
        
    case 'cluster'
        numClusters = 3;
        TrueSpectrum = floor(numClusters*rand(p,1))+1;
        
    case 'two-point'
        numClusters = 2;
        TrueSpectrum = ones(p,1)+diff*(floor(numClusters*rand(p,1)));
    case 'MP'
        
        X = 1/sqrt(n)*randn(n,p)*sqrt(Sigma);
        TrueSpectrum = svd(X,'econ').^2;
        
    case 'beta'
        TrueSpectrum = 1+10*betarnd(1,10,p,1);
        
    case 'AR1'
        %rho = 0.3;
        r = rho.^(0:1:p-1);
        Sigma = toeplitz(r);
        TrueSpectrum = eig(Sigma);
        
end
TrueSpectrum = sort(TrueSpectrum);

Sigma = diag(TrueSpectrum);
sqrtSigma = diag(sqrt(TrueSpectrum));
X = 1/sqrt(n)*randn(n,p)*sqrtSigma;
