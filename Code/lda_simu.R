library(expm)

gamma = 1
n = 400
n.test = 10000
alpha = 2
rho = 0.8
lambda = 10

p = floor(gamma * n)

Sigma = toeplitz(rho^(0:(p-1)))
rSigma = sqrtm(Sigma)

res = replicate(4, {
	mu = alpha / sqrt(p) * rnorm(p)
	X.raw = matrix(rnorm(n * p), n, p) %*% rSigma
	Y = 2 * rbinom(n, 1, 0.5) - 1
	X = matrix(Y, n, 1) %*% matrix(mu, 1, p) + X.raw
	
	X.sgn = X * ( matrix(Y, n, 1) %*% matrix(1, 1, p))
	mu.hat = colMeans(X.sgn)
	sigma.hat = var(X.sgn)
	beta.hat = solve(sigma.hat + lambda * diag(rep(1, p)), mu.hat)
	
	Y.test = 2 * rbinom(n.test, 1, 0.5) - 1
	X.raw.test = matrix(rnorm(n.test * p), n.test, p) %*% rSigma
	X.test = matrix(Y.test, n.test, 1) %*% matrix(mu, 1, p) + X.raw.test
	Y.hat.test = sign(X.test %*% beta.hat)
	mean(Y.test != Y.hat.test)
})

sd(res/10)/mean(res)
res
