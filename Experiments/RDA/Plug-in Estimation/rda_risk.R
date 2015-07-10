library(MASS)

get.risk = function(X) {
	
	M = X %*% t(X)
	xx = svd(M, nu = 0, nv = 0)$d
	gamma = ncol(X) / nrow(X)
	
	v = Vectorize(function(lambda) {
		mean(1 / (xx + lambda))
	})
	vp = Vectorize(function(lambda) {
		mean(1 / (xx + lambda)^2)
	})
	m = Vectorize(function(lambda) {
		(v(lambda) - 1/lambda) / gamma + 1/lambda
	})
	
	tau = Vectorize(function(lambda) {
		lambda * m(lambda) * v(lambda)
	})
	eta = Vectorize(function(lambda) {
		(v(lambda) - lambda * vp(lambda)) / gamma
	})
	xi = Vectorize(function(lambda) {
		vp(lambda) / v(lambda)^2 - 1
	})
	
	return(list(m=m, v=v, vp=vp, tau=tau, eta=eta, xi=xi))
}

p = 100
n = 100

rho = 0.9
Sigma = outer(1:p, 1:p, function(x, y) rho^abs(x - y))
X = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
ss = get.risk(X)

#plot(function(ll) log(ss$tau(ll^2)), 0, 100)
#plot(function(ll) log(ss$eta(ll^2)), 0, 100)
#plot(function(ll) log(ss$xi(ll^2)), 0, 100)

alpha = 5
plot(function(ll) alpha^2 * ss$tau(ll^2) / sqrt(alpha^2 * ss$eta(ll^2) + ss$xi(ll^2)), 0, 10)

Sigma2 = diag(qexp((1:p) / (p+1)))
X2 = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma2)
ss2 = get.risk(X2)

plot(function(ll) log(ss2$tau(ll^2)), 0, 100)
plot(function(ll) log(ss2$eta(ll^2)), 0, 100)
plot(function(ll) log(ss2$xi(ll^2)), 0, 100)

alpha = 1000
plot(function(ll) alpha^2 * ss2$tau(ll^2) / sqrt(alpha^2 * ss2$eta(ll^2) + ss2$xi(ll^2)), 0, 100)

plot(function(ll) log(ss$tau(ll^2)), 0, 100)
lines(0:100, log(ss2$tau((0:100)^2)), col = 2)

plot(function(ll) log(ss$eta(ll^2)), 0, 100)
lines(0:100, log(ss2$eta((0:100)^2)), col = 2)

plot(function(ll) log(ss$tau(ll^2)) - log(ss$eta(ll^2)) / 2, 0, 100, ylim = c(-4, -2.5))
lines(0:100, log(ss2$tau((0:100)^2)) - log(ss2$eta((0:100)^2)) / 2, col = 2)

plot(function(ll) ss$v(ll^2) * ll^2, 0, 100)
lines(0:100, ss2$v((0:100)^2) * (0:100)^2, col = 2)