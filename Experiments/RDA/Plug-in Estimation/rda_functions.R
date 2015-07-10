get.stieltjes.func = function(X, Y) {
	
	M = ((sum(Y == 0) - 1) * var(X[Y == 0,])
	      + (sum(Y == 1) - 1) * var(X[Y == 1,])) /
	          (length(Y) - 2)
	xx = svd(M, nu = 0, nv = 0)$d
	gamma = ncol(X) / nrow(X)
	
	m = Vectorize(function(lambda) {
		mean(1 / (xx + lambda))
	})
	mp = Vectorize(function(lambda) {
		mean(1 / (xx + lambda)^2)
	})
	v = Vectorize(function(lambda) {
		(m(lambda) - 1/lambda) * gamma + 1/lambda
	})
	vp = Vectorize(function(lambda) {
		(mp(lambda) - 1/lambda^2) * gamma + 1/lambda^2
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

rda = function(X, Y, lambdas) {
	mu0 = colMeans(X[Y==0,])
	mu1 = colMeans(X[Y==1,])
	mu = (mu0 + mu1) / 2
	delta = (mu1 - mu0) / 2
	
	sigma = ((sum(Y==0) - 1) * var(X[Y==0,]) + (sum(Y==1) - 1) * var(X[Y==0,])) / (length(Y) - 2)	
	sigma.svd = svd(sigma)
	
	w = sapply(lambdas, function(lambda) {
		sigma.svd$v %*% diag(1 / (sigma.svd$d + lambda)) %*% t(sigma.svd$u) %*% delta
	})
	
	return(list(mu=mu, delta=delta, w=w))
}
		
		
		
		