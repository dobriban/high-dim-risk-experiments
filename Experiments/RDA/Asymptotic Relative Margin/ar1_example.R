library(MASS)
library(colorRamps)
library(RColorBrewer)

rm(list = ls())

setwd("~/git_local/high-dim-risk/Code/R/")
source("rda_functions.R")

gamma = 1
alpha = 1
rho = 0.9

#
# Stability test
#

get.theta = function(n) {
	p = gamma * n
	Sigma = outer(1:p, 1:p, function(x, y) rho^abs(x - y))
	Z = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
	
	mu = replicate(2, sqrt(2) * alpha * rnorm(p) / sqrt(p))
	Y = c(rep(0, n/2), rep(1, n/2))
	X = t(mu)[Y + 1,] + Z
	
	stj = get.stieltjes.func(X, Y)
	
	theta = function(lambda) {
		alpha^2 * stj$tau(lambda) / sqrt(alpha^2 * stj$eta(lambda) + stj$xi(lambda))
	}
	return(theta)
}

theta = get.theta(1500)

pdf("ar1_plugin.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

plot(NA, NA, xlim = c(0, 1.5), ylim = c(1, 2.3), xlab = expression(sqrt(lambda)), ylab = expression(Theta))

for(iter in 1:10) {
	theta.s = get.theta(150)
	plot(function(ll)theta.s(ll^2), 0, 1.5, add = TRUE, col = 2, lwd = 2)
}

plot(function(ll)theta(ll^2), 0, 1.5, lwd = 4, add = TRUE)

par = pardef
dev.off()


#
# Effect of alpha
#

gamma = 1
alpha = 10
rho = 0.9
n = 1500
p = gamma * n
Sigma = outer(1:p, 1:p, function(x, y) rho^abs(x - y))
Z = mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)

mu = replicate(2, sqrt(2) * alpha * rnorm(p) / sqrt(p))
Y = c(rep(0, n/2), rep(1, n/2))
X = mu[,Y + 1] + Z

stj = get.stieltjes.func(X, Y)

theta = function(lambda, alpha) {
	alpha^2 * stj$tau(lambda) / sqrt(alpha^2 * stj$eta(lambda) + stj$xi(lambda))
}

Delta = (1 + rho^2) / (1 + rho) / (1 - rho)

pdf("ar1_arm.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)

plot(NA, NA, xlim = c(0, 1.5), ylim = c(0, 0.25), xlab = expression(sqrt(lambda)), ylab = "Relative Margin")

alphas = seq(0.1, 2, by = 0.05)
colfunc = colorRampPalette(brewer.pal(8, "YlOrRd")[c(2, 4, 6, 7, 8)])
#colfunc = colorRampPalette(c("yellow", "red"))
cols = colfunc(length(alphas))

for(iter in 1:length(alphas)) {
	alpha = alphas[iter]
	plot(function(ll)theta(ll^2, alpha) / alpha / Delta, 0, 1.5, add = TRUE, col = cols[iter], lwd = 2)
}

par = pardef
dev.off()
