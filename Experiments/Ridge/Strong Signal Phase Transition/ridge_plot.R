library(RColorBrewer)
setwd("~/git_local/high-dim-risk/Experiments/ridge_plot/")


m = function(z, gamma) {
	2 / (1 - gamma - z + sqrt((1 + gamma - z)^2 - 4 * gamma))
}

r = function(s, gamma) {
	delta = gamma / s
	1 + gamma / (1 - gamma + gamma * delta * m(-delta, gamma) + delta)
}


xx = 1.1^(-10:100)


cols = brewer.pal(3, "Set1")

pdf("ridge_plot.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(NA, NA, xlim = range(xx), ylim = c(r(min(xx), 0.25), max(xx)) - 1,log = "xy", xlab = expression("Signal Strength:" ~ alpha^2), ylab = expression("Excess Error Rate: R*" - 1))
for (gamma in c(0.25, 0.5, 0.8, 0.9, 1, 1.1, 1.3, 2, 4, 8)) {
	yy = sapply(xx, function(s) r(s, gamma))
	if (gamma < 1) col = cols[2]
	if (gamma > 1) col = cols[1]
	if (gamma == 1) col = cols[3]
	lines(xx, yy - 1, lwd = 3, col = col)
}
legend("topleft", c(expression(gamma > 1), expression(gamma == 1), expression(gamma < 1)), lwd = 3, col = cols[c(1, 3, 2)], lty = 1, cex = 1.5)
par = pardef
dev.off()