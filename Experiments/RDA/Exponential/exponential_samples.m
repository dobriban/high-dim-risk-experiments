l = 1/10;
p = 500;
mu = 1/l*ones(p,1);
R = exprnd(mu);

figure, hist(R)