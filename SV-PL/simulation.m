mu = 0.3;
lambda = 0.1;
kappa = 0.8;
muh = -1.2;
phih = 0.85;
omega2h = 0.1;
h = zeros(1000,1);
h(1) = normrnd(muh,sqrt(omega2h/(1-phih^2)));
for i = 2:1000
   h(i) = muh + phih*(h(i-1)-muh) + normrnd(0,sqrt(omega2h)); 
end
y = zeros(1000,1);
for i = 1:1000
   y(i) = mu + lambda*exp(kappa*h(i)/2) + normrnd(0,exp(h(i)/2)); 
end