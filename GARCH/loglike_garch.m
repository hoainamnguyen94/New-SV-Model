function [llike,sig2] = loglike_garch(e,gam,sig20)

T = length(e);
e2 = e.^2;
expgam = exp(gam);
sig2 = zeros(T,1);
sig2(1) = expgam(1) + expgam(3)*sig20;

for t=2:T
    sig2(t) = expgam(1) + expgam(2)*e2(t-1) + expgam(3)*sig2(t-1);
end

llike = -T/2*log(2*pi) -.5*sum(log(sig2)) -.5*sum(e2./sig2);

end 