function llike = intlike_sv_m(y,mu,lam,muh,phih,omegah2,M)

s = y-mu;
T = length(s);
Hphi = speye(T) - sparse(2:T,1:(T-1),phih*ones(1,T-1),T,T);
HiSH = Hphi'*spdiags([(1-phih^2)/omegah2; 1/omegah2*ones(T-1,1)],0,T,T)*Hphi;
deltah = Hphi\[muh; muh*(1-phih)*ones(T-1,1)];
HinvSHdeltah = HiSH*deltah;
errh = 1;
ht = zeros(T,1);

while errh > 10^(-3)
    expht = exp(ht);
    lamhtexpht = lam*ht./expht;
    fh = -.5 + .5*(s.^2)./expht + lam*s./expht - s.*lamhtexpht + ...
          .5*lam*ht.*lamhtexpht - lam*lamhtexpht;
    Gh = .5*(s.^2)./expht + 2*lam*s./expht - s.*lamhtexpht + ...
         .5*lam*ht.*lamhtexpht - 2*lam*lamhtexpht + (lam^2)./expht;
    Kh = HiSH + sparse(1:T,1:T,Gh);
    newht = Kh\(fh+Gh.*ht+HinvSHdeltah);
    errh = max(abs(newht-ht));
    ht = newht; 
end

cholHh = chol(Kh,'lower');
c = -T/2*log(2*pi) -.5*(T*log(omegah2) - log(1-phih^2)) - sum(log(diag(cholHh)));
store_llike = zeros(M,1);
bigh = repmat(ht,1,M) + cholHh'\randn(T,M);

for i=1:M    
    h = bigh(:,i);   
    store_llike(i) = c - .5*(h-deltah)'*HiSH*(h-deltah) - .5*sum(h) ...
                       - .5*exp(-h)'*(y-mu-lam*h).^2 + .5*(h-ht)'*Kh*(h-ht);
end

maxllike = max(store_llike);
llike = maxllike + log(mean(exp(store_llike-maxllike)));

end

