function llike = intlike_sv_pl(y,mu,lam,kap,muh,phih,omegah2,hhat,M)

s2 = (y-mu).^2;
T = length(s2);
Hphi = speye(T) - sparse(2:T,1:(T-1),phih*ones(1,T-1),T,T);
HiSH = Hphi'*spdiags([(1-phih^2)/omegah2; 1/omegah2*ones(T-1,1)],0,T,T)*Hphi;
deltah = Hphi\[muh; muh*(1-phih)*ones(T-1,1)];
HinvSHdeltah = HiSH*deltah;
errh = 1; ht = hhat;

while errh > 10^(-3)
    expht = exp(ht);
    sinvexpht = s2./expht;
    lam2expht = lam^2.*expht;
    fh = -.5 + .5*sinvexpht - .5*lam2expht;
    Gh = .5*sinvexpht + .5*lam2expht;
    Kh = HiSH + sparse(1:T,1:T,Gh);
    newht = Kh\(fh+Gh.*ht+HinvSHdeltah);
    errh = max(abs(newht-ht));
    ht = newht;          
end

cholHh = chol(Kh,'lower');
c = -T/2*log(2*pi) -.5*(T*log(omegah2) - log(1-phih^2)) - sum(log(diag(cholHh)));
store_llike = zeros(M,1);
bigh = repmat(ht,1,M) + cholHh'\randn(T,M);

for i = 1:M    
    h = bigh(:,i);
    store_llike(i) = c - .5*(h-deltah)'*HiSH*(h-deltah) - .5*sum(h) ...
                       - .5*exp(-h)'*(y-mu-lam*(exp(h).^(kap/2))).^2 + .5*(h-ht)'*Kh*(h-ht);
end

maxllike = max(store_llike);
llike = maxllike + log(mean(exp(store_llike-maxllike)));

end

