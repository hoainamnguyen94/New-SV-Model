function [thetahat,thetastd,hhat,ml,mlstd] = SV_PL(y,T,nloop,burnin,M,cp_ml)

phih0 = .97; Vphih = .1^2;
mu0 = 0; Vmu = 10; 
muh0 = 1; Vmuh = 10; 
lam0 = 0; Vlam = 100; 
nuh = 5; Sh = .2^2*(nuh-1);
phih_const = 1/(normcdf(1,phih0,sqrt(Vphih))-normcdf(-1,phih0,sqrt(Vphih)));
prior = @(m,a,mh,ph,oh) -.5*log(2*pi*Vmu)-.5*(m-mu0)^2/Vmu ... 
        -.5*log(2*pi*Vlam) -.5*(a-lam0)^2/Vlam ... 
        - log(2) ... 
        -.5*log(2*pi*Vmuh)-.5*(mh-muh0)^2/Vmuh ... 
        -.5*log(2*pi*Vphih)+log(phih_const)-.5*(ph-phih0)^2/Vphih ... 
        + nuh*log(Sh)-gammaln(nuh)-(nuh+1)*log(oh)-Sh/oh;

    %% initialize the Markov chain
lam = 0; kap = 0.5;
muh = log(var(y)); phih = .98; 
omegah2 = .2^2;
h = muh + sqrt(omegah2)*randn(T,1);
exph = exp(h);

    %% initialize for storage
store_theta = zeros((nloop - burnin),6); % [mu lam kap muh phih omegah2]
store_h = zeros((nloop - burnin),T);

    %% compute a few things outside the loop   
Hphi = speye(T) - sparse(2:T,1:(T-1),phih*ones(1,T-1),T,T);
newnuh = T/2 + nuh;
counth = 0; countphih = 0;
disp('Starting SV-PL.... ');
disp(' ' );

start_time = clock;

for loop = 1:nloop 
        %% sample mu
    invDmu = 1/Vmu + sum(1./exph);
    muhat = invDmu\sum((y-lam*(exph.^(kap/2)))./exph);   
    mu = muhat + chol(invDmu,'lower')'\randn;
    
        %% sample lambda (rejection sampling) 
    n = y-mu;
    m = exph.^(kap/2);
    S1 = sum(m.*n./exph); S2 = sum(m.*m./exph);
    lam_max = S1/S2;
    i = 0;
    while i == 0
        lam_temp = lam0+(Vlam^.5)*randn();
        log_ratio = (lam_temp-lam_max)*S1 - .5*((lam_temp^2)-(lam_max^2))*S2;
        ratio = exp(log_ratio);
        U = rand();
        if U < ratio
            lam = lam_temp;
            i = i+1;
        end
    end
    
        %% sample kappa (importance sampling)
    kap_vec = -1+2*rand(1,50);
    exph_kap = exph.^(kap_vec/2);
    mean_y = mu+lam*exph_kap;
    log_wei_vec = -.5*sum(((y-mean_y).^2)./exph);
    wei_vec = exp(log_wei_vec-mean(log_wei_vec));
    kap = randsample(kap_vec,1,true,wei_vec);
    
        %% sample h    
    HiSH = Hphi'*spdiags([(1-phih^2)/omegah2; 1/omegah2*ones(T-1,1)],0,T,T)*Hphi;
    deltah = Hphi\[muh; muh*(1-phih)*ones(T-1,1)];
    HinvSHdeltah = HiSH*deltah;
    s2 = (y-mu).^2;
    errh = 1; ht = h;
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
    % AR-step:     
    lph = @(x) -.5*(x-deltah)'*HiSH*(x-deltah) -.5*sum(x) ...
        -.5*exp(-x)'*(y-mu-lam*(exp(x).^(kap/2))).^2;
    logc = lph(ht) + log(3);   
    flag = 0;
    while flag == 0
        hc = ht + cholHh'\randn(T,1);        
        alpARc =  lph(hc) + .5*(hc-ht)'*Kh*(hc-ht) - logc;
        if alpARc > log(rand)
            flag = 1;
        end
    end        
    % MH-step    
    alpAR = lph(h) + .5*(h-ht)'*Kh*(h-ht) - logc;
    if alpAR < 0
        alpMH = 1;
    elseif alpARc < 0
        alpMH = - alpAR;
    else
        alpMH = alpARc - alpAR;
    end    
    if alpMH > log(rand) || loop == 1
        h = hc;
        exph = exp(h);
        counth = counth + 1;
    end
    
        %% sample omegah2
    errh = [(h(1)-muh)*sqrt(1-phih^2);  h(2:end)-phih*h(1:end-1)-muh*(1-phih)];    
    newSh = Sh + sum(errh.^2)/2;    
    omegah2 = 1/gamrnd(newnuh, 1./newSh);
    
        %% sample phih
    Xphi = h(1:end-1)-muh;
    yphi = h(2:end) - muh;
    Dphi = 1/(1/Vphih + Xphi'*Xphi/omegah2);
    phihat = Dphi*(phih0/Vphih + Xphi'*yphi/omegah2);
    phic = phihat + sqrt(Dphi)*randn;
    g = @(x) -.5*log(omegah2./(1-x.^2))-.5*(1-x.^2)/omegah2*(h(1)-muh)^2;
    if abs(phic)<.9999
        alpMH = exp(g(phic)-g(phih));
        if alpMH>rand
            phih = phic;
            countphih = countphih + 1;
            Hphi = speye(T) - sparse(2:T,1:(T-1),phih*ones(1,T-1),T,T);
        end
    end    
    
        %% sample muh    
    Dmuh = 1/(1/Vmuh + ((T-1)*(1-phih)^2 + (1-phih^2))/omegah2);
    muhat = Dmuh*(muh0/Vmuh + (1-phih^2)/omegah2*h(1) + (1-phih)/omegah2*sum(h(2:end)-phih*h(1:end-1)));
    muh = muhat + sqrt(Dmuh)*randn;   
    if loop>burnin
        i = loop-burnin;     
        store_h(i,:)  = h';         
        store_theta(i,:) = [mu lam kap muh phih omegah2];  
    end    
    if (mod(loop,2000)==0)
        disp([num2str(loop) ' loops... '])
    end       
end

disp( ['MCMC takes ' num2str(etime(clock,start_time)) ' seconds']);
disp(' ');

hhat = mean(exp(store_h/2))';
thetahat = mean(store_theta)';
thetastd = std(store_theta)';

%% compute ml
if cp_ml
    start_time = clock;
    disp('Computing the marginal likelihood.... ');    
    [ml,mlstd] = ml_sv_pl(y,store_theta,prior,hhat,M);    
    disp(['ML computation takes ' num2str( etime( clock, start_time) ) ' seconds']);    
end

end

