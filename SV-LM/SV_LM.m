function [thetahat,thetastd,hhat] = SV_LM(y,T,nloop,burnin)

phih0 = .97; Vphih = .1^2;
mu0 = 0; Vmu = 10;
muh0 = 1; Vmuh = 10;
lam0 = 0; Vlam = 100; 
nuh = 5; Sh = .2^2*(nuh-1);

%% initialize the Markov chain
muh = log(var(y)); phih = .98; 
omegah2 = .2^2;
h = muh + sqrt(omegah2)*randn(T,1);
exph = exp(h);

%% initialize for storeage
store_theta = zeros((nloop - burnin),5); % [mu lam muh phih omegah2]
store_h = zeros((nloop - burnin),T);

%% compute a few things outside the loop    
Hphi = speye(T) - sparse(2:T,1:(T-1),phih*ones(1,T-1),T,T);
newnuh = T/2 + nuh;
counth = 0; countphih = 0;
disp('Starting new_SV.... ');
disp(' ' );

start_time = clock;

for loop = 1:nloop 
    %% sample mu and lam
    X = [ones(T,1) h];
    invSy = sparse(1:T,1:T,1./exph);
    XinvSy = X'*invSy;
    invDbeta = sparse(1:2,1:2,[1/Vmu 1/Vlam]) + XinvSy*X;
    betahat = invDbeta\([mu0/Vmu; lam0/Vlam] + XinvSy*y);
    beta = betahat + chol(invDbeta,'lower')'\randn(2,1);
    mu = beta(1); 
    lam = beta(2);
    
    %% sample h    
    HiSH = Hphi'*spdiags([(1-phih^2)/omegah2; 1/omegah2*ones(T-1,1)],0,T,T)*Hphi;
    deltah = Hphi\[muh; muh*(1-phih)*ones(T-1,1)];
    HinvSHdeltah = HiSH*deltah;
    s = y-mu;
    errh = 1; ht = h;
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
    % AR-step:     
    lph = @(x) -.5*(x-deltah)'*HiSH*(x-deltah) -.5*sum(x) ...
        -.5*exp(-x)'*(y-mu-lam*x).^2;
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
    errh = [(h(1)-muh)*sqrt(1-phih^2);h(2:end)-phih*h(1:end-1)-muh*(1-phih)];    
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
        store_theta(i,:) = [mu lam muh phih omegah2];  
    end    
    if ( mod( loop, 5000 ) ==0 )
        disp(  [ num2str( loop ) ' loops... ' ] )
    end    
end

disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

hhat = mean(exp(store_h/2))'; 
thetahat = mean(store_theta)';
thetastd = std(store_theta)';

end






