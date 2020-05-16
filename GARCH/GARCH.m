function [thetahat,thetastd,sighat,ml,mlstd] = GARCH(y,T,nloop,burnin,M,cp_ml)
    
    %% prior
mu0 = 0; Vmu = 10;
gam0 = log([exp(1) .1 .8]'); Vgam = diag([10 1 1]);
tempgam = repmat(gam0',10000,1) + (chol(Vgam,'lower')*randn(3,10000))';
prob_s = sum(sum(exp(tempgam(:,2:end)),2)<1)/10000;
c_prior = -.5*log(2*pi*Vmu) -3/2*log(2*pi) -.5*log(det(Vgam)) -log(prob_s);
lpri_gam = @(x) -.5*(x-gam0)'*(Vgam\(x-gam0));
lpri_mu = @(x) -.5*(x-mu0)^2/Vmu;
prior = @(m,g) c_prior + lpri_mu(m) + lpri_gam(g);

    %% initialize for storeage
store_theta = zeros(nloop - burnin,4); % mu, gam
store_sig2 = zeros(nloop - burnin,T); 
store_llike = zeros(nloop-burnin,1);

    %% initialize the Markov chain
mu = mean(y);
e = y-mu;
e2 = e.^2;
tempb = [0;0];
sig20 = var(y);
gamt = fminsearch(@(x)-loglike_garch(e,x,sig20)-lpri_gam(x),[tempb;tempb(2)]);
expgamt = exp(gamt);
sig2 = zeros(T,1);
dsig2 = zeros(T,3);
dsig2(1,1) = expgamt(1);
dsig2(1,3) = expgamt(3)*sig20;
sig2(1) = expgamt(1) + expgamt(3)*sig20;
for t = 2:T
    sig2(t) = expgamt(1) + expgamt(2)*e2(t-1) + expgamt(3)*sig2(t-1);
    dsig2(t,1) = expgamt(1) + expgamt(3)*dsig2(t-1,1);
    dsig2(t,2) = expgamt(2)*e2(t-1) + expgamt(3)*dsig2(t-1,2);
    dsig2(t,3) = expgamt(3)*(sig2(t-1) + dsig2(t-1,3));
end
S = repmat(.5./sig2.*(e2./sig2-1),1,3) .* dsig2;
I = S'*S;    
lprop = @(x) -.5*(x-gamt)'*I*(x-gamt);
Cgam = chol(I,'lower');
gam = gamt + Cgam'\randn(3,1);    
while sum(exp(gam(2:3))) > .999 
    gam = gamt + Cgam'\randn(3,1);    
end
ybar = mean(y);
s2 = var(y)/T;
countgam = 0;
countmu = 0;

disp('Starting GARCH.... ');
disp(' ' );

randn('seed',sum(clock*97)); rand('seed',sum(clock*37)); %#ok<RAND>
start_time = clock;

for loop = 1:nloop
    
        %% sample mu
    muc = ybar + sqrt(s2)*randn;
    [llikec,sig2c] = loglike_garch(y-muc,gam,sig20); %#ok<ASGLU>
    [llike,sig2] = loglike_garch(y-mu,gam,sig20); %#ok<ASGLU>
    alpMH = llikec + lpri_mu(muc) + .5*(muc-ybar)^2/s2 ...
        - (llike + lpri_mu(mu) + .5*(mu-ybar)^2/s2);
    if alpMH > log(rand)
        mu = muc;
        countmu = countmu + 1;        
    end
    
        %% sample gam      
    e = y-mu;
    [llike,sig2] = loglike_garch(e,gam,sig20);    
    if loop == 1 || rand>.5  %% maximize with prob 0.5
        e2 = e.^2;        
        gamt = fminsearch(@(x)-loglike_garch(e,x,sig20)-lpri_gam(x),gamt);        
        expgamt = exp(gamt);
        sig2 = zeros(T,1);
        dsig2 = zeros(T,3);
        dsig2(1,1) = expgamt(1);
        dsig2(1,3) = expgamt(3)*sig20;
        sig2(1) = expgamt(1) + expgamt(3)*sig20;
        for t=2:T
            sig2(t) = expgamt(1) + expgamt(2)*e2(t-1) + expgamt(3)*sig2(t-1);
            dsig2(t,1) = expgamt(1) + expgamt(3)*dsig2(t-1,1);
            dsig2(t,2) = expgamt(2)*e2(t-1) + expgamt(3)*dsig2(t-1,2);
            dsig2(t,3) = expgamt(3)*(sig2(t-1) + dsig2(t-1,3));
        end
        S = repmat(.5./sig2.*(e2./sig2-1),1,3) .* dsig2;
        I = S'*S;
        lprop = @(x) -.5*(x-gamt)'*I*(x-gamt);
        Cgam = chol(I,'lower');
    end
    gamc = gamt + Cgam'\randn(3,1);
    if sum(exp(gamc(2:3))) < .999 %% impose stationarity
        [llikec,sig2c] = loglike_garch(e,gamc,sig20);
        alpMH = llikec + lpri_gam(gamc) - lprop(gamc)...
             - (llike + lpri_gam(gam) - lprop(gam));
        if alpMH > log(rand)
            gam = gamc;
            sig2 = sig2c;
            llike = llikec;
            countgam = countgam + 1;        
        end
    end
    
    if loop>burnin
        i = loop-burnin;        
        store_theta(i,:) = [mu exp(gam)'];
        store_sig2(i,:) = sig2';
        store_llike(i) = llike;        
    end
    
    if ( mod( loop, 5000 ) ==0 )
        disp(  [ num2str( loop ) ' loops... ' ] )
    end    
    
end

disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

sighat = mean(sqrt(store_sig2))';  % plot std. dev.
thetahat = mean(store_theta)';
thetastd = std(store_theta)';

%% compute ml
if cp_ml
    start_time = clock;
    disp('Computing the marginal likelihood.... ');    
    [ml,mlstd] = ml_garch(y,store_theta,prior,sig20,M);
    disp( ['ML computation takes '  num2str( etime( clock, start_time) ) ' seconds' ] ); 
end

end