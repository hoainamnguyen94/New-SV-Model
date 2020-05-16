function [ml,mlstd] = ml_sv_pl(y,store_theta,prior,hhat,M)

M = 20*ceil(M/20);
thetahat = mean(store_theta)';
thetavar = var(store_theta);
thetastd = sqrt(thetavar);
[ahat,bhat] = unifit(store_theta(:,3));
temp = gamfit(1./store_theta(:,end));
nuomegah2hat = temp(1); Somegah2hat = 1./temp(2);    

theta_IS = zeros(M,6);
theta_IS(:,1) = thetahat(1) + thetastd(1)*randn(M,1);
theta_IS(:,2) = thetahat(2) + thetastd(2)*randn(M,1);
theta_IS(:,3) = ahat + (bhat - ahat)*rand(M,1);
theta_IS(:,4) = thetahat(4) + thetastd(4)*randn(M,1);
theta_IS(:,5) = tnormrnd(thetahat(5),thetavar(5),-.999,.999,M);
theta_IS(:,6) = 1./gamrnd(nuomegah2hat,1./Somegah2hat,M,1);
store_w = zeros(M,1);

phih_const = 1/(normcdf(1,thetahat(5),thetastd(5))-normcdf(-1,thetahat(5),thetastd(5)));
gIS = @(m,a,mh,ph,oh) -.5*log(2*pi*thetavar(1))-.5*(m-thetahat(1))^2/thetavar(1) ...
      -.5*log(2*pi*thetavar(2)) -.5*(a-thetahat(2))^2/thetavar(2)...
      - log(bhat - ahat) ... 
      -.5*log(2*pi*thetavar(4))-.5*(mh-thetahat(4))^2/thetavar(4) ...
      -.5*log(2*pi*thetavar(5))+log(phih_const)-.5*(ph-thetahat(5))^2/thetavar(5) ...
      + nuomegah2hat*log(Somegah2hat)-gammaln(nuomegah2hat)-(nuomegah2hat+1)*log(oh)-Somegah2hat/oh; 

for loop = 1:M
    theta = theta_IS(loop,:)';
    mu = theta(1);
    lam = theta(2);
    kap = theta(3);
    muh = theta(4);
    phih = theta(5);
    omegah2 = theta(6);    
    llike = intlike_sv_pl(y,mu,lam,kap,muh,phih,omegah2,hhat,50);
    store_w(loop) = llike + prior(mu,lam,muh,phih,omegah2) - gIS(mu,lam,muh,phih,omegah2);    
end

shortw = reshape(store_w,M/20,20);
maxw = max(shortw);

bigml = log(mean(exp(shortw-repmat(maxw,M/20,1)),1)) + maxw;
ml = mean(bigml);
mlstd = std(bigml)/sqrt(20);

end