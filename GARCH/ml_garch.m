function [ml,mlstd] = ml_garch(y,store_theta,prior,sig20,M)

M = 20*ceil(M/20);
store_gam = log(store_theta(:,2:end));
muhat = mean(store_theta(:,1))';
mupre = 1/var(store_theta(:,1));
mustd = sqrt(var(store_theta(:,1)));
gamhat = mean(store_gam)';
gamstd = chol(cov(store_gam),'lower');
gampre = cov(store_gam)\speye(3);

theta_IS = zeros(M,4);
theta_IS(:,1) = muhat + mustd*randn(M,1); 
theta_IS(:,2:end) = repmat(gamhat',M,1) + (gamstd*randn(3,M))';

store_w = zeros(M,1);
gIS = @(m,g) -.5*log(2*pi/mupre) -.5*(m-muhat)^2*mupre ...
    -3/2*log(2*pi) - sum(log(diag(gamstd))) -.5*(g-gamhat)'*gampre*(g-gamhat);

for loop = 1:M
    theta = theta_IS(loop,:)';
    if sum(exp(theta(3:4)))<1  % check the draw is stationary
        mu = theta(1);
        gam = theta(2:end);
        llike = loglike_garch(y-mu,gam,sig20);    
        store_w(loop) = llike + prior(mu,gam) - gIS(mu,gam);
    else
        store_w(loop) = -inf;
    end
end

shortw = reshape(store_w,M/20,20);
maxw = max(shortw);

bigml = log(mean(exp(shortw-repmat(maxw,M/20,1)),1)) + maxw;
ml = mean(bigml);
mlstd = std(bigml)/sqrt(20);

end