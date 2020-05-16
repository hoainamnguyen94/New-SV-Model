function q_VAR = f_SV_LM(thetahat,w_T,part_T,N)

%thetahat = [mu lambda muh phih omega2_h]
alpha = thetahat(3)-thetahat(4)*thetahat(3);
beta = thetahat(4);
gamma = sqrt(thetahat(5));
mu = thetahat(1); lambda = thetahat(2);

part = [part_T;zeros(1,N)];
part(2,:) = alpha+beta*part(1,:)+gamma*randn(1,N);
resampidx = randsample(N,N,true,w_T); 
part(2,:) = part(2,resampidx); 
yhat_T1 = normrnd(mu+lambda*part(2,:),exp(part(2,:)/2));
q_VAR = quantile(yhat_T1,[0.01 0.025 0.05 0.99 0.975 0.95]);

end