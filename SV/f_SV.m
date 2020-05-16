function q_VAR = f_SV(thetahat,w_T,part_T,N)

%thetahat = [mu muh phih omega2_h]
alpha = thetahat(2)-thetahat(3)*thetahat(2);
beta = thetahat(3);
gamma = sqrt(thetahat(4));

part = [part_T;zeros(1,N)];
part(2,:) = alpha+beta*part(1,:)+gamma*randn(1,N);
resampidx = randsample(N,N,true,w_T); 
part(2,:) = part(2,resampidx); 
yhat_T1 = normrnd(thetahat(1)*ones(1,N),exp(part(2,:)/2));
q_VAR = quantile(yhat_T1,[0.01 0.025 0.05 0.99 0.975 0.95]);

end