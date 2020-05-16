function q_VAR = f_SV_L(y_end,thetahat,w_T,part_T,N)

%thetahat = [mu rho muh phih omega2_h]
alpha = thetahat(3)-thetahat(4)*thetahat(3);
beta = thetahat(4);
gamma = sqrt(thetahat(5));
rho = thetahat(2);

part = [part_T;zeros(1,N)];
yhat = [y_end*ones(1,N);zeros(1,N)];
part(2,:) = alpha+beta*part(1,:)+gamma* ...
            (rho*(yhat(1,:)-thetahat(1)).*exp(-0.5*part(1,:))+ ...
            sqrt(1-(rho^2))*randn(1,N));
resampidx = randsample(N,N,true,w_T); 
part(2,:) = part(2,resampidx);
yhat(2,:) = normrnd(thetahat(1)*ones(1,N),exp(part(2,:)/2));
q_VAR = quantile(yhat(2,:),[0.01 0.025 0.05 0.99 0.975 0.95]);

end