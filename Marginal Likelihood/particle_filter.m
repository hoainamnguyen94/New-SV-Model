function like = particle_filter(obs,param,N,Nth)
%%Generic particle filter.
  T = length(obs);
  part = zeros(T,N);
  w = zeros(T,N);
  Neff = zeros(1,T);
  like = zeros(1,T);
  part(1,:) = generate_initial(param,N);
  wei = exp(observation_log_pdf(obs(1),part(1,:),param));
  w(1,:) = wei/sum(wei);    
  Neff(1) = 1/sum(w(1,:).^2);
  like(1) = sum(wei);
  for t = 2:T
    if Neff(t-1) <= Nth
      resampidx = randsample(N,N,true,w(t-1,:));
      part(t,:) = generate_kernel(part(t-1,resampidx),obs(t-1),param);
      wei = exp(observation_log_pdf(obs(t),part(t,:),param));
    else
      part(t,:) = generate_kernel(part(t-1,:),obs(t-1),param);
      wei = w(t-1,:).*exp(observation_log_pdf(obs(t),part(t,:),param));
     end
    w(t,:) = wei/sum(wei);
    Neff(t) = 1/sum(w(t,:).^2);
    like(t) = sum(wei);
  end
end
