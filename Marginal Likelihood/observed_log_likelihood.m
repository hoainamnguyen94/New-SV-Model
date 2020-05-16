function likelihood = observed_log_likelihood(obs,param,N,Nth)
  filter = particle_filter(obs,param,N,Nth);
  likelihood = sum(log(filter));
end