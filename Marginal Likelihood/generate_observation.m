function obs = generate_observation(h,param,N)
%%Generate observation.
  obs = param(1)+exp(h/2)*normrnd(0,1,[1,N]);
end
