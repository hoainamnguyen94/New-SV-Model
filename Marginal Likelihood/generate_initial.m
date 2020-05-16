function initial = generate_initial(param,N)
%%Generate initial state.
  initial = param(4)/sqrt(1-param(3)^2)*normrnd(0,1,[1,N]);
end