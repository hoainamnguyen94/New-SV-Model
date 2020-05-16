function kernel = generate_kernel(hprev,yprev,param)
%%Generate Markov process from linear Gaussian kernel.
  N = length(hprev);
  kernel = param(2)+param(3)*(hprev-param(2))+param(4)*normrnd(0,1,[1,N]);
end
