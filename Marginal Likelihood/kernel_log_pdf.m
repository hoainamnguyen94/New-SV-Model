function kernel_l_pdf = kernel_log_pdf(hprev,h,yprev,param) %#ok<*INUSL>
%%Pdf of Markov kernel.
  one.row = matrix(1, 1, length(h));
  m = param(2)+param(3)*(hprev-param(2));
  s2 = param(4)^2;
  kernel_l_pdf = -0.5*log(s2)-0.5*((t(one.row)*h-m*one.row)^2)/s2;
end
