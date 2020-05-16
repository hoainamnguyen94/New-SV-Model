function initial_l_pdf = initial_log_pdf(h,param) %#ok<*DEFNU>
%%Pdf of initial state.
  m = 0;
  s2 = (param(4)/sqrt(1-param(3)^2))^2;
  initial_l_pdf = -0.5*log(s2)-0.5*((h - m)^2)./s2;
end
