function obs_l_pdf = observation_log_pdf(y,h,param)
%%Pdf of observation distribution.
  m = param(1);
  s2 = exp(h);
  obs_l_pdf = -0.5*log(s2)-0.5*((y - m)^2)./s2;
end
