function [wsmooth,w2smooth] = particle_smoother(obs,wfilt,part,param)
%%Particle smoother (run after particle filter).
  T = length(part(:,1));
  N = length(part(1,:));
  one.row = ones(1,N);
  w2smooth = zeros(T-1,N,N);
  wsmooth = zeros(T,N);
  wsmooth(T,:) = wfilt(T,:);
  for t = (T-1):-1:1
    kernel = exp(kernel_log_pdf(part(t,:),part(t+1,:),obs(t),param));
    wsmooth(t,:) = wfilt(t,:)*t(kernel*t(wsmooth(t+1,:)/(wfilt(t,:)*kernel)));
    w2smooth(t,:) = (wfilt(t,:)*t(wsmooth(t+1,:))*kernel)/(t(one.row)*(wfilt(t,:)*kernel));
  end
end


