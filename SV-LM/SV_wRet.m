clear; clc;
energy = csvread('WeeklyReturns.csv',1,1);
nloop = 11000;
burnin = 1000;
N = 1000;
Nth = 667;
T_total = length(energy);
iter = fix((T_total-52)/4);
q_VAR = zeros(T_total-52,6);
result = zeros(T_total-52,54);

%select a dataset 
%1.WTI 2.Brent 3.Gasoline NYH 4.Gasoline US Gulf 
%5.NYHNo2HeatingOil 6.LA Diesel 7.US Gulf Kerosine
%8.Propane 9.HenryHubNG

for j = 1:9

data = energy(:,j);

parfor i = 1:iter
   
   start_T = 4*(i-1)+1;
   end_T = start_T+51;
   train = data(start_T:end_T); %#ok<*PFBNS>
   thetahat = SV_LM(train,52,nloop,burnin);
   
   [w_T,part_T] = particle_filter_SV_LM(data(start_T:end_T),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N); 
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+1):(end_T+1)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+2):(end_T+2)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+3):(end_T+3)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
      
end

delete(gcp('nocreate'));

i = iter+1;
start_T = 4*(i-1)+1;
end_T = start_T+51;
train = data(start_T:end_T);
thetahat = SV_LM(train,52,nloop,burnin);
iter_k = T_total-end_T;

for k = 1:iter_k
    
    new_start = start_T+k-1;
    new_end = end_t+k-1;
    [w_T,part_T] = particle_filter_SV_LM(data(new_start:new_end),thetahat,N,Nth); 
    q_VAR(new_start,:) = f_SV_LM(thetahat,w_T,part_T,N);
    
end

start_point = 6*(j-1)+1;
end_point = 6*j;
result(:,start_point:end_point) = q_VAR;

end

header = {'1.00%','2.50%','5.00%','99.00%','97.50%','95.00%', ...
          '1.00%','2.50%','5.00%','99.00%','97.50%','95.00%', ...
          '1.00%','2.50%','5.00%','99.00%','97.50%','95.00%', ...
          '1.00%','2.50%','5.00%','99.00%','97.50%','95.00%', ...
          '1.00%','2.50%','5.00%','99.00%','97.50%','95.00%', ...
          '1.00%','2.50%','5.00%','99.00%','97.50%','95.00%', ...
          '1.00%','2.50%','5.00%','99.00%','97.50%','95.00%', ...
          '1.00%','2.50%','5.00%','99.00%','97.50%','95.00%', ...
          '1.00%','2.50%','5.00%','99.00%','97.50%','95.00%'};
csvwrite_with_headers('SV_wRet.csv',result,header);

