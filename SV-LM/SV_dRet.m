clear; clc; 
energy = csvread('DailyReturns.csv',1,1);
nloop = 11000;
burnin = 1000;
N = 1000;
Nth = 667;
T_total = length(energy);
iter = fix((T_total-1095)/20);
q_VAR = zeros(T_total-1095,6);
result = zeros(T_total-1095,54);

%1.WTI 2.Brent 3.Gasoline NYH 4.Gasoline US Gulf 
%5.NYHNo2HeatingOil 6.LA Diesel 7.US Gulf Kerosine
%8.Propane 9.HenryHubNG

for j = 1:9

data = energy(:,j);

parfor i = 1:iter
   
   start_T = 20*(i-1)+1;
   end_T = start_T+1094;
   train = data(start_T:end_T); %#ok<*PFBNS>
   thetahat = SV_LM(train,1095,nloop,burnin);
      
   [w_T,part_T] = particle_filter_SV_LM(data(start_T:end_T),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N); 
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+1):(end_T+1)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+2):(end_T+2)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+3):(end_T+3)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+4):(end_T+4)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+5):(end_T+5)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+6):(end_T+6)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+7):(end_T+7)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+8):(end_T+8)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N); 
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+9):(end_T+9)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+10):(end_T+10)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+11):(end_T+11)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+12):(end_T+12)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+13):(end_T+13)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+14):(end_T+14)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+15):(end_T+15)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+16):(end_T+16)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+17):(end_T+17)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+18):(end_T+18)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_LM(data((start_T+19):(end_T+19)),thetahat,N,Nth); 
   q_VAR(i,:) = f_SV_LM(thetahat,w_T,part_T,N);
   
end

delete(gcp('nocreate'));

i = iter+1;
start_T = 20*(i-1)+1;
end_T = start_T+1094;
train = data(start_T:end_T);
thetahat = SV_LM(train,1095,nloop,burnin);
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
csvwrite_with_headers('SV_dRet.csv',result,header);

