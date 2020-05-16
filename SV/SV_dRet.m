clear; clc; 
energy = csvread('DailyReturns.csv',1,1);
nloop = 11000;
burnin = 1000;
N = 1000;
Nth = 667;
T_total = length(energy);
iter = fix((T_total-1095)/20);
final_result = zeros(20*(iter+1),54);    

%1.WTI 2.Brent 3.Gasoline NYH 4.Gasoline US Gulf 
%5.NYHNo2HeatingOil 6.LA Diesel 7.US Gulf Kerosine
%8.Propane 9.HenryHubNG

for j = 1:9

result = zeros(iter+1,120);    
data = energy(:,j);

parfor i = 1:iter
   
   start_T = 20*(i-1)+1;
   end_T = start_T+1094;
   train = data(start_T:end_T); %#ok<*PFBNS>
   thetahat = SV(train,1095,nloop,burnin);
      
   [w_T,part_T] = particle_filter_SV(data(start_T:end_T),thetahat,N,Nth); 
   q_VAR_1 = f_SV(thetahat,w_T,part_T,N); 
   [w_T,part_T] = particle_filter_SV(data((start_T+1):(end_T+1)),thetahat,N,Nth); 
   q_VAR_2 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+2):(end_T+2)),thetahat,N,Nth); 
   q_VAR_3 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+3):(end_T+3)),thetahat,N,Nth); 
   q_VAR_4 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+4):(end_T+4)),thetahat,N,Nth); 
   q_VAR_5 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+5):(end_T+5)),thetahat,N,Nth); 
   q_VAR_6 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+6):(end_T+6)),thetahat,N,Nth); 
   q_VAR_7 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+7):(end_T+7)),thetahat,N,Nth); 
   q_VAR_8 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+8):(end_T+8)),thetahat,N,Nth); 
   q_VAR_9 = f_SV(thetahat,w_T,part_T,N); 
   [w_T,part_T] = particle_filter_SV(data((start_T+9):(end_T+9)),thetahat,N,Nth); 
   q_VAR_10 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+10):(end_T+10)),thetahat,N,Nth); 
   q_VAR_11 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+11):(end_T+11)),thetahat,N,Nth); 
   q_VAR_12 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+12):(end_T+12)),thetahat,N,Nth); 
   q_VAR_13 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+13):(end_T+13)),thetahat,N,Nth); 
   q_VAR_14 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+14):(end_T+14)),thetahat,N,Nth); 
   q_VAR_15 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+15):(end_T+15)),thetahat,N,Nth); 
   q_VAR_16 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+16):(end_T+16)),thetahat,N,Nth); 
   q_VAR_17 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+17):(end_T+17)),thetahat,N,Nth); 
   q_VAR_18 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+18):(end_T+18)),thetahat,N,Nth); 
   q_VAR_19 = f_SV(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV(data((start_T+19):(end_T+19)),thetahat,N,Nth); 
   q_VAR_20 = f_SV(thetahat,w_T,part_T,N);
   
   result(i,:) = [q_VAR_1 q_VAR_2 q_VAR_3 q_VAR_4 q_VAR_5 q_VAR_6 q_VAR_7 ...
                  q_VAR_8 q_VAR_9 q_VAR_10 q_VAR_11 q_VAR_12 q_VAR_13 q_VAR_14 ...
                  q_VAR_15 q_VAR_16 q_VAR_17 q_VAR_18 q_VAR_19 q_VAR_20];
   
end

delete(gcp('nocreate'));

i = iter+1;
start_T = 20*(i-1)+1;
end_T = start_T+1094;
train = data(start_T:end_T);
thetahat = SV(train,1095,nloop,burnin);
iter_k = T_total-end_T;

for k = 1:iter_k
    
    new_start_T = start_T+k-1;
    new_end_T = end_T+k-1;
    start_point = 6*(k-1)+1;
    end_point = 6*k;
    [w_T,part_T] = particle_filter_SV(data(new_start_T:new_end_T),thetahat,N,Nth); 
    result(i,start_point:end_point) = f_SV(thetahat,w_T,part_T,N);
    
end

reshaped_mat = reshaping(result,iter);
start_col = 6*(j-1)+1;
end_col = 6*j;
final_result(:,start_col:end_col) = reshaped_mat;

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
csvwrite_with_headers('SV_dRet.csv',final_result,header);

