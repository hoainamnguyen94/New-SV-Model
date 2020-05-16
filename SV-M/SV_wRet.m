clear; clc;
energy = csvread('WeeklyReturns.csv',1,1);
nloop = 11000;
burnin = 1000;
N = 1000;
Nth = 667;
T_total = length(energy);
iter = fix((T_total-156)/4);
final_result = zeros(4*(iter+1),54);

%select a dataset 
%1.WTI 2.Brent 3.Gasoline NYH 4.Gasoline US Gulf 
%5.NYHNo2HeatingOil 6.LA Diesel 7.US Gulf Kerosine
%8.Propane 9.HenryHubNG

for j = 1:9

result = zeros(iter+1,24); 
data = energy(:,j);

parfor i = 1:iter
   
   start_T = 4*(i-1)+1;
   end_T = start_T+155;
   train = data(start_T:end_T); %#ok<*PFBNS>
   thetahat = SV_M(train,156,nloop,burnin);
   
   [w_T,part_T] = particle_filter_SV_M(data(start_T:end_T),thetahat,N,Nth); 
   q_VAR_1 = f_SV_M(thetahat,w_T,part_T,N); 
   [w_T,part_T] = particle_filter_SV_M(data((start_T+1):(end_T+1)),thetahat,N,Nth); 
   q_VAR_2 = f_SV_M(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_M(data((start_T+2):(end_T+2)),thetahat,N,Nth); 
   q_VAR_3 = f_SV_M(thetahat,w_T,part_T,N);
   [w_T,part_T] = particle_filter_SV_M(data((start_T+3):(end_T+3)),thetahat,N,Nth); 
   q_VAR_4 = f_SV_M(thetahat,w_T,part_T,N);
   
   result(i,:) = [q_VAR_1 q_VAR_2 q_VAR_3 q_VAR_4];
              
end

delete(gcp('nocreate'));

i = iter+1;
start_T = 4*(i-1)+1;
end_T = start_T+155;
train = data(start_T:end_T);
thetahat = SV_M(train,156,nloop,burnin);
iter_k = T_total-end_T;

for k = 1:iter_k
    
    new_start_T = start_T+k-1;
    new_end_T = end_T+k-1;
    start_point = 6*(k-1)+1;
    end_point = 6*k;
    [w_T,part_T] = particle_filter_SV_M(data(new_start_T:new_end_T),thetahat,N,Nth); 
    result(i,start_point:end_point) = f_SV_M(thetahat,w_T,part_T,N);
    
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
csvwrite_with_headers('SV_wRet.csv',final_result,header);

