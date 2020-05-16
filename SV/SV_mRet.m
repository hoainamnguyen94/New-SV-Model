clear; clc;
energy = csvread('MonthlyReturns.csv',1,1);
nloop = 11000;
burnin = 1000;
N = 1000;
Nth = 667;
T_total = length(energy);
q_VAR = zeros(T_total-36,6);
result = zeros(T_total-36,54);

%select a dataset 
%1.WTI 2.Brent 3.Gasoline NYH 4.Gasoline US Gulf 
%5.NYHNo2HeatingOil 6.LA Diesel 7.US Gulf Kerosine
%8.Propane 9.HenryHubNG

for j = 1:9

data = energy(:,j);

parfor i = 1:(T_total-36)
    
   train = data(i:(i+35)); %#ok<*PFBNS>
   thetahat = SV(train,36,nloop,burnin); 
   [w_T,part_T] = particle_filter_SV(train,thetahat,N,Nth); 
   q_VAR(i,:) = f_SV(thetahat,w_T,part_T,N);
   
end

delete(gcp('nocreate'));
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
csvwrite_with_headers('SV_mRet.csv',result,header);

