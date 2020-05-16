clear; clc;     
nloop = 11000;
burnin = 1000;
cp_ml = 1;
M = 10000;
 
load 'plsv_dataset.mat';

y = eurusd_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat,ml,mlstd] = GARCH_GJR(y,T,nloop,burnin,M,cp_ml); %#ok<*ASGLU>
dlmwrite('para.csv',[thetahat,thetastd]);

y = gold_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat,ml,mlstd] = GARCH_GJR(y,T,nloop,burnin,M,cp_ml);
dlmwrite('para.csv',[thetahat,thetastd],'-append');

y = sp500_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat,ml,mlstd] = GARCH_GJR(y,T,nloop,burnin,M,cp_ml);
dlmwrite('para.csv',[thetahat,thetastd],'-append');

y = wti_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat,ml,mlstd] = GARCH_GJR(y,T,nloop,burnin,M,cp_ml);
dlmwrite('para.csv',[thetahat,thetastd],'-append');