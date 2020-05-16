clear; clc;     
nloop = 11000;
burnin = 1000;
cp_ml = 1;
M = 10000;
 
load 'plsv_dataset.mat';

y = eurusd_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat,ml,mlstd] = GARCH(y,T,nloop,burnin,M,cp_ml); %#ok<*ASGLU>

y = gold_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat,ml,mlstd] = GARCH(y,T,nloop,burnin,M,cp_ml);

y = sp500_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat,ml,mlstd] = GARCH(y,T,nloop,burnin,M,cp_ml);

y = wti_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat,ml,mlstd] = GARCH(y,T,nloop,burnin,M,cp_ml);