clear; clc;     
nloop = 11000;
burnin = 1000;       

load 'plsv_dataset.mat';

y = eurusd_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat] = SV_M(y,T,nloop,burnin); %#ok<*ASGLU>

y = gold_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat] = SV_M(y,T,nloop,burnin);

y = sp500_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat] = SV_M(y,T,nloop,burnin);

y = wti_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat] = SV_M(y,T,nloop,burnin);