clear; clc;
cp_ml = 0;     
nloop = 11000;
burnin = 1000;
M = 10000;          

load('plsv_dataset.mat');  

%% run the estimation and calculate the marginal likelihood for 4 sub-intervals

filename = 'Results_PL_SV.xlsx';

tlower = [datetime(2003,1,1);datetime(2006,1,1); ...
          datetime(2010,1,1);datetime(2014,1,1)];
tupper = [datetime(2005,12,31);datetime(2009,12,31); ...
          datetime(2013,12,31);datetime(2017,12,31)];
      
% run the code for the eurusd_ret data set 

sheet = 'eurusd_ret';
col_headers = {'All','2003-2005','2006-2009','2010-2013','2014-2017'};
row_headers = {'mu';'lambda';'kap';'mu_h';'phi_h';'omega2_h';'log marginal likelihood'};
xlswrite(filename,col_headers,sheet,'B2:F2');
xlswrite(filename,col_headers,sheet,'G2:K2');
xlswrite(filename,row_headers,sheet,'A3');

for i = 1:4
    y = eurusd_ret(isbetween(ret_date,tlower(i),tupper(i)));
    y = y(y~=0);
    T = length(y);    
    [thetahat,thetastd,hhat] = Copy_of_SV_PL(y,T,nloop,burnin,M,cp_ml);
    switch i
        case 1
            xlRange = 'C3';
            xlRange_sd = 'H3';
        case 2
            xlRange = 'D3';
            xlRange_sd = 'I3';
        case 3
            xlRange = 'E3';
            xlRange_sd = 'J3';
        case 4
            xlRange = 'F3';
            xlRange_sd = 'K3';
    end
    xlswrite(filename,[thetahat;ml],sheet,xlRange);
    xlswrite(filename,[thetastd;mlstd],sheet,xlRange_sd);
end

% run the code for the gold_ret data set 

sheet = 'gold_ret';
col_headers = {'All','2003-2005','2006-2009','2010-2013','2014-2017'};
row_headers = {'mu';'lambda';'kap';'mu_h';'phi_h';'omega2_h';'log marginal likelihood'};
xlswrite(filename,col_headers,sheet,'B2:F2');
xlswrite(filename,col_headers,sheet,'G2:K2');
xlswrite(filename,row_headers,sheet,'A3');

for i = 1:4
    y = gold_ret(isbetween(ret_date,tlower(i),tupper(i)));
    y = y(y~=0);
    T = length(y);    
    [thetahat,thetastd,hhat,ml,mlstd] = SV_PL(y,T,nloop,burnin,M,cp_ml); %#ok<*ASGLU>
    switch i
        case 1
            xlRange = 'C3';
            xlRange_sd = 'H3';
        case 2
            xlRange = 'D3';
            xlRange_sd = 'I3';
        case 3
            xlRange = 'E3';
            xlRange_sd = 'J3';
        case 4
            xlRange = 'F3';
            xlRange_sd = 'K3';
    end
    xlswrite(filename,[thetahat;ml],sheet,xlRange);
    xlswrite(filename,[thetastd;mlstd],sheet,xlRange_sd);
end        

% run the code for the sp500_ret data set 

sheet = 'sp500_ret';
col_headers = {'All','2003-2005','2006-2009','2010-2013','2014-2017'};
row_headers = {'mu';'lambda';'kap';'mu_h';'phi_h';'omega2_h';'log marginal likelihood'};
xlswrite(filename,col_headers,sheet,'B2:F2');
xlswrite(filename,col_headers,sheet,'G2:K2');
xlswrite(filename,row_headers,sheet,'A3');

for i = 1:4
    y = sp500_ret(isbetween(ret_date,tlower(i),tupper(i)));
    y = y(y~=0);
    T = length(y);    
    [thetahat,thetastd,hhat,ml,mlstd] = SV_PL(y,T,nloop,burnin,M,cp_ml);
    switch i
        case 1
            xlRange = 'C3';
            xlRange_sd = 'H3';
        case 2
            xlRange = 'D3';
            xlRange_sd = 'I3';
        case 3
            xlRange = 'E3';
            xlRange_sd = 'J3';
        case 4
            xlRange = 'F3';
            xlRange_sd = 'K3';
    end
    xlswrite(filename,[thetahat;ml],sheet,xlRange);
    xlswrite(filename,[thetastd;mlstd],sheet,xlRange_sd);
end        

% run the code for the wti_ret data set 

sheet = 'wti_ret';
col_headers = {'All','2003-2005','2006-2009','2010-2013','2014-2017'};
row_headers = {'mu';'lambda';'kap';'mu_h';'phi_h';'omega2_h';'log marginal likelihood'};
xlswrite(filename,col_headers,sheet,'B2:F2');
xlswrite(filename,col_headers,sheet,'G2:K2');
xlswrite(filename,row_headers,sheet,'A3');

for i = 1:4
    y = wti_ret(isbetween(ret_date,tlower(i),tupper(i)));
    y = y(y~=0);
    T = length(y);    
    [thetahat,thetastd,hhat,ml,mlstd] = SV_PL(y,T,nloop,burnin,M,cp_ml);
    switch i
        case 1
            xlRange = 'C3';
            xlRange_sd = 'H3';
        case 2
            xlRange = 'D3';
            xlRange_sd = 'I3';
        case 3
            xlRange = 'E3';
            xlRange_sd = 'J3';
        case 4
            xlRange = 'F3';
            xlRange_sd = 'K3';
    end
    xlswrite(filename,[thetahat;ml],sheet,xlRange);
    xlswrite(filename,[thetastd;mlstd],sheet,xlRange_sd);
end    

%% run the estimation and calculate the marginal likelihood for the whole data set

y = eurusd_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat,ml,mlstd] = SV_PL(y,T,nloop,burnin,M,cp_ml);
xlswrite(filename,[thetahat;ml],'eurusd_ret','B3');
xlswrite(filename,[thetastd;mlstd],'eurusd_ret','G3');

y = gold_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat,ml,mlstd] = SV_PL(y,T,nloop,burnin,M,cp_ml);
xlswrite(filename,[thetahat;ml],'gold_ret','B3');
xlswrite(filename,[thetastd;mlstd],'gold_ret','G3');        

y = sp500_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat,ml,mlstd] = SV_PL(y,T,nloop,burnin,M,cp_ml);
xlswrite(filename,[thetahat;ml],'sp500_ret','B3');
xlswrite(filename,[thetastd;mlstd],'sp500_ret','G3');      

y = wti_ret;
y = y(y~=0);
T = length(y);    
[thetahat,thetastd,hhat,ml,mlstd] = SV_PL(y,T,nloop,burnin,M,cp_ml);
xlswrite(filename,[thetahat;ml],'wti_ret','B3');
xlswrite(filename,[thetastd;mlstd],'wti_ret','G3');     