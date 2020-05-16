N = 20000;
Nth = 12000;

load('plsv_dataset.mat');

%%SV

param.SP500 = [0.080413514, -0.4013264304, 0.9831013794, sqrt(0.0380012659)];
param.Gold = [0.0440828572, 0.0003960336, 0.9823195436, sqrt(0.0185030461)];
param.EUR_USD = [-0.008995041, -1.1753016167, 0.9885813952, sqrt(0.0097347788)];
param.WTI = [0.0543043629, 1.3919383723, 0.9873095217, sqrt(0.0171572449)];

like.SV.SP500 = observed_log_likelihood(sp500_ret,param.SP500,N,Nth);
like.SV.Gold = observed_log_likelihood(gold_ret,param.Gold,N,Nth);
like.SV.EUR_USD = observed_log_likelihood(eurusd_ret,param.EUR_USD,N,Nth);
like.SV.WTI = observed_log_likelihood(wti_ret,param.WTI,N,Nth);

%%SV-L

source(paste(path, 'model_SV-L_Chan.R', sep=''))

param.SP500 = c(0.041356572, -0.457887952, 0.9709170328, sqrt(0.0715171372), -0.7409014934)
param.Gold = c(0.0535496, 0.0151439377, 0.984919829, sqrt(0.0172128649), 0.2620463213)
param.EUR_USD = c(-0.008098302, -1.1807590235, 0.9881346836, sqrt(0.0102442947), 0.0487043241)
param.WTI = c(0.0148771394,  1.3987825071, 0.988953773, sqrt(0.0150123074), -0.4601780304)

like.SVL.SP500 = observed.log.likelihood(SP500, param.SP500, N, Nth)
like.SVL.Gold = observed.log.likelihood(Gold, param.Gold, N, Nth)
like.SVL.EeurUR_USD = observed.log.likelihood(EUR_USD, param.EUR_USD, N, Nth)
like.SVL.WTI = observed.log.likelihood(WTI, param.WTI, N, Nth)

%%SV-M

source(paste(path, 'model_SV-M_Chan.R', sep=''))

param.SP500 = c(0.1049405937, -0.0503712106, -0.3952940375, 0.9808390164, sqrt(0.0438097516))
param.Gold = c(0.063708202, -0.0212184514, -0.0027243754, 0.9809495372, sqrt(0.019875346))
param.EUR_USD = c(-0.022481912, 0.0496130677, -1.1734698901, 0.9891856791, sqrt(0.0092361053))
param.WTI = c(0.102770065, -0.0145845099, 1.3873692491, 0.9872292542, sqrt(0.0174591279))

like.SVM.SP500 = observed.log.likelihood(SP500, param.SP500, N, Nth)
like.SVM.Gold = observed.log.likelihood(Gold, param.Gold, N, Nth)
like.SVM.EUR_USD = observed.log.likelihood(EUR_USD, param.EUR_USD, N, Nth)
like.SVM.WTI = observed.log.likelihood(WTI, param.WTI, N, Nth)

%%SV-LM (new SV)

source(paste(path, 'model_SV-LM_Chan.R', sep=''))

param.SP500 = c(0.0071677443, -0.1309236107, -0.1717601166, 0.9747810533, sqrt(0.0387617283))
param.Gold = c(0.0416149372, -0.0065975583, -0.001775019, 0.9831729308, sqrt(0.0174682997))
param.EUR_USD = c(-0.0006804484, 0.0055880468, -1.1798950506, 0.9886292105, sqrt(0.009645321))
param.WTI = c(0.0071677443, -0.1309236107, -0.1717601166, 0.9747810533, sqrt(0.0387617283))

like.SVLM.SP500 = observed.log.likelihood(SP500, param.SP500, N, Nth)
like.SVLM.Gold = observed.log.likelihood(Gold, param.Gold, N, Nth)
like.SVLM.EUR_USD = observed.log.likelihood(EUR_USD, param.EUR_USD, N, Nth)
like.SVLM.WTI = observed.log.likelihood(WTI, param.WTI, N, Nth)

%%result table

res = matrix(c(like.SV.SP500, like.SV.Gold, like.SV.EUR_USD, like.SV.WTI,
               like.SVL.SP500, like.SVL.Gold, like.SVL.EUR_USD, like.SVL.WTI,
               like.SVM.SP500, like.SVM.Gold, like.SVM.EUR_USD, like.SVM.WTI,
               like.SVLM.SP500, like.SVLM.Gold, like.SVLM.EUR_USD, like.SVLM.WTI), ncol=4, byrow=TRUE)
colnames(res) = c('S&P500', 'Gold', 'EUR/USD', 'WTI')
rownames(res) = c('SV', 'SV-L', 'SV-M', 'SV-LM')
likes = as.table(res)
print(likes)