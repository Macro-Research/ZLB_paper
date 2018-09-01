clear;clc;close all;
dynare NKPC_REE_Estimation;
dataset=oo_.endo_simul;
clearvars -except dataset;

yy_lag=dataset(1,1:end-1)';
yy=dataset(1,2:end)';
pinf_lag=dataset(2,1:end-1)';
pinf=dataset(2,2:end)';
cc=ones(length(yy_lag),1);

regr=([cc yy_lag]'*[cc yy_lag])^(-1)*([cc yy_lag]'*yy);
alpha_y=regr(1);beta_y=regr(2);
regr=([cc pinf_lag]'*[cc pinf_lag])^(-1)*([cc pinf_lag]'*pinf);
alpha_pinf=regr(1);beta_pinf=regr(2);

rr_y=([cc yy_lag]'*[cc yy_lag])/length(yy_lag);
rr_pinf=([cc pinf_lag]'*[cc pinf_lag])/length(pinf_lag);

save AR1_initial_beliefs.mat alpha_y beta_y rr_y...
                             alpha_pinf beta_pinf rr_pinf;