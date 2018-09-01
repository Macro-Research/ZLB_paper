clear;clc;close all;
dynare NKPC_REE_Estimation;
dataset=oo_.endo_simul;
clearvars -except dataset;

yy_lag=dataset(1,1:end-1)';
yy=dataset(1,2:end)';
pinf_lag=dataset(2,1:end-1)';
pinf=dataset(2,2:end)';
cc=ones(length(yy_lag),1);
rr=dataset(3,2:end)';
rr_lag=dataset(3,1:end-1)';
uu_y=dataset(4,2:end)';
uu_pinf=dataset(5,2:end)';

regressor=[cc yy_lag pinf_lag rr_lag];
regressand=[yy pinf rr];

regression=(regressor'*regressor)^(-1)*(regressor'*regressand);
rr_init=(regressor'*regressor)/length(yy);

alpha_init=regression(1,:)';
beta_init=regression(2:4,:)';

save VAR1_initial_beliefs.mat alpha_init beta_init rr_init;

