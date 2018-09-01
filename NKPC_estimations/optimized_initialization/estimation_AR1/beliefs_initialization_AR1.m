clear;clc;close all;
dynare NKPC_REE_Estimation;
dataset=oo_.endo_simul;
alpha_init=[oo_.posterior_mode.parameters.alpha_y oo_.posterior_mode.parameters.alpha_pinf];
beta_init=[oo_.posterior_mode.parameters.beta_y oo_.posterior_mode.parameters.beta_pinf];

yy_lag=dataset(1,1:end-1)';
yy=dataset(1,2:end)';
pinf_lag=dataset(2,1:end-1)';
pinf=dataset(2,2:end)';
cc=ones(length(yy_lag),1);

rr_y=([cc yy_lag]'*[cc yy_lag])/length(yy_lag);
rr_pinf=([cc pinf_lag]'*[cc pinf_lag])/length(pinf_lag);

save AR1_initial_beliefs.mat alpha_init beta_init rr_y rr_pinf;

