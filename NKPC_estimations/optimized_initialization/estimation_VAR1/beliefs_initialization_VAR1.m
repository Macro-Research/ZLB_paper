clear;clc;close all;
dynare NKPC_REE_Estimation;
dataset=oo_.endo_simul;

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

alpha_init=[oo_.posterior_mode.parameters.a_1;
    oo_.posterior_mode.parameters.a_2;
    oo_.posterior_mode.parameters.a_3];

beta_init=[oo_.posterior_mode.parameters.b_11,oo_.posterior_mode.parameters.b_12,oo_.posterior_mode.parameters.b_13;
    oo_.posterior_mode.parameters.b_21,oo_.posterior_mode.parameters.b_22,oo_.posterior_mode.parameters.b_23;
    oo_.posterior_mode.parameters.b_31,oo_.posterior_mode.parameters.b_32,oo_.posterior_mode.parameters.b_33];

save VAR1_initial_beliefs.mat alpha_init beta_init rr_init;

