clear;clc;close all;
%dynare SW_Estimation_REE;
%save('initial_ree_estimation.mat');
load('initial_ree_estimation.mat');
%dataset=oo_.endo_simul;
% clearvars -except dataset;
% forward_indices=[3 5 6 7 9 10 11];
% shock_indices=[14 15 16 17 18 19 20];
% backward_indices=[6 7 8 10 12 13];
% lagged=dataset(backward_indices,1:end-1)';
% contemp=dataset(forward_indices,2:end)';
% shocks=dataset(shock_indices,2:end)';
burn_in=50000;
lagged=[ c inve y pinf w r kp];
lagged=lagged(1+burn_in:end-1,:);
contemp=[rk q c inve lab pinf w];
contemp=contemp(2+burn_in:end,:);
shocks=[eps_a eps_b eps_g eps_i eps_r eps_p eps_w];
shocks=shocks(2+burn_in:end,:);
const=ones(length(lagged),1);
clearvars -except lagged contemp shocks const;
regressor=[const lagged shocks];
regressand=contemp;

regression=(regressor'*regressor)^(-1)*(regressor'*regressand);
rr_init=(regressor'*regressor)/length(contemp);

alpha_init=regression(1,:)';
beta_init=regression(2:8,:)';
cc_init=regression(9:end,:)';

save initial_beliefs_msv alpha_init beta_init cc_init rr_init;


