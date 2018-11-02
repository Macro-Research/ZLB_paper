clear;clc;close all;
dynare SW_Estimation_REE;
dataset=oo_.endo_simul;
clearvars -except dataset;
forward_indices=[3 5 6 7 9 10 11];
lagged_1=dataset(forward_indices,2:end-1)';
lagged_2=dataset(forward_indices,1:end-2)';
contemp=dataset(forward_indices,3:end)';

const=ones(length(lagged_1),1);
beta_init=zeros(length(forward_indices),1);

for jj=1:length(forward_indices)
    regr=([const lagged_1(:,jj) lagged_2(:,jj)]'*...
        [const lagged_1(:,jj) lagged_2(:,jj)])^(-1)*...
        ([const lagged_1(:,jj) lagged_2(:,jj)]'*contemp(:,jj));
    beta1_init(jj)=regr(2);
    beta2_init(jj)=regr(3);
    rr_init(:,:,jj)=([const lagged_1(:,jj) lagged_2(:,jj)]'*...
        [const lagged_1(:,jj) lagged_2(:,jj)])/length(const);
end

save AR1_initial_beliefs.mat beta1_init beta2_init rr_init;
