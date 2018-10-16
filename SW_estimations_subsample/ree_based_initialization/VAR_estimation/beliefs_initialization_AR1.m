clear;clc;close all;
dynare SW_Estimation_REE;
dataset=oo_.endo_simul;
clearvars -except dataset;
forward_indices=[3 5 6 7 9 10 11];
lagged=dataset(forward_indices,1:end-1)';
contemp=dataset(forward_indices,2:end)';
const=ones(length(lagged),1);
beta_init=zeros(length(forward_indices),1);

for jj=1:length(forward_indices)
    regr=([const lagged(:,jj)]'*[const lagged(:,jj)])^(-1)*([const lagged(:,jj)]'*contemp(:,jj));
    beta_init(jj)=regr(2);
    rr_init(:,:,jj)=([const lagged(:,jj)]'*[const lagged(:,jj)])/length(const);
end

save AR1_initial_beliefs.mat beta_init rr_init;
